from dask_jobqueue import htcondor
import fnmatch
import shutil
import copy
from time import sleep, time
from hashlib import md5
import numpy as np
from itertools import repeat
from distributed.comm import get_address_host
from coffea.processor.executor import _compression_wrapper, _reduce, WorkItem
from coffea.processor.accumulator import AccumulatorABC
from socket import getfqdn
import lz4.frame as lz4f
import cloudpickle
from distributed import WorkerPlugin, get_worker, as_completed
from coffea.processor.dask import ColumnCache
from zict import LRU
from threading import Lock
from collections import defaultdict, Counter
from distributed.diagnostics.progressbar import TextProgressBar
from distributed.diagnostics.plugin import SchedulerPlugin
from distributed.scheduler import KilledWorker
from distributed.core import Status
from tornado.ioloop import PeriodicCallback
from operator import attrgetter
from dask.base import tokenize
from tqdm.auto import tqdm
from utils.util import iter_chunks
from functools import cached_property
from subprocess import check_output
import os


def make_target(base, values):
    target = {}
    for pat, val in values.items():
        for key in fnmatch.filter(base, pat):
            target[key] = val
    for key, val in target.items():
        if val is None:
            del target[key]
    return target


class Reschedule(SchedulerPlugin):
    sat_tresh = 60
    sat_fract = 0.5

    def __init__(self, scheduler):
        self.scheduler = scheduler

        self._pc = pc = PeriodicCallback(callback=self.balance, callback_time=5000)
        self.scheduler.periodic_callbacks["rescheduling"] = pc
        self.scheduler.plugins.append(self)
        self.scheduler.extensions["rescheduling"] = self

    def teardown(self):
        self._pc.stop()

    def balance(self):
        sched = self.scheduler

        if not sched.idle:
            return

        evict = {}
        for ws in sched.saturated:
            proc = ws.processing
            tot = ws.occupancy
            tresh = max(tot * self.sat_fract, self.sat_tresh)
            for ts in sorted(proc.keys(), key=attrgetter("priority"), reverse=True):
                tot -= proc[ts]
                if tresh < tot:
                    evict[ts] = tot
                else:
                    break

        evict = sorted(evict.keys(), key=evict.get, reverse=True)
        evict = self.filter_evictions(evict)

        if evict:
            sched.transitions({ts.key: "released" for ts in evict[::-1]})

    def filter_evictions(self, evict):
        return evict


class TaskWorkerDistance:
    hash = md5
    digest_size = hash().digest_size
    host_replica = 8

    def __init__(self, scheduler, factor={}, ignore_prefix=()):
        self.scheduler = scheduler
        self.factor = factor
        self.ignore_prefix = ignore_prefix
        self._lock = Lock()
        self._tkey2i = {}
        self._host2j = {}
        self._thash = self._thash0 = np.zeros((0, self.digest_size), dtype=np.int8)
        self._hhash = self._hhash0 = np.zeros(
            (0, self.host_replica, self.digest_size), dtype=np.int8
        )
        self._hfact = self._hfact0 = np.zeros((0,), dtype=float)
        self._dist = np.zeros((0, 0))

    def get(self, ts, ws):
        if ts.prefix_key in self.ignore_prefix:
            return 0
        tkey = ts.key
        host = self.ws2host(ws)
        i = self._tkey2i.get(tkey, -1)
        j = self._host2j.get(host, -1)
        if i < 0 or j < 0:
            self.update()
            i = self._tkey2i[tkey]
            j = self._host2j[host]
        return self._dist.item(i, j)

    def update(self):
        with self._lock:
            i, j = Mi, Mj = self._dist.shape

            tkey2i = {}
            thash = []
            for tkey, ts in self.scheduler.tasks.items():
                if tkey in self._tkey2i:
                    continue
                if ts.prefix_key in self.ignore_prefix:
                    continue
                h = self.ts2hash(ts)
                if len(h) != self.digest_size:
                    continue
                thash.append(h)
                tkey2i[tkey] = i
                i += 1
            thash = np.array(thash) if thash else self._thash0

            host2j = {}
            hhash = []
            hfact = []
            for host in set(map(self.ws2host, self.scheduler.workers.values())):
                if host in self._host2j:
                    continue
                hhash.append(self.host2hash(host))
                hfact.append(self.host2factor(host))
                host2j[host] = j
                j += 1
            hhash = np.array(hhash) if hhash else self._hhash0
            hfact = np.array(hfact) if hfact else self._hfact0

            if (Mi, Mj) == (i, j):
                return

            dist = np.empty((i, j))

            for si, t in (
                np.s_[0:Mi, self._thash],
                np.s_[Mi:i, thash],
            ):
                for sj, h, f in (
                    np.s_[0:Mj, self._hhash, self._hfact],
                    np.s_[Mj:j, hhash, hfact],
                ):
                    if (si, sj) == np.s_[0:Mi, 0:Mj]:
                        dist[si, sj] = self._dist
                    else:
                        # t.shape = (n_tasks, digest_size)
                        # h.shape = (n_hosts, host_replica, digest_size)
                        # f.shape = (n_hosts, )
                        # d.shape = (n_tasks, n_hosts, host_replica, digest_size)
                        d = t[:, None, None, :] - h[None, ...]
                        dist[si, sj] = np.abs(d).sum(axis=-1).min(axis=-1) * f[None, :]

            self._dist = dist
            self._thash = np.concatenate((self._thash, thash), axis=0)
            self._hhash = np.concatenate((self._hhash, hhash), axis=0)
            self._hfact = np.concatenate((self._hfact, hfact), axis=0)
            self._tkey2i.update(tkey2i)
            self._host2j.update(host2j)

    def ts2hash(self, ts):
        return np.frombuffer(bytes.fromhex(ts.key.rsplit("-", 1)[-1]), dtype=np.int8)
        # return np.frombuffer(self.hash(ts.key.encode("utf8")).digest(), dtype=np.int8)

    def host2hash(self, host):
        s = self.hash(host.encode("utf8"))
        h = np.empty((self.host_replica, s.digest_size), dtype=np.int8)
        for i in range(len(h)):
            h[i] = np.frombuffer(s.digest(), dtype=np.int8)
            if i:
                s.update(s.digest())
        return h

    def ws2host(self, ws):
        return get_address_host(ws.address)

    def host2factor(self, host):
        return (1 / self.factor.get(host, 1)) ** (1 / md5().digest_size)


class AffineRescheduling(Reschedule):
    sat_tresh = 30
    idle_tresh = 10
    fill_tresh = 3 * idle_tresh
    fill_max = 5
    fill_min = 2

    def __init__(self, scheduler, **kwargs):
        super().__init__(scheduler)
        self.distance = TaskWorkerDistance(scheduler, **kwargs)
        self.evict = ()

        self.orig_worker_objective = scheduler.worker_objective
        self.orig_check_idle_saturated = scheduler.check_idle_saturated
        self.scheduler.worker_objective = self.worker_objective
        self.scheduler.check_idle_saturated = self.check_idle_saturated

    def teardown(self):
        super().teardown()
        self.scheduler.worker_objective = self.orig_worker_objective
        self.scheduler.check_idle_saturated = self.orig_check_idle_saturated

    def check_idle_saturated(self, ws, occ=None):
        sched = self.scheduler
        ntot = sched.total_nthreads
        if ntot == 0 or ws.nthreads == 0 or ws.status == Status.closed or not ws.memory_limit:
            return
        if occ is None:
            occ = ws.occupancy
        nc = ws.nthreads
        p = len(ws.processing)

        if p <= nc or occ < self.idle_tresh:
            sched.idle.add(ws)
            sched.saturated.discard(ws)
        else:
            sched.idle.discard(ws)

            avg = sched.total_occupancy / ntot
            pending = occ * (p - nc) / p / nc

            if p > nc and pending > max(avg, self.sat_tresh):
                sched.saturated.add(ws)
            else:
                sched.saturated.discard(ws)

    def worker_objective(self, ts, ws):
        start_time = self.orig_worker_objective(ts, ws)[-2]

        dist = self.distance.get(ts, ws) + 1
        nproc = len(ws.processing)
        reloc = (
            ws.memory_limit  # ignore unbound workers
            and (  # need some form of filling?
                nproc < self.fill_min * ws.nthreads  # slot
                or ws.occupancy < self.fill_tresh  # occupancy
            )
            and nproc < self.fill_max * ws.nthreads  # dont over-fill slots
            and ts.key in self.evict  # task was actually marked for rescheduling
        )

        return (
            -len(ws.actors) if ts.actor else 0,
            -1 / dist if reloc else 0,
            start_time if reloc else 0,
            dist,
            start_time,
            ws.nbytes,
        )

    def filter_evictions(self, evict):
        # free: occupancy, minmum, slots
        ofree = sfree = mfree = 0
        for ws in self.scheduler.workers.values():
            if not ws.memory_limit:
                continue
            nproc = len(ws.processing)
            ofree += max(0, self.fill_tresh - ws.occupancy) * ws.nthreads
            sfree += max(0, self.fill_max * ws.nthreads - nproc)
            mfree += max(0, self.fill_min * ws.nthreads - nproc)
        evict = evict[: max(mfree, sfree)]
        for i, ts in enumerate(evict):
            ofree -= self.scheduler.get_task_duration(ts)
            if ofree < 0:
                evict = evict[: max(mfree, i + 1)]
                break
        self.evict = set(map(attrgetter("key"), evict))
        return evict


from distributed.scheduler import TaskPrefix

da_decay = 0.02


@property
def duration_average(self):
    return self._duration_average


@duration_average.setter
def duration_average(self, value):
    last = self._duration_average
    if value is not None and last is not None:
        # undo: avg_duration = 0.5 * old_duration + 0.5 * new_duration
        # new = 2 * value - last
        # value = (1 - da_decay) * last + da_decay * new
        value = (1 - 2 * da_decay) * last + da_decay * 2 * value
        # value = da_decay * (last + 2 * value)
    self._duration_average = value


TaskPrefix.duration_average = duration_average
TaskPrefix._duration_average = None


def wi_tokenize(self):
    return self.treename, self.entrystart, self.entrystop, self.fileuuid


WorkItem.__dask_tokenize__ = wi_tokenize


class HTCondorJob(htcondor.HTCondorJob):
    sg_target = "openports"
    executable = shutil.which("sg")

    def job_script(self):
        """ Construct a job submission script """
        quoted_arguments = htcondor.quote_arguments([self.sg_target, "-c", self._command_template])
        quoted_environment = htcondor.quote_environment(self.env_dict)
        job_header_lines = "\n".join("%s = %s" % (k, v) for k, v in self.job_header_dict.items())
        return self._script_template % {
            "shebang": self.shebang,
            "job_header": job_header_lines,
            "quoted_environment": quoted_environment,
            "quoted_arguments": quoted_arguments,
            "executable": self.executable,
        }


class HTCondorCluster(htcondor.HTCondorCluster):
    job_cls = HTCondorJob
    machine_sep = "+"

    def __init__(self, *args, maintain=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.pc_maintain = PeriodicCallback(self._correct_state, callback_time=30000)
        self.maintain = maintain

    @property
    def maintain(self):
        return self.pc_maintain.is_running()

    @maintain.setter
    def maintain(self, val):
        if bool(val) != self.maintain:
            if val:
                self.loop.add_callback(self.pc_maintain.start)
            else:
                self.pc_maintain.stop()

    @classmethod
    def worker2machine(cls, worker):
        if not isinstance(worker, str):
            worker = worker.name
        parts = worker.split(cls.machine_sep)
        return parts[0] if len(parts) == 2 else None

    def _machine_status(self, query="1"):
        ret = Counter()
        for line in (
            check_output(
                [
                    "condor_status",
                    "-constraint",
                    'SlotType!="Dynamic"',
                    "-af",
                    "Machine",
                    query,
                ]
            )
            .decode()
            .splitlines()
        ):
            key, value = line.split()
            ret[key] += int(value)
        return ret

    @property
    def machines_known(self):
        return self._machine_status("TotalSlotCpus")

    @property
    def machines_free(self):
        return self._machine_status("Cpus")

    @property
    def machines(self):
        ret = Counter(map(self.worker2machine, self.workers.keys()))
        ret.pop(None, None)
        return ret

    @property
    def machines_live(self):
        ret = Counter(map(self.worker2machine, self.scheduler.workers.values()))
        ret.pop(None, None)
        return ret

    @property
    def machines_pending(self):
        ml = self.machines_live
        return {m for m in self.machines.keys() if m not in ml}

    @property
    def workers_live(self):
        return set(w.name for w in self.scheduler.workers.values())

    @property
    def workers_pending(self):
        return set(self.workers.keys()) - self.workers_live

    @property
    def workers_machine_pending(self):
        mp = self.machines_pending
        return set(w for w in self.workers_pending if self.worker2machine(w) in mp)

    def status_table(self):
        from rich.table import Table, Column
        from rich import box

        info = (
            self.machines,
            self.machines_live,
            self.machines_free,
            self.machines_known,
        )

        b, l, f, t = [sum(i.values()) for i in info]

        def pct(a, b, digits=0):
            return "%.*f" % (
                max(0, digits),
                round(100 * a / (b or 1), digits),
            )

        tab = Table(
            Column("Machine \\ CPUs", "Total"),
            Column("live %", pct(l, b), justify="right"),
            Column("live", str(l), justify="right"),
            Column("booked", str(b), justify="right"),
            Column("total", str(t), justify="right"),
            Column("free", str(f), justify="right"),
            Column("free %", pct(f, t), justify="right"),
            show_footer=True,
            box=box.SIMPLE,
        )
        for m in sorted(info[-1]):
            b, l, f, t = [i.get(m, 0) for i in info]
            tab.add_row(m, pct(l, b), str(l), str(b), str(t), str(f), pct(f, t))

        return tab

    def scale_machines(
        self,
        machines,
        wait_machines=True,
        timeout=np.inf,
        poll=1,
        stop_pending=False,
    ):
        # remove excess
        max_target = make_target(self.machines, machines)
        for name in self.worker_spec.keys():
            if self.machine_sep not in str(name):
                continue
            machine, num = name.split(self.machine_sep)
            if machine not in max_target:
                continue
            if max_target[machine] <= num:
                del self.worker_spec[name]

        # spawn new
        name_new = []
        for machine, num in make_target(self.machines_known, machines).items():
            for i in range(num):
                name = self.machine_sep.join((machine, str(i)))
                assert name not in self.worker_spec, name
                self.worker_spec[name] = self.new_machine_spec(machine)
                name_new.append(name)

        # sync internal sate
        self.sync(self._correct_state)

        if not poll:
            return

        endtime = time() + timeout
        while time() < endtime and (
            self.machines_pending if wait_machines else self.workers_pending
        ):
            sleep(poll)

        if stop_pending:
            wp = self.workers_pending
            if stop_pending == "machine":
                wp = self.workers_machine_pending
            for w in wp:
                del self.worker_spec[w]

    def new_machine_spec(self, machine):
        spec = copy.deepcopy(self.new_spec)
        req = 'Machine=="%s"' % machine
        job_extra = spec.setdefault("options", {}).setdefault("job_extra", {})
        if "Requirements" in job_extra:
            req = "(%s)&&(%s)" % (req, job_extra["Requirements"])
        job_extra["Requirements"] = req
        return spec


class WorkerLookup:
    def __init__(self, client):
        self.client = client

    @cached_property
    def wkey2fqdn(self):
        return self.client.run(getfqdn)

    @cached_property
    def host2fqdn(self):
        return {get_address_host(wkey): fqdn for wkey, fqdn in self.wkey2fqdn.items()}

    def fqdn2hosts(self, *fqdns, skip=()):
        return set(
            host
            for host, fqdn in self.host2fqdn.items()
            if (fqdns == (all,) or any(fnmatch.fnmatch(fqdn, pat) for pat in fqdns))
            and not any(fnmatch.fnmatch(fqdn, spat) for spat in skip)
        )


def dask_executor(
    items,
    function,
    accumulator,
    client,
    treereduction=20,
    status=True,
    compression=1,
    function_name=None,
    direct_heavy=None,
    pure=True,
    retries=3,
    priority=0,
    allow_other_workers=False,
    workers=None,
    merge_workers=None,
    live_callback=None,
    **kwargs,
):
    """Execute using dask futures
    Parameters
    ----------
        items : list
            List of input arguments
        function : callable
            A function to be called on each input, which returns an accumulator instance
        accumulator : AccumulatorABC
            An accumulator to collect the output of the function
        client : distributed.client.Client
            A dask distributed client instance
        treereduction : int, optional
            Tree reduction factor for output accumulators (default: 20)
        status : bool, optional
            If true (default), enable progress bar
        compression : int, optional
            Compress accumulator outputs in flight with LZ4, at level specified (default 1)
            Set to ``None`` for no compression.
        priority : int, optional
            Task priority, default 0
        retries : int, optional
            Number of retries for failed tasks (default: 3)
        function_name : str, optional
            Name of the function being passed
    """
    if len(items) == 0:
        return accumulator
    reducer = _reduce()

    items_orig = items
    keys = ["Processor-%s" % tokenize(item) for item in items]
    if pure and len(cloudpickle.dumps(function)) > 1000:
        func = client.scatter(function, broadcast=True, hash=False, direct=direct_heavy)
        fkey = "-%s-" % tokenize(func.key)[:8]
        keys = [key.replace("-", fkey, 1) for key in keys]
        items = list(zip(repeat(func), items))
        function = lambda fa: fa[0](fa[1])
    assert len(set(keys)) == len(keys)
    key2item = dict(zip(keys, items_orig))

    if compression is not None:
        function = _compression_wrapper(compression, function, name=function_name)
        reducer = _compression_wrapper(compression, reducer)

    work = client.map(
        function,
        items,
        key=keys,
        pure=pure,
        retries=retries,
        priority=priority,
        workers=workers,
        allow_other_workers=allow_other_workers,
    )

    if callable(live_callback):
        ac = as_completed(work, with_results=True)
        total = sum(
            wi.entrystop - wi.entrystart for wi in items_orig
        )  # if this fails, use a different pre_executor ie futures_executor
        with tqdm(total=total, unit="event", unit_scale=True, disable=not status) as prog:
            for batch in ac.batches():
                for future, result in batch:
                    if future.status == "error":
                        typ, exc, tb = result
                        raise exc.with_traceback(tb)
                    accumulator += _maybe_decompress(result)
                prog.update(
                    accumulator["metrics"]["entries"].value - prog.n
                )  # if this fails, use savemetrics
                pf = live_callback(accumulator)
                if isinstance(pf, dict):
                    prog.set_postfix(pf)
        return accumulator

    if status is True:
        status = function_name or "Processing"

    pc = PeriodicCallback(lambda: client._replicate(func), callback_time=30000)
    client.loop.add_callback(pc.start)

    ac = as_completed(work)

    def logic(ac):
        seen = set()
        # operate in batches
        for batch in ac.batches():
            bad = []
            good = []
            for future in batch:
                # did this task kill it's worker?
                if (
                    future.status == "error"
                    and future.key not in seen
                    and isinstance(future.exception(), KilledWorker)
                ):
                    seen.add(future.key)
                    bad.append(future)
                else:
                    good.append(future)
            if good:
                # move output asap to merge_workers
                if merge_workers and not allow_other_workers:
                    moved = client.map(
                        lambda x: x,
                        good,
                        pure=True,
                        retries=retries,
                        priority=priority + 1,
                        workers=merge_workers,
                        allow_other_workers=False,
                    )
                    # release explicit retention
                    for future in good:
                        future.release()
                    good = moved
                yield from good
            if bad:
                # rescheduler killer tasks on merge_workers
                for f in bad:
                    item = key2item[f.key]
                    print(
                        "KilledWorker by key=%s, offset=%d, filename=%s"
                        % (f.key, item.entrystart, item.filename)
                    )
                client._send_to_scheduler(
                    {
                        "op": "update-graph",
                        "tasks": {},
                        "keys": [f.key for f in bad],
                        "restrictions": {f.key: list(merge_workers or ()) for f in bad},
                        "client": client.id,
                    }
                )
                client.retry(bad)
                ac.update(bad)

    tq = tqdm(
        logic(ac),
        total=ac.count(),
        unit="task",
        desc=status or "",
        disable=not status,
        smoothing=1 / 20,
    )
    for chunk in iter_chunks(tq, treereduction):
        for future in chunk:
            if future.status == "error":
                pc.stop()
                future.result()
        work = client.submit(
            reducer,
            chunk,
            pure=True,
            retries=retries,
            priority=priority + 2,
            workers=merge_workers,
            allow_other_workers=allow_other_workers if merge_workers else None,
        )
        # release explicit retention
        for future in chunk:
            future.release()
        if not ac.is_empty():
            ac.add(work)
            tq.total += 1
    tq.close()

    accumulator += _maybe_decompress(work.result())
    pc.stop()
    return accumulator


def _maybe_decompress(item):
    if isinstance(item, bytes):
        item = cloudpickle.loads(lz4f.decompress(item))
    if isinstance(item, AccumulatorABC):
        return item
    raise ValueError(
        "Executors can only reduce accumulators or LZ4-compressed pickled accumulators"
    )


class TQDMProgressBar(TextProgressBar):
    def __init__(self, *args, **kwargs):
        tqargs = dict(kwargs.pop("tqdm", ()))
        tqargs.setdefault("ncols", kwargs.get("width", None))
        tqargs.setdefault("unit", "task")
        self.tq = tqdm(**tqargs)
        kwargs["width"] = self.tq.ncols
        super().__init__(*args, **kwargs)

    def _draw_bar(self, remaining, all, **kwargs):
        self.tq.total = all
        self.tq.update(all - remaining - self.tq.n)

    def _draw_stop(self, **kwargs):
        self.tq.close()


class MEMCache(ColumnCache):
    def __init__(self, maxmem=5e8):
        self._maxmem = maxmem

    def setup(self, worker):
        self.cache = LRU(n=self._maxmem, d={}, weight=lambda k, v: v.nbytes)
        self.lock = Lock()
        self.hits = 0
        self.misses = 0


class ConfigureXRootD(WorkerPlugin):
    name = "user_proxy"

    def __init__(self, proxy_file=None):
        """
        If proxy_file is None, look for it in default location
        """
        if not proxy_file:
            file = os.environ.get("X509_USER_PROXY", "/tmp/x509up_u%d" % os.getuid())
        else:
            file = proxy_file
        self._proxy = open(file, "rb").read()

    def setup(self, worker):
        self._location = os.path.join(worker.local_directory, "userproxy")
        with open(self._location, "wb") as fout:
            fout.write(self._proxy)
        os.environ["X509_USER_PROXY"] = self._location
        os.environ["XRD_CONNECTIONWINDOW"] = "10"
        os.environ["XRD_STREAMTIMEOUT"] = "10"
        os.environ["XRD_TIMEOUTRESOLUTION"] = "2"
        os.environ["XRD_WORKERTHREADS"] = "4"
        os.environ["XRD_REQUESTTIMEOUT"] = "60"

    def teardown(self, worker):
        os.remove(self._location)
        del os.environ["X509_USER_PROXY"]


def register_plugins(client, add=defaultdict(dict)):
    """
    Usage:
        plugins = {
                    "MEMCache": {"maxmem": 5e8},
                    "ConfigureXRootD": {"proxy_file": None}
                  }
        register_plugins(client, add=plugins)
    """
    plugins = set()
    for p in client.run(lambda: set(get_worker().plugins)).values():
        plugins |= p
    for name, opts in add.items():
        plugin = globals()[name]
        if plugin.name not in plugins:
            client.register_worker_plugin(plugin(**opts))
