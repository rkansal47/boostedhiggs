from law.util import flatten
import os
from uproot_methods.classes.TH1 import from_numpy
from utils.uproot import bulk_write
from tqdm.auto import tqdm
from utils.shape import Shape
from utils.export import export1d_hist
from enum import IntEnum
from tools.evil import pin
from functools import cached_property
import uproot
import numpy as np
from abc import abstractmethod
from operator import itemgetter
from numbers import Number
from collections import OrderedDict
import fnmatch
from rich.console import Console


class ProcessesMixin:
    """
    @unique
    class Processes(IntEnum):

        # Signal Processes
        # They have to be <= 0
        ggHH_kl_0_kt_1_2B2VTo2L2Nu = -2
        ggHH_kl_1_kt_1_2B2VTo2L2Nu = -1
        ggHH_kl_5_kt_1_2B2VTo2L2Nu = 0

        # Background Processes
        # They have to be > 0
        ggH = 1
        VBFH = 2
        WH = 3
        ZH = 4
        dy = 5


    class StatModel(Datacard):
        processes = Processes

        def __init__(self, *args, **kwargs):
            print(self.signal_processes)
            >>> ["ggHH_kl_0_kt_1_2B2VTo2L2Nu", "ggHH_kl_1_kt_1_2B2VTo2L2Nu", "ggHH_kl_5_kt_1_2B2VTo2L2Nu"]
    """

    _processes = None

    @property
    def _sorted(self):
        return {k: v for k, v in sorted(self._processes.__members__.items(), key=lambda kv: kv[1])}

    @property
    def background_processes(self):
        return [n for n, v in self._sorted.items() if v > 0]

    @property
    def signal_processes(self):
        return [n for n, v in self._sorted.items() if v <= 0]

    @property
    def processes(self):
        return list(self._sorted.keys())

    @property
    def processes_ids(self):
        return list(self._sorted.values())


def auto_repr(cls):
    def repr(self):
        return f"{self.__class__.__name__}({', '.join([f'{k}={v}' for k, v in self.__dict__.items() if not k.startswith('_')])})"

    cls.__repr__ = repr
    return cls


class DIR(IntEnum):
    DOWN = 0
    UP = 1


@auto_repr
class Datacard(ProcessesMixin):
    """Example usage:

    @unique
    class Processes(IntEnum):
        # all processes can be accessed with `self.processes`

        # Signal Processes
        # They have to be <= 0
        # can be accessed with `self.signal_processes`
        ggHH_kl_0_kt_1_2B2VTo2L2Nu = -2
        ggHH_kl_1_kt_1_2B2VTo2L2Nu = -1
        ggHH_kl_5_kt_1_2B2VTo2L2Nu = 0

        # Background Processes
        # They have to be > 0
        # can be accessed with `self.background_processes`
        ggH = 1
        VBFH = 2
        WH = 3
        ZH = 4
        dy = 5
        rare = 6
        st = 7
        tt = 8
        ttV = 9
        ttVH = 10
        ttVV = 11
        vv = 12
        vvv = 13
        wjets = 14


    class StatModel(Datacard):
        _processes = Processes

        # add custom lines to the end of the datacard
        @property
        def custom_lines(self):
            return ["* autoMCStats 10 0 1"]

        def build_systematics(self):

            # add many systs with individual strengths
            self.add_systematics(
                names={
                    "CMS_lumi_13TeV_2017": 1.02,
                    # correlated lumi
                    "CMS_lumi_13TeV_XYFact": 1.008,
                    "CMS_lumi_13TeV_LScale": 1.003,
                    "CMS_lumi_13TeV_BBDefl": 1.004,
                    "CMS_lumi_13TeV_DynBeta": 1.005,
                    "CMS_lumi_13TeV_CurrCalib": 1.003,
                    "CMS_lumi_13TeV_Ghosts": 1.001,
                },
                type="lnN",
                processes=self.processes,
            )

            # or like this
            self.add_systematics(
                names=["lnNSys4", "lnNSys5"],
                type="lnN",
                strengths=[1.1, 1.2],
                processes=self.background_processes,
            )

            or many with the same strength
            self.add_systematics(
                names=["lnNSys6", "lnNSys7"],
                type="lnN",
                strengths=1.1,
                processes=self.background_processes,
            )

            or a single nuisance
            self.add_systematics(
                names="lnNSys8",
                type="lnN",
                strengths=1.1,
                processes=["DY"],
            )

            # also wildcards
            self.add_systematics(
                names="CMS_JES*",
                type="shape",
                strength=1.0,
                processes=self.processes,
            )

    Run by hand:

    card = StatModel(
        analysis="bbww_dl",
        category="mumu_diHiggs",
        variable="dnn_score_max",
        hists_path="/net/scratch/cms/dihiggs/store/bbww_dl/Run2_pp_13TeV_2017/Rebin/dev1/rebinned_hists.root",
        config_inst=...
    )
    card.dump("./")
    """

    sep = os.linesep
    dashes = 80 * "-"

    def __init__(self, analysis, category, variable, hists, config_inst, year):
        if not self._processes:
            raise Exception("Processes have to be defined!")
        if not issubclass(self._processes, IntEnum):
            raise Exception(f"{self.processes} has to be subclassed from <IntEnum>")
        # e.g.: category = "ee_2b"
        channel = category.split("_")[0]
        _nuisances_names = []
        _valid_processes = []
        _nuisances_renaming = {}
        _newhists = {}
        console = Console()
        pin(locals())

    @classmethod
    def requires(cls, task):
        return task.base_requires()

    @classmethod
    def rebin(cls, variable, h):
        return {c: h[:, c, :, :] for c in h.axes["category"]}

    def rate(self, process, systematic="nominal"):
        h = self.hists[
            process,
            systematic,
            :,
        ]
        arr = h.view(flow=True)["value"]
        # we remove negative bins later anyhow
        return np.sum(arr[arr >= 0.0])

    @property
    def observation(self):
        return -1

    @property
    def bin(self):
        return self.category

    @property
    def available_systematics(self):
        return [*self.hists.axes["systematic"]]

    @property
    def available_processes(self):
        return [*self.hists.axes["process"]]

    @cached_property
    def spaces(self):
        return (
            max(
                *[
                    len(p)
                    for p in self.processes
                    + [
                        "data_obs",
                    ]
                ],
                len(self.bin),
                *[len(" ".join(tpl)) for tpl in self._nuisances_names],
            )
            + 5
        )

    @property
    def nominal_pattern(self):
        return Shape.ch_template_nom.format(variable=str(self.variable)).replace("$BIN", self.bin)

    @property
    def sys_pattern(self):
        return Shape.ch_template_sys.format(variable=str(self.variable)).replace("$BIN", self.bin)

    @property
    def custom_lines(self):
        return []

    def nuisance_idx(self, name):
        return list(map(itemgetter(0), self._nuisances_names)).index(name)

    def add_systematics(self, names, type, processes, **kwargs):
        up, down = kwargs.pop("suffix", ("Up", "Down"))
        # remove non-valid processes
        if not (set(processes) <= set(self.available_processes)):
            self.console.log(
                f"Skip unknown process(es): {set(processes) - set(self.available_processes)} for nuisance(s) {names}."
            )

        processes = [p for p in processes if p in self.vps]
        if isinstance(names, str):
            rpl = lambda x: x.replace(up, "").replace(down, "")
            names = fnmatch.filter(set(map(rpl, self.available_systematics)), names)
        if isinstance(names, dict):
            names, strengths = list(names.keys()), list(names.values())
        strength = kwargs.pop("strength", None)
        if strength:
            if isinstance(strength, Number):
                strengths = len(names) * (strength,)
            else:
                strengths = strength
        # at this point names has to be a valid list
        assert isinstance(names, (list, tuple))
        assert len(names) == len(strengths)
        for name, strength in zip(names, strengths):
            if isinstance(strength, Number):
                strength = 2 * (strength,)
            assert len(strength) == 2
            # first check if systematic is empty, only for shape
            if type == "shape":
                name_up = name + up
                name_down = name + down
                # verify that the shape exists
                assert name_up in self.available_systematics, name_up
                assert name_down in self.available_systematics, name_down
                # track malformed shift namings for combine...
                if up != "Up":
                    self._nuisances_renaming[name_up] = name + "Up"
                if down != "Down":
                    self._nuisances_renaming[name_down] = name + "Down"
                to_be_removed = []
                for p in processes:
                    if (self.rate(p, name_up) <= 0.0) or (self.rate(p, name_down) <= 0.0):
                        self.console.log(
                            f"Remove Systematic (NULL Yield): {name} for process {p} (up: {self.rate(p, name_up)}, down: {self.rate(p, name_down)})."
                        )
                        to_be_removed.append(p)
                processes = list(set(processes) - set(to_be_removed))
            if not processes:
                continue
            self._nuisances_names.append((name, type))
            # get index:
            nuisance_idx = self.nuisance_idx(name)
            process_idx = [self.vps.index(p) for p in processes]
            self._nuisances_arr = np.concatenate(
                (self._nuisances_arr, np.zeros((1, self._nuisances_arr.shape[1], 2))), axis=0
            )
            self._nuisances_arr[nuisance_idx, process_idx, DIR.DOWN] = strength[DIR.DOWN]
            self._nuisances_arr[nuisance_idx, process_idx, DIR.UP] = strength[DIR.UP]

    @abstractmethod
    def build_systematics(self):
        pass

    @property
    def data_name(self):
        # corresponds to the process name of data in the hist.Hist
        return "data"

    def check_processes(self):
        i = 0
        for p, id in zip(self.processes, self.processes_ids):
            # id > 0, we don't want to remove signal
            if (id > 0) and ((r := self.rate(p)) <= 0.0):
                self.console.log(f"Remove: process {p} has negative or no contribution (rate: {r})")
                i += 1
            else:
                self._valid_processes.append((p, id - i))

    @property
    def vps(self):
        # valid processes
        return [tpl[0] for tpl in self._valid_processes]

    @property
    def vids(self):
        # valid processes ids
        return [tpl[1] for tpl in self._valid_processes]

    def fix_histograms(self):
        # relevant histograms:
        templates = []
        with tqdm(
            total=len(self.vps) * len(self.available_systematics),
            unit="hist",
            desc="Fix histograms",
        ) as pbar:
            for p in self.vps:
                for s in self.available_systematics:
                    # skip if the yield is <= 0?
                    syst = self._nuisances_renaming.get(s, s)
                    h = export1d_hist(self.hists[p, s, :].project(self.variable))
                    w = np.array(h)
                    w2 = h.allvariances
                    # fix negative bins
                    bad = w <= 0.0
                    w[bad] = 0.0
                    w2[bad] = 0.0
                    # make new hist
                    h[:] = w.tolist()
                    h._fSumw2 = w2
                    h._fTsumw2 = w2[1:-1].sum()
                    h._fTsumw = w[1:-1].sum()
                    n = Shape.template.format(
                        category=self.bin,
                        process=p,
                        variable=self.variable,
                        systematic=syst,
                    )
                    self._newhists[n] = h
                    pbar.update(1)

    def build(self):
        # check that processes are available in histograms
        assert set(self.processes) <= set(self.available_processes)

        # now check if there is any process with no contribution
        self.check_processes()

        if self.data_name != "data_obs":
            self.console.log(
                f"For combine you probably want to rename the data process in the future, I will do it for you now: [bold]{self.data_name} :right_arrow: 'data_obs'"
            )
        # just add data to output histograms
        self._newhists[
            Shape.template.format(
                category=self.category,
                process="data_obs",
                variable=self.variable,
                systematic="nominal",
            )
        ] = export1d_hist(self.hists[self.data_name, "nominal", :].project(self.variable))

        # now check if there is no data
        if (r := self.rate(self.data_name)) <= 0.0:
            self.console.log(f"No data, no fun! (rate: {r})")

        # axis = 0: nuisance name
        # axis = 1: valid (!) processess
        # axis = 2: up/down strength
        self._nuisances_arr = np.zeros((0, len(self.vps), 2))

        self.build_systematics()

        # in the end fix histograms, i.e. negative bins
        self.fix_histograms()

    def dump(self, directory, txt="datacard.txt", shapes="shapes.root"):
        txt_path = os.path.join(directory, txt)
        shapes_path = os.path.join(directory, shapes)
        card = open(txt_path, "w")
        # first write header
        card.write(f"# {repr(self)}{self.sep}")
        # build:
        #  - check histograms
        #  - fix histograms
        #  - build systematics
        self.build()
        card.write(f"imax 1 number of bins{self.sep}")
        card.write(f"jmax {len(self.vps)-1} number of processes minus 1{self.sep}")
        card.write(f"kmax * number of nuisance parameters{self.sep}")
        card.write(f"{self.dashes}{self.sep}")
        card.write(
            f"shapes * {self.bin} {shapes} {self.nominal_pattern} {self.sys_pattern}{self.sep}"
        )
        card.write(f"{self.dashes}{self.sep}")
        card.write(f"bin {self.bin}{self.sep}")
        card.write(f"observation {self.observation}{self.sep}")
        card.write(f"{self.dashes}{self.sep}")
        card.write("bin".ljust(self.spaces))  # category
        card.write("".join([f"{self.bin.ljust(self.spaces)}"] * len(self.vps)) + f"{self.sep}")
        card.write("process".ljust(self.spaces))  # process names
        card.write("".join([f"{p}".ljust(self.spaces) for p in self.vps]) + f"{self.sep}")
        card.write("process".ljust(self.spaces))  # process ids
        card.write("".join([f"{i}".ljust(self.spaces) for i in self.vids]) + f"{self.sep}")
        card.write("rate".ljust(self.spaces))
        card.write(
            "".join([f"{self.rate(p)}".ljust(self.spaces) for p in self.vps]) + f"{self.sep}"
        )
        card.write(f"{self.dashes}{self.sep}")
        if self._nuisances_names:
            for name, type in self._nuisances_names:
                nidx = self.nuisance_idx(name)
                card.write(f"{name} {type}".ljust(self.spaces))
                for pidx in range(len(self.vps)):
                    s = self._nuisances_arr[nidx, pidx, ...]
                    assert s.size == 2
                    if not np.any(s):
                        card.write("-".ljust(self.spaces))
                        continue
                    if np.any(s) and (s[DIR.DOWN] == s[DIR.UP]):
                        card.write(f"{round(s[DIR.DOWN], 3)}".ljust(self.spaces))
                        continue
                    else:
                        card.write(
                            f"{round(s[DIR.DOWN], 3)}/{round(s[DIR.UP], 3)}".ljust(self.spaces)
                        )
                card.write(f"{self.sep}")

        # at the end: add custom lines
        if self.custom_lines:
            card.write(self.sep.join(self.custom_lines))
        card.close()

        # finally dump new/fixed histograms
        with uproot.recreate(shapes_path) as root_file:
            bulk_write(
                root_file,
                tqdm(self._newhists.items(), total=len(self._newhists), desc="write to file"),
            )
