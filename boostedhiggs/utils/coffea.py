import numpy as np
from coffea.processor.accumulator import (
    dict_accumulator,
    defaultdict_accumulator,
    column_accumulator,
    set_accumulator,
    AccumulatorABC,
)
import hist
from hist import Hist
from coffea.processor import PackedSelection as _PS, ProcessorABC as _pABC
from coffea.lookup_tools import evaluator
from abc import abstractmethod
from operator import itemgetter, attrgetter
from time import time
from contextlib import contextmanager
from collections import ChainMap
import numbers


class PackedSelection(_PS):
    def add(self, name, selection):
        if name.startswith("!"):
            raise RuntimeError("selection name must not start with '!'")
        return super().add(name, selection)

    def all(self, *names):
        return self.require(
            **{
                (name if want else name[1:]): want
                for name, want in ((n, n.startswith("!")) for n in names)
            }
        )


class ProcessorABC(_pABC):
    output = "data.coffea"
    debug_dataset = "st"
    debug_uuids = None
    debug_paths = None

    @classmethod
    def requires(cls, task):
        return task.base_requires()

    @classmethod
    def live_callback(cls, accumulator):
        return {}

    @classmethod
    def save(cls, task, output, **kwargs):
        assert not callable(cls.output)
        target = task.output()
        target.parent.touch()
        if cls.output == "*.npy":
            output = {
                key: value.value if isinstance(value, column_accumulator) else value
                for key, value in output.items()
            }
            for name, value in output.items():
                if not isinstance(value, np.ndarray):
                    print(
                        "output[%r]=%r expected to be %r but is %r"
                        % (name, value, np.array, type(value))
                    )
                    del output[name]
            target.dump(output, **kwargs)
            target.touch()
        else:
            target.dump(output, **kwargs)

    def postprocess(self, accumulator):
        return accumulator


class BaseProcessor(ProcessorABC):
    individal_weights = False
    jes_shifts = False
    jec_objects = dict(Jet="AK4PFchs", FatJet="AK8PFPuppi", CorrT1METJet="AK4PFchs")
    dataset_shifts = False

    def __init__(self, task):
        self.publish_message = task.publish_message if task.debug else None
        self.config = task.config_inst
        self.year = task.year
        self.corrections = task.load_corrections()

        self.dataset_axis = hist.axis.StrCategory(
            [], name="dataset", label="Primary dataset", growth=True
        )
        self.dataset_shift_axis = hist.axis.StrCategory(
            [], name="dataset_shift", label="Dataset shift", growth=True
        )
        self.category_axis = hist.axis.StrCategory(
            [], name="category", label="Category selection", growth=True
        )
        self.syst_axis = hist.axis.StrCategory(
            [], name="systematic", label="Shift of systematic uncertainty", growth=True
        )

        self._accumulator = dict_accumulator(
            n_events=defaultdict_accumulator(int),
            sum_gen_weights=defaultdict_accumulator(float),
            object_cutflow=defaultdict_accumulator(int),
            cutflow=hist_accumulator(
                Hist(
                    self.dataset_axis,
                    self.category_axis,
                    hist.axis.Regular(10, 0, 10, name="cutflow", label="Cut index"),
                    storage=hist.storage.Weight(),
                )
            ),
        )

    @contextmanager
    def timeit(self, msg, publish=False):
        if "%" not in msg:
            msg += ": %.2fs"
        t = time()
        try:
            yield
        finally:
            msg %= time() - t
            if self.publish_message:
                if publish:
                    self.publish_message(msg)
                print(repr(self), msg)

    @property
    def accumulator(self):
        return self._accumulator

    def select(self, events, unc="nominal", shift=None):
        raise NotImplementedError("%r needs to be overridden!" % self.select)
        return locals()  # pass needed select_output on as mapping

    def get_dataset(self, events):
        return self.config.get_dataset(events.metadata["dataset"][0])

    def get_dataset_shift(self, events):
        return events.metadata["dataset"][1]

    def get_lfn(self, events, default=None):
        ds = self.get_dataset(events)
        fn = events.metadata["filename"].rsplit("/", 1)[-1]

        for lfn in ds.info[self.get_dataset_shift(events)].aux["lfns"]:
            if lfn.endswith(fn):
                return lfn
        else:
            if default is not None:
                return default
            else:
                raise RuntimeError(
                    "could not find original LFN for: %s" % events.metadata["filename"]
                )

    def get_pu_key(self, events):
        ds = self.get_dataset(events)
        if ds.is_data:
            return "data"
        else:
            lfn = self.get_lfn(events, default="")
            for name, hint in ds.campaign.aux.get("pileup_lfn_scenario_hint", {}).items():
                if hint in lfn:
                    return name
            else:
                return "MC"

    def select_with_jec(self, events):
        sources = ["nominal"]
        ds = self.get_dataset(events)
        do_shifts = ds.is_mc and self.jes_shifts and self.get_dataset_shift(events) == "nominal"
        if do_shifts:
            sources += self.config.aux["jes_sources"]
        if ds.is_mc:
            self.corrections["jet"].inject_nPtGen(events.Jet, events.GenJet, dr=0.4)
            self.corrections["jet"].inject_nPtGen(events.FatJet, events.GenJetAK8, dr=0.8)

        for ev, src, dir in self.corrections["jet"].generate(
            events=events,
            objs=self.jec_objects,
            sources=sources,
            run=ds.aux["run"] if ds.is_data else "mc",
            met="METFixEE2017" if self.year == "2017" else "MET",
        ):
            yield self.select(ev, src, dir)


class Histogramer(BaseProcessor):
    jes_shifts = True
    dataset_shifts = True

    @property
    def variables(self):
        return self.config.variables

    def __init__(self, task):
        super().__init__(task)

        self._accumulator["histograms"] = dict_accumulator(
            {
                variable.name: hist_accumulator(
                    Hist(
                        self.dataset_axis,
                        self.category_axis,
                        self.syst_axis,
                        hist.axis.Regular(
                            variable.binning[0],
                            variable.binning[1],
                            variable.binning[2],
                            name=variable.name,
                            label=variable.x_title,
                        ),
                        metadata={
                            "name": variable.name,
                            "x_title": variable.x_title,
                            "y_title": variable.y_title,
                        },
                        storage=hist.storage.Weight(),
                    )
                )
                for variable in self.variables
            }
        )

    def category_variables(self, category):
        return self.variables

    def variable_full_shifts(self, variable):
        return True

    def process(self, events):
        output = self.accumulator.identity()
        dsname, dsshift = events.metadata["dataset"]
        for select_output in self.select_with_jec(events):
            categories = select_output["categories"]
            selection = select_output["selection"]
            weights = select_output["weights"]

            # determine jes_shift
            jes_shift = str(select_output["unc"])
            if select_output["shift"] is not None:
                jes_shift += "_{}".format(select_output["shift"])

            # collect cutflow output
            if jes_shift == "nominal":
                output += select_output["output"]

            datasets = select_output.get("datasets", {})
            datasets.setdefault(None, ())
            for dsn, dscuts in datasets.items():
                if dsn is None:
                    dsn = dsname
                for category, cuts in categories.items():
                    cut = selection.all(*dscuts, *cuts)
                    if jes_shift == "nominal" == dsshift and not np.any(cut):
                        continue
                    weight = weights.weight()[cut]
                    variations = {jes_shift: None}
                    if jes_shift == "nominal":
                        if dsshift == "nominal":
                            variations.update({var: var for var in weights.variations})
                        else:
                            variations = {dsshift: None}
                    else:
                        assert dsshift == "nominal"  # select_with_jec assures this

                    for variable in self.category_variables(category):
                        values = {}
                        values["dataset"] = dsn
                        values["category"] = category
                        try:
                            values[variable.name] = eval(
                                variable.expression, globals(), ChainMap({}, select_output)
                            )[cut]
                        except Exception as e:
                            raise RuntimeError(
                                f"variable {variable.name} (expr: {variable.expression}) erred in {self.get_lfn(events)}"
                            ) from e

                        full_shifts = self.variable_full_shifts(variable)
                        for variation, modifier in variations.items():
                            if full_shifts or variation == "nominal":
                                weight = weights.weight(modifier=modifier)[cut]
                                values["weight"] = weight
                                values["systematic"] = variation
                                output["histograms"][variable.name].hist.fill(**values)

                    if jes_shift == "nominal" == dsshift:
                        # cutflow
                        cutflow_cuts = set()
                        output["cutflow"].hist.fill(
                            dataset=dsn,
                            category=category,
                            cutflow=np.array([0]),
                            weight=np.array([weights.weight().sum()]),
                        )
                        for i, cutflow_cut in enumerate(cuts):
                            cutflow_cuts.add(cutflow_cut)
                            cutflow_cut = selection.all(*cutflow_cuts)
                            output["cutflow"].hist.fill(
                                dataset=dsn,
                                category=category,
                                cutflow=np.array([i + 1]),
                                weight=np.array([weights.weight()[cutflow_cut].sum()]),
                            )

        return output


class ArrayExporter(BaseProcessor):
    output = "*.npy"
    dtype = None
    sep = "_"

    def __init__(self, task):
        super().__init__(task)

        self._accumulator["arrays"] = dict_accumulator()

    @abstractmethod
    def arrays(self, select_output):
        """
        select_output is the output of self.select
        this function should return an dict of numpy arrays, the "weight" key is reserved
        """
        pass

    def categories(self, select_output):
        selection = select_output.get("selection")
        categories = select_output.get("categories")
        return (
            {cat: selection.all(*cuts) for cat, cuts in categories.items()}
            if selection and categories
            else {"all": slice(None)}
        )

    def process(self, events):
        select_output = self.select(events, unc="nominal", shift=None)
        categories = self.categories(select_output)
        output = select_output["output"]

        if categories:
            arrays = self.arrays(ChainMap({}, select_output))
            dataset = self.get_dataset(events)

            weight = select_output["weights"].weight()
            arrays.setdefault(
                "weight", np.stack([np.full_like(weight, dataset.id), weight], axis=-1)
            )

            assert all(not a.dtype.hasobject for a in arrays.values())

            if self.dtype:
                arrays = {key: array.astype(self.dtype) for key, array in arrays.items()}
        else:
            arrays = {}

        output["arrays"] = dict_accumulator(
            {
                category: dict_accumulator(
                    {key: array_accumulator(array[cut, ...]) for key, array in arrays.items()}
                )
                for category, cut in categories.items()
            }
        )

        return output

    def postprocess(self, output):
        dsids, weights = np.array(
            sorted(
                (
                    (self.config.get_dataset(dsname).id, 1.0 / sum_gen_weight)
                    for (dsname, dsshift), sum_gen_weight in output["sum_gen_weights"].items()
                    if dsshift == "nominal"
                ),
                key=itemgetter(0),
            )
        ).T

        for arrays in output["arrays"].values():
            dsid, weight = arrays["weight"].value.T
            arrays["weight"] = array_accumulator(weight * weights[np.searchsorted(dsids, dsid)])

        out = {
            self.sep.join((category,) + (key if isinstance(key, tuple) else (key,))): aa
            for category, arrays in output["arrays"].items()
            for key, aa in arrays.items()
        }

        output.clear()
        output.update(out)

    @classmethod
    def live_callback(cls, accumulator):
        return {cat: len(arrays["weight"]) for cat, arrays in accumulator["out"]["arrays"].items()}


class _Preheater(_pABC):
    prefix = "preheat_"

    def __init__(self, proc):
        assert isinstance(proc, BaseProcessor)
        self.proc = proc

    @property
    def accumulator(self):
        return dict_accumulator({})

    def process(self, events):
        from distributed import worker_client, Variable, Lock

        assert isinstance(self.proc, BaseProcessor)
        assert not isinstance(self.proc, _Preheater)

        s = self.proc.get_dataset(events).data_source
        d = self.prefix + s

        with worker_client(separate_thread=False) as c:
            v = Variable(d, c)
            l = Lock(d, c)

            if l.acquire(blocking=False):
                self.proc.process(events)

                cols = set()
                for col in events.materialized:
                    col = col.replace("_", ".", 1)
                    try:
                        attrgetter(col)(events)
                    except AttributeError:
                        pass
                    else:
                        cols.add(col)
                cols = sorted(cols)
                v.set(cols)
                return dict_accumulator({s: set_accumulator(cols)})
            else:
                cols = v.get()

        for ag in map(attrgetter, cols):
            data = ag(events)
            data = getattr(data, "content", data)
            if callable(getattr(data, "materialize")):
                data.materialize()
        return dict_accumulator({})

    def postprocess(self, accumulator):
        return accumulator


class array_accumulator(column_accumulator):
    """ column_accumulator with delayed concatenate """

    def __init__(self, value):
        self._empty = value[:0]
        self._value = [value]

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.value)

    def identity(self):
        return self.__class__(self._empty)

    def add(self, other):
        assert self._empty.shape == other._empty.shape
        assert self._empty.dtype == other._empty.dtype
        self._value.extend(v for v in other._value if len(v))

    @property
    def value(self):
        if len(self._value) > 1:
            self._value = [np.concatenate(self._value)]
        return self._value[0]

    def __len__(self):
        return sum(map(len, self._value))


class hist_accumulator(AccumulatorABC):
    """ accumulator for hist """

    def __init__(self, hist):
        assert isinstance(hist, Hist)
        self.hist = hist

    def identity(self):
        return self.__class__(self.hist.copy().reset())

    def add(self, other):
        if isinstance(other, self.__class__):
            other = other.hist
        self.hist += other


# https://github.com/CoffeaTeam/coffea/blob/c940f4118bb0fbbd7e2ae4fe27a7ddcd80634c2b/coffea/hist/hist_tools.py#L1244-L1280
def scale(self, factor, axis=None):
    assert isinstance(self, hist.Hist)
    """Scale histogram in-place by factor

    Parameters
    ----------
        factor : float or dict
            A number or mapping of identifier to number
        axis : optional
            Which (sparse) axis the dict applies to

    Examples
    --------
    This function is useful to quickly reweight according to some
    weight mapping along a sparse axis, such as the ``species`` axis
    in the `Hist` example:

    >>> h.scale({'ducks': 0.3, 'geese': 1.2}, axis='species')
    """
    if self._sumw2 is None:
        self._init_sumw2()
    if isinstance(factor, numbers.Number) and axis is None:
        for key in self._sumw.keys():
            self._sumw[key] *= factor
            self._sumw2[key] *= factor ** 2
    elif isinstance(factor, dict):
        if not isinstance(axis, tuple):
            axis = (axis,)
        axis = tuple(map(self.axis, axis))
        kget = itemgetter(*map(self._isparse, axis))
        factor = {tuple(a.index(e) for a, e in zip(axis, k)): v for k, v in factor.items()}
        for key in self._sumw.keys():
            fkey = kget(key)
            if fkey in factor:
                self._sumw[key] *= factor[fkey]
                self._sumw2[key] *= factor[fkey] ** 2
    elif isinstance(factor, np.ndarray):
        axis = self.axis(axis)
        raise NotImplementedError("Scale dense dimension by a factor")
    else:
        raise TypeError("Could not interpret scale factor")


def sf(correction, particle):
    return correction(particle.eta, particle.pt).prod()


def mk_dense_evaluator(values, edges=None):
    if edges is not None:
        if not isinstance(edges, dict):
            edges = {k: edges for k in values.keys()}
        values = {k: (v, edges[k]) for k, v in values.items()}

    return evaluator(
        names={k: k for k in values.keys()},
        types={k: "dense_lookup" for k in values.keys()},
        primitives=values,
    )
