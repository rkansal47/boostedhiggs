import numpy as np
from scipy.stats import chi2
from scipy.special import erf
from scipy.interpolate import RectBivariateSpline
from tools.evil import pin
from collections import ChainMap
from collections.abc import Mapping, Sequence
from bisect import bisect_right

try:
    from functools import cached_property
except ImportError:
    from backports.cached_property import cached_property
from operator import methodcaller
from pathlib import Path
from warnings import warn
import yaml

CL1 = erf(1 / np.sqrt(2.0))


def meanstd(values, ddof=1):
    s0, s1, s2 = 0, 0, 0
    for v in values:
        s0 += 1
        s1 += v
        s2 += v * v
    s1 /= s0
    s2 /= s0
    sv = s2 - s1 * s1
    sv *= s0 / (s0 - ddof) if s0 > ddof else 0.0
    return s1, np.sqrt(np.maximum(sv, 0))


def parse_config(path, config=None, first=False):
    with open(path, "r") as f:
        c = yaml.safe_load(f)
    if isinstance(config, ChainMap):
        return config.new_child(c)
    if isinstance(config, Mapping):
        return ChainMap(c, config)
    else:
        return c


class coprop(cached_property):
    def __set__(self, obj, value):
        obj.__dict__[self.attrname] = value


class _indirect_tuple(tuple):
    def __new__(self, base, ids):
        ret = tuple.__new__(_indirect_tuple, ids)
        ret.base = base
        return ret

    def __getitem__(self, idx):
        id = super().__getitem__(idx)
        if isinstance(idx, slice):
            return type(self)(self.base, id)
        else:
            return self.base[id]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class AutoDict(dict):
    def __init__(self, func):
        assert callable(func)
        pin(locals())

    def __missing__(self, key):
        return self.func(key)


class CacheDict(AutoDict):
    def __missing__(self, key):
        return self.setdefault(key, super().__missing__(key))


class StateSansCache:
    def _clean_cache(self, __dict__=None):
        cls = self.__class__
        if __dict__ is None:
            __dict__ = self.__dict__
        for key in list(__dict__.keys()):
            val = getattr(cls, key, None)
            if isinstance(val, cached_property) and not isinstance(val, coprop):
                del __dict__[key]
        return __dict__

    def __getstate__(self):
        return self._clean_cache(self.__dict__.copy())


flavors = dict(
    enumerate(
        [
            "down",
            "up",
            "strange",
            "charm",
            "bottom",
            "top",
        ],
        start=1,
    )
)


# actual classes
class LHAPDF(StateSansCache):
    def __init__(self, path, config=None):
        path = Path(path)
        assert path.is_dir()
        config = parse_config(path.joinpath("lhapdf.conf"), config)
        pin(locals())

    @cached_property
    def index(self):
        with self.path.joinpath("pdfsets.index").open("r") as file:
            return {int(id): name for id, name, _ in map(methodcaller("split"), file)}

    @cached_property
    def _index(self):
        return tuple(sorted(self.index.keys()))

    def id2name(self, id):
        assert isinstance(id, int)
        if id in self.index:
            return self.index[id]
        idx = bisect_right(self._index, id)
        if idx < 1:
            raise ValueError("id %r not in index" % id)
        return self.index[self._index[idx - 1]]

    def getSet(self, set):
        if isinstance(set, int):
            set = self.id2name(set)
        return PDFSet(self.path.joinpath(set), self.config)


class PDFSet(StateSansCache, Sequence):
    def __init__(self, path, config=None):
        path = Path(path)
        if path.suffix == ".info" and path.is_file():
            config = parse_config(path, config)
            path = path.parent
        elif path.is_dir():
            config = parse_config(path.joinpath(path.name).with_suffix(".info"), config)
        else:
            raise ValueError("path not understood: %r" % path)
        _pdfs = CacheDict(self._make_pdf)
        pin(locals())

    def _make_pdf(self, idx):
        if not 0 <= idx < len(self):
            raise IndexError("index %r out of range [0,%d)" % (idx, len(self)))
        ret = PDF(
            self.path.joinpath("%s_%04d.dat" % (self.path.name, idx)),
            self.config,
        )
        ret.idx = idx
        return ret

    def __len__(self):
        return self.config["NumMembers"]

    def __getitem__(self, num):
        return self._pdfs[num]

    @property
    def name(self):
        return self.path.name

    @property
    def description(self):
        return self.config["SetDesc"]

    @property
    def lhapdfID(self):
        return int(self.config.get("SetIndex", -1))

    @property
    def dataversion(self):
        return int(self.config.get("DataVersion", -1))

    @property
    def errorType(self):
        return self.config.get("ErrorType", "UNKNOWN")

    @property
    def hasReplicas(self):
        return self.errorType.startswith("replicas") if self.errorType else None

    @property
    def errorConfLevel(self):
        return self.config.get(
            "ErrorConfLevel",
            -1 if self.hasReplicas else 100 * CL1,
        )

    def central(self, inp, alternative=False):
        if self.hasReplicas:
            res = self.uncertainty(inp, cl=-1, alternative=alternative)
        else:
            res = PDFResult(self, inp)
        return res.central

    def uncertainty(self, inp, cl=100 * CL1, alternative=False):
        setCL = CL1 if self.hasReplicas else self.errorConfLevel / 100
        assert 0 < setCL < 1
        reqCL = cl / 100 if cl >= 0 else setCL
        assert 0 < reqCL < 1

        ret = PDFResult(self, inp)
        if not ret.members:
            raise RuntimeError("PDFSet must contain members (besides variations)")

        ret.cl = cl
        if setCL != reqCL:
            ret.scale = np.sqrt(chi2.ppf(reqCL, 1) / chi2.ppf(setCL, 1))

        if alternative:
            assert ret.replicas
            values = np.stack(list(ret.replicas), axis=0)
            ret.central = np.median(values, axis=0, overwrite_input=True)
            minus, plus = np.percentile(
                values,
                50 * (1 + np.array([-reqCL, reqCL])),
                axis=0,
                overwrite_input=True,
                interpolation="nearest",
            )
            ret.errminus_pdf = ret.central - minus
            ret.errplus_pdf = plus - ret.central

        elif ret.replicas:
            ret.central, ret.errsymm_pdf = meanstd(ret.replicas)
        elif ret.symmhessian:
            ret.errsymm_pdf = np.sqrt(sum(np.square(v - ret.central) for v in ret.symmhessian))
        elif ret.hessian:
            errplus = errminus = errsymm = 0
            for down, up in zip(*ret.hessian):
                a, b = down - self.central, up - self.central
                errplus += np.square(np.maximum(np.maximum(a, b), 0.0))
                errminus += np.square(np.minimum(np.minimum(a, b), 0.0))
                errsymm += np.square(down - up)
            ret.errsymm_pdf = 0.5 * np.sqrt(errsymm)
            ret.errplus_pdf = np.sqrt(errplus) * ret.scale
            ret.errminus_pdf = np.sqrt(errminus) * ret.scale
        else:
            raise RuntimeError("ErrorType: %r not supported by PDFSet.uncertainty" % self.errorType)

        if not alternative:
            ret.errsymm_pdf *= ret.scale

        if ret.variations:
            err_par = sum(np.square(a - b) for a, b in ret.variations.values())
            ret.err_par = ret.scale * 0.5 * np.sqrt(err_par)

        return ret

    def correlation(self, a, b):
        assert isinstance(a, PDFResult)
        assert self is a.pdfset
        assert isinstance(b, PDFResult)
        assert self is b.pdfset
        raise NotImplementedError()

    def randomValueFromHessian(self, res, random, symmetrise):
        assert isinstance(res, PDFResult)
        assert self is res.pdfset
        random = np.asarray(random)
        assert random.shape == (len(self),) + res.central.shape
        raise NotImplementedError()

    def __repr__(self):
        return "<%s(name=%r, lhapdfID=%d)>" % (type(self).__name__, self.name, self.lhapdfID)

    def __eq__(self, other):
        return isinstance(other, PDFSet) and self.lhapdfID == other.lhapdfID != -1

    def __ne__(self, other):
        return not self == other


class PDFResult:
    def __init__(self, pdfset, values):
        assert isinstance(pdfset, PDFSet)
        if isinstance(values, PDFInput):
            inp = values
            values = AutoDict(self._from_input)
        else:
            assert isinstance(values, (Sequence, np.ndarray))
            assert len(pdfset) == len(values)
        pin(locals())

    def _from_input(self, id):
        return self.pdfset[id].xfx(self.inp)

    @cached_property
    def variations(self):
        return {
            k: _indirect_tuple(self.values, range(len(self))[-2 * i :][:2])
            for i, k in enumerate(self.errorType.split("+")[1:][::-1], start=1)
        }

    @cached_property
    def all_values(self):
        return _indirect_tuple(self.values, range(len(self)))

    @cached_property
    def members(self):
        return _indirect_tuple(self.values, range(1, len(self) - len(self.variations) * 2))

    @cached_property
    def replicas(self):
        return self.members if self.errorType.startswith("replicas") else ()

    @cached_property
    def hessian(self):
        if not self.errorType.startswith("hessian"):
            return ()
        assert len(self.members) & 1 == 0
        return self.members[0::2], self.members[1::2]

    @cached_property
    def symmhessian(self):
        return self.members if self.errorType.startswith("symmhessian") else ()

    @coprop
    def central(self):
        return self.values[0]

    @coprop
    def scale(self):
        return 1.0

    @coprop
    def errsymm_pdf(self):
        return 0.5 * (self.errplus_pdf + self.errminus_pdf)

    @coprop
    def errplus_pdf(self):
        return self.errsymm_pdf

    @coprop
    def errminus_pdf(self):
        return self.errsymm_pdf

    @coprop
    def err_par(self):
        return 0.0

    @coprop
    def errsymm(self):
        return np.sqrt(np.square(self.errsymm_pdf) + np.square(self.err_par))

    @coprop
    def errplus(self):
        return np.sqrt(np.square(self.errplus_pdf) + np.square(self.err_par))

    @coprop
    def errminus(self):
        return np.sqrt(np.square(self.errminus_pdf) + np.square(self.err_par))

    @property
    def errs(self):
        return [self.errminus, self.errplus]

    @property
    def errs_pdf(self):
        return [self.errminus_pdf, self.errplus_pdf]

    @property
    def errmin(self):
        return self.central - self.errminus

    @property
    def errmax(self):
        return self.central + self.errplus

    @property
    def errrange(self):
        return self.errmin, self.errmax

    @property
    def errmin_pdf(self):
        return self.central - self.errminus_pdf

    @property
    def errmax_pdf(self):
        return self.central + self.errplus_pdf

    @property
    def errrange_pdf(self):
        return self.errmin_pdf, self.errmax_pdf

    def __len__(self):
        return len(self.pdfset)

    def __getattr__(self, name):
        return getattr(self.pdfset, name)


class PDFInput(StateSansCache):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __getitem__(self, mask):
        return PDFInputMasked(_base=self, _mask=mask)

    def __call__(self, **kwargs):
        kwargs.setdefault("_base", self)
        return PDFInputBased(**kwargs)

    def _reduce(self, array):
        axis = getattr(self, "_axis", ())
        return array if axis == () else array.prod(axis=axis)

    def _own_attrib(self, *names):
        pass

    @coprop
    def q2(self):
        self._own_attrib("q")
        return self.q * self.q

    @coprop
    def logq2(self):
        self._own_attrib("q", "q2")
        return np.log(self.q2)

    @coprop
    def ix(self):
        self._own_attrib("x")
        return 1 / self.x

    @coprop
    def logx(self):
        self._own_attrib("x")
        return np.log(self.x)

    @property
    def xq2(self):
        return self.x, self.q2

    @property
    def logxq2(self):
        return self.logx, self.logq2

    def broadcast_call(self, keys, func, *args, **kwargs):
        keys = tuple(keys)
        assert keys
        vals = np.broadcast_arrays(*(getattr(self, k) for k in keys))
        inp = self(**{k: v.reshape(-1) for k, v in zip(keys, vals)})
        return func(inp, *args, **kwargs).reshape(vals[0].shape)


class PDFInputBased(PDFInput):
    def __init__(self, **kwargs):
        assert isinstance(kwargs.get("_base", None), PDFInput)
        super().__init__(**kwargs)

    def _own_attrib(self, *names):
        for name in names:
            if name in self.__dict__:
                return
        raise AttributeError("%r has no attribute %r" % (self, name))

    def __getattr__(self, name):
        return getattr(self._base, name)


class PDFInputMasked(PDFInputBased):
    def __init__(self, **kwargs):
        assert "_mask" in kwargs
        super().__init__(**kwargs)

    def __getattr__(self, name):
        val = getattr(self._base, name)
        if not np.asarray(val).ndim:
            return val
        if isinstance(self._mask, np.ndarray) and self._mask.dtype == bool:
            val = np.broadcast_to(val, self._mask.shape)
        return val[self._mask]


class PDF(StateSansCache):
    dtype = float

    def _parse_config(self, path, config=None):
        with open(path, "r") as f:
            l = yaml.SafeLoader(f)
            c = l.get_data()
            h = l.line + 1
            l.dispose()
        if isinstance(config, ChainMap):
            c = config.new_child(c)
        elif isinstance(config, Mapping):
            c = ChainMap(c, config)
        return c, h

    @staticmethod
    def _parse_line(line):
        return tuple(map(float, line.strip().split()))

    def _parse_blocks(self, path, skip):
        with open(path, "rt") as file:
            lines = iter(line for line in file if not line.startswith("#"))
            # read yaml header
            for i, line in zip(range(skip), lines):
                pass
            assert line.strip() == "---"
            # read blocks
            while True:
                axes = tuple(
                    np.asarray(self._parse_line(l), dtype=t)
                    for t, l in zip((self.dtype, self.dtype, int), lines)
                )
                if not axes:
                    break
                x, q, pdgId = axes
                n = len(x) * len(q)
                xfx = np.empty((n, len(pdgId)), dtype=self.dtype)
                for i, l in zip(range(n), lines):
                    xfx[i] = self._parse_line(l)
                line = next(lines)
                assert line.strip() == "---"
                xfx = xfx.reshape(x.shape + q.shape + pdgId.shape)
                xfx = np.moveaxis(xfx, -1, 0)
                eq = np.zeros_like(q, dtype=bool)
                eq[0] = eq[-1] = True
                yield PDFBlock(pdf=self, x=x, q=q, pdgId=pdgId.tolist(), xfx=xfx, eq=eq)

    @staticmethod
    def _merge_blocks(blocks):
        A = None
        for V in blocks:
            if A is None:
                A = V
            else:
                try:
                    A = A.merge(A, V)
                except AssertionError:
                    yield A
                    A = V
        if A is not None:
            yield A

    def __init__(self, path, config=None):
        self.name = path.name
        config, skip = self._parse_config(path, config)
        assert config["Format"] == "lhagrid1"
        blocks = list(self._merge_blocks(self._parse_blocks(path, skip)))
        pin(locals())

    def fx(self, inp, check=True):
        return inp._reduce(self._xfx(inp, check=check) * inp.ix)

    def xfx(self, inp, check=True):
        return inp._reduce(self._xfx(inp, check=check))

    def _xfx(self, inp, check=True):
        x, q2 = inp.xq2
        if check:
            assert np.all(0 <= x)
            assert np.all(x <= 1)
            assert np.all(0 <= q2)
        if len(self.blocks) == 1:
            return np.maximum(self.blocks[0](inp), self.xfxMin)
        # multiple blocks
        x, q2, pdgId = np.broadcast_arrays(x, q2, inp.pdgId)
        out = np.empty_like(x)
        idx = np.searchsorted(self._block_q2edges, q2, side="right")
        for i, block in enumerate(self.blocks):
            mask = idx == i
            out[mask] = block(inp[mask])
        return np.maximum(out, self.xfxMin)

    def __repr__(self):
        return "<%s(name=%r)>" % (type(self).__name__, self.name)

    @property
    def xfxMin(self):
        return {0: -np.inf, 1: 0, 2: 1e-10}[self.config["ForcePositive"]]

    @cached_property
    def _block_q2edges(self):
        assert all(l.q2[-1] == h.q2[0] for l, h in zip(self.blocks, self.blocks[1:]))
        return np.array([b.q2[-1] for b in self.blocks[:-1]])

    @cached_property
    def pdgId(self):
        return np.sort(np.concatenate([block.pdgId for block in self.blocks], axis=0))

    @cached_property
    def x(self):
        return np.sort(np.concatenate([block.x for block in self.blocks], axis=0))

    @cached_property
    def q(self):
        return np.sort(np.concatenate([block.q for block in self.blocks], axis=0))

    @cached_property
    def q2(self):
        return np.square(self.q)

    @cached_property
    def _alphaS(self):
        return AlphaS.make(self.config)

    def alphaS(self, inp):
        return self._alphaS(inp)

    def _alphaScheck(self, other, inp, tolerance):
        if np.all(tolerance < 0):
            return True
        a = self.alphaS(inp)
        b = other.alphaS(inp)
        bad = (2 * np.abs(a - b) / (np.abs(a) + np.abs(b)) > tolerance) & (tolerance >= 0)
        if np.any(bad):
            warn(
                "WARNING: %d (%.1f%%) alphaS mismatch(es) at q2=%s: %s=%s vs %s=%s"
                % (
                    bad.sum(),
                    100 * bad.mean(),
                    np.broadcast_to(inp.q2, bad.shape)[bad],
                    self,
                    a[bad],
                    other,
                    b[bad],
                )
            )
            return False
        return True

    def reweight(self, target, inp, aschk=5e-2):
        self._alphaScheck(target, inp, tolerance=aschk)
        return target.xfx(inp) / self.xfx(inp)


class PDFBlock:
    @classmethod
    def merge(cls, A, V):
        if A.q[-1] == V.q[0]:
            assert np.all(A.x == V.x)
            assert np.all(A.pdgId == V.pdgId)
            assert np.all(A.xfx[..., -1] == V.xfx[..., 0])
            return cls(
                pdf=A.pdf,
                x=A.x,
                q=np.concatenate((A.q[:-1], V.q), axis=0),
                pdgId=A.pdgId,
                xfx=np.concatenate((A.xfx[..., :-1], V.xfx), axis=-1),
                eq=np.concatenate((A.eq[:-1], V.eq), axis=0),
            )
        assert A.q[0] == V.q[-1]
        return cls.merge(V, A)

    def __init__(self, pdf, x, q, pdgId, xfx, eq):
        assert xfx.ndim == 3
        assert xfx.shape == (len(pdgId),) + x.shape + q.shape
        xfx = np.ascontiguousarray(xfx)
        pin(locals())

    @cached_property
    def q2(self):
        return np.square(self.q)

    @cached_property
    def pdgId2XFX(self):
        return CacheDict(self._pdgId2XFX)

    def _pdgId2XFX(self, pdgId):
        if pdgId in self.pdgId:
            return XFX(self, self.xfx[self.pdgId.index(pdgId)])
        warn("pdgId %d not available in %r" % (pdgId, self.pdf))
        return XFX.unavailable

    def __call__(self, inp):
        x, q2 = inp.xq2
        pdgId = np.asarray(inp.pdgId, dtype=int)
        if 0 in pdgId:
            pdgId = np.where(pdgId == 0, 21, pdgId)
        pdgIds = np.unique(pdgId)
        if len(pdgIds) == 1:
            pdgId = pdgId.flat[0]
        if pdgId.ndim:
            x, q2, pdgId = np.broadcast_arrays(x, q2, pdgId)
            xfx = np.empty_like(x)
            for pId in pdgIds:
                mask = pdgId == pId
                xfx[mask] = self.pdgId2XFX[pId](inp[mask])
        else:
            xfx = inp.broadcast_call(("x", "q", "q2"), self.pdgId2XFX[pdgId])
        return xfx


class AlphaS(StateSansCache):
    @classmethod
    def make(cls, config):
        typ = config["AlphaS_Type"]
        name = "%s_%s" % (cls.__name__, typ)
        for c in cls.__subclasses__():
            if c.__name__ == name:
                return c(config)
        else:
            raise NotImplementedError("unsupported AlphaS type: %r" % typ)

    def threshold(self, key):
        key = flavors.get(key, key)
        key = "Threshold%s" % key.capitalize()
        return self.config.get("AlphaS_%s" % key, self.config.get(key, None))

    def mass(self, key):
        key = flavors.get(key, key)
        key = "M%s" % key.capitalize()
        return self.config.get("AlphaS_%s" % key, self.config.get(key, None))

    def __init__(self, config):
        pin(locals())
        # assert None not in set(map(self.threshold, flavors))
        assert None not in set(map(self.mass, flavors))


class AlphaS_ipol(AlphaS):
    @cached_property
    def blocks(self):
        q = np.asarray(self.config["AlphaS_Qs"])
        v = np.asarray(self.config["AlphaS_Vals"])
        assert q.shape == v.shape
        assert q.ndim == v.ndim == 1

        dx = np.diff(q)
        assert np.all(dx >= 0)

        # TODO: embed edges in x where dy == 0

        disc = np.where(dx == 0)[0] + 1

        return tuple(
            PDFInput(q=q[b:e], v=v[b:e]) for b, e in zip(np.r_[0, disc], np.r_[disc, len(q)])
        )

    @cached_property
    def edges(self):
        return tuple(b.logq2[-1] for b in self.blocks[:-1])

    @cached_property
    def ipols(self):
        return tuple(CubicHermiteSpline(x=b.logq2, y=b.v) for b in self.blocks)

    @cached_property
    def c(self):
        b = self.blocks[0]
        return np.log10(b.v[1] / b.v[0]) / np.log10(b.q2[1] / b.q2[0])

    @cached_property
    def a(self):
        return self.blocks[0].v[0] * self.blocks[0].q2[0] ** -self.c

    def __call__(self, inp):
        x = np.minimum(inp.logq2, self.blocks[-1].logq2[-1])

        if len(self.blocks) == 1:
            ret = self.ipols[0](x)
        else:
            ret = np.empty_like(x)
            # import ipdb; ipdb.set_trace()
            idx = np.searchsorted(self.edges, x, side="right")
            for i, ipol in enumerate(self.ipols):
                mask = i == idx
                ret[mask] = ipol(x[mask])

        ret[np.abs(ret) >= 2] = np.finfo(ret.dtype).max

        l = x < self.blocks[0].logq2[0]
        ret[l] = self.a * np.exp(x[l] * self.c)
        # ret[l] = self.ref.v[0] * (inp.q2 / self.logq2.v[0])**self.c

        return ret


class XFX(StateSansCache):
    @classmethod
    def unavailable(cls, inp):
        return np.zeros_like(inp.x)

    def __init__(self, block, xfx):
        pin(locals())

    @property
    def ipol(self):
        ipol = self.block.pdf.config["Interpolator"]
        assert ipol in ("linear", "cubic", "log", "logcubic")
        return ipol

    @property
    def epol(self):
        epol = self.block.pdf.config["Extrapolator"]
        assert epol in ("nearest", "error", "continuation")
        return epol

    @property
    def ipol_degree(self):
        return 3 if "cubic" in self.ipol else 1

    def ipol_xq2(self, inp):
        if "log" in self.ipol:
            return inp.logxq2
        return inp.xq2

    @property
    def xmin(self):
        return self.block.x[0]

    @property
    def xmax(self):
        return self.block.x[-1]

    @property
    def q2min(self):
        return self.block.q2[0]

    @property
    def q2max(self):
        return self.block.q2[-1]

    @property
    def pxmin(self):
        return self.block.x[:2]

    @property
    def pq2min(self):
        return [self.q2min, 1.01 * self.q2min]

    @property
    def pq2max(self):
        return self.block.q2[-2:]

    @cached_property
    def interpolator(self):
        x, q2 = self.ipol_xq2(PDFInput(x=self.block.x, q2=self.block.q2))
        if "cubic" in self.ipol:
            if len(x) < 4:
                raise RuntimeError("too few knots in x for %r interpolation" % self.ipol)
            if len(q2) < 2:
                raise RuntimeError("too few knots in q2 for %r interpolation" % self.ipol)
            if len(q2) >= 4:
                return BicubicHermiteSpline(x=x, y=q2, z=self.xfx, ey=self.block.eq)
        return RectBivariateSpline(x=x, y=q2, z=self.xfx, kx=1, ky=1)

    def interpolate(self, inp, agrid=True):
        x, q2 = map(np.asarray, self.ipol_xq2(inp))
        if agrid:
            x = x.reshape(x.shape + (1,) * q2.ndim)

        return self.interpolator.ev(x, q2)

    def extrapolate_IO(self, inp):
        v = self.interpolate(inp(q2=self.pq2max))
        return _extrapolateLinear(inp.logq2, *np.log(self.pq2max), *v.T)

    def extrapolate_UO(self, inp):
        v = self.interpolate(PDFInput(x=self.pxmin, q2=self.pq2max))
        w = _extrapolateLinear(
            x=inp.logq2[..., None],
            xl=np.log(self.pq2max[0, None, None]),
            xh=np.log(self.pq2max[1, None, None]),
            yl=v[None, :, 0],
            yh=v[None, :, 1],
        )
        return _extrapolateLinear(inp.logx, *np.log(self.pxmin), *w.T)

    def extrapolate_UI(self, inp):
        v = self.interpolate(inp(x=self.pxmin))
        return _extrapolateLinear(inp.logx, *np.log(self.pxmin), *v)

    def extrapolate_UU(self, inp):
        v = self.interpolate(PDFInput(x=self.pxmin, q2=self.pq2min))
        w = _extrapolateLinear(
            x=inp.logx[..., None],
            xl=np.log(self.pxmin[0, None, None]),
            xh=np.log(self.pxmin[1, None, None]),
            yl=v[None, 0, :],
            yh=v[None, 1, :],
        )
        return self._extrapolate__U(inp.q2, *w.T)

    def extrapolate_IU(self, inp):
        v = self.interpolate(inp(q2=self.pq2min))
        return self._extrapolate__U(inp.q2, *v.T)

    # https://gitlab.com/hepcedar/lhapdf/-/blob/master/src/ContinuationExtrapolator.cc#L95-113
    def _extrapolate__U(self, q2, fq2Min, fq2Min1):
        assert q2.shape == fq2Min.shape == fq2Min1.shape
        anom = np.maximum(
            -2.5,
            (fq2Min1 - fq2Min) / fq2Min / 0.01,
            where=np.abs(fq2Min) >= 1e-5,
            out=np.ones_like(q2),
        )
        q2r = q2 / self.q2min
        return fq2Min * q2r ** (anom * q2r + 1.0 - q2r)

    def extrapolate_fm(self, inp):
        x, q2 = inp.xq2
        xU = x < self.xmin
        qU = q2 < self.q2min
        qO = self.q2max < q2
        yield self.extrapolate_IO, ~xU & qO
        yield self.extrapolate_UO, xU & qO
        yield self.extrapolate_UI, xU & ~(qO | qU)
        yield self.extrapolate_UU, xU & qU
        yield self.extrapolate_IU, ~xU & qU

    def __call__(self, inp):
        x, q2 = inp.xq2
        assert max(x.ndim, q2.ndim) <= 1
        if self.epol == "error":
            assert np.all(self.xmin <= x)
            assert np.all(x <= self.xmax)
            assert np.all(self.q2min <= q2)
            assert np.all(q2 <= self.q2max)
        xfx = self.interpolate(inp, agrid=False)
        if self.epol == "continuation":
            assert np.all(x <= self.xmax)
            for func, mask in self.extrapolate_fm(inp):
                if np.any(mask):
                    xfx[mask] = func(inp[mask])
        return xfx


class CubicHermiteSpline(StateSansCache):
    def __init__(self, x, y, e=None):
        x = Coord(x, e)
        del e
        pin(locals())

    @cached_property
    def vdLH(self):
        vd = np.diff(self.y, axis=0) / self.x.dv
        vdp = np.r_[vd[0], vd, vd[-1]]
        vdc = (vd[:-1] + vd[1:]) * 0.5
        vdc = np.r_[np.nan, vdc, np.nan]
        dv = self.x.dv
        return (
            np.r_[dv, dv[-1]] * np.where(self.x.e, vdp[1:], vdc),
            np.r_[dv[0], dv] * np.where(self.x.e, vdp[:-1], vdc),
        )

    def __call__(self, x):
        ix, i1, dx, tx = self.x.prep(x)

        return _interpolateCubic(
            T=tx,
            VL=self.y[ix],
            VH=self.y[i1],
            VDL=self.vdLH[0][ix],
            VDH=self.vdLH[1][i1],
        )


class BicubicHermiteSpline(StateSansCache):
    dtype = float

    def __init__(self, x, y, z, ex=None, ey=None):
        assert ex is None
        x = Coord(x, ex)
        y = Coord(y, ey)
        z = np.asarray(z, dtype=self.dtype)
        assert x.v.shape + y.v.shape == z.shape
        pin(locals())

    @cached_property
    def zdx(self):
        v = np.diff(self.z, axis=0) * self.x.idv[:, None]
        v0 = np.zeros_like(self.y.v)[None, :]
        vx = np.concatenate((v0, v, v0), axis=0)
        vs = vx[1:] + vx[:-1]
        vs[1:-1] *= 0.5
        return vs

    def ev(self, x, y):
        # TODO: consider linear fall through for y
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        x, y = np.broadcast_arrays(x, y)
        shape = x.shape
        if max(x.ndim, y.ndim) > 1:
            x, y = map(np.ravel, (x, y))

        oy = np.arange(-1, 3, dtype=int)
        self.zdx

        ix, i1, dx, tx = self.x.prep(x)
        iy, i3, dy, ty = self.y.prep(y)

        i2 = iclip(iy[:, None] + oy[None, 0:], self.y, low=True)

        v = _interpolateCubic(
            T=tx[:, None],
            VL=self.z[ix[:, None], i2],
            VH=self.z[i1[:, None], i2],
            VDL=dx[:, None] * self.zdx[ix[:, None], i2],
            VDH=dx[:, None] * self.zdx[i1[:, None], i2],
        )

        idy = self.y.idv
        vd = np.diff(v, axis=-1)
        vd *= idy[iclip(iy[:, None] + oy[None, :-1], idy, low=True)]
        vdc = (vd[:, :-1] + vd[:, 1:]) * 0.5

        ey = self.y.e
        ret = _interpolateCubic(
            T=ty,
            VL=v[:, 1],
            VH=v[:, 2],
            VDL=dy * np.where(ey[iy], vd[:, 1], vdc[:, 0]),
            VDH=dy * np.where(ey[i3], vd[:, 1], vdc[:, 1]),
        )

        return ret.reshape(shape)


class Coord(StateSansCache):
    def __init__(self, v, e=None):
        v = np.asarray(v)
        assert v.ndim == 1
        if e is None:
            e = np.zeros_like(v, dtype=bool)
        else:
            e = np.asarray(e, dtype=bool)
            assert e.shape == v.shape
        e[0] = e[-1] = True
        pin(locals())

    @cached_property
    def dv(self):
        return np.diff(self.v)

    @cached_property
    def idv(self):
        return 1 / self.dv

    def prep(self, v):
        i = iclip(np.searchsorted(self.v, v, side="right") - 1, self.v, low=True)
        j = iclip(i, self.dv)

        dv = self.dv[j]
        tv = (v - self.v[i]) * self.idv[j]

        return i, iclip(i + 1, self.v), dv, tv

    @property
    def shape(self):
        return self.v.shape


def iclip(idx, ref, axis=0, low=False, high=True):
    if high:
        l = ref.shape[axis] - 1
        if low:
            return np.clip(idx, 0, l)
        else:
            return np.minimum(idx, l)
    elif low:
        return np.maximum(idx, 0)
    else:
        return idx


# https://gitlab.com/hepcedar/lhapdf/-/blob/master/src/BicubicInterpolator.cc#L21-36
def _interpolateCubic(T, VL, VDL, VH, VDH):
    # Pre-calculate powers of T
    t2 = T * T
    t3 = t2 * T

    # Calculate left point
    p0 = (2 * t3 - 3 * t2 + 1) * VL
    m0 = (t3 - 2 * t2 + T) * VDL

    # Calculate right point
    p1 = (-2 * t3 + 3 * t2) * VH
    m1 = (t3 - t2) * VDH

    return p0 + m0 + p1 + m1


# https://gitlab.com/hepcedar/lhapdf/-/blob/master/src/ContinuationExtrapolator.cc#L14-24
def _extrapolateLinear(x, xl, xh, yl, yh):
    r = (x - xl) / (xh - xl)
    l = (yl > 1e-3) & (yh > 1e-3)
    a = np.any(l)
    if a:
        yl = np.log(yl, out=yl.copy(), where=l)
        yh = np.log(yh, out=yh.copy(), where=l)
        # yl[l] = np.log(yl[l])
        # yh[l] = np.log(yh[l])
    r = yl + r * (yh - yl)
    if a:
        r = np.exp(r, out=r, where=l)
        # r[l] = np.exp(r[l])
    return r
