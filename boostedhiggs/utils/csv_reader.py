import numpy as np
from collections import defaultdict
from string import ascii_lowercase
import csv
import re
from operator import attrgetter
from enum import Enum

from tqdm import tqdm
import sympy


class IntEnum(int, Enum):
    @classmethod
    def _missing_(cls, value):
        return cls._member_map_.get(value, None)


class WP(IntEnum):
    loose = 0
    medium = 1
    tight = 2
    reshape = 3


class FID(IntEnum):
    b = 0
    c = 1
    l = 2


class BTagSFCSVReader:
    type = None
    cols = [
        "OP|OperatingPoint",
        "measurement(Type|)",
        "sysType",
        "jetFlavou?r",
        "etaMin",
        "etaMax",
        "ptMin",
        "ptMax",
        "discrMin",
        "discrMax",
        "formula",
    ]

    # https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
    def __init__(self, fn):
        raw = defaultdict(list)
        with open(fn, "r") as f:
            for i, row in enumerate(csv.reader(f)):
                if i == 0 and ";" in row[0]:
                    self.type, row[0] = row[0].split(";")
                    row = [r.strip() for r in row]
                    assert len(row) == len(self.cols)
                    for col, val in zip(self.cols, row):
                        assert re.match(f"({col})$", val), (col, val)
                    self.cols = row
                    continue
                row = [r.strip().strip('"') for r in row]
                raw[tuple(row[:4])].append(row[4:])

        self.rules = {}
        for (wp, meth, syst, fid), entries in tqdm(raw.items(), total=len(raw), unit="rule"):
            wp, fid = WP(int(wp)), FID(int(fid))

            self.rules.setdefault((wp, fid, meth), {})[syst] = BTagSFRule(
                entries, isShape=wp == WP.reshape
            )

    def _get(self, wp, fid, meth):
        if meth is None:
            if wp == WP.reshape:
                meth = "iterativefit"
            elif fid == FID.l:
                meth = "incl"
            else:
                avail = self.methods(wp, fid)
                if "comb" in avail:
                    meth = "comb"
                else:
                    meth, = avail

        return self.rules[(wp, fid, meth)] if meth else []

    def get(self, wp, b=None, c=None, l=None):
        wp = WP(wp)

        return BTagSFGroup({
            fid: self._get(wp, fid, v)
            for fid, v in {FID.b: b, FID.c: c, FID.l: l}.items()
            if v is not False
        })

    def fids(self, wp):
        wp = WP(wp)
        return set(k[1] for k in self.rules.keys() if k[0] == wp)

    def methods(self, wp, fid):
        key = WP(wp), FID(fid)
        return set(k[2] for k in self.rules.keys() if k[:2] == key)


class BTagSFGroup:
    def __init__(self, rules):
        assert all(WP(wp) == wp for wp in rules.keys())
        assert all("central" in val for val in rules.values())
        assert all(all(isinstance(v, BTagSFRule) for v in val.values()) for val in rules.values())
        self.rules = rules

    @property
    def shifts(self):
        return set().union(*(val.keys() for val in self.rules.values()))

    def reduce(self, shifts=(), updown=("up_%s", "down_%s")):
        updown = updown or ("%s",)
        shifts = sum([[f % p for f in updown] for p in shifts], ["central"])
        return type(self)(
            {
                wp: {syst: rule for syst, rule in systs.items() if syst in shifts}
                for wp, systs in self.rules.items()
            }
        )

    def eval_awk(self, obj, discr, pt="pt", eta="eta", flavor="hadronFlavour"):
        def devirt(arr):
            return getattr(arr, "array", arr)
        obj = devirt(obj)
        args = tuple(getattr(obj, a) for a in (flavor, eta, pt, discr))
        ref = args[-1]
        args = tuple(devirt(a.content) for a in args)
        assert all(isinstance(a, np.ndarray) for a in args)
        return {syst: ref.copy(content=sf) for syst, sf in self.eval(*args).items()}

    def eval(self, flavor, eta, pt, discr):
        systs = self.shifts - {"central"}
        systs = ["central"] + sorted(systs)

        fl = np.abs(flavor)
        fl[fl > 5] = 3  # gluons are light
        fl[fl < 3] = 3
        mask = {fid: fl == (5 - fid) for fid in self.rules.keys()}
        del fl

        args = {
            fid: dict(
                eta=eta[m],
                pt=pt[m],
                discr=discr[m],
            )
            for fid, m in mask.items()
        }

        ret = {}
        oob = np.empty_like(pt, dtype=bool)
        for syst in systs:

            sf = np.ones_like(pt)
            oob[:] = False
            for fid, rules in self.rules.items():
                if syst in rules:
                    m = mask[fid]
                    sf[m], oob[m] = rules[syst].eval_oob(**args[fid])

            if syst != "central" and np.any(oob):
                sf[oob] *= 2
                sf[oob] -= ret["central"][oob]

            ret[syst] = sf

        return ret


class BTagSFRule(object):
    pt_bias = 1e-4
    dtype = np.float64

    def __init__(self, entries, isShape):
        self.isShape = bool(isShape)

        edge1d = np.array([e[:-1] for e in entries], dtype=self.dtype)
        edge1d.shape = (-1, 3, 2)

        self.edge3d = [np.unique(edge1d[:, i]) for i in range(3)]
        Eeta, Ept, Ediscr = Etotal = tuple(e.shape[0] - 1 for e in self.edge3d)
        self.bounds_eta = np.full((Ediscr, 2), [np.inf, -np.inf], dtype=self.dtype)
        self.bounds_pt = np.full((Eeta, Ediscr, 2), [np.inf, -np.inf], dtype=self.dtype)
        self.abseta = np.all(self.edge3d[0] >= 0).item()

        self.expr, self.func, coeff = self.parse(e[-1] for e in entries)

        # if coeff is None:
        #     coeff = np.zeros(edge1d.shape[:1] + (0,), dtype=self.dtype)
        self.coeff = np.zeros(Etotal + coeff.shape[1:])

        for c, bounds in zip(coeff, edge1d):
            # Seta, Spt, Sdiscr = Stotal = tuple(slice(l, h) for l, h in zip(
            #     self.idx(*bounds[..., 0], side="right"),
            #     self.idx(*bounds[..., 1], side="left"),
            # ))
            Beta, Bpt, C = self.get(*bounds, slice=True)
            C[...] = c[None, None, None, ...]

            for bound, val in (
                (Beta, bounds[0, :]),
                (Bpt, bounds[1, :]),
            ):
                for idx, func in enumerate([np.minimum, np.maximum]):
                    bound[..., idx] = func(bound[..., idx], val[idx])

        self.bounds_eta = self._bounds_reduce(self.bounds_eta, np.s_[:1])
        self.bounds_pt = self._bounds_reduce(self.bounds_pt, np.s_[:1])
        self.bounds_pt = self._bounds_reduce(self.bounds_pt, np.s_[:, :1])

    def _bounds_reduce(self, bounds, sli):
        ref = bounds[sli]
        return ref if np.all(bounds == ref) else bounds

    def idx(self, eta, pt, discr=0, side="right"):
        return tuple(
            np.clip(np.searchsorted(e, x, side=side) - 1, 0, len(e) - 2)
            for x, e in zip(np.broadcast_arrays(eta, pt, discr), self.edge3d)
        )

    def sli(self, eta, pt, discr=0):
        return tuple(
            slice(l, h + 1)
            for l, h in zip(
                self.idx(eta[..., 0], pt[..., 0], discr[..., 0], side="right"),
                self.idx(eta[..., 1], pt[..., 1], discr[..., 1], side="left"),
            )
        )

    def get(self, eta, pt, discr=0, slice=False):
        Xeta, Xpt, Xdiscr = Xtotal = (self.sli if slice else self.idx)(eta, pt, discr)
        Seta = Xdiscr
        if not slice:
            if self.bounds_eta.shape[0] == 1:
                Seta = 0
            if self.bounds_pt.shape[0] == 1:
                Xeta = 0
                if self.bounds_pt.shape[1] == 1:
                    Xdiscr = 0
            elif self.bounds_pt.shape[1] == 1:
                Xdiscr = np.s_[:1]
        return self.bounds_eta[Seta], self.bounds_pt[Xeta, Xdiscr], self.coeff[Xtotal]

    def eval(self, eta, pt, discr=0):
        return self.eval_oob(eta, pt, discr)[0]

    def eval_oob(self, eta, pt, discr=0):
        if self.abseta:
            eta = np.abs(eta)
        Beta, Bpt, C = self.get(eta, pt, discr)

        # https://github.com/cms-sw/cmssw/blob/master/CondTools/BTau/src/BTagCalibrationReader.cc#L142-L152
        pt_low = pt <= Bpt[..., 0]
        pt_high = pt > Bpt[..., 1]
        oob = pt_low | pt_high
        if np.any(oob):
            pt = np.copy(pt)
            pt[pt_low] = Bpt[... if Bpt.ndim == 1 else pt_low, 0] + self.pt_bias
            pt[pt_high] = Bpt[... if Bpt.ndim == 1 else pt_high, 1] - self.pt_bias

        ret = self.func(discr if self.isShape else pt, C.T)

        # https://github.com/cms-sw/cmssw/blob/master/CondTools/BTau/src/BTagCalibrationReader.cc#L124-L140
        ret[(eta <= Beta[..., 0]) | (eta > Beta[..., 1])] = 1

        return ret, oob

    class _spce(list):
        x = "x"
        vs = ascii_lowercase.replace(x, "")

        def sub(self, val):
            idx = len(self)
            self.append(val)
            return sympy.Symbol(self.vs[idx])

        @classmethod
        def extract(cls, expr):
            ret = cls()
            return expr.replace(attrgetter("is_Number"), ret.sub), ret

        @classmethod
        def symbols(cls, num):
            return sympy.symbols(" ".join(cls.vs[:num]), seq=True)

    @classmethod
    def parse(cls, raw):
        expr, coeffs = cls._parse(raw)
        x = sympy.Symbol(cls._spce.x)
        c = sorted(expr.free_symbols - {x}, key=attrgetter("name"))
        func = sympy.lambdify([x, tuple(c)], expr)
        return expr, func, coeffs

    @classmethod
    def _horner_test(cls, expr):
        p = expr.as_poly(sympy.Symbol(cls._spce.x))
        return p.degree() > 1 if p else False

    @classmethod
    def _parse(cls, raw):
        raw = list(raw)
        x = sympy.Symbol(cls._spce.x)
        exprs = list(map(sympy.parse_expr, raw))

        # attemps polynom extraction
        try:
            # polys = [sympy.poly(f, x, domain="RR") for f in exprs]
            polys = sympy.parallel_poly_from_expr(exprs, x, domain="RR")[0]
        except sympy.PolynomialError:
            pass
        else:
            ncoeff = max(len(p.all_coeffs()) for p in polys)
            coeffs = np.zeros((len(polys), ncoeff), dtype=cls.dtype)
            coeffs[:, 0] = 1
            for i, p in enumerate(polys):
                c = p.all_coeffs()[::-1]
                coeffs[i, : len(c)] = c

            expr = sympy.horner(sympy.Poly.from_list(cls._spce.symbols(ncoeff)[::-1], x))
            return expr, coeffs

        # extraction
        for formulas in list(map(sympy.simplify, exprs)), exprs:

            expr = set()
            coeffs = []
            for e, c in map(cls._spce.extract, formulas):
                expr.add(e)
                if len(expr) > 1:
                    break
                coeffs.append(c)
            else:
                (expr,) = list(expr)
                coeff0 = np.array(coeffs[0], dtype=object)
                coeffs = np.array(coeffs, dtype=cls.dtype)
                const = np.all(coeffs[:1] == coeffs, axis=0)
                symbols = np.array(cls._spce.symbols(len(const)))
                # put back constants
                expr = expr.xreplace(dict(zip(symbols[const], coeff0[const])))
                # shift symbols
                for old, new in zip(symbols[~const], symbols):
                    if old != new:
                        assert new not in expr.free_symbols
                        expr.replace(old, new)
                # shift coeffs
                coeffs = coeffs[:, ~const]

                # horner-ify
                expr = expr.replace(cls._horner_test, sympy.horner)

                return expr, coeffs

        return RuntimeError("\t\n".join(["can't unify expression:"] + raw))
