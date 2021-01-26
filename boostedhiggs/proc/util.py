# coding: utf-8

import functools
from functools import reduce
from operator import and_, or_, attrgetter
import logging

import numpy as np
import awkward
from uproot_methods import TLorentzVectorArray
from coffea.analysis_objects import JaggedCandidateArray
from coffea.nanoaod.nanoevents import NanoCollection

from boostedhiggs.utils.util import parametrized
from boostedhiggs.utils.order import NanoAOD_version


logger = logging.getLogger(__name__)


def df_object_overlap(toclean, cleanagainst, dr=0.4):
    particle_pair = toclean["p4"].cross(cleanagainst["p4"], nested=True)
    return (particle_pair.i0.delta_r(particle_pair.i1)).min() > dr


def nano_object_overlap(toclean, cleanagainst, dr=0.4):
    particle_pair = toclean.cross(cleanagainst, nested=True)
    return (particle_pair.i0.delta_r(particle_pair.i1)).min() > dr


def df_mask_or(df, masks):
    decision = reduce(or_, (df[mask] for mask in masks))
    return decision


def df_mask_and(df, masks):
    decision = reduce(and_, (df[mask] for mask in masks))
    return decision


def nano_mask_or(events, masks, skipable=()):
    masks = set(map(str, masks))
    if skipable is True:
        skipable = masks
    else:
        skipable = set(map(str, skipable))
    return reduce(
        or_,
        (getattr(events, mask) for mask in masks if mask not in skipable or hasattr(events, mask)),
    )


def nano_mask_and(events, masks):
    decision = reduce(and_, (getattr(events, str(mask)) for mask in masks))
    return decision


def nano_cut(what, *cuts):
    return what[reduce(and_, cuts)]


def reduce_and(*what):
    return reduce(and_, what)


def reduce_or(*what):
    return reduce(or_, what)


def objUnion(objs, sortBy=None):
    objs = tuple(objs)
    assert objs
    assert all(isinstance(obj, awkward.JaggedArray) for obj in objs)
    assert all(obj.shape == objs[0].shape for obj in objs)
    tags = awkward.JaggedArray.concatenate(
        [
            awkward.JaggedArray(obj.starts, obj.stops, np.full(obj.shape, i, dtype=np.intp))
            for i, obj in enumerate(objs)
        ],
        axis=1,
    )
    index = awkward.JaggedArray.concatenate(
        [
            awkward.JaggedArray(obj.starts, obj.stops, np.arange(len(obj.content), dtype=np.intp))
            for obj in objs
        ],
        axis=1,
    )
    if sortBy:
        if not callable(sortBy):
            sortBy = attrgetter(sortBy)
        shuf = awkward.JaggedArray.concatenate(list(map(sortBy, objs)), axis=1).argsort()
        tags = tags[shuf]
        index = index[shuf]
    # index = tags.copy(content=index.content) # reduces memory footprint
    return tags.copy(
        content=awkward.UnionArray(tags.content, index.content, [obj.content for obj in objs])
    )


def lazyConcat1(arrays, methods=None, name=None, cols_ignore=(), ref_col=None, check_types=True):
    """
    lazy concatenation on axis=1 for jagged arrays of tables

    ref_col: is used for output start&stops calculation, can be explicitly given to
        avoid materializing an unneded column (is optional)

    check_types: assert the same type in all concatenateable arrays. Skip (check_types=False) if you know what you do...
    """
    arrays = tuple(arrays)
    assert all(isinstance(a, awkward.JaggedArray) for a in arrays)
    assert all(isinstance(a.content, awkward.Table) for a in arrays)

    JaggedArray = NanoCollection._get_mixin(methods, awkward.JaggedArray)
    Table = NanoCollection._get_mixin(methods, awkward.Table)
    data = Table() if name is None else Table.named(name)
    types = arrays[0].content.type.to

    cols = reduce(and_, (set(a.columns) for a in arrays)) - set(cols_ignore)
    cols = set(col for col in cols if isinstance(types[col], np.dtype))

    if ref_col is None:
        ref_col = list(cols)[0]
    ref = awkward.concatenate([a[ref_col] for a in arrays], axis=1)
    tot = ref.stops.max()

    def _do(col):
        res = awkward.concatenate([a[col] for a in arrays], axis=1)
        assert (ref.starts == res.starts).all() and (ref.stops == res.stops).all()
        return res.content

    # data.type.takes = tot
    for col in cols:
        if check_types:
            assert len(set(a[col].type for a in arrays)) == 1
        data.contents[col] = awkward.VirtualArray(
            _do, args=(col,), type=awkward.type.ArrayType(tot, types[col])
        )
    data.contents[ref_col] = ref.content  # avoid recreating

    return JaggedArray(ref.starts, ref.stops, data)


def apply_PSWeight(events, weights, dataset, prefix="PSWeight"):

    w = events.PSWeight.regular()

    if NanoAOD_version(dataset) < 7:
        try:
            w *= events.LHEWeight.originalXWGTUP.regular()[:, None]
        except AttributeError:
            return

    # PS weight shifts: ISR & FSR
    ones = np.ones(events.size)
    ones.flags.writeable = False
    if w.shape[1:] == (4,):
        weights.add(
            "%s_ISR" % prefix,
            ones,
            weightUp=w[..., 2],
            weightDown=w[..., 0],
        )
        weights.add(
            "%s_FSR" % prefix,
            ones,
            weightUp=w[..., 3],
            weightDown=w[..., 1],
        )
    else:
        assert w.shape[1:] == (1,)
        weights.add("%s_ISR" % prefix, ones, ones, ones)
        weights.add("%s_FSR" % prefix, ones, ones, ones)


def top_pT_sf_formula(pt):
    return np.exp(
        -2.02274e-01 + 1.09734e-04 * pt + -1.30088e-07 * pt ** 2 + (5.83494e01 / (pt + 1.96252e02))
    )


def top_pT_reweighting(gen):
    """
    Apply this SF only to TTbar datasets!

    Documentation:
        - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        - https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
    """
    top = gen[(gen.pdgId == 6) & gen.hasFlags(["isLastCopy"])]
    anti_top = gen[(gen.pdgId == -6) & gen.hasFlags(["isLastCopy"])]
    return np.sqrt(top_pT_sf_formula(top.pt) * top_pT_sf_formula(anti_top.pt))


@parametrized
def padflat(func, n_particles=1):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        return res.pad(n_particles).fillna(np.nan).flatten()

    return wrapper


def bregcorr(jets):
    return TLorentzVectorArray.from_ptetaphim(
        jets.pt * jets.bRegCorr, jets.eta, jets.phi, jets.mass * jets.bRegCorr
    )


@padflat(1)
def m_bb(jets):
    lead_jets = jets[..., :2].distincts()
    bb = lead_jets.i0 + lead_jets.i1
    m_bb = bb.mass
    return m_bb


def get_ht(jets):
    return jets.pt.sum()


def min_dr_part1_part2(part1, part2, getn=0, fill=np.nan):
    """
    For each particle in part1 returns the minimum
    delta_r between this and each particle in part2
    """
    a, b = part1.cross(part2, nested=True).unzip()
    r = a.delta_r(b).min()
    if 0 < getn:
        r = r.pad(getn, clip=True).fillna(np.nan)
        return tuple(r[:, i] for i in range(getn))
    else:
        return r


def get_metp4(met):
    return TLorentzVectorArray.from_ptetaphim(met.pt, met.pt * 0, met.phi, met.pt * 0)


def invariant_fourvec(*particles):
    return reduce((lambda x, y: x + y), [p.sum() for p in particles])


def invariant_mass(*particles):
    return invariant_fourvec(*particles).mass


def min_dr(particles):
    di_particles = particles.distincts()
    return di_particles.i0.delta_r(di_particles.i1).min()


def min_dphi(particles):
    di_particles = particles.distincts()
    return abs(di_particles.i0.delta_phi(di_particles.i1)).min()


def get_met_ld(jets, leps, met, met_coef=0.6, mht_coef=0.4):
    mht = jets.sum() + leps.sum()
    return met_coef * met.pt + mht_coef * mht.pt


def get_cone_pt(part, n=1):
    padded_part = part.pad(n)
    return [padded_part[..., i].cone_pt.fillna(np.nan) for i in range(n)]


def pairs(objs):
    pairs = objs.distincts()
    return pairs.i0 + pairs.i1


def lead_diobj(objs):
    two = objs[..., :2]
    ret = two.sum()
    a, b = two.distincts().unzip()
    ret["deltaR"] = a.delta_r(b).min()
    ret["deltaPhi"] = a.delta_phi(b).min()
    return ret


def where_mul(cond, true, false):
    return cond * true + ~cond * false


class chunked:
    def __init__(self, func, chunksize=10000):
        self.func = func
        self.chunksize = chunksize

    def __call__(self, *args, **kwargs):
        lens = set(map(len, args))
        if len(lens) != 1:
            raise ValueError("inconsistent *args len")
        return awkward.concatenate(
            [
                self.func(*(a[off : off + self.chunksize] for a in args), **kwargs)
                for off in range(0, max(lens), self.chunksize)
            ]
        )


# https://github.com/HEP-KBFI/hh-bbww/commit/739fd138e89013ad625ec03c50557bbb4c892337
class Trigger:
    def __init__(self, HLT, trigger, config, dataset):
        assert not isinstance(dataset.aux.get("channels", None), str)
        self.HLT = HLT
        self.trigger = trigger
        self.config = config
        self.dataset = dataset

    @property
    def run(self):
        return self.dataset.aux["run"] if self.dataset.is_data else None

    @property
    def channels(self):
        return self.dataset.aux["channels"] if self.dataset.is_data else None

    def _has_channel(self, channel):
        return self.dataset.is_mc or channel.name in self.channels

    @staticmethod
    def _in_range(value, range):
        if value is None or range is all:
            return True
        if "-" in range:
            if isinstance(range, str):
                value = value.lower()
                range = range.lower()
            low, high = range.split("-")
            return low <= value <= high
        else:
            return value == range

    def _get(self, channel):
        ch = self.config.get_channel(channel)
        tr = self.trigger[ch]
        if isinstance(tr, dict):
            tr = [name for name, range in tr.items() if self._in_range(self.run, range)]
        return nano_mask_or(self.HLT, tr, skipable=True)

    def get(self, *channels):
        ret = False
        for i, ch in enumerate(map(self.config.get_channel, channels)):
            if not self._has_channel(ch):
                continue
            if self.dataset.is_data and i:
                vetos = set(channels[:i]) - set(self.channels)
            else:
                vetos = ()
            good = self._get(ch)
            for veto in vetos:
                good = good & ~self._get(veto)
            ret = ret | good
        return ret
