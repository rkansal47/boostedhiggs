# coding: utf-8

from tools.evil import pin, ccall
import functools
from coffea import hist
import numpy as np


def _short_repr(cls):
    def repr(self):
        return f"{self.__class__.__name__}({', '.join([k for k in self.__dict__.keys() if not k.startswith('_')])})"

    cls.__repr__ = repr
    return cls


class Variable:
    def __init__(
        self,
        label=None,
        xtitle=None,
        ytitle="Entries",
        binning=[],
        unit=None,
        metadata={},
        register=True,
    ):
        pin(locals())

    def __call__(self, func, *args, **kwargs):
        @functools.wraps(func)
        def variable(*args, **kwargs):
            return func(*args, **kwargs)

        variable.label = self.label
        variable.xtitle = self.xtitle
        variable.ytitle = self.ytitle
        variable.binning = self.binning
        variable.unit = self.unit
        variable.metadata = self.metadata
        return variable

    @staticmethod
    def to_dense(var):
        assert var.label and var.binning
        return hist.Bin(var.__name__, var.label, var.binning[0], var.binning[1], var.binning[2])


@_short_repr
class Group:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @property
    def variables(self):
        return [k for k in self.__dict__.keys() if not k.startswith("_")]


class Vars(ccall(Group)):
    pred = np.random.random((10, 1000))

    class HighLevel(ccall(Group)):
        @Variable(
            label=r"mll",
            binning=[20, 0, 200],
            xtitle=r"$m_{ll}$",
            ytitle="Entries",
            unit="GeV",
        )
        def mll(l1, l2):
            return (l1 + l2).mass()

        @Variable(
            binning=[40, 0, 1],
            xtitle=r"DNN score",
            ytitle="Entries",
        )
        def dnn_score_max():
            return np.max(pred, axis=-1)


class LowLevel(ccall(Group)):
    # ak4 jets
    @Variable(
        label=r"jet1_pt",
        binning=[40, 0, 200],
        xtitle=r"Jet 1 $p_T$",
        unit="GeV",
    )
    def jet1_pt(jets):
        return jets.pt.pad(4, clip=True).fillna(np.nan)[:, 0].flatten()

    @Variable(
        label=r"jet2_pt",
        binning=[40, 0, 200],
        xtitle=r"Jet 2 $p_T$",
        unit="GeV",
    )
    def jet2_pt(jets):
        return jets.pt.pad(4, clip=True).fillna(np.nan)[:, 1].flatten()

    @Variable(
        label=r"jet1_btag",
        binning=[25, 0, 1],
        xtitle=r"Jet 1 $Tag^B_{{DeepFlavour}}$",
    )
    def jet1_btag(jets):
        return jets.btagDeepFlavB.pad(4, clip=True).fillna(np.nan)[:, 0].flatten()

    @Variable(
        label=r"jet2_btag",
        binning=[25, 0, 1],
        xtitle=r"Jet 2 $Tag^B_{{DeepFlavour}}$",
    )
    def jet2_btag(jets):
        return jets.btagDeepFlavB.pad(4, clip=True).fillna(np.nan)[:, 1].flatten()
