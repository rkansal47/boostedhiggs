import uproot4
import awkward1 as ak

# %matplotlib inline
import matplotlib.pyplot as plt
import mplhep as hep
import coffea.hist as hist

import numpy as np
from cycler import cycler

# from skhep.math.vectors import LorentzVector
from coffea.nanoevents.methods.vector import PtEtaPhiELorentzVector

from tqdm import tqdm
from os import listdir

from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)


plt.style.use(hep.style.ROOT)
def init_hists(vars, labels, binrs, fatJet=True, name=None):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    lpf = ["Leading Jet", "Sub-Leading Jet", "3rd Leading Jet"] if fatJet else ["Gen Higgs 1", "Gen Higgs 2"]
    nbins = 3 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        if type(binrs[0]) is not list:
            temp = binrs
            binrs = []
            for i in range(len(vars)):
                binrs.append(temp)

        for i in range(len(vars)):
            for j in range(nbins):
                bin = hist.Bin("bin", r"{} {}".format(lpf[j], labels[i]), *binrs[i])
                hists[name + str(j) + vars[i]] = hist.Hist("Events", hist.Cat("sample", "Sample"), bin)
    else:
        for j in range(nbins):
            bin = hist.Bin("bin", r"{} {}".format(lpf[j], labels[i]), *binrs[i])
            hists[name + str(j) + vars] = hist.Hist("Events", hist.Cat("sample", "Sample"), bin)





def init_hists(vars, labels, binrs, fatJet=True, name=None):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    lpf = ["Leading Jet", "Sub-Leading Jet", "3rd Leading Jet"] if fatJet else ["Gen Higgs 1", "Gen Higgs 2"]
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        if type(binrs[0]) is not list:
            temp = binrs
            binrs = []
            for i in range(len(vars)):
                binrs.append(temp)

        for i in range(len(vars)):
            bins = []
            for j in range(nbins):
                bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels[i]), *binrs[i]))

            hists[name + vars[i]] = hist.Hist("Events", hist.Cat("sample", "Sample"), *bins)
    else:
        bins = []
        for j in range(nbins):
            bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels), *binrs))

        hists[name + vars] = hist.Hist("Events", hist.Cat("sample", "Sample"), *bins)


weights = {
    "HH4V": None,
    "HHbbWWqq": evtDict["HHbbWWqq"]['weight'],
    "QCD": evtDict["QCD"]['weight'],
    "tt": evtDict["tt"]['weight'],
}


scale_factor = {
        'HH4V': 1 / len(evtDict["HH4V"]["weight"]),
        'HHbbWWqq': 1 / np.sum(evtDict["HHbbWWqq"]["weight"]),
        'HHbbWWqq - Hbb': 2 / np.sum(evtDict["HHbbWWqq"]["weight"]),
        'HHbbWWqq - HWW': 2 / np.sum(evtDict["HHbbWWqq"]["weight"]),
        'QCD': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
        'tt': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
        }


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        for i in range(len(vars)):
            if not type(vars[i]) is dict:
                temp = vars[i]
                vars[i] = {}
                for s in evtDict.keys():
                    vars[i][s] = temp

            for s, evts in evtDict.items():
                if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                    if not fatJet or s != "HHbbWWqq":
                        kwargs = {}
                        for j in range(nbins):
                            kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[i][s])]
                        hists[name + vars[i]["HHbbWWqq"]].fill(sample=s, weight=weights[s], **kwargs)
                    else:
                        # kwargs = {}
                        # kwargs["jet1"] = evts["fatJet1{}".format(vars[i][s])][fatJet1W]
                        # kwargs["jet2"] = evts["fatJet2{}".format(vars[i][s])][fatJet2W]
                        hists[name + vars[i]["HHbbWWqq"]].fill(sample=s + " - HWW", weight=weights[s][fatJet1W], jet1=evts["fatJet1{}".format(vars[i][s])][fatJet1W])
                        hists[name + vars[i]["HHbbWWqq"]].fill(sample=s + " - HWW", weight=weights[s][fatJet2W], jet1=evts["fatJet2{}".format(vars[i][s])][fatJet2W])

                        # kwargs["jet1"] = evts["fatJet1{}".format(vars[i][s])][~fatJet1W]
                        # kwargs["jet2"] = evts["fatJet2{}".format(vars[i][s])][~fatJet2W]
                        # hists[name + vars["HHbbWWqq"]].fill(sample=s + " - Hbb", weight=weights[s], **kwargs)

                        hists[name + vars[i]["HHbbWWqq"]].fill(sample=s + " - Hbb", weight=weights[s][~fatJet1W], jet1=evts["fatJet1{}".format(vars[i][s])][~fatJet1W])
                        hists[name + vars[i]["HHbbWWqq"]].fill(sample=s + " - Hbb", weight=weights[s][~fatJet2W], jet1=evts["fatJet2{}".format(vars[i][s])][~fatJet2W])

            if scale: hists[name + vars[i]["HHbbWWqq"]].scale(scale_factor, axis='sample')
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        for s, evts in evtDict.items():
            hists[name + vars["HHbbWWqq"]].fill(

            if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                if not fatJet or s != "HHbbWWqq":
                    kwargs = {}
                    for j in range(nbins):
                        kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[s])]
                    hists[name + vars["HHbbWWqq"]].fill(sample=s, weight=weights[s], **kwargs)
                else:
                    kwargs = {}
                    kwargs["jet1"] = evts["fatJet1{}".format(vars[s])] * fatJet1W
                    kwargs["jet2"] = evts["fatJet2{}".format(vars[s])] * fatJet2W
                    # kwargs["jet3"] = ak.zeros_like(fatJet1W)
                    hists[name + vars["HHbbWWqq"]].fill(sample=s + " - HWW", weight=weights[s], **kwargs)

                    kwargs = {}
                    kwargs["jet1"] = evts["fatJet1{}".format(vars[s])] * ~fatJet1W
                    kwargs["jet2"] = evts["fatJet2{}".format(vars[s])] * ~fatJet2W
                    # kwargs["jet3"] = ak.zeros_like(fatJet1W)
                    hists[name + vars["HHbbWWqq"]].fill(sample=s + " - Hbb", weight=weights[s], **kwargs)


        if scale: hists[name + vars["HHbbWWqq"]].scale(scale_factor, axis='sample')


fatJet1W
fatJet2W

bb = evtDict["HHbbWWqq"]["fatJet1Pt"] * ~fatJet1W
ww = evtDict["HHbbWWqq"]["fatJet1Pt"] * fatJet1W
bb[np.array(bb).nonzero()]
ww[np.array(ww).nonzero()]

bb = evtDict["HHbbWWqq"]["fatJet2Pt"] * ~fatJet2W
ww = evtDict["HHbbWWqq"]["fatJet2Pt"] * fatJet2W
bb[np.array(bb).nonzero()]
ww[np.array(ww).nonzero()]


colors_cycle = ['#e31a1c', '#33a02c', '#a6cee3', '#1f78b4', '#b2df8a', '#fb9a99']

fill_opts = {
    'edgecolor': (0, 0, 0, 0.3),
    'linewidth': 2,
    'alpha': 0.8
}

line_opts = {
    'linestyle': '-',
    'linewidth': 3,
    'alpha': 0.8
}

hists["fatJetPt"].identifiers(axis='sample')

def plot_hists(vars, fname, fatJet=True, name=None):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        fig, axs = plt.subplots(len(vars), nbins, figsize=(nbins * 9, len(vars) * 9))
        for i in range(len(vars)):
            var = vars[i]
            for j in range(nbins):
                ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                hist.plot1d(hists[name + var][:-2].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
                axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                axs[i, j].set_ylim(0, ylim)

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()
    else:
        fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))
        for j in range(nbins):
            ylim = np.max(list(hists[name + vars].project("sample", bpf + str(j + 1)).values().values())) * 1.1
            axs[j].set_prop_cycle(cycler(color=colors_cycle))
            hist.plot1d(hists[name + vars][:-2].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
            axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
            ax = hist.plot1d(hists[name + vars][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
            axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
            axs[j].set_ylim(0, ylim)

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()
