import uproot4
import awkward1 as ak

%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import mplhep as hep
from matplotlib.patches import Rectangle

import coffea.hist as hist

import numpy as np
from cycler import cycler

# from skhep.math.vectors import LorentzVector
from coffea.nanoevents.methods.vector import PtEtaPhiELorentzVector

from tqdm import tqdm
import pandas

import re

from os import listdir
from copy import copy, deepcopy
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
plt.style.use(hep.style.ROOT)



samples = {
    "HH4V": "data/new/pnet_h4q/GluGluToHHTo4V_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8_1pb_weighted.root",
    "HHbbVV4q": 'data/new/pnet_h4q/HHToBBVVToBBQQQQ_cHHH1_1pb_weighted.root',
    "QCD": "data/new/pnet_h4q/QCD_HT*.root",
    "tt": "data/new/pnet_h4q/TTTo*.root",
    # "HH4b": "data/weighted/GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8_1pb_weighted.root"
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")

# del evtDict['HH4b']

# for i in range(1, 55):
nevents = uproot4.open(samples["HHbbVV4q"])['NEvents'].values()[1]
postevents = len(evtDict["HHbbVV4q"]['weight'])

nevents
postevents


data17 = uproot4.concatenate('data/new/pnet_h4q/JetHT*.root:tree')
evtDict["data"] = data17

# trigger_effs = uproot4.open('triggers/JetHTTriggerEfficiency_2017.root')
# trigger_effs.keys()
# trigger_effs['efficiency_ptmass_Xbb0p98To1p0'].to_numpy()


triggers17 = ['HLT_PFJet500',
                'HLT_AK8PFJet500',
                'HLT_AK8PFJet360_TrimMass30',
                'HLT_AK8PFJet380_TrimMass30',
                'HLT_AK8PFJet400_TrimMass30',
                'HLT_AK8PFHT800_TrimMass50',
                'HLT_AK8PFJet330_PFAK8BTagCSV_p17']

evtDict["HHbbVV4q"].fields

# Weights, all values in pb
RUN2LUMI = 137000
LUMI17 = 40000
XSECHHBBVV4Q = 31.05e-3 * 0.5807 * (0.2154 * 0.676 ** 2 + 0.02643 * 0.692 ** 2) * 2
XSECHH4V = 31.05e-3 * (0.2154 + 0.02643) ** 2
XSECHHBBBB = 31.05e-3 * 0.5807**2
ACCEPTANCE = 0.08239891267414204 * 0.27748453608247425
BBVV_ACC = 0.10435006608924793
HHBBWW4Q_NEVENTS = 185417
HH4B_NEVENTS = 807198
DATA_MC_SF = 0.7709694064237897
QCD_MC_SF = 0.7677139264268172

XSECHH4V
XSECHHBBVV4Q

LUMI = LUMI17

weights = {}

weights["data"] = np.ones(len(data17))
weights["HH4V"] = evtDict["HH4V"]["totalWeight"] * LUMI * XSECHH4V * ACCEPTANCE
weights["QCD"] = evtDict["QCD"]["totalWeight"] * LUMI * QCD_MC_SF
weights["tt"] = evtDict["tt"]["totalWeight"] * LUMI
weights["HHbbVV4q"] = evtDict["HHbbVV4q"]["totalWeight"] * LUMI * XSECHHBBVV4Q
# weights["HH4b"] = evtDict["HH4b"]["totalWeight"] * LUMI  # * ACCEPTANCE
# weights["fullHH4b"] = 0.0299 * np.ones(len(full_hh4b_samples['genWeight'])) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])
# weights["fullHH4b"] = ((full_hh4b_samples['genWeight'] > 0) * 2 - 1) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])

ak.sum(weights["QCD"])
# ak.sum(weights["HH4b"])
ak.sum(weights["HHbbVV4q"])
ak.sum(weights["HH4V"])
ak.sum(weights["tt"])

ak.sum(weights["tt"]) + ak.sum(weights["QCD"])
ak.sum(weights["data"])




# Jet Matching By Geometry

keys = ["genHiggs1W1", "genHiggs1W2", "genHiggs2W1", "genHiggs2W2"]
fvecs = {}

for key in keys:
    fvecs[key] = ak.zip({
        "pt": evtDict['HHbbVV4q'][key + 'Pt'],
        "eta": evtDict['HHbbVV4q'][key + 'Eta'],
        "phi": evtDict['HHbbVV4q'][key + 'Phi'],
        "mass": evtDict['HHbbVV4q'][key + 'M'],
    }, with_name="PtEtaPhiMLorentzVector")


genHiggs1 = ak.zip({"pt": evtDict['HHbbVV4q']['genHiggs1Pt'], "eta": evtDict['HHbbVV4q']['genHiggs1Eta'], "phi": evtDict['HHbbVV4q']['genHiggs1Phi'], "mass": ak.full_like(evtDict['HHbbVV4q']['genHiggs1Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")
genHiggs2 = ak.zip({"pt": evtDict['HHbbVV4q']['genHiggs2Pt'], "eta": evtDict['HHbbVV4q']['genHiggs2Eta'], "phi": evtDict['HHbbVV4q']['genHiggs2Phi'], "mass": ak.full_like(evtDict['HHbbVV4q']['genHiggs2Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")

fatJet1 = ak.zip({"pt": evtDict['HHbbVV4q']['fatJet1Pt'], "eta": evtDict['HHbbVV4q']['fatJet1Eta'], "phi": evtDict['HHbbVV4q']['fatJet1Phi'], "mass": evtDict['HHbbVV4q']['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet2 = ak.zip({"pt": evtDict['HHbbVV4q']['fatJet2Pt'], "eta": evtDict['HHbbVV4q']['fatJet2Eta'], "phi": evtDict['HHbbVV4q']['fatJet2Phi'], "mass": evtDict['HHbbVV4q']['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet3 = ak.zip({"pt": evtDict['HHbbVV4q']['fatJet3Pt'], "eta": evtDict['HHbbVV4q']['fatJet3Eta'], "phi": evtDict['HHbbVV4q']['fatJet3Phi'], "mass": evtDict['HHbbVV4q']['fatJet3Mass']}, with_name="PtEtaPhiMLorentzVector")


ak15fatJet1 = ak.zip({"pt": evtDict['HHbbVV4q']['ak15fatJet1Pt'], "eta": evtDict['HHbbVV4q']['ak15fatJet1Eta'], "phi": evtDict['HHbbVV4q']['ak15fatJet1Phi'], "mass": evtDict['HHbbVV4q']['ak15fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
ak15fatJet2 = ak.zip({"pt": evtDict['HHbbVV4q']['ak15fatJet2Pt'], "eta": evtDict['HHbbVV4q']['ak15fatJet2Eta'], "phi": evtDict['HHbbVV4q']['ak15fatJet2Phi'], "mass": evtDict['HHbbVV4q']['ak15fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")


# Javier's method for matching fatJets with gen WW or bb jets

dR = 0.8

fatJet1H1 = fatJet1.delta_r(genHiggs1) < dR
fatJet1H1W = fatJet1H1 * (fatJet1.delta_r(fvecs["genHiggs1W1"]) < dR) * (fatJet1.delta_r(fvecs["genHiggs1W2"]) < dR)
fatJet1H2 = fatJet1.delta_r(genHiggs2) < dR
fatJet1H2W = fatJet1H2 * (fatJet1.delta_r(fvecs["genHiggs2W1"]) < dR) * (fatJet1.delta_r(fvecs["genHiggs2W2"]) < dR)
fatJet1W = fatJet1H1W + fatJet1H2W

fatJet2H1 = fatJet2.delta_r(genHiggs1) < dR
fatJet2H1W = fatJet2H1 * (fatJet2.delta_r(fvecs["genHiggs1W1"]) < dR) * (fatJet2.delta_r(fvecs["genHiggs1W2"]) < dR)
fatJet2H2 = fatJet2.delta_r(genHiggs2) < dR
fatJet2H2W = fatJet2H2 * (fatJet2.delta_r(fvecs["genHiggs2W1"]) < dR) * (fatJet2.delta_r(fvecs["genHiggs2W2"]) < dR)
fatJet2W = fatJet2H1W + fatJet2H2W

fatJet3H1 = fatJet3.delta_r(genHiggs1) < dR
fatJet3H1W = fatJet3H1 * (fatJet3.delta_r(fvecs["genHiggs1W1"]) < dR) * (fatJet3.delta_r(fvecs["genHiggs1W2"]) < dR)
fatJet3H2 = fatJet3.delta_r(genHiggs2) < dR
fatJet3H2W = fatJet3H2 * (fatJet3.delta_r(fvecs["genHiggs2W1"]) < dR) * (fatJet3.delta_r(fvecs["genHiggs2W2"]) < dR)
fatJet3W = fatJet3H1W + fatJet3H2W

fatJet1b = (fatJet1H1 + fatJet1H2) * ~(fatJet1W)
fatJet2b = (fatJet2H1 + fatJet2H2) * ~(fatJet2W)


dR = 0.8

ak15fatJet1H1 = ak15fatJet1.delta_r(genHiggs1) < dR
ak15fatJet1H1W = ak15fatJet1H1 * (ak15fatJet1.delta_r(fvecs["genHiggs1W1"]) < dR) * (ak15fatJet1.delta_r(fvecs["genHiggs1W2"]) < dR)
ak15fatJet1H2 = ak15fatJet1.delta_r(genHiggs2) < dR
ak15fatJet1H2W = ak15fatJet1H2 * (ak15fatJet1.delta_r(fvecs["genHiggs2W1"]) < dR) * (ak15fatJet1.delta_r(fvecs["genHiggs2W2"]) < dR)
ak15fatJet1W = ak15fatJet1H1W + ak15fatJet1H2W

ak15fatJet2H1 = ak15fatJet2.delta_r(genHiggs1) < dR
ak15fatJet2H1W = ak15fatJet2H1 * (ak15fatJet2.delta_r(fvecs["genHiggs1W1"]) < dR) * (ak15fatJet2.delta_r(fvecs["genHiggs1W2"]) < dR)
ak15fatJet2H2 = ak15fatJet2.delta_r(genHiggs2) < dR
ak15fatJet2H2W = ak15fatJet2H2 * (ak15fatJet2.delta_r(fvecs["genHiggs2W1"]) < dR) * (ak15fatJet2.delta_r(fvecs["genHiggs2W2"]) < dR)
ak15fatJet2W = ak15fatJet2H1W + ak15fatJet2H2W

ak15fatJet1b = (ak15fatJet1H1 + ak15fatJet1H2) * ~(ak15fatJet1W)
ak15fatJet2b = (ak15fatJet2H1 + ak15fatJet2H2) * ~(ak15fatJet2W)



# Matching ak15 fatjets with ak8 fatjets

dR = 1.0

ak15fj1ak8fj1 = ak15fatJet1.delta_r(fatJet1) < dR
ak15fj1ak8fj2 = ak15fatJet1.delta_r(fatJet2) < dR

ak15fj2ak8fj1 = ak15fatJet2.delta_r(fatJet1) < dR
ak15fj2ak8fj2 = ak15fatJet2.delta_r(fatJet2) < dR

ak.sum(ak15fj1ak8fj1)
ak.sum(ak15fj1ak8fj2)
ak.sum(ak15fj1ak8fj1 * ak15fj1ak8fj2)
ak.sum(ak15fj1ak8fj1) + ak.sum(ak15fj1ak8fj2) - (ak.sum(ak15fj1ak8fj1 * ak15fj1ak8fj2))

ak.sum(ak15fj2ak8fj1)
ak.sum(ak15fj2ak8fj2)
ak.sum(ak15fj2ak8fj1 * ak15fj2ak8fj2)
ak.sum(ak15fj2ak8fj1) + ak.sum(ak15fj2ak8fj2) - (ak.sum(ak15fj2ak8fj1 * ak15fj2ak8fj2))

ak.sum(ak15fj1ak8fj1 * ak15fj2ak8fj1)
ak.sum(ak15fj1ak8fj2 * ak15fj2ak8fj2)

ak15fj1both = ak15fj1ak8fj1 * ak15fj1ak8fj2
ak15fj1ak8fj1 = ak15fj1ak8fj1 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) <= ak15fatJet1.delta_r(fatJet2)) * ak15fj1both
ak15fj1ak8fj2 = ak15fj1ak8fj2 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) > ak15fatJet1.delta_r(fatJet2)) * ak15fj1both

ak15fj2both = ak15fj2ak8fj1 * ak15fj2ak8fj2
ak15fj2ak8fj1 = ak15fj2ak8fj1 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) <= ak15fatJet2.delta_r(fatJet2)) * ak15fj2both
ak15fj2ak8fj2 = ak15fj2ak8fj2 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) > ak15fatJet2.delta_r(fatJet2)) * ak15fj2both

ak.sum(ak15fj1ak8fj1)
ak.sum(ak15fj1ak8fj2)
ak.sum(ak15fj1ak8fj1 * ak15fj1ak8fj2)
ak.sum(ak15fj1ak8fj1) + ak.sum(ak15fj1ak8fj2) - (ak.sum(ak15fj1ak8fj1 * ak15fj1ak8fj2))

ak.sum(ak15fj2ak8fj1)
ak.sum(ak15fj2ak8fj2)
ak.sum(ak15fj2ak8fj1 * ak15fj2ak8fj2)
ak.sum(ak15fj2ak8fj1) + ak.sum(ak15fj2ak8fj2) - (ak.sum(ak15fj2ak8fj1 * ak15fj2ak8fj2))

ak8fj1match = ak15fj1ak8fj1 + ak15fj2ak8fj1
ak8fj2match = ak15fj1ak8fj2 + ak15fj2ak8fj2

ak.sum(ak8fj1match)
ak.sum(ak8fj2match)





hists = {}


def init_hists(vars, labels, binrs, fatJet=True, name=None, scale=True, bbsorted=False, wwsorted=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    lpf = ["Leading Jet", "Sub-Leading Jet", "3rd Leading Jet"] if fatJet else ["Gen Higgs 1", "Gen Higgs 2"]
    if bbsorted: lpf = ["Hbb Candidate Jet", "HWW Candidate Jet"]
    if wwsorted: lpf = ["Leading H4q Candidate Jet", "Sub-Leading H4q Candidate Jet"]
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    ename = "Events (Unit Normalized)" if scale else "Events"

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

            hists[name + vars[i]] = hist.Hist(ename, hist.Cat("sample", "Sample"), *bins)
    else:
        bins = []
        for j in range(nbins):
            bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels), *binrs))

        hists[name + vars] = hist.Hist(ename, hist.Cat("sample", "Sample"), *bins)


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True, pevtDict=None, useEDWeights=False, data=False, ak15=False, jet2ak15=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf
    if pevtDict is None: pevtDict = evtDict
    if useEDWeights:
        pweights = {}
        for s, evts in pevtDict.items():
            pweights[s] = evts.weight
    else: pweights = weights

    if type(vars) is list:
        for i in range(len(vars)):
            if not type(vars[i]) is dict:
                temp = vars[i]
                vars[i] = {}
                for s in pevtDict.keys():
                    vars[i][s] = temp

            for s, evts in pevtDict.items():
                if (fatJet or "HH" in s) and (s != "HH4V" or hh4v) and (s != 'data' or data):
                    kwargs = {}
                    for j in range(nbins):
                        if fatJet and (ak15 or (j == 1 and jet2ak15)): rnpf = "ak15" + npf
                        else: rnpf = npf
                        kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(rnpf, j + 1, vars[i][s])]
                    hists[name + vars[i]["HHbbVV4q"]].fill(sample=s, weight=pweights[s], **kwargs)

            if scale: hists[name + vars[i]["HHbbVV4q"]].scale(scale_factor, axis='sample')
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in pevtDict.keys():
                vars[s] = temp

        for s, evts in pevtDict.items():
            if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                if fatJet and ak15: rnpf = "ak15" + npf
                else: rnpf = npf
                kwargs = {}
                for j in range(nbins):
                    kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(rnpf, j + 1, vars[s])]
                hists[name + vars["QCD"]].fill(sample=s, weight=pweights[s], **kwargs)

        if scale: hists[name + vars["QCD"]].scale(scale_factor, axis='sample')


colors_cycle = ['#e31a1c', '#33a02c', '#1f78b4', '#a6cee3', '#b2df8a', '#fb9a99', '#f1005f', '#d09d09']

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

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}

regbg = re.compile('^(QCD|tt)$')
regsig = re.compile('^(HH4V.*|HHbbVV4q)$')
regnotdata = re.compile('^(?!data).*')
regmc = re.compile('^(QCD|tt|HHbbVV4q)$')
reghh4v = re.compile('HH4V.*')


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True, hh4b=False, sepsig=True, stackall=False, log=False, lumilabel=None, data=False, ratio=True, rat_ylim=1, ak15=False, blinding=None):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        fig, axs = plt.subplots(len(vars), nbins, figsize=(nbins * 9, len(vars) * 9))
        for i in range(len(vars)):
            if not type(vars[i]) is dict:
                temp = vars[i]
                vars[i] = {}
                for s in evtDict.keys():
                    vars[i][s] = temp

            var = vars[i]["HHbbVV4q"]
            if type(binrs[0]) is not list:
                temp = binrs
                binrs = []
                for k in range(len(vars)):
                    binrs.append(temp)

            if fatJet:
                if sepsig:
                    rnpf = "ak15" + npf if ak15 else npf
                    if ak15:
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    else:
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)


                    for j in range(nbins):
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        if hh4v: hist.plot1d(hists[name + var][reghh4v].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)

                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                        axs[i, j].set_ylim(0, ylim)
                else:
                    for j in range(nbins):
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][regsig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                        axs[i, j].set_ylim(0, ylim)

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        var = vars["HHbbVV4q"]
        if data and ratio: fig, (axs, rax) = plt.subplots(2, nbins, figsize=(nbins * 9, 12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        else: fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))

        if fatJet:
            if sepsig:
                rnpf = "ak15" + npf if ak15 else npf
                if ak15:
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                else:
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)

                for j in range(nbins):
                    ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                    axs[j].set_prop_cycle(cycler(color=colors_cycle))
                    if hh4v: hist.plot1d(hists[name + var][reghh4v].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)

                    axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                    hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                    axs[j].set_ylim(0, ylim)
            else:
                for j in range(nbins):
                    if data:
                        # hist.plot1d(hists[name + var][regbg].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['tt', 'QCD'])
                        hist.plot1d(hists[name + var]['QCD'].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['QCD'])
                        hist.plot1d(hists[name + var]['data'].project("sample", bpf + str(j + 1)), ax=axs[j], clear=False, error_opts=data_err_opts)
                        ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_ylim(0, ylim)

                        if blinding is not None:
                            axs[j].add_patch(Rectangle((blinding[0] - 1, 0), blinding[1] - blinding[0] + 2, 1e12, facecolor="white", zorder=100))

                        if ratio:
                            axs[j].set_xlabel(None)
                            hist.plotratio(num=hists[name + var]['data'].project("sample", bpf + str(j + 1)).sum("sample"), denom=hists[name + var][regbg].project("sample", bpf + str(j + 1)).sum("sample"), ax=rax[j], error_opts=data_err_opts, unc='num')
                            rax[j].set_ylabel('data/MC')
                            rax[j].set_ylim(0, rat_ylim)
                            rax[j].grid()
                            if blinding is not None:
                                rax[j].add_patch(Rectangle((blinding[0] - 1, 0), blinding[1] - blinding[0] + 2, 1e12, facecolor="white", zorder=100))
                    elif stackall:
                        hist.plot1d(hists[name + var][regnotdata].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbVV4q', 'tt', 'QCD'])
                    else:
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][regsig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
                        axs[j].set_prop_cycle(cycler(color=colors_cycle[3::-1]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts, order=['tt', 'QCD'])
                        axs[j].set_ylim(0, ylim)

                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

                    if log:
                        axs[j].set_ylim(1)
                        axs[j].set_yscale('log')

                    if lumilabel is not None: hep.label.lumitext(str(lumilabel) + "fb$^{-1}$", ax=axs[j])

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()


def ak815_comp_plot(vars, varsl, bins, ak8events, ak15events, name, qcd=True, lumilabel=40):
    for i in range(len(vars)):
        histname = "ak8_ak15_{}_comp".format(vars[i])
        hists[histname] = hist.Hist("Events",
                                        hist.Cat("sample", "Sample"),
                                        hist.Bin("jet1", r"Jet {}".format(varsl[i]), *bins[i]),
                                        hist.Bin("jet2", r"Jet {}".format(varsl[i]), *bins[i]),
                                        )

        for s in ["HHbbVV4q", "QCD"]:
            hists[histname].fill(sample=s + " AK8",
                                 jet1=ak8events[s]["fatJet1" + vars[i]],
                                 jet2=ak8events[s]["fatJet2" + vars[i]],
                                 weight=ak8events[s]["weight"]
                                 )
            hists[histname].fill(sample=s + " AK15",
                                 jet1=ak15events[s]["ak15fatJet1" + vars[i]],
                                 jet2=ak15events[s]["ak15fatJet2" + vars[i]],
                                 weight=ak15events[s]["weight"]
                                 )


    fig, axs = plt.subplots(len(vars), 2, figsize=(len(vars) * 9, 18))

    xlabs = ['Hbb', 'HWW']

    loc = [7, 2] if qcd else [1, 1]

    for i in range(2):
        for j in range(len(vars)):
            histname = "ak8_ak15_{}_comp".format(vars[j])
            # print(str(i) + ", " + str(j))
            hist.plot1d(hists[histname]["HHbbVV4q AK15"].project("sample", "jet" + str(i + 1)), clear=False, line_opts=line_opts, ax=axs[j, i])
            hist.plot1d(hists[histname]["HHbbVV4q AK8"].project("sample", "jet" + str(i + 1)), clear=False, line_opts=line_opts, ax=axs[j, i])
            # axs[j, i].hist(events_ak15bb_sorted[sample]["ak15fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='red', weights=events_ak15bb_sorted[sample]['weight'], label='AK15', linewidth=2)
            if qcd:
                axqcd = axs[j, i].twinx()
                axqcd.set_prop_cycle(cycler(color=colors_cycle[0:]))
                hist.plot1d(hists[histname]["QCD AK15"].project("sample", "jet" + str(i + 1)), clear=False, line_opts={"linestyle": "dotted", "linewidth": 5}, ax=axqcd)
                hist.plot1d(hists[histname]["QCD AK8"].project("sample", "jet" + str(i + 1)), clear=False, line_opts={"linestyle": "dotted", "linewidth": 5}, ax=axqcd)
                axqcd.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16}, loc=1)
                axqcd.set_ylabel("QCD Events")
            # axqcd.hist(events_bb_sorted['QCD']["fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='green', weights=events_bb_sorted['QCD']['weight'], label='AK8', linewidth=3, linestyle='dashed')
            # axqcd.hist(events_ak15bb_sorted['QCD']["ak15fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='red', weights=events_ak15bb_sorted['QCD']['weight'], label='AK15', linewidth=3, linestyle='dashed')
            axs[j, i].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16}, loc=loc[j])
            axs[j, i].set_ylabel("HHbbVV4q Events")
            axs[j, i].set_xlabel("{} Candidate Jet {}".format(xlabs[i], varsl[j]))
            axs[j, i].set_xlim([bins[j][1], bins[j][2]])
            # hep.label.lumitext(, ax=axs[j, i])

            axs[j, i].text(
                x=0.15,
                y=1.005,
                s=str(lumilabel) + "fb$^{-1}$",
                transform=axs[j, i].transAxes,
                ha="right",
                va="bottom",
                fontsize=plt.rcParams["font.size"] * 0.95,
                fontweight="normal",
                # fontname=fontname,
            )

    plt.tight_layout(pad=2.0)
    plt.savefig("figs/{}.pdf".format(name))
    plt.show()


def tagger2d_plots(histname, ak8events, ak15events, vmaxhh=0.01, vmaxqcd=250, name=None, lumilabel=40):
    if name is None: name = histname
    titles = ['HHbbVV4q bb1 vs VV2 tagger', 'QCD bb1 vs VV2 tagger', 'QCD bb1 vs bb2 tagger']

    hists[histname] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin("jet1bb", r"Jet1 Tagger Score", 20, 0.9, 1),
                                    hist.Bin("jet2WW", r"Jet2 Tagger Score", 20, 0.9, 1),
                                    )


    s = "HHbbVV4q"
    hists[histname].fill(sample="AK8 " + titles[0],
                         jet1bb = ak8events[s]["fatJet1PNetXbb"],
                         jet2WW = ak8events[s]["fatJet2PNetHqqqq"],
                         weight = ak8events[s]["weight"],
                        )

    s = 'QCD'
    hists[histname].fill(sample="AK8 " + titles[1],
                                 jet1bb = ak8events[s]["fatJet1PNetXbb"],
                                 jet2WW = ak8events[s]["fatJet2PNetHqqqq"],
                                 weight = ak8events[s]["weight"],
                                )

    hists[histname].fill(sample="AK8 " + titles[2],
                                 jet1bb = ak8events[s]["fatJet1PNetXbb"],
                                 jet2WW = ak8events[s]["fatJet2PNetXbb"],
                                 weight = ak8events[s]["weight"],
                                )


    s = "HHbbVV4q"
    hists[histname].fill(sample="AK15 " + titles[0],
                         jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                         jet2WW = ak15events[s]["ak15fatJet2PNetHqqqq"],
                         weight = ak15events[s]["weight"],
                        )

    s = 'QCD'
    hists[histname].fill(sample="AK15 " + titles[1],
                                 jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                                 jet2WW = ak15events[s]["ak15fatJet2PNetHqqqq"],
                                 weight = ak15events[s]["weight"],
                                )

    hists[histname].fill(sample="AK15 " + titles[2],
                                 jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                                 jet2WW = ak15events[s]["ak15fatJet2PNetMDXbb"],
                                 weight = ak15events[s]["weight"],
                                )

    patch_opts = {
        'cmap': 'jet',
        'vmin': 0,
        'vmax': vmaxhh,
    }

    samps = ["AK8 ", "AK15 "]

    fig, axs = plt.subplots(2, 3, figsize=(3*9, 2*8))
    for j in range(3):
        if j == 1:
            patch_opts['vmax'] = vmaxqcd
            # patch_opts['vmin'] = 1
            # patch_opts['norm'] = mpl.colors.LogNorm()
        for i in range(2):
            hist.plot2d(hists[histname][samps[i] + titles[j]].sum('sample'), 'jet1bb', ax=axs[i][j], patch_opts=patch_opts)
            axs[i][j].set_title(samps[i] + titles[j], size=24)

            axs[i][j].text(
                x=1.2,
                y=1.01,
                s=str(lumilabel) + "fb$^{-1}$",
                transform=axs[i][j].transAxes,
                ha="right",
                va="bottom",
                fontsize=plt.rcParams["font.size"] * 0.95,
                fontweight="normal",
                # fontname=fontname,
            )

    plt.tight_layout(pad=1)
    plt.savefig("figs/{}.pdf".format(name), bbox_inches='tight')
    plt.show()


plt.rcParams.update({'font.size': 16})
plt.style.use(hep.style.CMS)


init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=False)

tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


tot_events


scale_factor = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbVV4q': 1 / tot_events['HHbbVV4q'],
        # 'HH4b': 1 / tot_events['HH4b'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}



init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=True, hh4v=False)
plot_hists("Pt", "kin_tests", [40, 200, 2000], hh4v=False)


# init_hists("MassSD", "Soft Drop Mass (GeV)", [50, 1, 400], name="data", scale=False)
# fill_hists("MassSD", scale=False, data=True, name='data')
# plot_hists("MassSD", "data_masssd_rat_pre_tf_qcd_only", [50, 1, 400], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)
#
#
# init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], name="data", scale=False)
# fill_hists("Pt", scale=False, data=True, name='data')
# plot_hists("Pt", "data_pt_rat_post_tf_test", [40, 200, 2000], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)


# vars = ["Pt", "Mass", "MassSD"]
# varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
# bins = [[40, 200, 2000], [50, 1, 400], [50, 1, 400]]
# init_hists(vars, varsl, bins, fatJet=True)
# fill_hists(copy(vars), fatJet=True, scale=True)
# plot_hists(vars, "jet_kin", bins)
#
#
vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [30, 50, 200], [30, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="fine")
fill_hists(copy(vars), fatJet=True, scale=True, name="fine")
plot_hists(vars, "jet_kin_fine_new", bins, name="fine")


hists['fatJetPt'].values()

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 1000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True, ak15=True, hh4v=False)
plot_hists(vars, "jet_kin_ak15", bins, ak15=True, hh4v=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 1000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True, ak15=False, hh4v=False)
plot_hists(vars, "jet_kin_new", bins, ak15=False, hh4v=False)


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 1000], [50, 1, 250], [50, 1, 250]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True, ak15=False, hh4v=False)
plot_hists(vars, "jet_kin_new_nosep", bins, ak15=False, hh4v=False, sepsig=False)





vars = {
    # "HH4V": "PNetXbb_alt",
    "HHbbVV4q": "PNetXbb_alt",
    # "HH4b": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
    # "data": "PNetXbb",
}
disc_bin = [100, 0, 1]
init_hists("PNetXbb", "Particle Net Xbb", [100, 0, 1], fatJet=True)
fill_hists(vars, fatJet=True, scale=True, hh4v=False)
plot_hists(vars, "pnetxbb_new", [100, 0, 1], hh4v=False)


vars = "PNetMDXbb"
disc_bin = [100, 0, 1]
init_hists(vars, "Particle Net Xbb", [100, 0, 1], fatJet=True)
fill_hists(vars, fatJet=True, scale=True, hh4v=False, ak15=True)
plot_hists(vars, "pnetxbb_ak15", [100, 0, 1], hh4v=False, ak15=True)



init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, hh4v=False)
plot_hists("PNetHqqqq", "deepak8mdh4q", [100, 0, 1], hh4v=False)






evtDict["HHbbVV4q"].fields


events_bb_sorted = {}
fatjet_vars = ["Pt", "MassSD", "Mass", "PNetMDXbb", "PNetHqqqq"]

for s, evts in evtDict.items():
    if s != "HH4V" and s != 'HH4b':
        pnet_key = "PNetXbb_alt"  # if s == "HHbbVV4q" else "PNetXbb"
        jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

        ak_dict =   {
                        "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading),
                        "fatJet2PNetXbb": (evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading),
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            if var != "PNetMDXbb":
                ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)
                ak_dict["fatJet2" + var] = (evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading)
            if var != "PNetHqqqq":
                ak_dict["ak15fatJet1" + var] = (evts["ak15fatJet1" + var] * jet1_bb_leading + evts["ak15fatJet2" + var] * ~jet1_bb_leading)
                ak_dict["ak15fatJet2" + var] = (evts["ak15fatJet1" + var] * ~jet1_bb_leading + evts["ak15fatJet2" + var] * jet1_bb_leading)

        events_bb_sorted[s] = ak.zip(ak_dict)

events_bb_sorted




hists['ddthistak8'] = hist.Hist("Events",
                    hist.Cat("sample", "Sample"),
                    hist.Bin("pt", r"$p_T$ (GeV)", 25, 250, 1000),
                    hist.Bin("rho", r"$\rho$", 25, -4, -0),
                    hist.Bin("pnethqqqq", r"pnethqqqq", 1000, 0, 1),
                    )

hists['ddthistak8'].fill(sample='QCD',
                     pt = events_bb_sorted['QCD']['fatJet2Pt'],
                     rho = np.log(events_bb_sorted['QCD']['fatJet2Mass'] ** 2 / events_bb_sorted['QCD']['fatJet2Pt'] ** 2),
                     pnethqqqq = events_bb_sorted['QCD']['fatJet2PNetHqqqq'],
                     weight = events_bb_sorted['QCD']["weight"]
                     )

val_QCD = hists['ddthistak8']['QCD'].values(overflow='allnan')[('QCD',)]

qcd_maxval_temp = np.cumsum(val_QCD, axis=2)
qcd_maxval = qcd_maxval_temp[:,:,-1]
norma = qcd_maxval_temp / np.maximum(1e-10,qcd_maxval[:,:,np.newaxis])

hist_y_QCD = deepcopy(hists['ddthistak8'])
template = hist_y_QCD.sum('pnethqqqq',) #pt, rho base
hist_y_QCD.clear()
hist_y_QCD._sumw = {():norma}

res = np.apply_along_axis(lambda norma: norma.searchsorted(0.01), axis = 2, arr = norma)
res[res>1000]=1
qmap = (res - 1) / 1000
qmap

template.clear()
template._sumw = {():qmap}
template.label = 'cut at ' + str(int(1)) + '%'

hist.plot2d(template.sum('sample'), xaxis = "rho", patch_opts={'norm': mpl.colors.LogNorm(), 'vmax': 0.5})

# save(template, '../boostedhiggs/ddtmap_%s.coffea'%postfix)
plot(template, 'ddt_%i_%s'%(int(100*percentile),postfix))



plt.clf()

#ax = hist.plot2d(template, xaxis = "jet_rho",  patch_opts={"vmin":0.5, "vmax":0.99})#,xoverflow='all',yoverflow='all')
ax = hist.plot2d(template, xaxis = "rho",  patch_opts={"vmin":0.5, "vmax":1.})#,xoverflow='all',yoverflow='all')
cmstext = plt.text(0.0, 1., "CMS",fontsize=12,horizontalalignment='left',verticalalignment='bottom', fontweight='bold',transform=ax.transAxes)
addtext = plt.text(0.11, 1., "Simulation Preliminary",fontsize=10,horizontalalignment='left',verticalalignment='bottom', style='italic', transform=ax.transAxes)

plt.ylim(200,1200)
plt.xlim(-5.5,-2)

# events_bb_sorted_triggered = {}
# fatjet_vars = ["Pt", "MassSD", "Mass", "PNetMDXbb", "PNetHqqqq"]
#
# for s, evts in evtDict.items():
#     if s != "HH4V" and s != 'HH4b':
#         triggered = evts[triggers17[0]]
#         for i in range(1, len(triggers17)): triggered = triggered + evts[triggers17[i]]
#
#         pnet_key = "PNetXbb_alt"  # if s == "HHbbVV4q" else "PNetXbb"
#         jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]
#
#         ak_dict =   {
#                         "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading)[triggered],
#                         "fatJet2PNetXbb": (evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading)[triggered],
#                         "weight": weights[s][triggered],
#                     }
#
#         for var in fatjet_vars:
#             if var != "PNetMDXbb":
#                 ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)[triggered]
#                 ak_dict["fatJet2" + var] = (evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading)[triggered]
#             if var != "PNetHqqqq":
#                 ak_dict["ak15fatJet1" + var] = (evts["ak15fatJet1" + var] * jet1_bb_leading + evts["ak15fatJet2" + var] * ~jet1_bb_leading)[triggered]
#                 ak_dict["ak15fatJet2" + var] = (evts["ak15fatJet1" + var] * ~jet1_bb_leading + evts["ak15fatJet2" + var] * jet1_bb_leading)[triggered]
#
#         events_bb_sorted[s] = ak.zip(ak_dict)
#
# events_bb_sorted_triggered


events_ak15bb_sorted = {}

for s, evts in evtDict.items():
    if s != "HH4V" and s != 'HH4b':
        triggered = evts[triggers17[0]]
        for i in range(1, len(triggers17)): triggered = triggered + evts[triggers17[i]]

        pnet_key = "PNetXbb_alt"  # if s == "HHbbVV4q" else "PNetXbb"
        jet1_bb_leading = evts["ak15fatJet1PNetMDXbb"] > evts["ak15fatJet2PNetMDXbb"]

        ak_dict =   {
                        "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading),
                        "fatJet2PNetXbb": (evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading),
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            if var != "PNetMDXbb":
                ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)
                ak_dict["fatJet2" + var] = (evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading)
            if var != "PNetHqqqq":
                ak_dict["ak15fatJet1" + var] = (evts["ak15fatJet1" + var] * jet1_bb_leading + evts["ak15fatJet2" + var] * ~jet1_bb_leading)
                ak_dict["ak15fatJet2" + var] = (evts["ak15fatJet1" + var] * ~jet1_bb_leading + evts["ak15fatJet2" + var] * jet1_bb_leading)

        events_ak15bb_sorted[s] = ak.zip(ak_dict)

events_ak15bb_sorted


events_8bb_15VV_sorted = {}
fatjet_vars = ["Pt", "MassSD"]

for s, evts in evtDict.items():
    if s != "HH4V" and s != 'HH4b':
        triggered = evts[triggers17[0]]
        for i in range(1, len(triggers17)): triggered = triggered + evts[triggers17[i]]

        fatJet1 = ak.zip({"pt": evts['fatJet1Pt'], "eta": evts['fatJet1Eta'], "phi": evts['fatJet1Phi'], "mass": evts['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
        fatJet2 = ak.zip({"pt": evts['fatJet2Pt'], "eta": evts['fatJet2Eta'], "phi": evts['fatJet2Phi'], "mass": evts['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")

        ak15fatJet1 = ak.zip({"pt": evts['ak15fatJet1Pt'], "eta": evts['ak15fatJet1Eta'], "phi": evts['ak15fatJet1Phi'], "mass": evts['ak15fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
        ak15fatJet2 = ak.zip({"pt": evts['ak15fatJet2Pt'], "eta": evts['ak15fatJet2Eta'], "phi": evts['ak15fatJet2Phi'], "mass": evts['ak15fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")

        dR = 0.8

        ak15fj1ak8fj1 = ak15fatJet1.delta_r(fatJet1) < dR
        ak15fj1ak8fj2 = ak15fatJet1.delta_r(fatJet2) < dR

        ak15fj2ak8fj1 = ak15fatJet2.delta_r(fatJet1) < dR
        ak15fj2ak8fj2 = ak15fatJet2.delta_r(fatJet2) < dR

        ak15fj1both = ak15fj1ak8fj1 * ak15fj1ak8fj2
        ak15fj1ak8fj1 = ak15fj1ak8fj1 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) <= ak15fatJet1.delta_r(fatJet2)) * ak15fj1both
        ak15fj1ak8fj2 = ak15fj1ak8fj2 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) > ak15fatJet1.delta_r(fatJet2)) * ak15fj1both

        ak15fj2both = ak15fj2ak8fj1 * ak15fj2ak8fj2
        ak15fj2ak8fj1 = ak15fj2ak8fj1 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) <= ak15fatJet2.delta_r(fatJet2)) * ak15fj2both
        ak15fj2ak8fj2 = ak15fj2ak8fj2 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) > ak15fatJet2.delta_r(fatJet2)) * ak15fj2both

        # overlap when bb fatjet matches geometrically with VV fatjet:
        # fatjet1 is bb, ak15fatjet1 is VV, match
        # fatjet1 is bb, ak15fatjet2 is VV, match
        # fatjet2 is bb, ak15fatjet1 is VV, match
        # fatjet2 is bb, ak15fatjet2 is VV, match

        jet1_bb_leading = evts["fatJet1PNetXbb_alt"] > evts["fatJet2PNetXbb_alt"]
        jet1_VV_leading = evts["ak15fatJet1PNetHqqqq"] > evts["ak15fatJet2PNetHqqqq"]

        overlap =   (jet1_bb_leading * jet1_VV_leading * ak15fj1ak8fj1) + \
                    (jet1_bb_leading * ~jet1_VV_leading * ak15fj2ak8fj1) + \
                    (~jet1_bb_leading * jet1_VV_leading * ak15fj1ak8fj2) + \
                    (~jet1_bb_leading * ~jet1_VV_leading * ak15fj2ak8fj2)

        jet1_VV_cand = (jet1_VV_leading * ~overlap) + (~jet1_VV_leading * overlap)


        ak_dict =   {
                        "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading),
                        "ak15fatJet2PNetHqqqq": (evts["ak15fatJet1PNetHqqqq"] * jet1_VV_cand + evts["ak15fatJet2PNetHqqqq"] * ~jet1_VV_cand),
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            if var != "PNetMDXbb" and var != "PNetHqqqq":
                ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)

            ak_dict["ak15fatJet2" + var] = (evts["ak15fatJet1" + var] * jet1_VV_cand + evts["ak15fatJet2" + var] * ~jet1_VV_cand)

        events_8bb_15VV_sorted[s] = ak.zip(ak_dict)

events_8bb_15VV_sorted


events_8bb_15VV_sorted = {}
fatjet_vars = ["Pt", "MassSD"]

for s, evts in evtDict.items():
    if s != "HH4V" and s != 'HH4b':
        fatJet1 = ak.zip({"pt": evts['fatJet1Pt'], "eta": evts['fatJet1Eta'], "phi": evts['fatJet1Phi'], "mass": evts['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
        fatJet2 = ak.zip({"pt": evts['fatJet2Pt'], "eta": evts['fatJet2Eta'], "phi": evts['fatJet2Phi'], "mass": evts['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")

        ak15fatJet1 = ak.zip({"pt": evts['ak15fatJet1Pt'], "eta": evts['ak15fatJet1Eta'], "phi": evts['ak15fatJet1Phi'], "mass": evts['ak15fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
        ak15fatJet2 = ak.zip({"pt": evts['ak15fatJet2Pt'], "eta": evts['ak15fatJet2Eta'], "phi": evts['ak15fatJet2Phi'], "mass": evts['ak15fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")

        dR = 0.8

        ak15fj1ak8fj1 = ak15fatJet1.delta_r(fatJet1) < dR
        ak15fj1ak8fj2 = ak15fatJet1.delta_r(fatJet2) < dR

        ak15fj2ak8fj1 = ak15fatJet2.delta_r(fatJet1) < dR
        ak15fj2ak8fj2 = ak15fatJet2.delta_r(fatJet2) < dR

        ak15fj1both = ak15fj1ak8fj1 * ak15fj1ak8fj2
        ak15fj1ak8fj1 = ak15fj1ak8fj1 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) <= ak15fatJet1.delta_r(fatJet2)) * ak15fj1both
        ak15fj1ak8fj2 = ak15fj1ak8fj2 * ~ak15fj1both + (ak15fatJet1.delta_r(fatJet1) > ak15fatJet1.delta_r(fatJet2)) * ak15fj1both

        ak15fj2both = ak15fj2ak8fj1 * ak15fj2ak8fj2
        ak15fj2ak8fj1 = ak15fj2ak8fj1 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) <= ak15fatJet2.delta_r(fatJet2)) * ak15fj2both
        ak15fj2ak8fj2 = ak15fj2ak8fj2 * ~ak15fj2both + (ak15fatJet2.delta_r(fatJet1) > ak15fatJet2.delta_r(fatJet2)) * ak15fj2both

        # overlap when bb fatjet matches geometrically with VV fatjet:
        # fatjet1 is bb, ak15fatjet1 is VV, match
        # fatjet1 is bb, ak15fatjet2 is VV, match
        # fatjet2 is bb, ak15fatjet1 is VV, match
        # fatjet2 is bb, ak15fatjet2 is VV, match

        jet1_bb_leading = evts["fatJet1PNetXbb_alt"] > evts["fatJet2PNetXbb_alt"]
        jet1_VV_leading = evts["ak15fatJet1PNetHqqqq"] > evts["ak15fatJet2PNetHqqqq"]

        overlap =   (jet1_bb_leading * jet1_VV_leading * ak15fj1ak8fj1) + \
                    (jet1_bb_leading * ~jet1_VV_leading * ak15fj2ak8fj1) + \
                    (~jet1_bb_leading * jet1_VV_leading * ak15fj1ak8fj2) + \
                    (~jet1_bb_leading * ~jet1_VV_leading * ak15fj2ak8fj2)

        jet1_VV_cand = (jet1_VV_leading * ~overlap) + (~jet1_VV_leading * overlap)


        ak_dict =   {
                        "fatJet1PNetXbb": evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading,
                        "ak15fatJet2PNetHqqqq": evts["ak15fatJet1PNetHqqqq"] * jet1_VV_cand + evts["ak15fatJet2PNetHqqqq"] * ~jet1_VV_cand,
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            if var != "PNetMDXbb" and var != "PNetHqqqq":
                ak_dict["fatJet1" + var] = evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading

            ak_dict["ak15fatJet2" + var] = evts["ak15fatJet1" + var] * jet1_VV_cand + evts["ak15fatJet2" + var] * ~jet1_VV_cand

        events_8bb_15VV_sorted[s] = ak.zip(ak_dict)

events_8bb_15VV_sorted




var = "PNetMDXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists(var, varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(var, "jet_pnet_ak15", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False, ak15=True)


var = "PNetMDXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists(var, varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bb_sorted, useEDWeights=True, ak15=True)
plot_hists(var, "jet_pnet_ak15_ak8_sorted", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False, ak15=True)


init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak15bb_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_ak15bb_leading", [100, 0, 1], name="bb_sorted_cut", hh4v=False, sepsig=False)


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bb_sorted, useEDWeights=True)
plot_hists(vars, "jet_kin_bb_sorted", bins, name="bb_sorted", sepsig=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_bb_sorted_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bb_sorted, useEDWeights=True)
plot_hists(vars, "jet_kin_bb_sorted_fine", bins, name="bb_sorted", sepsig=False)

vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_bb_sorted_fine_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
ak815_comp_plot(vars, varsl, bins, events_bb_sorted, events_ak15bb_sorted, "ak815_or_cut")


var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists("PNetXbb", varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bb_sorted, useEDWeights=True)
plot_hists(var, "jet_cut_pnet2", bins, name="bb_sorted_cut", hh4v=False, sepsig=False)



vars = ["Pt", "Mass", "MassSD"]
bins = [[50, 200, 600], [50, 1, 250], [50, 1, 250]]

fig, axs = plt.subplots(2, len(vars), figsize=(len(vars) * 10, 18))

for i in range(len(vars)):
    for j in range(2):
        if j == 0:
            ak8hists = events_bb_sorted["HHbbVV4q"]["fatJet1{}".format(vars[i])] * ak8fj1match
            ak15hists = events_bb_sorted["HHbbVV4q"]["ak15fatJet1{}".format(vars[i])] * ak15fj1ak8fj1 + events_bb_sorted["HHbbVV4q"]["ak15fatJet2{}".format(vars[i])] * ak15fj2ak8fj1
        else:
            ak8hists = events_bb_sorted["HHbbVV4q"]["fatJet2{}".format(vars[i])] * ak8fj2match
            ak15hists = events_bb_sorted["HHbbVV4q"]["ak15fatJet1{}".format(vars[i])] * ak15fj1ak8fj2 + events_bb_sorted["HHbbVV4q"]["ak15fatJet2{}".format(vars[i])] * ak15fj2ak8fj2

        h2d = axs[j, i].hist2d(ak.to_numpy(ak8hists), ak.to_numpy(ak15hists), bins=bins[i][0], range=[[bins[i][1], bins[i][2]], [bins[i][1], bins[i][2]]])
        ax_divider = make_axes_locatable(axs[j, i])
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        fig.colorbar(h2d[3], cax=cax)
        axs[j, i].set_xlabel("FatJet{} {} (GeV)".format(j + 1, vars[i]))
        axs[j, i].set_ylabel("AK15FatJet {} {} (GeV)".format(j + 1, vars[i]))

fig.tight_layout(pad=2)
plt.savefig("figs/ak15ak82dplots.pdf")
plt.show()


plt.hist2d(np.concatenate([ak.to_numpy(events_bb_sorted["HHbbVV4q"]["fatJet1PNetXbb"] * ak8fj1match), ak.to_numpy(events_bb_sorted["HHbbVV4q"]["fatJet2PNetXbb"] * ak8fj2match)]),
            np.concatenate([ak.to_numpy(events_bb_sorted["HHbbVV4q"]["ak15fatJet1PNetMDXbb"] * ak15fj1ak8fj1 + events_bb_sorted["HHbbVV4q"]["ak15fatJet2PNetMDXbb"] * ak15fj2ak8fj1),
                            ak.to_numpy(events_bb_sorted["HHbbVV4q"]["ak15fatJet1PNetMDXbb"] * ak15fj1ak8fj2 + events_bb_sorted["HHbbVV4q"]["ak15fatJet2PNetMDXbb"] * ak15fj2ak8fj2)]),
            bins=25, range=[[0, 1], [0, 1]], norm=mpl.colors.LogNorm())
plt.colorbar()
plt.xlabel("AK8 PNet Score")
plt.ylabel("AK15 PNet Score")
plt.savefig("figs/pnet_2d.pdf")


vars = ["PNet"]
bins = [[50, 200, 600], [50, 1, 250], [50, 1, 250]]

fig, axs = plt.subplots(1, 2, figsize=(len(vars) * 10, 18))

for j in range(2):
    h2d = axs[j, i].hist2d(ak.to_numpy(events_bb_sorted["HHbbVV4q"]["fatJet{}{}".format(j + 1, vars[i])]), ak.to_numpy(events_bb_sorted["HHbbVV4q"]["ak15fatJet{}{}".format(j + 1, vars[i])]), bins=bins[i][0], range=[[bins[i][1], bins[i][2]], [bins[i][1], bins[i][2]]])
    ax_divider = make_axes_locatable(axs[j, i])
    cax = ax_divider.append_axes("right", size="7%", pad="2%")
    fig.colorbar(h2d[3], cax=cax)
    axs[j, i].set_xlabel("FatJet{} {} (GeV)".format(j + 1, vars[i]))
    axs[j, i].set_ylabel("AK15FatJet {} {} (GeV)".format(j + 1, vars[i]))

fig.tight_layout(pad=2)
plt.savefig("figs/ak15ak82dplots.pdf")
plt.show()



triggers17



events_bb_sorted_triggered = {}
fatjet_vars = ["Pt", "MassSD", "PNetHqqqq"]

for s, evts in evtDict.items():
    if s != "HH4V":
        triggered = evts[triggers17[0]]
        for i in range(1, len(triggers17)): triggered = triggered + evts[triggers17[i]]

        pnet_key = "PNetXbb_alt"  # if s == "HHbbVV4q" else "PNetXbb"
        jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

        ak_dict =   {
                        "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading),
                        "fatJet2PNetXbb": (evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading),
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)
            ak_dict["fatJet2" + var] = (evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading)

        events_bb_sorted_triggered[s] = ak.zip(ak_dict)



events_bb_sorted_triggered




# events_WW_sorted = {}
# fatjet_vars = ["Pt", "MassSD", "PNetHqqqq", 'Tau2OverTau1', 'Tau4OverTau3', 'DeepAK8_H']
#
# for s, evts in evtDict.items():
#     print(s)
#     if s != "HH4b" and s != "HHbbVV4q":
#         if s == 'HH4V':
#             pnet_key = "PNetHqqqq"
#             jet1_WW_leading = evts["fatJet1" + pnet_key][hh8q] > evts["fatJet2" + pnet_key][hh8q]
#
#             ak_dict =   {
#                             "weight": weights[s][hh8q],
#                         }
#
#             fatjet_vars.append("4q")
#
#             for var in fatjet_vars:
#                 ak_dict["fatJet1" + var] = evts["fatJet1" + var][hh8q] * jet1_WW_leading + evts["fatJet2" + var][hh8q] * ~jet1_WW_leading
#                 ak_dict["fatJet2" + var] = evts["fatJet1" + var][hh8q] * ~jet1_WW_leading + evts["fatJet2" + var][hh8q] * jet1_WW_leading
#
#             fatjet_vars = fatjet_vars[:-1]
#
#             events_WW_sorted['HH4V8q'] = ak.zip(ak_dict)
#         else:
#             pnet_key = "PNetHqqqq"
#             jet1_WW_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]
#
#             ak_dict =   {
#                             "weight": weights[s],
#                         }
#
#             for var in fatjet_vars:
#                 ak_dict["fatJet1" + var] = evts["fatJet1" + var] * jet1_WW_leading + evts["fatJet2" + var] * ~jet1_WW_leading
#                 ak_dict["fatJet2" + var] = evts["fatJet1" + var] * ~jet1_WW_leading + evts["fatJet2" + var] * jet1_WW_leading
#
#             events_WW_sorted[s] = ak.zip(ak_dict)
#
#
# events_WW_sorted

# vars = ["Pt", "MassSD"]
# varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
# bins = [[50, 250, 600], [50, 30, 200]]
# init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
# fill_hists(copy(vars), fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
# plot_hists(vars, "jet_kin_bb_leading", bins, name="bb_sorted", hh4v=False, sepsig=False)
#
# init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD Score", [100, 0, 1], fatJet=True, name="bb_sorted", bbsorted=True)
# fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
# plot_hists("PNetHqqqq", "jet_deepak8mdh4q_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)
#
# init_hists("PNetXbb", "ParticleNet Xbb Tagger Score", [100, 0, 1], fatJet=True, name="bb_sorted", bbsorted=True)
# fill_hists("PNetXbb", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
# plot_hists("PNetXbb", "jet_pnetxbb_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)



def cutflow_func(var_cuts, events=events_bb_sorted):
    cutflow = {}
    events_cuts = {}
    for s, evts in events.items():
        cuts = []
        for var, brange in var_cuts.items():
            if '+' in var:
                vars = var.split('+')
                cut1 = evts[vars[0]] > brange[0]
                cut2 = evts[vars[0]] < brange[1]
                for tvars in vars[1:]:
                    cut1 = cut1 + (evts[tvars] > brange[0])
                    cut2 = cut2 + (evts[tvars] < brange[1])

                cuts.append(cut1)
                cuts.append(cut2)
            else:
                cuts.append(evts[var] > brange[0])
                cuts.append(evts[var] < brange[1])

        cut = cuts[0]
        cutflow[s] = []
        for i in np.arange(1, len(cuts)):
            cutflow[s].append(np.sum(evts[cut].weight))
            cut = cut * cuts[i]

        cutflow[s].append(np.sum(evts[cut].weight))
        events_cuts[s] = evts[cut]

    return cutflow, events_cuts


def cftable_func(cutflow, var_cuts, cut_idx=None):
    cut_labels = []

    for var, cuts in var_cuts.items():
        varname = var.split('fat')[-1]
        if cuts[0] > 0: cut_labels.append("{} > {}".format(varname, cuts[0]))
        if cuts[1] < 9999: cut_labels.append("{} < {}".format(varname, cuts[1]))

    if cut_idx is None:
        i = 0
        cut_idx = []
        for brange in var_cuts.values():
            for j in brange:
                if j != 9999 and j != -9999:
                    cut_idx.append(i)
                i += 1

    return pandas.DataFrame(np.round(np.array(list(cutflow.values()))[:, cut_idx], 3), list(cutflow.keys()), cut_labels)



# ak8 kin ntuple cuts

var_cuts = {
    "fatJet1Pt": [250, 9999],
    "fatJet2Pt": [250, 9999],
    "fatJet1MassSD": [20, 9999],
    "fatJet2MassSD": [20, 9999],
}

hhbbVV_cutflow, events_bbs_ak8_cuts = cutflow_func(var_cuts)

cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak8_cuts, useEDWeights=True)
plot_hists(vars, "jet_kin_ak8_cuts_bb_sorted", bins, name="bb_sorted", sepsig=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak8_cuts, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_ak8_cuts_bb_sorted_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)



# ak15 kin ntuple cuts

var_cuts = {
    "ak15fatJet1Pt": [250, 9999],
    "ak15fatJet2Pt": [250, 9999],
    "ak15fatJet1MassSD": [20, 9999],
    "ak15fatJet2MassSD": [20, 9999],
}

hhbbVV_cutflow, events_bbs_ak15_cuts = cutflow_func(var_cuts, events_ak15bb_sorted)

cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
ak815_comp_plot(vars, varsl, bins, events_bbs_ak8_cuts, events_bbs_ak15_cuts, "ak815_resp_cut")


tagger2d_plots('tagger2d_ntuple_cuts', events_bbs_ak8_cuts, events_bbs_ak15_cuts, vmaxhh=0.01, vmaxqcd=250)

ak.sum(events_bbs_ak8_cuts['data']['weight'])
ak.sum(events_bbs_ak15_cuts['data']['weight'])

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_ak8_cuts, useEDWeights=True)
plot_hists(var, "jetak8_ntuple_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_bbs_ak15_cuts, useEDWeights=True)
plot_hists(var, "jetak15_ntuple_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])




# pt cuts only

kin_cuts = {
    "ak15fatJet1Pt": [250, 9999],
    "ak15fatJet2Pt": [280, 9999],
    "ak15fatJet1MassSD": [115, 145],
}

cutflow, events_bbs_pt_cuts = cutflow_func(kin_cuts)

cftable = cftable_func(cutflow, kin_cuts)
cftable


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted_ak15", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted_ak15", pevtDict=events_bbs_pt_cuts, useEDWeights=True, ak15=True)
plot_hists(var, "jet_pt_cuts_masssd_ak15", bins, name="bb_sorted_ak15", sepsig=False, ak15=True)





# kin cuts only ak15

kin_cuts = {
    "ak15fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "ak15fatJet1MassSD": [110, 140],
    "ak15fatJet2MassSD": [115, 145],
}

cutflow, events_ak15bbs_kin_cuts = cutflow_func(kin_cuts, events_ak15bb_sorted)


# var = "PNetMDXbb"
# varl = "PNetXbb Score"
# bins =  [100, 0, 1]
# init_hists(var, varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
# fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True, ak15=True)
# plot_hists(var, "jet_cut_pnet2_ak15", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False, ak15=True)
#
#
# init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
# fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
# plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_bb_leading2_ak15", [100, 0, 1], name="bb_sorted_cut", hh4v=False, sepsig=False)




# kin cuts only ak8

kin_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet2MassSD": [115, 145],
}

cutflow, events_bbs_kin_cuts = cutflow_func(kin_cuts)


tagger2d_plots('tagger2d_kin_cuts', events_bbs_kin_cuts, events_ak15bbs_kin_cuts, vmaxhh=0.005, vmaxqcd=10)



var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jetak8_kin_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_ak15bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jetak15_kin_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_ak8_cuts_masssd", bins, name="bb_sorted", sepsig=False)

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted_ak15", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted_ak15", pevtDict=events_bbs_kin_cuts, useEDWeights=True, ak15=True)
plot_hists(var, "jet_ak8_cuts_masssd_ak15", bins, name="bb_sorted_ak15", sepsig=False, ak15=True)


var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists("PNetXbb", varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_cut_pnet2", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False)

init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_bb_leading2", [100, 0, 1], name="bb_sorted_cut", hh4v=False, sepsig=False)


yields = {"sig": [], "bg": []}
for pnetcutoff in np.arange(0.99, 0.995, 0.001)[:-1]:
    sigy = []
    bgy = []
    for dakcutoff in np.arange(0.9, 1, 0.01):
        cuts = {}
        for s, evts in events_bbs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetXbb"] > pnetcutoff) * (evts["fatJet2PNetHqqqq"] > dakcutoff)

        sig = np.sum(events_bbs_kin_cuts["HHbbVV4q"][cuts["HHbbVV4q"]].weight)
        bg = np.sum(events_bbs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_bbs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)






# tagger cuts only ak15

tagger_cuts = {
    "ak15fatJet1PNetMDXbb": [0, 9999],
    "ak15fatJet2PNetHqqqq": [0.1, 9999],
}

cutflow, events_ak15bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak15_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

1.125/170000

# tagger cuts only ak8

tagger_cuts = {
    "fatJet1PNetXbb": [0, 9999],
    "fatJet2PNetHqqqq": [0.1, 9999],
}

cutflow, events_ak8bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak8_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

0.971 / 109000

# ak815_comp_plot(["Pt", "MassSD"], ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"], [[20, 250, 500], [25, 50, 200]], events_ak8bbs_tagger_cuts, events_ak15bbs_tagger_cuts, "ak815_tagger_cuts", qcd=True)
# ak815_comp_plot(["Pt", "MassSD"], ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"], [[20, 250, 500], [25, 50, 200]], events_ak8bbs_tagger_cuts, events_ak15bbs_tagger_cuts, "ak815_tagger_cuts_no_qcd", qcd=False)

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_tagger_cut", pevtDict=events_ak8bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak8_loose_tagger_cuts_{}_{}_masssd_rat".format(tagger_cuts["fatJet1PNetXbb"][0], tagger_cuts["fatJet2PNetHqqqq"][0]), bins, name="bb_sorted_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=True, rat_ylim=2.5, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak15_loose_tagger_cuts_{}_{}_masssd_rat".format(tagger_cuts["fatJet1PNetXbb"][0], tagger_cuts["fatJet2PNetHqqqq"][0]), bins, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=True, rat_ylim=2.5, blinding=[115, 145])



# tagger cuts only ak8 + ak15

tagger_cuts = {
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}

cutflow, events_ak815bbs_tagger_cuts = cutflow_func(tagger_cuts, events_8bb_15VV_sorted)
cftable = cftable_func(cutflow, tagger_cuts)
cftable


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[20, 250, 500], [25, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted815", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted815", pevtDict=events_ak815bbs_tagger_cuts, useEDWeights=True, jet2ak15=True)
plot_hists(vars, "jet_kin_tagger_cuts_ak815", bins, name="bb_sorted815", sepsig=False, ak15=True)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_tagger_cut", pevtDict=events_ak8bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak8_tagger_cuts_masssd", bins, name="bb_sorted_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak15_tagger_cuts_masssd", bins, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])

# kin + tagger cuts

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2PNetHqqqq": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts)
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [10, 50, 150]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_bbs_cuts, useEDWeights=True)
plot_hists(var, "jet_cuts_masssd", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, lumilabel=True, stackall=True)


init_hists("MassSD", "Soft Drop Mass (GeV)", [50, 1, 400], name="data", scale=False)
fill_hists("MassSD", scale=False, data=True, name='data')
plot_hists("MassSD", "data_masssd_rat_pre_tf", [50, 1, 400], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=0.5)







# AK15

# ak15 all cuts

var_cuts = {
    "ak15fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "ak15fatJet1MassSD": [110, 140],
    "ak15fatJet2MassSD": [115, 145],
    "ak15fatJet1PNetMDXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_ak15bb_sorted)[0]

cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

cftable.to_csv('hhbbVV_cutflow_ak15.csv')


# bg estimatation using mass2 sidebands

var_cuts = {
    "ak15fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "ak15fatJet1MassSD": [110, 140],
    # "ak15fatJet2MassSD": [115, 145],
    "ak15fatJet1PNetMDXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}

cutflow, events_ak15bbs_nomass2_cuts = cutflow_func(var_cuts, events_ak15bb_sorted)

cftable = cftable_func(cutflow, var_cuts)
cftable


# mass sidebands

mass_sb1 = {
    "ak15fatJet2MassSD": [100, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_ak15bbs_nomass2_cuts)[0]
msb1_cf


mass_sb2 = {
    "ak15fatJet2MassSD": [145, 160],
}

msb2_cf = cutflow_func(mass_sb2, events_ak15bbs_nomass2_cuts)[0]
msb2_cf






# AK8

# ak8 all cuts

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2PNetHqqqq": [0.96, 9999],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_bbs_ak8_cuts)[0]

cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

cftable.to_csv('hhbbVV_cutflow_triggered.csv')

0.029/11

0.023/7

# bg estimatation using mass2 sidebands

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    # "fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2PNetHqqqq": [0.96, 9999],
}

cutflow, events_bbs_nomass2_cuts = cutflow_func(var_cuts, events_bbs_ak8_cuts)
cftable = cftable_func(cutflow, var_cuts)
cftable

# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [100, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_nomass2_cuts)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [145, 160],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_nomass2_cuts)[0]
msb2_cf



ak815_comp_plot(["Pt", "MassSD"], ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"], [[20, 250, 500], [25, 50, 200]], events_bbs_cuts, events_ak15bbs_cuts, "ak815_nomassd_cut", qcd=False)

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_nomass2_cuts, useEDWeights=True)
plot_hists(var, "jetak8_all-mass2_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_ak15bbs_nomass2_cuts, useEDWeights=True)
plot_hists(var, "jetak15_all-mass2_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])



#
# # AK15
#
# # ak15 all cuts
#
# var_cuts = {
#     "fatJet1Pt": [275, 9999],
#     "ak15fatJet2Pt": [300, 9999],
#     "fatJet1MassSD": [110, 140],
#     "ak15fatJet2MassSD": [115, 145],
#     "fatJet1PNetXbb": [0.99, 9999],
#     "ak15fatJet2PNetHqqqq": [0.935, 9999],
# }
#
# hhbbVV_cutflow = cutflow_func(var_cuts, events_8bb_15VV_sorted)[0]
#
# cftable = cftable_func(hhbbVV_cutflow, var_cuts)
# cftable
#
# cftable.to_csv('hhbbVV_cutflow_ak815.csv')
#
#
# # bg estimatation using mass2 sidebands
#
# var_cuts = {
#     "fatJet1Pt": [275, 9999],
#     "ak15fatJet2Pt": [300, 9999],
#     "fatJet1MassSD": [110, 140],
#     # "ak15fatJet2MassSD": [115, 145],
#     "fatJet1PNetXbb": [0.99, 9999],
#     "ak15fatJet2PNetHqqqq": [0.935, 9999],
# }
#
# cutflow, events_ak815bbs_cuts = cutflow_func(var_cuts, events_8bb_15VV_sorted)
#
# cftable = cftable_func(cutflow, var_cuts)
# cftable
#
#
# # mass sidebands
#
# mass_sb1 = {
#     "ak15fatJet2MassSD": [100, 115],
# }
#
# msb1_cf = cutflow_func(mass_sb1, events_ak815bbs_cuts)[0]
# msb1_cf
#
#
# mass_sb2 = {
#     "ak15fatJet2MassSD": [145, 160],
# }
#
# msb2_cf = cutflow_func(mass_sb2, events_ak815bbs_cuts)[0]
# msb2_cf
#
#
#
#
#





# AK8 for bb, AK15 for WW


# + mass2

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [115, 145],
    "ak15fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_8bb_15VV_sorted)[0]

# del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

cftable.to_csv('hhbbVV_ak8bb15VV_cutflow.csv')


# kin + bbVV tagger cuts - masssd

var_cuts = {
    "fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [115, 145],
    # "ak15fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}


hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts, events_8bb_15VV_sorted)
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

# mass sidebands

mass_sb1 = {
    "ak15fatJet2MassSD": [100, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts)[0]
msb1_cf


mass_sb2 = {
    "ak15fatJet2MassSD": [145, 160],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_cuts)[0]
msb2_cf












# 4b kin + bbVV tagger cuts + triggers

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2PNetHqqqq": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts_triggered = cutflow_func(var_cuts, events_bb_sorted_triggered)

cut_labels = [
                'Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
             ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable


# 4b kin + bbVV tagger cuts + mass2 + triggers

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2PNetHqqqq": [0.94, 9999],
    "fatJet2MassSD": [115, 135],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_bb_sorted_triggered)[0]

cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
                'Jet2 MassSD > 115',
                'Jet2 MassSD < 135',
            ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable

cftable.to_csv('hhbbVV_cutflow_triggered.csv')



# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [105, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts_triggered)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 145],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_cuts_triggered)[0]
msb2_cf





# 4b cuts - jet2 MassSD

hh4b_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    # "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
}

hh4bcf, hh4b_cuts_evts = cutflow_func(hh4b_var_cuts)
hh4bcf


cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.985',
                'Jet2 PNetXbb > 0.985',
            ]

del(hh4bcf['HHbbVV4q'])

cftable = cftable_func(hh4bcf, hh4b_var_cuts, cut_labels)
cftable


cftable.to_csv('hh4b_cutflow_2.csv')


# 4b cut

hh4b_mass2_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
    "fatJet2MassSD": [105, 135],
}

hh4bm2cf, hh4b_mass2_cuts_evts = cutflow_func(hh4b_mass2_cuts)


cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.985',
                'Jet2 PNetXbb > 0.985',
                'Jet2 MassSD > 105',
                'Jet2 MassSD < 135',
            ]

del(hh4bm2cf['HHbbVV4q'])

cftable = cftable_func(hh4bm2cf, hh4b_mass2_cuts, cut_labels)
cftable

cftable.to_csv('hh4b_cutflow_40.csv')


# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [90, 105],
}

hh4b_msb1_cf = cutflow_func(mass_sb1, hh4b_cuts_evts)[0]
hh4b_msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 150],
}

hh4b_msb2_cf = cutflow_func(mass_sb2, hh4b_cuts_evts)[0]
hh4b_msb2_cf




evtDict["HH4V"].fields

init_hists("PNetHqqqq", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("DeepAK8_H", "DeepAK8MD H", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("DeepAK8_H", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("DeepAK8_H", "jet_deepak8h_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau2OverTau1", "tau2/tau1", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau2OverTau1", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau2OverTau1", "jet_tau2otau1_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau4OverTau3", "tau4/tau3", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau4OverTau3", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau4OverTau3", "jet_tau4otau3_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("PNetHqqqq", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)



init_hists("Pt", "$p_T$", [100, 250, 1200], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Pt", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Pt", "jet_pt_WW_leading", [100, 250, 1200], name="WW_sorted", hh4v=True, sepsig=False)


# WW kin cuts only

kin_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [75, 140],
    "fatJet2MassSD": [75, 140],
}

cutflow, events_WWs_kin_cuts = cutflow_func(kin_cuts, events_WW_sorted)

cutflow

init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="WW_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted_cut", pevtDict=events_WWs_kin_cuts, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted_cut", hh4v=True, sepsig=False)

sig

hists['WW_sorted_cutPNetHqqqq']['HH4V']


yields = {"sig": [], "bg": []}
for pnetcutoff1 in np.arange(0.8, 0.85, 0.01):
    sigy = []
    bgy = []
    for pnetcutoff2 in np.arange(0.8, 0.85, 0.01):
        cuts = {}
        for s, evts in events_WWs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetHqqqq"] > pnetcutoff1) * (evts["fatJet2PNetHqqqq"] > pnetcutoff2)

        sig = np.sum(events_WWs_kin_cuts["HH4V"][cuts["HH4V"]].weight)
        bg = np.sum(events_WWs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_WWs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)


yields
















# Transfer factor study

# Cutting above the data trigger turn-ons
var_cuts = {
    "fatJet1Pt": [500, 9999],
    "fatJet2Pt": [500, 9999],
    # "fatJet1MassSD": [200, 9999],
    # "fatJet2MassSD": [200, 9999],
}

cutflow_control, events_tf_cuts = cutflow_func(var_cuts)

cftable = cftable_func(cutflow_control, var_cuts)
cftable

cutflow_control['data'][-1] / (cutflow_control['tt'][-1] + cutflow_control['QCD'][-1])
(cutflow_control['data'][-1] - cutflow_control['tt'][-1]) / cutflow_control['QCD'][-1]

ak.sum(weights['data'])

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 0, 400]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_tf_cuts, useEDWeights=True)
plot_hists(var, "jet_cuts_masssd_rat_post_pt_cut", bins, data=True, name="bb_sorted_cut", stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)

var = "Pt"
varl = "$p_T$ (GeV)"
bins =  [38, 500, 2000]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_tf_cuts, useEDWeights=True)
plot_hists(var, "jet_cuts_pt_rat_post_tf", bins, data=True, name="bb_sorted_cut", stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)












# Tagger Plots/ROC Curves

init_hists("PNetXbb_alt", "Particle Net Xbb", disc_bin, fatJet=True)
fill_hists(vars, fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetXbb_alt"][fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetXbb_alt"][fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetXbb_alt"][fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetXbb_alt"][fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "PNetXbb_alt"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['bg'].append(fatJet1WXbb)
disc_vals['jet2']['bg'].append(fatJet2WXbb)
disc_vals['jet1']['sig'].append(fatJet1bXbb)
disc_vals['jet2']['sig'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('ParticleNetXbb ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_bb_sep_pxbb_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('ParticleNetXbb Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_pxbb_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()









# Tagger Plots/ROC Curves AK15

disc_bin = [100, 0, 1]

init_hists("PNetMDXbb", "Particle Net Xbb", disc_bin, fatJet=True)
fill_hists("PNetMDXbb", fatJet=True, scale=False, hh4v=False, ak15=True)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet1PNetMDXbb"][ak15fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet1PNetMDXbb"][ak15fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet2PNetMDXbb"][ak15fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet2PNetMDXbb"][ak15fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "PNetMDXbb"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['bg'].append(fatJet1WXbb)
disc_vals['jet2']['bg'].append(fatJet2WXbb)
disc_vals['jet1']['sig'].append(fatJet1bXbb)
disc_vals['jet2']['sig'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('ParticleNetMDXbb ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/ak15jet_bb_sep_pxbb_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('ParticleNetXbb Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/ak15jet_pxbb_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()

evtDict["HHbbVV4q"].fields



#PNetH4q ak15

disc_bin = [10000, 0, 1]

init_hists("PNetHqqqq", "ParticleNet H4q", disc_bin, fatJet=True)
fill_hists("PNetHqqqq", fatJet=True, scale=False, hh4v=False, ak15=True)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet1PNetHqqqq"][ak15fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet1PNetHqqqq"][ak15fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet2PNetHqqqq"][ak15fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["ak15fatJet2PNetHqqqq"][ak15fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][ak15fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "PNetHqqqq"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['sig'].append(fatJet1WXbb)
disc_vals['jet2']['sig'].append(fatJet2WXbb)
disc_vals['jet1']['bg'].append(fatJet1bXbb)
disc_vals['jet2']['bg'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetH4q ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/ak15jet_bb_sep_pneth4q_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetH4q Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/ak15jet_pneth4q_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()


# PNetH4q ak8

disc_bin = [100, 0, 1]

init_hists("PNetHqqqq", "PNetHqqqq", disc_bin, fatJet=True)
fill_hists("PNetHqqqq", fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetHqqqq"][fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetHqqqq"][fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetHqqqq"][fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetHqqqq"][fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "PNetHqqqq"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['sig'].append(fatJet1WXbb)
disc_vals['jet2']['sig'].append(fatJet2WXbb)
disc_vals['jet1']['bg'].append(fatJet1bXbb)
disc_vals['jet2']['bg'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetH4q ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_bb_sep_pneth4q_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetH4q Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_pneth4q_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()





# PNetHqqqq

disc_bin = [10000, 0, 1]

init_hists("PNetHqqqq", "PNetHqqqq", disc_bin, fatJet=True)
fill_hists("PNetHqqqq", fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetHqqqq"][fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1PNetHqqqq"][fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetHqqqq"][fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2PNetHqqqq"][fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "PNetHqqqq"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['sig'].append(fatJet1WXbb)
disc_vals['jet2']['sig'].append(fatJet2WXbb)
disc_vals['jet1']['bg'].append(fatJet1bXbb)
disc_vals['jet2']['bg'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetHqqqq ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_bb_sep_deepak8h4q_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('PNetHqqqq Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_deepak8h4q_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()






evtDict["HHbbVV4q"].fields

# Deep AK8H
disc_bin = [100, 0, 1]

init_hists("DeepAK8_H", "Deep AK8H", disc_bin, fatJet=True)
fill_hists("DeepAK8_H", fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1DeepAK8_H"][fatJet1W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet1DeepAK8_H"][fatJet1b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2DeepAK8_H"][fatJet2W]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbVV4q"]["fatJet2DeepAK8_H"][fatJet2b]), weights=np.array(evtDict["HHbbVV4q"]["weight"][fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

var = "DeepAK8_H"
disc_vals = {}
for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk] = {'sig': [], 'bg': []}
    for s in evtDict.keys():
        if 'HH' not in s:
            disc_vals[jk]['bg'].append(hists['fatJet' + var].project("sample", jk).values()[(s, )])

disc_vals['jet1']['bg'].append(fatJet1WXbb)
disc_vals['jet2']['bg'].append(fatJet2WXbb)
disc_vals['jet1']['sig'].append(fatJet1bXbb)
disc_vals['jet2']['sig'].append(fatJet2bXbb)

for i in range(2):
    jk = 'jet' + str(i + 1)
    disc_vals[jk]['sig'] = np.sum(np.array(disc_vals[jk]['sig']), axis=0)
    disc_vals[jk]['bg'] = np.sum(np.array(disc_vals[jk]['bg']), axis=0)
    disc_vals[jk]['tpr'] = np.cumsum(disc_vals[jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['sig']))
    disc_vals[jk]['fpr'] = np.cumsum(disc_vals[jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[jk]['bg']))
    disc_vals[jk]['tpr'] = np.append(disc_vals[jk]['tpr'], 0)
    disc_vals[jk]['fpr'] = np.append(disc_vals[jk]['fpr'], 0)


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('DeepAK8H ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_bb_sep_deepak8h_roc.pdf", bbox_inches='tight')
plt.show()


for j in range(2):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[jk]['tpr'][:-1] + disc_vals[jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[jk]['fpr'][:-1] - disc_vals[jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[jk]['fpr'], disc_vals[jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title('DeepAK8H Semilog ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_deepak8h_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()











# BSM check comparing bbVV with 4b

hists['genHiggsPt'] = hist.Hist("Events (Unit Normalized)",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Leading Higgs $p_T$ ", 51, 0, 1000),
                             )

hists['genHiggsPt'].fill(sample='HH4b', h1=evtDict['HH4b']['genHiggs1Pt'], weight=weights['HH4b'])
hists['genHiggsPt'].fill(sample='HHbbVV4q', h1=evtDict['HHbbVV4q']['genHiggs1Pt'], weight=weights['HHbbVV4q'])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])

hists['genHiggsPt'].scale(scale_factor, axis='sample')

hist.plot1d(hists['genHiggsPt'], line_opts=line_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.savefig('figs/hhbbVV_vs_hh4b_norm_new.pdf', bbox_inches='tight')
plt.show()


fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['genHiggsPt'], line_opts=line_opts, ax=axs, density=True)
# axs.set_yscale('log')
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
hist.plotratio(hists['genHiggsPt']['HHbbVV4q'].sum('sample'), hists['genHiggsPt']['HH4b'].sum('sample'), ax=rax, error_opts=data_err_opts, unc='num')
rax.set_ylabel("bbVV/4b")
rax.set_ylim(0, 2)
plt.savefig('figs/hhbbVV_vs_hh4b_norm.pdf', bbox_inches='tight')
plt.show()




hists['genHiggsMass'] = hist.Hist("Events (Unit Normalized)",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"HH Mass", 51, 200, 2500),
                             )


hists['genHiggsMass'].fill(sample='HH4b', h1=evtDict['HH4b']['genHH_mass'], weight=weights['HH4b'])
hists['genHiggsMass'].fill(sample='HHbbVV4q', h1=evtDict['HHbbVV4q']['genHH_mass'], weight=weights['HHbbVV4q'])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])
hists['genHiggsMass'].scale(scale_factor, axis='sample')

hist.plot1d(hists['genHiggsMass'], line_opts=line_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.savefig('figs/hhbbVV_vs_hh4b_mass_norm.pdf', bbox_inches='tight')
plt.show()


fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['genHiggsMass'], line_opts=line_opts, ax=axs, density=True)
# axs.set_yscale('log')
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
hist.plotratio(hists['genHiggsMass']['HHbbVV4q'].sum('sample'), hists['genHiggsMass']['HH4b'].sum('sample'), ax=rax, error_opts=data_err_opts, unc='num')
rax.set_ylabel("bbVV/4b")
rax.set_ylim(0, 5)
plt.savefig('figs/hhbbVV_vs_hh4b_mass_norm_new.pdf', bbox_inches='tight')
plt.show()













# Invariant Mass Calc


masses = {}
for s, evts in evtDict.items():
    masses[s] = (ak.zip({"pt": evtDict[s]['fatJet1Pt'], "eta": evtDict[s]['fatJet1Eta'], "phi": evtDict[s]['fatJet1Phi'], "mass": evtDict[s]['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
                + ak.zip({"pt": evtDict[s]['fatJet2Pt'], "eta": evtDict[s]['fatJet2Eta'], "phi": evtDict[s]['fatJet2Phi'], "mass": evtDict[s]['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")).mass


hists["jets_mass"] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet12", r"Leading 2 Jets' Invariant Mass (GeV)", 100, 400, 4000),
                                )

for s in evtDict.keys():
    if s == 'HH4V':
        hists["jets_mass"].fill(sample=s,
                             jet12 = masses[s],
                             )
    else:
        hists["jets_mass"].fill(sample=s,
                             jet12 = masses[s],
                             weight = evtDict[s]["weight"]
                             )

hists['jets_mass'].scale(scale_factor, axis='sample')

plt.gca().set_prop_cycle(cycler(color=colors_cycle))
ylim = np.max(list(hists["jets_mass"].project("sample", "jet12").values().values())) * 1.1
hist.plot1d(hists["jets_mass"][:2].project("sample", "jet12"), overlay='sample', clear=False, line_opts=line_opts)
plt.gca().set_prop_cycle(cycler(color=colors_cycle[2:]))
hist.plot1d(hists["jets_mass"][2:].project("sample", "jet12"), overlay='sample', stack=True, clear=False, fill_opts=fill_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.ylim(0, ylim)
plt.savefig("figs/jets_inv_mass.pdf", bbox_inches='tight')
