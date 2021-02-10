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
from copy import copy
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)


plt.style.use(hep.style.ROOT)


evts = uproot4.open("./data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root")["tree"]
evts = uproot4.open("./data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root")["tree"]
evts.keys()


evts.arrays()
evts['genLeptonId']


samples = {
    "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWWqq": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")

# sample_names = {
#     "HH4V": "HHToVVVV",
#     "HHbbWWqq": "HHToBBVV",
#     "QCD": "QCD",
#     "TT": "TT",
# }
#
# fnames = listdir('data/weighted')
# fnames
#
# evtDict = {}
# for name in fnames:
#     fname = 'data/weighted/' + name
#     print(fname)
#     for sample, fn in sample_names.items():
#         if fn in name:
#             if sample in evtDict:
#                 evtDict[sample] = ak.concatenate((evtDict[sample], uproot4.open(fname)["tree"].arrays()))
#             else:
#                 evtDict[sample] = uproot4.open(fname)["tree"].arrays()




# evtDict['HHbbWWqq']['genHiggs1W2Decay']
# evtDict['HHbbWWqq']['genHiggs1W2Decay']
# evtDict['HHbbWWqq']['genHiggs1W1Decay']
# evtDict['HHbbWWqq']['genHiggs2W1Decay']
# evtDict['HHbbWWqq']['genHiggs2W2Decay']
#
# np.unique(np.array((evtDict['HHbbWWqq']['genHiggs1W2Decay'] + evtDict['HHbbWWqq']['genHiggs2W1Decay'])))
#
# evtDict['HHbbWWqq']['genHiggs1Pt']
# evtDict['HHbbWWqq']['genHiggs1W1Pt']
# evtDict['HHbbWWqq']['fatJet1Pt']
# evtDict['HHbbWWqq']['fatJet2Pt']
# evtDict['HHbbWWqq']['genHiggs2Pt']
#
#
# evtDict['HHbbWWqq']['genHiggs1Eta'][:5]
# evtDict['HHbbWWqq']['fatJet1Eta'][:5]
# evtDict['HHbbWWqq']['fatJet2Eta'][:5]
# evtDict['HHbbWWqq']['genHiggs2Eta'][:5]


genHiggs1 = ak.zip({"pt": evtDict['HHbbWWqq']['genHiggs1Pt'], "eta": evtDict['HHbbWWqq']['genHiggs1Eta'], "phi": evtDict['HHbbWWqq']['genHiggs1Phi'], "mass": ak.full_like(evtDict['HHbbWWqq']['genHiggs1Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")
genHiggs2 = ak.zip({"pt": evtDict['HHbbWWqq']['genHiggs2Pt'], "eta": evtDict['HHbbWWqq']['genHiggs2Eta'], "phi": evtDict['HHbbWWqq']['genHiggs2Phi'], "mass": ak.full_like(evtDict['HHbbWWqq']['genHiggs2Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")

# genHiggs12 = ak.zip({"pt": [evtDict['HHbbWWqq']['genHiggs1Pt'], evtDict['HHbbWWqq']['genHiggs2Pt']],
#                     "eta": [evtDict['HHbbWWqq']['genHiggs1Eta'], evtDict['HHbbWWqq']['genHiggs2Eta']],
#                     "phi": [evtDict['HHbbWWqq']['genHiggs1Phi'], evtDict['HHbbWWqq']['genHiggs2Phi']],
#                     "mass": ak.full_like([evtDict['HHbbWWqq']['genHiggs1Pt'], evtDict['HHbbWWqq']['genHiggs2Pt']], 125.1)})#, with_name="PtEtaPhiMLorentzVector")

fatJet1 = ak.zip({"pt": evtDict['HHbbWWqq']['fatJet1Pt'], "eta": evtDict['HHbbWWqq']['fatJet1Eta'], "phi": evtDict['HHbbWWqq']['fatJet1Phi'], "mass": evtDict['HHbbWWqq']['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet2 = ak.zip({"pt": evtDict['HHbbWWqq']['fatJet2Pt'], "eta": evtDict['HHbbWWqq']['fatJet2Eta'], "phi": evtDict['HHbbWWqq']['fatJet2Phi'], "mass": evtDict['HHbbWWqq']['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet3 = ak.zip({"pt": evtDict['HHbbWWqq']['fatJet3Pt'], "eta": evtDict['HHbbWWqq']['fatJet3Eta'], "phi": evtDict['HHbbWWqq']['fatJet3Phi'], "mass": evtDict['HHbbWWqq']['fatJet3Mass']}, with_name="PtEtaPhiMLorentzVector")


# print(genHiggs1.delta_r(fatJet1)[:10])
# print(genHiggs1.delta_r(fatJet2)[:10])
# print(genHiggs1.delta_r(fatJet3)[:10])
#
# print(genHiggs2.delta_r(fatJet1)[:10])
# print(genHiggs2.delta_r(fatJet2)[:10])
# print(genHiggs2.delta_r(fatJet3)[:10])
#
# print(fatJet1.delta_r(genHiggs1)[:10])
# print(fatJet1.delta_r(genHiggs2)[:10])

H1ltH2 = fatJet1.delta_r(genHiggs1) < fatJet1.delta_r(genHiggs2)
fatJet1H1 = fatJet1.delta_r(genHiggs1) < 0.4
fatJet1H2 = fatJet1.delta_r(genHiggs2) < 0.4
fatJet1H12 = fatJet1H1 * fatJet1H2

fatJet1H1 = fatJet1H1 * (~fatJet1H12 + fatJet1H12 * H1ltH2)
fatJet1H2 = fatJet1H1 * (~fatJet1H12 + fatJet1H12 * ~H1ltH2)


H1ltH2 = fatJet2.delta_r(genHiggs1) < fatJet2.delta_r(genHiggs2)
fatJet2H1 = fatJet2.delta_r(genHiggs1) < 0.4
fatJet2H2 = fatJet2.delta_r(genHiggs2) < 0.4
fatJet2H12 = fatJet2H1 * fatJet2H2

fatJet2H1 = fatJet2H1 * (~fatJet2H12 + fatJet2H12 * H1ltH2)
fatJet2H2 = fatJet2H1 * (~fatJet2H12 + fatJet2H12 * ~H1ltH2)

H1W = evtDict['HHbbWWqq']['genHiggs1W1Decay'] == 0
H2W = evtDict['HHbbWWqq']['genHiggs2W1Decay'] == 0

fatJet1W = fatJet1H1 * H1W + fatJet1H2 * H2W
fatJet2W = fatJet2H1 * H1W + fatJet2H2 * H2W



# np.sum(np.array((fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)))
np.sum(np.array((fatJet1.delta_r(genHiggs1) < 0.4) * (fatJet1.delta_r(genHiggs2) < 0.4)))
np.sum(np.array((fatJet2.delta_r(genHiggs1) < 0.4) * (fatJet2.delta_r(genHiggs2) < 0.4)))
np.sum(np.array((fatJet3.delta_r(genHiggs1) < 0.4) + (fatJet3.delta_r(genHiggs2) < 0.4)))

# sum(genHiggs1.pt < 300)
#
# fatjet1_nohiggs = np.array((fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4))
# fatjet2_nohiggs = np.array((fatJet2.delta_r(genHiggs1) > 0.4) * (fatJet2.delta_r(genHiggs2) > 0.4))
#
# fatjet12_nohiggs = fatjet1_nohiggs * fatjet2_nohiggs
#
# np.sum(np.array(fatjet12_nohiggs))
#
# plt.hist(evtDict['HHbbWWqq']['genHiggs1W1Pt'], histtype='step', bins=np.linspace(0, 1200, 51))
# plt.xlabel("Gen Higgs W1 $p_T$ (GeV)")
# plt.ylabel("Events")
# # plt.title("Events where fat jets 1, 2 match with neither gen higgs")
# plt.show()
#
#
# plt.hist(genHiggs1[fatjet12_nohiggs].pt, histtype='step', bins=np.linspace(0, 1200, 51))
# plt.xlabel("Gen Higgs 1 $p_T$ (GeV)")
# plt.ylabel("Events")
# plt.title("Events where fat jets 1, 2 match with neither gen higgs")
# plt.show()
#
# np.sum(np.array((fatJet2.delta_r(genHiggs1) > 0.4) * (fatJet2.delta_r(genHiggs2) > 0.4)))
#
# np.sum(np.array((fatJet3.delta_r(genHiggs1) < 90)))
#
# print(fatJet2.delta_r(genHiggs1)[(fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)][:10])
# print(fatJet2.delta_r(genHiggs2)[(fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)][:10])
#
# print(fatJet3.delta_r(genHiggs1)[(fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)][:10])
# print(fatJet3.delta_r(genHiggs2)[(fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)][:10])
#
#
# print(fatJet2.delta_r(genHiggs1)[:10])
# print(fatJet2.delta_r(genHiggs2)[:10])
#
# print(fatJet3.delta_r(genHiggs1)[:10])
# print(fatJet3.delta_r(genHiggs2)[:10])


# plt.hist(evtDict["HHbbWWqq"]['fatJet1PNetXbb'], bins=np.linspace(-1, 1, 100))
# plt.title("HHbbWWqq")
# plt.xlabel("Fat Jet 1 PNet Score")
#
# plt.hist(evtDict["HHbbWWqq"]['fatJet1DeepAK8H'], bins=np.linspace(-1, 1, 100))
# plt.title("HHbbWWqq")
# plt.xlabel("Fat Jet 1 Deep AK8H Score")



hists = {}


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


# weights = {
#     "HH4V": None,
#     "HHbbWWqq": evtDict["HHbbWWqq"]['weight'],
#     "QCD": evtDict["QCD"]['weight'],
#     "tt": evtDict["tt"]['weight'],
# }


# scale_factor = {
#         'HH4V': 1 / len(evtDict["HH4V"]["weight"]),
#         'HHbbWWqq': 1 / np.sum(evtDict["HHbbWWqq"]["weight"]),
#         'HHbbWWqq - Hbb': 2 / np.sum(evtDict["HHbbWWqq"]["weight"]),
#         'HHbbWWqq - HWW': 2 / np.sum(evtDict["HHbbWWqq"]["weight"]),
#         'QCD': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
#         'tt': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
#         }


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
                    kwargs = {}
                    for j in range(nbins):
                        kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[i][s])]
                    hists[name + vars[i]["HHbbWWqq"]].fill(sample=s, weight=weights[s], **kwargs)

            if scale: hists[name + vars[i]["HHbbWWqq"]].scale(scale_factor, axis='sample')
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        for s, evts in evtDict.items():


            if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                kwargs = {}
                for j in range(nbins):
                    kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[s])]
                hists[name + vars["HHbbWWqq"]].fill(sample=s, weight=weights[s], **kwargs)

        if scale: hists[name + vars["HHbbWWqq"]].scale(scale_factor, axis='sample')


colors_cycle = ['#e31a1c', '#33a02c', '#a6cee3', '#1f78b4', '#b2df8a', '#fb9a99', '#f1005f', '#d09d09']

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


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True):
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

            var = vars[i]["HHbbWWqq"]
            if type(binrs[0]) is not list:
                temp = binrs
                binrs = []
                for k in range(len(vars)):
                    binrs.append(temp)

            for j in range(nbins):
                if fatJet:
                    ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                    axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                    if hh4v: hist.plot1d(hists[name + var][:1].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)

                    _ = axs[i, j].hist(evtDict["HHbbWWqq"]["fatJet{}{}".format(j + 1, var)] * ~fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbWWqq - Hbb', weights=evtDict["HHbbWWqq"]["weight"] * scale_factor["HHbbWWqq"], histtype='step', **line_opts)
                    _ = axs[i, j].hist(evtDict["HHbbWWqq"]["fatJet{}{}".format(j + 1, var)] * fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbWWqq - HWW', weights=evtDict["HHbbWWqq"]["weight"] * scale_factor["HHbbWWqq"], histtype='step', **line_opts)

                    axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                    ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
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

        var = vars["HHbbWWqq"]
        fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))
        for j in range(nbins):
            if fatJet:
                ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                axs[j].set_prop_cycle(cycler(color=colors_cycle))
                if hh4v: hist.plot1d(hists[name + var][:1].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)

                _ = axs[j].hist(evtDict["HHbbWWqq"]["fatJet{}{}".format(j + 1, var)] * ~fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbWWqq - Hbb', weights=evtDict["HHbbWWqq"]["weight"] * scale_factor["HHbbWWqq"], histtype='step', **line_opts)
                _ = axs[j].hist(evtDict["HHbbWWqq"]["fatJet{}{}".format(j + 1, var)] * fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[7], label='HHbbWWqq - HWW', weights=evtDict["HHbbWWqq"]["weight"] * scale_factor["HHbbWWqq"], histtype='step', **line_opts)

                axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
                axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                axs[j].set_ylim(0, ylim)

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()


init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=False)

tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])

scale_factor = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbWWqq': 1 / tot_events['HHbbWWqq'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}

init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=True)
plot_hists("Pt", "kin_tests", [40, 200, 2000])

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 2000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True)
plot_hists(vars, "jet_kin", bins)


vars = {
    "HH4V": "PNetXbb_alt",
    "HHbbWWqq": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
}
disc_bin = [100, 0, 1]
init_hists("PNetXbb_alt", "Particle Net Xbb", [100, 0, 1], fatJet=True)
fill_hists(vars, fatJet=True, scale=True, hh4v=False)
plot_hists(vars, "pnetxbb", [100, 0, 1], hh4v=False)



hists['fatJet' + var].project("sample", jk).values()
hists['fatJet' + var].project("sample", jk).values()


evtDict["QCD"]["weight"]
evtDict["tt"]["weight"]


evtDict["HHbbWWqq"]["weight"] * len(evtDict["HHbbWWqq"]["weight"])

np.unique(np.array(evtDict["HHbbWWqq"]["weight"]))
plt.hist(evtDict["HHbbWWqq"]["weight"])

hists


init_hists("PNetXbb_alt", "Particle Net Xbb", disc_bin, fatJet=True)
fill_hists(vars, fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbWWqq"]["fatJet1PNetXbb_alt"][fatJet1W]), weights=np.array(evtDict["HHbbWWqq"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbWWqq"]["fatJet1PNetXbb_alt"][~fatJet1W]), weights=np.array(evtDict["HHbbWWqq"]["weight"][~fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbWWqq"]["fatJet2PNetXbb_alt"][fatJet2W]), weights=np.array(evtDict["HHbbWWqq"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbWWqq"]["fatJet2PNetXbb_alt"][~fatJet2W]), weights=np.array(evtDict["HHbbWWqq"]["weight"][~fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

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

plt.title(var + ' ROC Curve')
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

plt.title(var + ' ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_pxbb_sep_roc_semilog.pdf", bbox_inches='tight')
plt.show()



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






var = "PNetXbb_alt"

xbb = ak.zip({"fatJet1PNetXbb": evtDict["HHbbWWqq"]["fatJet1PNetXbb_alt"],
              "fatJet2PNetXbb": evtDict["HHbbWWqq"]["fatJet2PNetXbb_alt"],
              "weight": evtDict["HHbbWWqq"]["weight"],
              })



pnetxbbs = {}
for s, evts in evtDict.items():
    if s != "HH4V":
        key = "PNetXbb_alt" if "HH" in s else "PNetXbb"
        pnetxbbs[s] = ak.zip({
                                "fatJet1PNetXbb": evts["fatJet1" + key],
                                "fatJet2PNetXbb": evts["fatJet2" + key],
                                "weight": evts["weight"],
                             })

pnetxbb_cutoff = 0.99

yields = {}
yields["jet1xbb"] = []
yields["jet2xbb"] = []
yields["jet12xbb"] = []

for cutoff in np.arange(0.9, 1, 0.005):
    sig = np.sum(pnetxbbs["HHbbWWqq"][pnetxbbs["HHbbWWqq"]["fatJet1PNetXbb"] > cutoff].weight)
    bg = np.sum(pnetxbbs["QCD"][pnetxbbs["QCD"]["fatJet1PNetXbb"] > cutoff].weight) + np.sum(pnetxbbs["tt"][pnetxbbs["tt"]["fatJet1PNetXbb"] > cutoff].weight)
    yields["jet1xbb"].append(sig / np.sqrt(sig ** 2 + bg ** 2))

    sig = np.sum(pnetxbbs["HHbbWWqq"][pnetxbbs["HHbbWWqq"]["fatJet2PNetXbb"] > cutoff].weight)
    bg = np.sum(pnetxbbs["QCD"][pnetxbbs["QCD"]["fatJet2PNetXbb"] > cutoff].weight) + np.sum(pnetxbbs["tt"][pnetxbbs["tt"]["fatJet2PNetXbb"] > cutoff].weight)
    yields["jet2xbb"].append(sig / np.sqrt(sig ** 2 + bg ** 2))

    sig = np.sum(pnetxbbs["HHbbWWqq"][(pnetxbbs["HHbbWWqq"]["fatJet1PNetXbb"] > cutoff) + (pnetxbbs["HHbbWWqq"]["fatJet2PNetXbb"] > cutoff)].weight)
    bg = np.sum(pnetxbbs["QCD"][(pnetxbbs["QCD"]["fatJet1PNetXbb"] > cutoff) + (pnetxbbs["QCD"]["fatJet2PNetXbb"] > cutoff)].weight) + np.sum(pnetxbbs["tt"][(pnetxbbs["tt"]["fatJet1PNetXbb"] > cutoff) + (pnetxbbs["tt"]["fatJet2PNetXbb"] > cutoff)].weight)
    yields["jet12xbb"].append(sig / np.sqrt(sig ** 2 + bg ** 2))

yields

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













disc_var = 'DeepAK8MD_H4qvsQCD'
disc_var_label = 'DeepAK8MD_H4qvsQCD'
disc_var_bins = [100, -1, 1]



var = 'PNetXbb'
hists["jet" + var] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet1", r"Leading Jet " + var, 100, 0, 1),
                                hist.Bin("jet2", r"Sub-Leading Jet " + var, 100, 0, 1),
                                hist.Bin("jet3", r"3rd-Leading Jet " + var, 100, 0, 1)
                                )

for s, evts in evtDict.items():
    if s != 'HH4V':
        if 'HH' in s:
            hists["jet" + var].fill(sample=s,
                                 jet1 = evts["fatJet1PNetXbb_alt"],
                                 jet2 = evts["fatJet2PNetXbb_alt"],
                                 jet3 = evts["fatJet3PNetXbb_alt"],
                                 weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                                )
        else:
            hists["jet" + var].fill(sample=s,
                                 jet1 = evts["fatJet1" + var],
                                 jet2 = evts["fatJet2" + var],
                                 jet3 = evts["fatJet3" + var],
                                 weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                                )

hists['jet' + var].scale(scale, axis='sample')

fig, axs = plt.subplots(1, 3, figsize=(27, 9))

for j in range(3):
    ylim = np.max(list(hists["jet" + var].project("sample", "jet" + str(j + 1)).values().values())) * 1.1
    axs[j].set_prop_cycle(cycler(color=colors_cycle))
    hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
    axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
    ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
    axs[j].set_ylim(0, ylim)
    # axs[i, j].ticklabel_format(axis='x', style='sci')
plt.tight_layout(0.5)

plt.tight_layout(0.5)
plt.savefig("figs/jet_pxbb.pdf", bbox_inches='tight')
plt.show()


#
# disc_vars = ['DeepAK8H', 'DeepAK8H4qMD', 'PNetXbb']
# disc_vars_labels = ['DeepAK8H', 'DeepAK8H4qMD', 'PNetXbb']
# disc_bins = [[100, -0.1, 0.1], [100, -0.1, 0.1], [100, -1, 1]]


disc_vars = ['DeepAK8_H', 'DeepAK8MD_H4qvsQCD']
disc_vars_labels = ['DeepAK8_H', 'DeepAK8MD_H4qvsQCD']
disc_bins = [[100, 0, 1], [100, 0, 1]]


for i in range(len(disc_vars)):
    var = disc_vars[i]
    varl = disc_vars_labels[i]
    hists["jet" + var] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin("jet1", r"Leading Jet " + varl, disc_bins[i][0], disc_bins[i][1], disc_bins[i][2]),
                                    hist.Bin("jet2", r"Sub-Leading Jet " + varl, disc_bins[i][0], disc_bins[i][1], disc_bins[i][2]),
                                    hist.Bin("jet3", r"3rd-Leading Jet " + varl, disc_bins[i][0], disc_bins[i][1], disc_bins[i][2])
                                    )

for s, evts in evtDict.items():
    # if s == 'HH4V':
    #     for i in range(len(disc_vars)):
    #         var = disc_vars[i]
    #         hists["jet" + var].fill(sample=s,
    #                              jet1 = evts["fatJet1" + var],
    #                              jet2 = evts["fatJet2" + var],
    #                              jet3 = evts["fatJet3" + var],
    #                              # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
    #                             )
    # else:
    if s != 'HH4V':
        # print(evts["totalWeight"][:10])
        for i in range(len(disc_vars)):
            var = disc_vars[i]
            hists["jet" + var].fill(sample=s,
                                 jet1 = evts["fatJet1" + var],
                                 jet2 = evts["fatJet2" + var],
                                 jet3 = evts["fatJet3" + var],
                                 weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                                )


for var in disc_vars:
    hists['jet' + var].scale(scale, axis='sample')




for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        ylim = np.max(list(hists["jet" + var].project("sample", "jet" + str(j + 1)).values().values())) * 1.1
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
        hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
        ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
        axs[i, j].set_ylim(0, ylim)

plt.tight_layout(0.5)
plt.savefig("figs/jet_disc_vars_weighted.pdf", bbox_inches='tight')
plt.show()

# fig, axs = plt.subplots(1, 3, figsize=(28, len(disc_vars) * 9))
# plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
# # axs.set_prop_cycle(cycler(color=colors))
# # for i in range(len(disc_vars)):
# var = disc_vars[i]
# for j in range(3):
#     axs[j].set_prop_cycle(cycler(color=colors_cycle))
#     hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
#     axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
#     ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
#     axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
#     # axs[i, j].ticklabel_format(axis='x', style='sci')
# plt.tight_layout(0.5)
#
# plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
# plt.savefig("figs/jet_disc_vars_weighted.pdf", bbox_inches='tight')
# plt.show()

colors1 = ['#e31a1c', '#a6cee3', '#1f78b4']
colors2 = ['#a6cee3', '#1f78b4']

fig, axs = plt.subplots(len(disc_vars), 3, figsize=(28, len(disc_vars) * 9))
for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        ylim = np.max(list(hists["jet" + var].project("sample", "jet" + str(j + 1)).values().values())) * 1.1
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
        hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
        ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
        axs[i, j].set_ylim(0, ylim)

plt.tight_layout(0.5)
plt.savefig("figs/jet_disc_vars_weighted.pdf", bbox_inches='tight')
plt.show()

disc_vals = {}
disc_vars = ["PNetXbb"]

for var in disc_vars:
    disc_vals[var] = {}
    for i in range(3):
        jk = 'jet' + str(i + 1)
        disc_vals[var][jk] = {'sig': [], 'bg': []}
        for s in evtDict.keys():
            if s != 'HH4V':
                if 'HH' in s:
                    disc_vals[var][jk]['sig'].append(hists['jet' + var].project("sample", jk).values()[(s, )])
                else:
                    disc_vals[var][jk]['bg'].append(hists['jet' + var].project("sample", jk).values()[(s, )])

        disc_vals[var][jk]['sig'] = np.sum(np.array(disc_vals[var][jk]['sig']), axis=0)
        disc_vals[var][jk]['bg'] = np.sum(np.array(disc_vals[var][jk]['bg']), axis=0)
        disc_vals[var][jk]['tpr'] = np.cumsum(disc_vals[var][jk]['sig'][::-1])[::-1] / np.sum(np.array(disc_vals[var][jk]['sig']))
        disc_vals[var][jk]['fpr'] = np.cumsum(disc_vals[var][jk]['bg'][::-1])[::-1] / np.sum(np.array(disc_vals[var][jk]['bg']))
        disc_vals[var][jk]['tpr'] = np.append(disc_vals[var][jk]['tpr'], 0)
        disc_vals[var][jk]['fpr'] = np.append(disc_vals[var][jk]['fpr'], 0)


disc_vals[var][jk]['tpr'][:-1]  disc_vals[var][jk]['tpr'][1:]
(disc_vals[var][jk]['fpr'][:-1] - disc_vals[var][jk]['fpr'][1:]) / 2




for j in range(3):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[var][jk]['tpr'][:-1] + disc_vals[var][jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[var][jk]['fpr'][:-1] - disc_vals[var][jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.plot(disc_vals[var][jk]['fpr'], disc_vals[var][jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title(var + ' ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_pxbb_roc.pdf", bbox_inches='tight')
plt.show()



for j in range(3):
    jk = 'jet' + str(j + 1)
    tpr_aves = (disc_vals[var][jk]['tpr'][:-1] + disc_vals[var][jk]['tpr'][1:]) / 2
    fpr_diffs = (disc_vals[var][jk]['fpr'][:-1] - disc_vals[var][jk]['fpr'][1:])
    auc = np.sum(tpr_aves * fpr_diffs)
    plt.semilogx(disc_vals[var][jk]['fpr'], disc_vals[var][jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

plt.title(var + ' ROC Curve')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_pxbb_roc_semilog.pdf", bbox_inches='tight')
plt.show()



fig, axs = plt.subplots(len(disc_vars), 3, figsize=(28, len(disc_vars) * 9))

for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        jk = 'jet' + str(i + 1)
        axs[i, j].plot(disc_vals[var][jk]['fpr'], disc_vals[var][jk]['tpr'])
        axs[i, j].set_title('Fat Jet {} {}'.format(j + 1, var))
        axs[i, j].set_xlabel('FPR')
        axs[i, j].set_ylabel('TPR')

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_disc_vars_roc.pdf", bbox_inches='tight')
plt.show()


fig, axs = plt.subplots(1, len(disc_vars), figsize=(len(disc_vars) * 9, 9))

for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        jk = 'jet' + str(j + 1)
        auc = np.sum(disc_vals[var][jk]['tpr']) / 100
        axs[i].plot(disc_vals[var][jk]['fpr'], disc_vals[var][jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

    axs[i].set_title(var + ' ROC Curve')
    axs[i].set_xlabel('FPR')
    axs[i].set_ylabel('TPR')
    axs[i].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_disc_vars_roc2.pdf", bbox_inches='tight')
plt.show()




fig, axs = plt.subplots(1, len(disc_vars), figsize=(len(disc_vars) * 9, 9))

for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        jk = 'jet' + str(j + 1)
        auc = np.sum(disc_vals[var][jk]['tpr']) / 100
        axs[i].semilogx(disc_vals[var][jk]['fpr'], disc_vals[var][jk]['tpr'], label='Fat Jet {} AUC = {:.2f}'.format(j + 1, auc))

    axs[i].set_title(var + ' ROC Curve')
    axs[i].set_xlabel('FPR')
    axs[i].set_ylabel('TPR')
    axs[i].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_disc_vars_logx.pdf", bbox_inches='tight')
plt.show()





jk = 'jet1'
for i in range(len(disc_vars)):
    var = disc_vars[i]
    auc = np.sum(vals[var][jk]['tpr']) / 100
    plt.plot(vals[var][jk]['fpr'], vals[var][jk]['tpr'], label='{} AUC = {:.2f}'.format(var, auc))

plt.title('ROC Curves')
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.savefig("figs/jet_disc_vars_roc3.pdf", bbox_inches='tight')
plt.show()


vals[var][jk]['tpr']


vals[var][jk]['fpr']

np.cumsum(vals[var][jk]['sig'][::-1])[::-1]

np.sum(vals[var][jk]['sig'])

len(vals[var][jk]['sig'])



roc_vals = {}

for var in disc_vars:
    roc_vals[var] = {}
    for i in range(3):
        jk = 'jet' + str(i + 1)
        roc_vals[var][jk] = {'tpr': [], 'fpr': []}
        cumsum = np.cumsum(vals[var][jk]['sig'])
        tot = np.sum(vals[var][jk]['sig'])

        roc_vals[var][jk]['tpr'] = cumsum / tot

        cumsum = np.cumsum(vals[var][jk]['bg'])
        cumsumr = vals[var][jk]['bg'][::-1].cumsum()[::-1]

        for j in range(len(vals[var][jk]['sig'])):

            vals[var][jk]['tpr'].append(vals[var][jk]['sig'] / (vals[var][jk]['sig'] + np.sum(vals[var][jk]['sig'][j]))

        vals[var][jk]['tpr'] = vals[var][jk]['sig'] / (vals[var][jk]['sig'] + vals[var][jk]['bg'])

            vals[var][jk]['sig'].append(hists['jet' + var].project("sample", jk).values()[(s, )])
        else:
            vals[var][jk]['bg'].append(hists['jet' + var].project("sample", jk).values()[(s, )])

        vals[var][jk]['sig'] = np.sum(np.array(vals[var][jk]['sig']), axis=0)
        vals[var][jk]['bg'] = np.sum(np.array(vals[var][jk]['bg']), axis=0)




np.sum(np.array(vals[disc_vars[0]]['jet3']['sig']), axis=0)


vals = hists['jetDeepAK8_H'].project("sample", "jet1").values()
vals



massfns = listdir('data/mass/')
massfns

for i in range(len(massfns)):
    massfns[i] = 'data/mass/' + massfns[i]

mass_dict = {
    "HH4V": np.load('data/mass/HHToVVVV_node_SM_Pt300_1pb_weighted.root_inv_mass.npy'),
    "HHbbWWqq": np.load('data/mass/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root_inv_mass.npy'),
}


for fn in massfns:
    if "QCD" in fn:
        if "QCD" not in mass_dict: mass_dict["QCD"] = np.load(fn)
        else: mass_dict["QCD"] = np.concatenate((mass_dict["QCD"], np.load(fn)), axis=1)
    elif "TT" in fn:
        if "tt" not in mass_dict: mass_dict["tt"] = np.load(fn)
        else: mass_dict["tt"] = np.concatenate((mass_dict["tt"], np.load(fn)), axis=1)


hists["jets_mass"] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet12", r"Leading 2 Jets' Invariant Mass (GeV)", 100, 0, 5000),
                                )

for s in evtDict.keys():
    if s == 'HH4V':
        hists["jets_mass"].fill(sample=s,
                             jet12 = mass_dict[s][0],
                             )
    else:
        hists["jets_mass"].fill(sample=s,
                             jet12 = mass_dict[s][0],
                             weight = evtDict[s]["weight"]
                             )

len(mass_dict["QCD"][0])
evtDict["QCD"]["fatJet1Pt"]


tot_events_mass = {}
vals = hists['jets_mass'].project("sample", "jet12").values()
for s in evtDict.keys():
    tot_events_mass[s] = np.sum(vals[(s,)])

scale_mass = {'HH4V': 1 / tot_events_mass['HH4V'],
        'HHbbWWqq': 1 / tot_events_mass['HHbbWWqq'],
        'QCD': 1 / (tot_events_mass['QCD'] + tot_events_mass['tt']),
        'tt': 1 / (tot_events_mass['QCD'] + tot_events_mass['tt'])}



hists['jets_mass'].scale(scale_mass, axis='sample')

plt.gca().set_prop_cycle(cycler(color=colors_cycle))
ylim = np.max(list(hists["jets_mass"].project("sample", "jet12").values().values())) * 1.1
hist.plot1d(hists["jets_mass"][:2].project("sample", "jet12"), overlay='sample', clear=False, line_opts=line_opts)
plt.gca().set_prop_cycle(cycler(color=colors_cycle[2:]))
hist.plot1d(hists["jets_mass"][2:].project("sample", "jet12"), overlay='sample', stack=True, clear=False, fill_opts=fill_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.ylim(0, ylim)
plt.savefig("figs/jets_mass.pdf", bbox_inches='tight')





for i in range(len(kin_vars)):
    var = kin_vars[i]
    varl = kin_vars_labels[i]
    hists["jet" + var] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin("jet1", r"Leading Jet " + varl, kin_bins[i][0], kin_bins[i][1], kin_bins[i][2]),
                                    hist.Bin("jet2", r"Sub-Leading Jet " + varl, kin_bins[i][0], kin_bins[i][1], kin_bins[i][2]),
                                    hist.Bin("jet3", r"3rd-Leading Jet " + varl, kin_bins[i][0], kin_bins[i][1], kin_bins[i][2])
                                    )

for s, evts in evtDict.items():
    if s == 'HH4V':
        for i in range(len(kin_vars)):
            var = kin_vars[i]
            hists["jet" + var].fill(sample=s,
                                 jet1 = evts["fatJet1" + var],
                                 jet2 = evts["fatJet2" + var],
                                 jet3 = evts["fatJet3" + var],
                                 # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                                )
    else:
        for i in range(len(kin_vars)):
            var = kin_vars[i]
            hists["jet" + var].fill(sample=s,
                                 jet1 = evts["fatJet1" + var],
                                 jet2 = evts["fatJet2" + var],
                                 jet3 = evts["fatJet3" + var],
                                 weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                                )




tot_events = {}
vals = hists['jetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])

tot_events


hists["met"] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("met", r"MET (GeV)", 100, 0, 500),
                                )


for s, evts in evtDict.items():
    if s == 'HH4V':
        hists["met"].fill(sample=s,
                             met = evts["MET"],
                             # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                            )
    else:
        hists["met"].fill(sample=s,
                             met = evts["MET"],
                             weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                            )


tot_events_met = {}
vals = hists['met'].project("sample", "met").values()
for s in evtDict.keys():
    tot_events_met[s] = np.sum(vals[(s,)])


scale_met = {'HH4V': 1 / tot_events_met['HH4V'],
    'HHbbWWqq': 1 / tot_events_met['HHbbWWqq'],
    'QCD': 1 / (tot_events_met['QCD'] + tot_events_met['tt']),
    'tt': 1 / (tot_events_met['QCD'] + tot_events_met['tt'])}


hists['met'].scale(scale_met, axis='sample')

plt.gca().set_prop_cycle(cycler(color=colors_cycle))
ylim = np.max(list(hists["met"].project("sample", "met").values().values())) * 1.1
hist.plot1d(hists["met"][:2].project("sample", "met"), overlay='sample', clear=False, line_opts=line_opts)
plt.gca().set_prop_cycle(cycler(color=colors_cycle[2:]))
hist.plot1d(hists["met"][2:].project("sample", "met"), overlay='sample', stack=True, clear=False, fill_opts=fill_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.ylim(0, ylim)
plt.savefig("figs/met.pdf", bbox_inches='tight')



hists["genhiggs"] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("genHiggs1", r"Gen Higgs 1 $p_T$ (GeV)", 100, 0, 1200),
                                hist.Bin("genHiggs2", r"Gen Higgs 2 $p_T$ (GeV)", 100, 0, 1200),
                                )


for s, evts in evtDict.items():
    if s == 'HH4V':
        hists["genhiggs"].fill(sample=s,
                             genHiggs1 = evts["genHiggs1Pt"],
                             genHiggs2 = evts["genHiggs2Pt"],
                             # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                            )
    elif 'HH' in s:
        hists["genhiggs"].fill(sample=s,
                             genHiggs1 = evts["genHiggs1Pt"],
                             genHiggs2 = evts["genHiggs2Pt"],
                             weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                            )


hists['genhiggs'].scale(scale, axis='sample')


fig, axs = plt.subplots(1, 2, figsize=(18, 9))
# axs.set_prop_cycle(cycler(color=colors))

axs[0].set_prop_cycle(cycler(color=colors_cycle))
hist.plot1d(hists["genhiggs"].project("sample", "genHiggs1"), overlay='sample', ax=axs[0], clear=False, line_opts=line_opts)
axs[0].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

axs[1].set_prop_cycle(cycler(color=colors_cycle))
hist.plot1d(hists["genhiggs"].project("sample", "genHiggs2"), overlay='sample', ax=axs[1], clear=False, line_opts=line_opts)
axs[1].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
    # axs[i].ticklabel_format(axis='x', style='sci')
plt.tight_layout(0.5)
plt.savefig("figs/genhiggspt.pdf", bbox_inches='tight')



plt.gca().set_prop_cycle(cycler(color=colors_cycle))
ylim = np.max(list(hists["met"].project("sample", "met").values().values())) * 1.1
hist.plot1d(hists["met"][:2].project("sample", "met"), overlay='sample', clear=False, line_opts=line_opts)
plt.gca().set_prop_cycle(cycler(color=colors_cycle[2:]))
hist.plot1d(hists["met"][2:].project("sample", "met"), overlay='sample', stack=True, clear=False, fill_opts=fill_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.ylim(0, ylim)
plt.savefig("figs/met.pdf", bbox_inches='tight')







genhh12mass = {}
for s, evts in evtDict.items():
    genhh12mass[s] = []
    if "HH" in s:
        for i in tqdm(range(len(evtDict[s]["fatJet1Pt"]))):
            higgs1 = LorentzVector()
            higgs1.setptetaphim(evts["genHiggs1Pt"][i], evts["genHiggs1Eta"][i], evts["genHiggs1Phi"][i], 125.1)
            higgs2 = LorentzVector()
            higgs2.setptetaphim(evts["genHiggs2Pt"][i], evts["genHiggs2Eta"][i], evts["genHiggs2Phi"][i], 125.1)

            genhh12mass[s].append((higgs1 + higgs2).mass)




plt.hist(genhh12mass["HH4V"], bins=np.linspace(40, 4000, 101))


genhh_vars = ['pt', 'mass']
genhh_vars_labels = [r'$p_T$', 'Mass']
bins = [[100, 0, 2000], [100, 200, 4000]]


for i in range(len(genhh_vars)):
    var = genhh_vars[i]
    varl = genhh_vars_labels[i]
    hists["hh" + var] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin('hh', r"Gen HH " + varl, bins[i][0], bins[i][1], bins[i][2]),
                                    )


for s, evts in evtDict.items():
    if "HH" in s:
        for i in range(len(genhh_vars)):
            var = genhh_vars[i]
            hists["hh" + var].fill(sample=s,
                                 hh = evts["genHH_" + var],
                                 )

for var in genhh_vars:
    hists['hh' + var].scale(scale, axis='sample')

fig, axs = plt.subplots(1, len(genhh_vars), figsize=(len(genhh_vars) * 9, 9))
# axs.set_prop_cycle(cycler(color=colors))
for i in range(len(genhh_vars)):
    var = genhh_vars[i]
    axs[i].set_prop_cycle(cycler(color=colors_cycle))
    hist.plot1d(hists["hh" + var], overlay='sample', ax=axs[i], clear=False, line_opts=line_opts)
    axs[i].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
    # axs[i].ticklabel_format(axis='x', style='sci')
plt.tight_layout(0.5)

# plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/genhh_vars.pdf", bbox_inches='tight')
plt.show()
