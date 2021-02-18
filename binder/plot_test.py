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
import pandas

from os import listdir
from copy import copy
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
plt.style.use(hep.style.ROOT)


hhbbwwevts = uproot4.open("./data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root")["tree"]
# evts = uproot4.open("./data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root")["tree"]
hhbbwwevts.keys()



samples = {
    "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWW4q": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")



evts2 = uproot4.open("data/hh4b/nano_1.root")


evts2['Events'][]

evts2['Events'].keys()

evts2['Events']['FatJet_pt'].arrays()[3]
evts2['Events']['genWeight'].arrays()

ak.array(evts2['Events']['FatJet_pt'])

hh4b_samples = uproot4.concatenate("data/hh4b/*.root:Events")

hh4bfjpt = ak.pad_none(hh4b_samples[:]['FatJet_pt'], 2, axis=1)[:, :2]
hh4bfjmassd = ak.pad_none(hh4b_samples[:]['FatJet_msoftdrop'], 2, axis=1)[:, :2]
hh4bfjpnetxbb = ak.pad_none(hh4b_samples[:]['FatJet_ParticleNetMD_probXbb'], 2, axis=1)[:, :2]
hh4bfjpnetxbb = ak.pad_none(hh4b_samples[:]['FatJet_ParticleNetMD_probXbb'], 2, axis=1)[:, :2]


hh4bjpt = ak.pad_none(hh4b_samples[:]['GenJetAK8_pt'], 3, axis=1)

fjptnon = hh4bfjpt[~ak.is_none(hh4bfjpt[:, 1])]
len(fjptnon)

hh4bfjpt[~ak.is_none(hh4bfjpt[:, 1])]

jptnon = hh4bjpt[~ak.is_none(hh4bjpt[:, 0])]
jptnon


genpart_pt = hh4b_samples['GenPart_pt']
genpart_status = hh4b_samples['GenPart_status']

hh4b_samples['GenPart_status']
hh4b_samples['GenPart_statusFlags']

higgses = (hh4b_samples['GenPart_pdgId'] == 25)
np.unique(hh4b_samples['GenPart_status'][higgses])
np.unique(ak.to_numpy(ak.pad_none(hh4b_samples['GenPart_status'][higgses], 40, axis=1)))

ak.pad_none(hh4b_samples['GenPart_status'][higgses], 40, axis=1)


fjptgencut = hh4bfjpt[ak.sort(higgs_pt, axis=1)[:, -1] > 300]
fjmsdgencut = hh4bfjmassd[ak.sort(higgs_pt, axis=1)[:, -1] > 300]

ak.sum()
ak.sum((fjptgencut[:, 0] > 250) * (fjptgencut[:, 1] > 250)  * (fjmsdgencut[:, 0] > 20) * (fjmsdgencut[:, 1] > 20))

23452 / 94068



higgs_pt = genpart_pt[(hh4b_samples['GenPart_pdgId'] == 25)]

ak.sort(higgs_pt, axis=1)[0]

ak.to_numpy(higgs_pt)

higgs_pt[2]

ak.sum(higgs_pt[:, 0] > 300) / 372000

ak.sum(ak.sort(higgs_pt, axis=1)[:, -1] > 300) / 372000

plt.hist(ak.sort(higgs_pt, axis=1)[:, -1], bins=np.linspace(0, 1000, 101))

(hh4b_samples['GenPart_pdgId'] == 25)[0]
hh4b_samples['GenPart_pdgId'][0]

ak.to_numpy(hh4b_samples['GenPart_mass'][0])

hh4b_samples['GenPart_pdgId'] == 25

genpart_pt[0]

ak.sum((fjptnon[:, 0] > 250) * (fjptnon[:, 1] > 250))

ak.sum(jptnon[:, 0] > 300) / 372000

plt.hist(jptnon[:, 0], bins=np.linspace(0, 1000, 101))
plt.xlabel("Gen Jet 1 pT")


plt.hist(fjptnon[:, 0], bins=np.linspace(0, 1000, 101))
plt.xlabel("Fat Jet 1 pT")



ak.to_numpy(ak.pad_none(hh4bfjpt, 3, axis=1))

ak.pad_none(hh4bfjpt, 3, axis=1)[:, 1]

ak.firsts(hh4bfjpt)

del evtDict["HH4b"]
















# in fb
RUN2LUMI = 137
XSECHHBBWWQQ = 1.82
XSECHHBBBB = 31.05 * 0.58**2
ACCEPTANCE = 0.16928763440860214

weights = {}

weights["HH4V"] = None
weights["QCD"] = evtDict["QCD"]["weight"] * RUN2LUMI
weights["tt"] = evtDict["tt"]["weight"] * RUN2LUMI
weights["HHbbWW4q"] = evtDict["HHbbWW4q"]["totalWeight"] * RUN2LUMI * XSECHHBBWWQQ * ACCEPTANCE

evtDict["HHbbWW4q"]["weight"]
evtDict["HHbbWW4q"]["totalWeight"]
evtDict["HHbbWW4q"]["pileupWeight"]
evtDict["HHbbWW4q"]["pileupWeight"] * evtDict["HHbbWW4q"]["weight"]


plt.hist(evtDict["HHbbWW4q"]["pileupWeight"])
plt.xlabel('Pileup Weight')

plt.hist(evtDict["HHbbWW4q"]["totalWeight"])
plt.xlabel('Total Weight', loc='center')

RUN2LUMI * XSECHHBBWWQQ
RUN2LUMI * XSECHHBBWWQQ * ACCEPTANCE



# sample_names = {
#     "HH4V": "HHToVVVV",
#     "HHbbWW4q": "HHToBBVV",
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




# evtDict['HHbbWW4q']['genHiggs1W2Decay']
# evtDict['HHbbWW4q']['genHiggs1W2Decay']
# evtDict['HHbbWW4q']['genHiggs1W1Decay']
# evtDict['HHbbWW4q']['genHiggs2W1Decay']
# evtDict['HHbbWW4q']['genHiggs2W2Decay']
#
# np.unique(np.array((evtDict['HHbbWW4q']['genHiggs1W2Decay'] + evtDict['HHbbWW4q']['genHiggs2W1Decay'])))
#
# evtDict['HHbbWW4q']['genHiggs1Pt']
# evtDict['HHbbWW4q']['genHiggs1W1Pt']
# evtDict['HHbbWW4q']['fatJet1Pt']
# evtDict['HHbbWW4q']['fatJet2Pt']
# evtDict['HHbbWW4q']['genHiggs2Pt']
#
#
# evtDict['HHbbWW4q']['genHiggs1Eta'][:5]
# evtDict['HHbbWW4q']['fatJet1Eta'][:5]
# evtDict['HHbbWW4q']['fatJet2Eta'][:5]
# evtDict['HHbbWW4q']['genHiggs2Eta'][:5]


keys = ["genHiggs1W1", "genHiggs1W2", "genHiggs2W1", "genHiggs2W2"]
fvecs = {}

for key in keys:
    fvecs[key] = ak.zip({
        "pt": evtDict['HHbbWW4q'][key + 'Pt'],
        "eta": evtDict['HHbbWW4q'][key + 'Eta'],
        "phi": evtDict['HHbbWW4q'][key + 'Phi'],
        "mass": evtDict['HHbbWW4q'][key + 'M'],
    }, with_name="PtEtaPhiMLorentzVector")


genHiggs1 = ak.zip({"pt": evtDict['HHbbWW4q']['genHiggs1Pt'], "eta": evtDict['HHbbWW4q']['genHiggs1Eta'], "phi": evtDict['HHbbWW4q']['genHiggs1Phi'], "mass": ak.full_like(evtDict['HHbbWW4q']['genHiggs1Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")
genHiggs2 = ak.zip({"pt": evtDict['HHbbWW4q']['genHiggs2Pt'], "eta": evtDict['HHbbWW4q']['genHiggs2Eta'], "phi": evtDict['HHbbWW4q']['genHiggs2Phi'], "mass": ak.full_like(evtDict['HHbbWW4q']['genHiggs2Pt'], 125.1)}, with_name="PtEtaPhiMLorentzVector")

# genHiggs12 = ak.zip({"pt": [evtDict['HHbbWW4q']['genHiggs1Pt'], evtDict['HHbbWW4q']['genHiggs2Pt']],
#                     "eta": [evtDict['HHbbWW4q']['genHiggs1Eta'], evtDict['HHbbWW4q']['genHiggs2Eta']],
#                     "phi": [evtDict['HHbbWW4q']['genHiggs1Phi'], evtDict['HHbbWW4q']['genHiggs2Phi']],
#                     "mass": ak.full_like([evtDict['HHbbWW4q']['genHiggs1Pt'], evtDict['HHbbWW4q']['genHiggs2Pt']], 125.1)})#, with_name="PtEtaPhiMLorentzVector")

fatJet1 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet1Pt'], "eta": evtDict['HHbbWW4q']['fatJet1Eta'], "phi": evtDict['HHbbWW4q']['fatJet1Phi'], "mass": evtDict['HHbbWW4q']['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet2 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet2Pt'], "eta": evtDict['HHbbWW4q']['fatJet2Eta'], "phi": evtDict['HHbbWW4q']['fatJet2Phi'], "mass": evtDict['HHbbWW4q']['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet3 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet3Pt'], "eta": evtDict['HHbbWW4q']['fatJet3Eta'], "phi": evtDict['HHbbWW4q']['fatJet3Phi'], "mass": evtDict['HHbbWW4q']['fatJet3Mass']}, with_name="PtEtaPhiMLorentzVector")


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



# Old method

# H1ltH2 = fatJet1.delta_r(genHiggs1) < fatJet1.delta_r(genHiggs2)
# fatJet1H1 = fatJet1.delta_r(genHiggs1) < 0.4
# fatJet1H2 = fatJet1.delta_r(genHiggs2) < 0.4
# fatJet1H12 = fatJet1H1 * fatJet1H2
#
# fatJet1H1 = fatJet1H1 * (~fatJet1H12 + fatJet1H12 * H1ltH2)
# fatJet1H2 = fatJet1H1 * (~fatJet1H12 + fatJet1H12 * ~H1ltH2)
#
#
# H1ltH2 = fatJet2.delta_r(genHiggs1) < fatJet2.delta_r(genHiggs2)
# fatJet2H1 = fatJet2.delta_r(genHiggs1) < 0.4
# fatJet2H2 = fatJet2.delta_r(genHiggs2) < 0.4
# fatJet2H12 = fatJet2H1 * fatJet2H2
#
# fatJet2H1 = fatJet2H1 * (~fatJet2H12 + fatJet2H12 * H1ltH2)
# fatJet2H2 = fatJet2H1 * (~fatJet2H12 + fatJet2H12 * ~H1ltH2)
#
# H1W = evtDict['HHbbWW4q']['genHiggs1W1Decay'] == 0
# H2W = evtDict['HHbbWW4q']['genHiggs2W1Decay'] == 0
#
# fatJet1W_old = fatJet1H1 * H1W + fatJet1H2 * H2W
# fatJet2W_old = fatJet2H1 * H1W + fatJet2H2 * H2W


# New Javier's method

fatJet1H1 = fatJet1.delta_r(genHiggs1) < 0.8
fatJet1H1W = fatJet1H1 * (fatJet1.delta_r(fvecs["genHiggs1W1"]) < 0.8) * (fatJet1.delta_r(fvecs["genHiggs1W2"]) < 0.8)
fatJet1H2 = fatJet1.delta_r(genHiggs2) < 0.8
fatJet1H2W = fatJet1H2 * (fatJet1.delta_r(fvecs["genHiggs2W1"]) < 0.8) * (fatJet1.delta_r(fvecs["genHiggs2W2"]) < 0.8)
fatJet1W = fatJet1H1W + fatJet1H2W

fatJet2H1 = fatJet2.delta_r(genHiggs1) < 0.8
fatJet2H1W = fatJet2H1 * (fatJet2.delta_r(fvecs["genHiggs1W1"]) < 0.8) * (fatJet2.delta_r(fvecs["genHiggs1W2"]) < 0.8)
fatJet2H2 = fatJet2.delta_r(genHiggs2) < 0.8
fatJet2H2W = fatJet2H2 * (fatJet2.delta_r(fvecs["genHiggs2W1"]) < 0.8) * (fatJet2.delta_r(fvecs["genHiggs2W2"]) < 0.8)
fatJet2W = fatJet2H1W + fatJet2H2W

fatJet3H1 = fatJet3.delta_r(genHiggs1) < 0.8
fatJet3H1W = fatJet3H1 * (fatJet3.delta_r(fvecs["genHiggs1W1"]) < 0.8) * (fatJet3.delta_r(fvecs["genHiggs1W2"]) < 0.8)
fatJet3H2 = fatJet3.delta_r(genHiggs2) < 0.8
fatJet3H2W = fatJet3H2 * (fatJet3.delta_r(fvecs["genHiggs2W1"]) < 0.8) * (fatJet3.delta_r(fvecs["genHiggs2W2"]) < 0.8)
fatJet3W = fatJet3H1W + fatJet3H2W

fatJet1b = (fatJet1H1 + fatJet1H2) * ~(fatJet1W)
fatJet2b = (fatJet2H1 + fatJet2H2) * ~(fatJet2W)

#
#
# np.sum(fatJet1W * fatJet1W_old)
#
# missed_H1Ws = H1W * ~(fatJet1H1W + fatJet2H1W)
# missed_H2Ws = H2W * ~(fatJet1H2W + fatJet2H2W)
#
#
# np.sum(missed_H1Ws)
# np.sum(missed_H2Ws)
# np.sum(missed_H1Ws * (fatJet1H1 + fatJet2H1))
# np.sum(missed_H2Ws * (fatJet1H2 + fatJet2H2))
#
# np.sum(missed_H1Ws * (fatJet1H1 + fatJet2H1))
#
# np.sum(missed_H1Ws * fatJet2H2)
#
# fatJet1[missed_H1Ws]
# genHiggs1W
#
#
# np.sum(fatJet1H1W) + np.sum(fatJet2H1W) + np.sum(fatJet3H1W)
# np.sum(fatJet1H2W) + np.sum(fatJet2H2W) + np.sum(fatJet3H2W)
#
# np.sum(fatJet1H1W) + np.sum(fatJet2H1W)
# np.sum(fatJet1H2W) + np.sum(fatJet2H2W)
#
#
#
# np.sum(fatJet1H1W * ~H1W)
# np.sum(fatJet2H1W * ~H1W)
# np.sum(fatJet3H1W * ~H1W)
# np.sum(fatJet1H2W * ~H2W)
# np.sum(fatJet2H2W * ~H2W)
# np.sum(fatJet3H2W * ~H2W)
#
# # fatjet3s are matching with -90's of w's ==> only ~2000 fatjet3 matches
#
# np.sum(H1W)
# np.sum(H2W)
#
#
# np.sum(fatJet1H1W)
# np.sum(fatJet1H2W)
# np.sum(fatJet2H1W)
# np.sum(fatJet2H2W)
# np.sum(fatJet3H1W)
# np.sum(fatJet3H2W)
#
# np.sum(fatJet1W)
# np.sum(fatJet2W)
# np.sum(fatJet1b)
# np.sum(fatJet2b)
#
#
#
# np.sum(fatJet1W) + np.sum(fatJet2W)
# np.sum(fatJet1b) + np.sum(fatJet2b)
#
# np.sum(fatJet3W)
#
# np.sum(fatJet1W * fatJet1W_old)
# np.sum(fatJet2W * fatJet2W_old)
#
# np.sum(fatJet1W_old)
# np.sum(fatJet2W_old)
#
# # np.sum(np.array((fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4)))
# np.sum(np.array((fatJet1.delta_r(genHiggs1) < 0.4) * (fatJet1.delta_r(genHiggs2) < 0.4)))
# np.sum(np.array((fatJet2.delta_r(genHiggs1) < 0.4) * (fatJet2.delta_r(genHiggs2) < 0.4)))
# np.sum(np.array((fatJet3.delta_r(genHiggs1) < 0.4) + (fatJet3.delta_r(genHiggs2) < 0.4)))

# sum(genHiggs1.pt < 300)
#
# fatjet1_nohiggs = np.array((fatJet1.delta_r(genHiggs1) > 0.4) * (fatJet1.delta_r(genHiggs2) > 0.4))
# fatjet2_nohiggs = np.array((fatJet2.delta_r(genHiggs1) > 0.4) * (fatJet2.delta_r(genHiggs2) > 0.4))
#
# fatjet12_nohiggs = fatjet1_nohiggs * fatjet2_nohiggs
#
# np.sum(np.array(fatjet12_nohiggs))
#
# plt.hist(evtDict['HHbbWW4q']['genHiggs1W1Pt'], histtype='step', bins=np.linspace(0, 1200, 51))
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


# plt.hist(evtDict["HHbbWW4q"]['fatJet1PNetXbb'], bins=np.linspace(-1, 1, 100))
# plt.title("HHbbWW4q")
# plt.xlabel("Fat Jet 1 PNet Score")
#
# plt.hist(evtDict["HHbbWW4q"]['fatJet1DeepAK8H'], bins=np.linspace(-1, 1, 100))
# plt.title("HHbbWW4q")
# plt.xlabel("Fat Jet 1 Deep AK8H Score")



hists = {}


def init_hists(vars, labels, binrs, fatJet=True, name=None, scale=True, bbsorted=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    lpf = ["Leading Jet", "Sub-Leading Jet", "3rd Leading Jet"] if fatJet else ["Gen Higgs 1", "Gen Higgs 2"]
    if bbsorted: lpf = ["Hbb Candidate Jet", "HWW Candidate Jet"]
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


# weights = {
#     "HH4V": None,
#     "HHbbWW4q": evtDict["HHbbWW4q"]['weight'],
#     "QCD": evtDict["QCD"]['weight'],
#     "tt": evtDict["tt"]['weight'],
# }


# scale_factor = {
#         'HH4V': 1 / len(evtDict["HH4V"]["weight"]),
#         'HHbbWW4q': 1 / np.sum(evtDict["HHbbWW4q"]["weight"]),
#         'HHbbWW4q - Hbb': 2 / np.sum(evtDict["HHbbWW4q"]["weight"]),
#         'HHbbWW4q - HWW': 2 / np.sum(evtDict["HHbbWW4q"]["weight"]),
#         'QCD': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
#         'tt': 1 / (np.sum(evtDict["QCD"]["weight"]) + np.sum(evtDict["tt"]["weight"])),
#         }


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True, pevtDict=None, useEDWeights=False):
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
                for s in evtDict.keys():
                    vars[i][s] = temp

            for s, evts in pevtDict.items():
                if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                    kwargs = {}
                    for j in range(nbins):
                        kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[i][s])]
                    hists[name + vars[i]["HHbbWW4q"]].fill(sample=s, weight=pweights[s], **kwargs)

            if scale: hists[name + vars[i]["HHbbWW4q"]].scale(scale_factor, axis='sample')
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        for s, evts in pevtDict.items():
            if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                kwargs = {}
                for j in range(nbins):
                    kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(npf, j + 1, vars[s])]
                hists[name + vars["HHbbWW4q"]].fill(sample=s, weight=pweights[s], **kwargs)

        if scale: hists[name + vars["HHbbWW4q"]].scale(scale_factor, axis='sample')


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


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True, sepsig=True, stackall=False, log=False, lumilabel=False):
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

            var = vars[i]["HHbbWW4q"]
            if type(binrs[0]) is not list:
                temp = binrs
                binrs = []
                for k in range(len(vars)):
                    binrs.append(temp)

            if fatJet:
                if sepsig:
                    _ = axs[i, 0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                    _ = axs[i, 0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                    _ = axs[i, 1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                    _ = axs[i, 1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)

                    for j in range(nbins):
                        ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        if hh4v: hist.plot1d(hists[name + var][:1].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)

                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                        ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                        axs[i, j].set_ylim(0, ylim)
                else:
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
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        var = vars["HHbbWW4q"]
        fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))
        if fatJet:
            if sepsig:
                _ = axs[0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)

                for j in range(nbins):
                    ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                    axs[j].set_prop_cycle(cycler(color=colors_cycle))
                    if hh4v: hist.plot1d(hists[name + var][:1].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)

                    axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                    ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                    axs[j].set_ylim(0, ylim)
            else:
                for j in range(nbins):
                    if stackall:
                        # hist.plot1d(hists[name + var].sum("sample", bpf + str((j * -1) + 2)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbWW4q', 'tt', 'QCD'])
                        hist.plot1d(hists[name + var].sum(bpf + str((j * -1) + 2)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbWW4q', 'tt', 'QCD'])
                    else:
                        ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][:-2].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
                        axs[j].set_prop_cycle(cycler(color=colors_cycle[3::-1]))
                        ax = hist.plot1d(hists[name + var][-2:].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts, order=['tt', 'QCD'])
                        axs[j].set_ylim(0, ylim)

                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                    if log:
                        axs[j].set_ylim(1)
                        axs[j].set_yscale('log')

                    if lumilabel: hep.label.lumitext("137fb$^{-1}$", ax=axs[j])

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()

hists

init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=False)

tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


tot_events

tot_events['HHbbWW4q'] / 0.173

scale_factor = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbWW4q': 1 / tot_events['HHbbWW4q'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}

init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=True)
plot_hists("Pt", "kin_tests", [40, 200, 2000])

hist.plot1d(hists['fatJetPt'][-1:].project('sample', 'jet1'))

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 2000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True)
plot_hists(vars, "jet_kin", bins)


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [30, 50, 200], [30, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="fine")
fill_hists(copy(vars), fatJet=True, scale=True, name="fine")
plot_hists(vars, "jet_kin_fine", bins, name="fine")


vars = {
    "HH4V": "PNetXbb_alt",
    "HHbbWW4q": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
}
disc_bin = [100, 0, 1]
init_hists("PNetXbb_alt", "Particle Net Xbb", [100, 0, 1], fatJet=True)
fill_hists(vars, fatJet=True, scale=True, hh4v=False)
plot_hists(vars, "pnetxbb", [100, 0, 1], hh4v=False)

init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, hh4v=False)
plot_hists("DeepAK8MD_H4qvsQCD", "deepak8mdh4q", [100, 0, 1], hh4v=False)



init_hists("PNetXbb_alt", "Particle Net Xbb", disc_bin, fatJet=True)
fill_hists(vars, fatJet=True, scale=False, hh4v=False)
fatJet1WXbb = np.histogram(np.array(evtDict["HHbbWW4q"]["fatJet1PNetXbb_alt"][fatJet1W]), weights=np.array(evtDict["HHbbWW4q"]["weight"][fatJet1W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet1bXbb = np.histogram(np.array(evtDict["HHbbWW4q"]["fatJet1PNetXbb_alt"][fatJet1b]), weights=np.array(evtDict["HHbbWW4q"]["weight"][fatJet1b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2WXbb = np.histogram(np.array(evtDict["HHbbWW4q"]["fatJet2PNetXbb_alt"][fatJet2W]), weights=np.array(evtDict["HHbbWW4q"]["weight"][fatJet2W]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]
fatJet2bXbb = np.histogram(np.array(evtDict["HHbbWW4q"]["fatJet2PNetXbb_alt"][fatJet2b]), weights=np.array(evtDict["HHbbWW4q"]["weight"][fatJet2b]), bins=np.linspace(disc_bin[1], disc_bin[2], disc_bin[0] + 1))[0]

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



events_bb_sorted = {}
fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD"]

for s, evts in evtDict.items():
    if s != "HH4V":
        pnet_key = "PNetXbb_alt" if "HH" in s else "PNetXbb"
        jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

        ak_dict =   {
                        "fatJet1PNetXbb": evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading,
                        "fatJet2PNetXbb": evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading,
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            ak_dict["fatJet1" + var] = evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading
            ak_dict["fatJet2" + var] = evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading

        events_bb_sorted[s] = ak.zip(ak_dict)


events_bb_sorted['HHbbWW4q']['fatJet1PNetXbb']

vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 600], [50, 30, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(copy(vars), fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
plot_hists(vars, "jet_kin_bb_leading", bins, name="bb_sorted", hh4v=False, sepsig=False)


init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD Score", [100, 0, 1], fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_deepak8mdh4q_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)

init_hists("PNetXbb", "ParticleNet Xbb Tagger Score", [100, 0, 1], fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists("PNetXbb", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
plot_hists("PNetXbb", "jet_pnetxbb_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)

kin_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
}

events_bbs_kin_cuts = {}
cutflow = {}

for s, evts in events_bb_sorted.items():
    cuts = []
    for var, brange in kin_cuts.items():
        cuts.append(evts[var] > brange[0])
        cuts.append(evts[var] < brange[1])

    cut = cuts[0]
    cutflow[s] = []
    for i in np.arange(1, len(cuts)):
        cutflow[s].append(np.sum(evts[cut].weight))
        cut = cut * cuts[i]

    cutflow[s].append(np.sum(evts[cut].weight))

    events_bbs_kin_cuts[s] = evts[cut]

cutflow

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.9, 9999],
}

events_bbs_cuts = {}
cutflow = {}

for s, evts in events_bb_sorted.items():
    cuts = []
    for var, brange in var_cuts.items():
        cuts.append(evts[var] > brange[0])
        cuts.append(evts[var] < brange[1])

    cut = cuts[0]
    cutflow[s] = []
    for i in np.arange(1, len(cuts)):
        cutflow[s].append(np.sum(evts[cut].weight))
        cut = cut * cuts[i]

    cutflow[s].append(np.sum(evts[cut].weight))

    events_bbs_cuts[s] = evts[cut]


cutflow

events_bbs_cuts['HHbbWW4q']['fatJet2MassSD']

hists['bb_sorted_cutMassSD'].project('sample', 'jet2').values()

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [10, 50, 150]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_bbs_cuts, useEDWeights=True)
plot_hists(var, "jet_cuts_masssd", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, lumilabel=True, stackall=True)




var = {
    "HHbbWW4q": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
}
var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [50, 0.9, 1]
init_hists("PNetXbb", varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_cut_pnet", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=True, log=False)


init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_cut_deepak8mdh4q_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)

events_bbs_kin_cuts

np.sum(events_bbs_kin_cuts["HHbbWW4q"].weight)
np.sum(events_bbs_kin_cuts["QCD"].weight) + np.sum(events_bbs_kin_cuts["tt"].weight)


yields = {"sig": [], "bg": []}
for pnetcutoff in np.arange(0.95, 1, 0.01):
    sig = np.sum(events_bbs_kin_cuts["HHbbWW4q"][(events_bbs_kin_cuts["HHbbWW4q"]["fatJet1PNetXbb"] > pnetcutoff)].weight)
    bg = np.sum(events_bbs_kin_cuts["QCD"][(events_bbs_kin_cuts["QCD"]["fatJet1PNetXbb"] > pnetcutoff)].weight) + np.sum(events_bbs_kin_cuts["tt"][(events_bbs_kin_cuts["tt"]["fatJet1PNetXbb"] > pnetcutoff)].weight)
    yields["sig"].append(sig)
    yields["bg"].append(bg)


yields

yields = {"sig": [], "bg": []}
for pnetcutoff in np.arange(0.95, 1, 0.01)[:-1]:
    sigy = []
    bgy = []
    for dakcutoff in np.arange(0.6, 1, 0.05):
        cuts = {}
        for s, evts in events_bbs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetXbb"] > pnetcutoff) * (evts["fatJet2DeepAK8MD_H4qvsQCD"] > dakcutoff)

        sig = np.sum(events_bbs_kin_cuts["HHbbWW4q"][cuts["HHbbWW4q"]].weight)
        bg = np.sum(events_bbs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_bbs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)



np.array(yields["sig"])
np.array(yields["bg"])

plt.plot(np.array(yields["sig"]))


non_mass_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [40, 9999],
    "fatJet2MassSD": [40, 9999],
    "fatJet1PNetXbb": [0.99, 9999],
}

events_bbs_nm_cuts = {}

for s, evts in events_bb_sorted.items():
    cuts = []
    for var, brange in kin_cuts.items():
        if var == "fatJet1PNetXbb" and "HH" in s: var = "fatJet1PNetXbb_alt"
        cuts.append(evts[var] > brange[0])
        cuts.append(evts[var] < brange[1])

    cut = cuts[0]
    for i in np.arange(1, len(cuts)):
        cut = cut * cuts[i]

    events_bbs_nm_cuts[s] = evts[cut]



del range
var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 40, 500]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_nm_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_nm_cut", pevtDict=events_bbs_nm_cuts)
plot_hists(var, "jet_nm_cut_no_stack", bins, name="bb_sorted_nm_cut", hh4v=False, sepsig=False, log=True, lumilabel=True)




ak.pad_none(hh4b_samples[:]['weight'], 1, axis=1)[:, :2]
ak.sum(hh4b_samples['weight'])


XSECHHBBBB
RUN2LUMI * XSECHHBBBB

len(pnetxbb) * RUN2LUMI * XSECHHBBBB / 372000
len(pnetxbb) / 372000


fatjet_vars = ["Pt", "MassSD"]
var_names = ["pt", "msoftdrop"]


pnetxbb = ak.pad_none(hh4b_samples[:]['FatJet_ParticleNetMD_probXbb'], 2, axis=1)[:, :2][~ak.is_none(hh4bfjpt[:, 1])]
jet1_bb_leading = pnetxbb[:, 0] > pnetxbb[:, 1]

ak_dict =   {
                "fatJet1PNetXbb": pnetxbb[:, 0] * jet1_bb_leading + pnetxbb[:, 1] * ~jet1_bb_leading,
                "fatJet2PNetXbb": pnetxbb[:, 0] * ~jet1_bb_leading + pnetxbb[:, 1] * jet1_bb_leading,
                "weight": np.ones(len(pnetxbb)) * RUN2LUMI * XSECHHBBBB / 372000,
            }

for var, name in zip(fatjet_vars, var_names):
    evts = ak.pad_none(hh4b_samples[:]['FatJet_' + name], 2, axis=1)[:, :2][~ak.is_none(hh4bfjpt[:, 1])]
    ak_dict["fatJet1" + var] = evts[:, 0] * jet1_bb_leading + evts[:, 1] * ~jet1_bb_leading
    ak_dict["fatJet2" + var] = evts[:, 0] * ~jet1_bb_leading + evts[:, 1] * jet1_bb_leading

ak_dict
hh4b_events_bb_sorted = ak.zip(ak_dict)
hh4b_events_bb_sorted

hh4b_kin_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1MassSD": [105, 135],
    # "fatJet2MassSD": [50, 150],
    # "fatJet1PNetXbb": [0.985, 9999],
    # "fatJet1PNetXbb": [0.985, 9999],
    # "fatJet2DeepAK8MD_H4qvsQCD": [0.9, 9999],
}

events_hh4b_kin_cuts = {}
hh4b_kincf = {}

for s, evts in events_bb_sorted.items():
    cuts = []
    for var, brange in hh4b_kin_var_cuts.items():
        cuts.append(evts[var] > brange[0])
        cuts.append(evts[var] < brange[1])

    cuts.append((evts['fatJet1Pt'] > 350) + (evts['fatJet2Pt'] > 350))

    cut = cuts[0]
    hh4b_kincf[s] = []
    for i in np.arange(1, len(cuts)):
        hh4b_kincf[s].append(np.sum(evts[cut].weight))
        cut = cut * cuts[i]

    hh4b_kincf[s].append(np.sum(evts[cut].weight))

    events_hh4b_kin_cuts[s] = evts[cut]

s = 'hh4b'
cuts = []
for var, brange in hh4b_kin_var_cuts.items():
    cuts.append(hh4b_events_bb_sorted[var] > brange[0])
    cuts.append(hh4b_events_bb_sorted[var] < brange[1])

cuts.append((hh4b_events_bb_sorted['fatJet1Pt'] > 350) + (hh4b_events_bb_sorted['fatJet2Pt'] > 350))
cut = cuts[0]
# hh4bcf[s] = []
for i in np.arange(1, len(cuts)):
    # hh4bcf[s].append(np.sum(evts[cut].weight))
    cut = cut * cuts[i]

# hh4bcf[s].append(np.sum(evts[cut].weight))

events_hh4b_kin_cuts[s] = hh4b_events_bb_sorted[cut]

hh4b_kincf


hh4b_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1MassSD": [105, 135],
    # "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
    # "fatJet2DeepAK8MD_H4qvsQCD": [0.9, 9999],
}

events_hh4b_cuts = {}
hh4bcf = {}

for s, evts in events_bb_sorted.items():
    cuts = []
    for var, brange in hh4b_var_cuts.items():
        cuts.append(evts[var] > brange[0])
        cuts.append(evts[var] < brange[1])
        if var == "fatJet2Pt":
            cuts.append((evts['fatJet1Pt'] > 350) + (evts['fatJet2Pt'] > 350))



    cut = cuts[0]
    hh4bcf[s] = []
    for i in np.arange(1, len(cuts)):
        hh4bcf[s].append(np.sum(evts[cut].weight))
        cut = cut * cuts[i]

    hh4bcf[s].append(np.sum(evts[cut].weight))

    events_hh4b_cuts[s] = evts[cut]

s = 'hh4b'
cuts = []
for var, brange in hh4b_var_cuts.items():
    cuts.append(hh4b_events_bb_sorted[var] > brange[0])
    cuts.append(hh4b_events_bb_sorted[var] < brange[1])
    if var == "fatJet2Pt":
        cuts.append((hh4b_events_bb_sorted['fatJet1Pt'] > 350) + (hh4b_events_bb_sorted['fatJet2Pt'] > 350))

cut = cuts[0]
hh4bcf[s] = []
for i in np.arange(1, len(cuts)):
    hh4bcf[s].append(np.sum(hh4b_events_bb_sorted[cut].weight))
    cut = cut * cuts[i]

hh4bcf[s].append(np.sum(hh4b_events_bb_sorted[cut].weight))

events_hh4b_cuts[s] = hh4b_events_bb_sorted[cut]


cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.985',
                'Jet2 PNetXbb > 0.985',
            ]

del(hh4bcf['HHbbWW4q'])
cut_idx = [0, 2, 4, 5, 6, 7, 9]

np.round(np.array(list(hh4bcf.values()))[:, cut_idx], 1)

cftable = pandas.DataFrame(np.round(np.array(list(hh4bcf.values()))[:, cut_idx], 1), list(hh4bcf.keys()), cut_labels)
cftable
cftable.to_csv('hh4b_cutflow.csv')


np.sum(events_hh4b_cuts['hh4b']['weight'])
np.sum(events_hh4b_kin_cuts['hh4b']['weight'])



NOT_NONE = 94068
NOT_NONE_EVTS = 94068 * (RUN2LUMI * XSECHHBBBB) / 372000

XSECHHBBBB

NOT_NONE_EVTS

np.sum(events_hh4b_cuts['hh4b']['weight'])
np.sum(hh4b_events_bb_sorted['weight'])
(np.sum(events_hh4b_cuts['QCD']['weight']) + np.sum(events_hh4b_cuts['tt']['weight']))
(np.sum(events_hh4b_kin_cuts['QCD']['weight']) + np.sum(events_hh4b_kin_cuts['tt']['weight']))




events_bb_sorted['HHbbWW4q']['fatJet1PNetXbb']







events_bbs_kin_cuts

var = "tagger2d_cut"
hists[var] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet1bb", r"Jet1 Tagger Score", 20, 0.8, 1),
                                hist.Bin("jet2WW", r"Jet2 Tagger Score", 20, 0.8, 1),
                                )




for s, evts in evtDict.items():
    if s != 'QCD' and s != 'HH4V':
        hists[var].fill(sample=s,
                             jet1bb = events_bbs_kin_cuts[s]["fatJet1PNetXbb"],
                             jet2WW = events_bbs_kin_cuts[s]["fatJet2DeepAK8MD_H4qvsQCD"],
                             weight = events_bbs_kin_cuts[s]["weight"],
                            )

s = 'QCD'
hists[var].fill(sample='QCD bb1 vs WW2 tagger',
                             jet1bb = events_bbs_kin_cuts[s]["fatJet1PNetXbb"],
                             jet2WW = events_bbs_kin_cuts[s]["fatJet2DeepAK8MD_H4qvsQCD"],
                             weight = events_bbs_kin_cuts[s]["weight"],
                            )

hists[var].fill(sample='QCD bb1 vs bb2 tagger',
                             jet1bb = events_bbs_kin_cuts[s]["fatJet1PNetXbb"],
                             jet2WW = events_bbs_kin_cuts[s]["fatJet2PNetXbb"],
                             weight = events_bbs_kin_cuts[s]["weight"],
                            )


events_hh4b_cuts['hh4b']["weight"]
s = 'hh4b'
hists[var].fill(sample='HH4b',
                 jet1bb = ak.to_numpy(events_hh4b_cuts[s]["fatJet1PNetXbb"]),
                 jet2WW = ak.to_numpy(events_hh4b_cuts[s]["fatJet2PNetXbb"]),
                 weight = ak.to_numpy(events_hh4b_cuts[s]["weight"]),
                )

hists[var][2:3].project('sample').values().keys()

patch_opts = {
    'cmap': 'jet',
    'vmin': 0,
    # 'vmax': 5,
}

# titles = ['HH4b', 'HHbbWW4q', 'QCD', 'tt']
titles = ['HHbbWW4q', 'QCD bb1 vs WW2 tagger', 'QCD bb1 vs bb2 tagger', 'tt']
fig, axs = plt.subplots(1, 4, figsize=(4 * 9, 9))
for j in range(4):
    hist.plot2d(hists[var][j:j + 1].sum('sample'), 'jet1bb', ax=axs[j], patch_opts=patch_opts)
    axs[j].set_title(titles[j])

plt.tight_layout(pad=3)
plt.savefig("figs/tagger_2d_cut_q.pdf", bbox_inches='tight')
plt.show()





d

























































np.sum(weights["HHbbWW4q"])
np.sum(weights["QCD"])
np.sum(weights["tt"])
np.sum(weights["QCD"]) + np.sum(weights["tt"])



events = {}
for s, evts in evtDict.items():
    if s != "HH4V":
        key = "PNetXbb_alt" if "HH" in s else "PNetXbb"
        events[s] = ak.zip({
                                "fatJet1Pt": evts["fatJet1Pt"],
                                "fatJet2Pt": evts["fatJet2Pt"],
                                "fatJet1MassSD": evts["fatJet1MassSD"],
                                "fatJet2MassSD": evts["fatJet2MassSD"],
                                "fatJet1PNetXbb": evts["fatJet1" + key],
                                "fatJet2PNetXbb": evts["fatJet2" + key],
                                "weight": weights[s],
                             })


jet1ptmincut = 350
jet2ptmincut = 300
jet1masssdmincut = 100
jet1masssdmaxcut = 150
jet2masssdmincut = 90
# jet2masssdmaxcut = 150

events_kin_cuts = {}

events_kin_cuts

for s, evts in events.items():
    cut = (evts["fatJet1Pt"] > jet1ptmincut) * \
          (evts["fatJet2Pt"] > jet2ptmincut) * \
          (evts["fatJet1MassSD"] > jet1masssdmincut) * \
          (evts["fatJet1MassSD"] < jet1masssdmaxcut) * \
          (evts["fatJet2MassSD"] > jet1masssdmincut)

    events_kin_cuts[s] = evts[cut]

yields = {"sig": [], "bg": []}
for cutoff in np.arange(0, 1, 0.5):
    sig = np.sum(events_kin_cuts["HHbbWW4q"][(events_kin_cuts["HHbbWW4q"]["fatJet1PNetXbb"] > cutoff) + (events_kin_cuts["HHbbWW4q"]["fatJet2PNetXbb"] > cutoff)].weight)
    bg = np.sum(events_kin_cuts["QCD"][(events_kin_cuts["QCD"]["fatJet1PNetXbb"] > cutoff) + (events_kin_cuts["QCD"]["fatJet2PNetXbb"] > cutoff)].weight) + np.sum(events_kin_cuts["tt"][(events_kin_cuts["tt"]["fatJet1PNetXbb"] > cutoff) + (events_kin_cuts["tt"]["fatJet2PNetXbb"] > cutoff)].weight)
    yields["sig"].append(sig)
    yields["bg"].append(bg)

np.sum(weights["HHbbWW4q"])
np.sum(weights["QCD"])
np.sum(weights["tt"])
np.sum(weights["QCD"]) + np.sum(weights["tt"])


yields


yields




for cutoff in np.arange(0.9, 1, 0.005):
    sig = np.sum(pnetxbbs["HHbbWW4q"][pnetxbbs["HHbbWW4q"]["fatJet1PNetXbb"] > cutoff ].weight)
    bg = np.sum(pnetxbbs["QCD"][pnetxbbs["QCD"]["fatJet1PNetXbb"] > cutoff].weight) + np.sum(pnetxbbs["tt"][pnetxbbs["tt"]["fatJet1PNetXbb"] > cutoff].weight)
    yields["jet1xbb"].append(sig / np.sqrt(sig ** 2 + bg ** 2))

    sig = np.sum(pnetxbbs["HHbbWW4q"][pnetxbbs["HHbbWW4q"]["fatJet2PNetXbb"] > cutoff].weight)
    bg = np.sum(pnetxbbs["QCD"][pnetxbbs["QCD"]["fatJet2PNetXbb"] > cutoff].weight) + np.sum(pnetxbbs["tt"][pnetxbbs["tt"]["fatJet2PNetXbb"] > cutoff].weight)
    yields["jet2xbb"].append(sig / np.sqrt(sig ** 2 + bg ** 2))

    sig = np.sum(pnetxbbs["HHbbWW4q"][(pnetxbbs["HHbbWW4q"]["fatJet1PNetXbb"] > cutoff) + (pnetxbbs["HHbbWW4q"]["fatJet2PNetXbb"] > cutoff)].weight)
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
    "HHbbWW4q": np.load('data/mass/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root_inv_mass.npy'),
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
        'HHbbWW4q': 1 / tot_events_mass['HHbbWW4q'],
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
    'HHbbWW4q': 1 / tot_events_met['HHbbWW4q'],
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
