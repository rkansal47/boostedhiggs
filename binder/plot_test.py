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

import re

from os import listdir
from copy import copy
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
plt.style.use(hep.style.ROOT)


samples = {
    "HH4V": "data/new/weighted/GluGluToHHTo4V_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8_1pb_weighted.root",
    "HHbbVV4q": 'data/new/weighted/HHToBBVVToBBQQQQ_cHHH1_1pb_weighted.root',
    "QCD": "data/new/weighted/QCD_HT*.root",
    "tt": "data/new/weighted/TTTo*.root",
    "HH4b": "data/weighted/GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8_1pb_weighted.root"
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")

# evtDict["HHbbVV4q"] = uproot4.concatenate('data/new/HHToBBVVToBBQQQQ_cHHH1/*.root')
#
# evtDict["HHbbVV4q"] = evtDict["HHbbVV4q"]

nevents = 0

uproot4.open('data/new/HHToBBVVToBBQQQQ_cHHH1/HHToBBWWNtuple_option1_ak15_{}.root'.format(i))['NEvents'].values()[1]

# for i in range(1, 55):
uproot4.open(samples["HHbbVV4q"])['NEvents'].values()[1]

nevents = 2368010.0


bbWW = uproot4.open('data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root')

bbWW['NEvents'].values()
len(bbWW['tree']['weight'].arrays())

86577 / 185417

247102 / 2368010



evtDict["HHbbVV4q"]


2368010


len(evtDict["HHbbVV4q"]) / nevents




Ωdata17 = uproot4.concatenate('data/weighted/JetHT*.root:tree')
evtDict["data"] = data17

trigger_effs = uproot4.open('triggers/JetHTTriggerEfficiency_2017.root')

trigger_effs.keys()


trigger_effs['efficiency_ptmass_Xbb0p98To1p0'].to_numpy()

evtDict["data"].fields


evtDict["HH4V"].fields


evtDict["HH4V"]['lep1Pt']
evtDict["HH4V"]['lep2Pt']
evtDict["HH4V"]['lep2Id']


triggers = data17['HLT_PFJet500'] + \
data17['HLT_AK8PFJet500'] + \
data17['HLT_AK8PFJet360_TrimMass30'] + \
data17['HLT_AK8PFJet380_TrimMass30'] + \
data17['HLT_AK8PFJet400_TrimMass30'] + \
data17['HLT_AK8PFHT800_TrimMass50'] + \
data17['HLT_AK8PFJet330_PFAK8BTagCSV_p17']

triggers17 = ['HLT_PFJet500',
                'HLT_AK8PFJet500',
                'HLT_AK8PFJet360_TrimMass30',
                'HLT_AK8PFJet380_TrimMass30',
                'HLT_AK8PFJet400_TrimMass30',
                'HLT_AK8PFHT800_TrimMass50',
                'HLT_AK8PFJet330_PFAK8BTagCSV_p17']

ak.sum(data17['HLT_AK8PFJet360_TrimMass30'])
ak.sum(data17['HLT_AK8PFJet380_TrimMass30'])
ak.sum(data17['HLT_AK8PFJet380_TrimMass30'] * data17['HLT_AK8PFJet360_TrimMass30'])
ak.sum(data17['HLT_AK8PFJet380_TrimMass30'] * ~data17['HLT_AK8PFJet360_TrimMass30'])


evtDict["HHbbVV4q"].fields




# acceptance calculation

full_hh4b_samples = uproot4.concatenate("data/hh4b/nano_*.root:Events")

genpart_pt = full_hh4b_samples['GenPart_pt']
genpart_status = full_hh4b_samples['GenPart_status']

full_hh4b_samples['GenPart_status']
full_hh4b_samples['GenPart_statusFlags']

higgses = (full_hh4b_samples['GenPart_pdgId'] == 25)
full_hh4b_samples['GenPart_status'][higgses]

np.unique(ak.to_numpy(ak.pad_none(full_hh4b_samples['GenPart_status'][higgses], 40, axis=1)))

ak.pad_none(hh4b_samples['GenPart_status'][higgses], 40, axis=1)
higgs_pt = genpart_pt[(full_hh4b_samples['GenPart_pdgId'] == 25)]

gt_300 = ak.sort(higgs_pt, axis=1)[:, -1] > 300


num_gt_300 = ak.sum(ak.sort(higgs_pt, axis=1)[:, -1] > 300)
num_gt_300 / len(higgs_pt)

fhiggs_pt = genpart_pt[(full_hh4b_samples['GenPart_pdgId'] == 25)][full_hh4b_samples['GenPart_status'][higgses] == 22]

jet_pt = ak.pad_none(full_hh4b_samples['FatJet_pt'], 2, axis=1)[:, :2][gt_300]
jet_msd = ak.pad_none(full_hh4b_samples['FatJet_msoftdrop'], 2, axis=1)[:, :2][gt_300]

ak.sum((jet_pt[:, 0] > 250) * (jet_pt[:, 1] > 250) * (jet_msd[:, 0] > 20) * (jet_msd[:, 1] > 20)) / len(jet_pt)

higgs_pt
jet_pt


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

weights["HH4V"] = evtDict["HH4V"]["totalWeight"] * LUMI * XSECHH4V * ACCEPTANCE
weights["QCD"] = evtDict["QCD"]["totalWeight"] * LUMI * QCD_MC_SF
# weights["QCD"] = evtDict["QCD"]["totalWeight"] * LUMI17 # * QCD_MC_SF
weights["tt"] = evtDict["tt"]["totalWeight"] * LUMI  # * DATA_MC_TF
weights["HHbbVV4q"] = evtDict["HHbbVV4q"]["totalWeight"] * LUMI * XSECHHBBVV4Q
weights["HH4b"] = evtDict["HH4b"]["totalWeight"] * LUMI  # * ACCEPTANCE
# weights["data"] = np.ones(len(data17))
# weights["fullHH4b"] = 0.0299 * np.ones(len(full_hh4b_samples['genWeight'])) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])
# weights["fullHH4b"] = ((full_hh4b_samples['genWeight'] > 0) * 2 - 1) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])

evtDict["QCD"]["weight"]


ak.sum(weights["QCD"])
ak.sum(weights["HH4b"])
ak.sum(weights["HHbbVV4q"])
ak.sum(weights["HH4V"])
ak.sum(weights["tt"])

ak.sum(weights["tt"]) + ak.sum(weights["QCD"])
ak.sum(weights["data"])

ak.sum(weights["HHbbVV4q"]) / ACCEPTANCE

h1W1q = (evtDict["HH4V"]['genHiggs1W1Decay'] == 1)
h1W2q = (evtDict["HH4V"]['genHiggs1W2Decay'] == 1)
h2W1q = (evtDict["HH4V"]['genHiggs2W1Decay'] == 1)
h2W2q = (evtDict["HH4V"]['genHiggs2W2Decay'] == 1)

h1W1l = (evtDict["HH4V"]['genHiggs1W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W1Decay'] < 5)
h1W2l = (evtDict["HH4V"]['genHiggs1W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W2Decay'] < 5)
h2W1l = (evtDict["HH4V"]['genHiggs2W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W1Decay'] < 5)
h2W2l = (evtDict["HH4V"]['genHiggs2W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W2Decay'] < 5)

hh8q = h1W1q * h1W2q * h2W1q * h2W2q
hh4q4l = (h1W1q * h1W2q * h2W1l * h2W2l) + \
        (h2W1q * h2W2q * h1W1l * h1W2l) + \
        (h1W1l * h1W2l * h2W1l * h2W2l) + \
        (h2W1q * h2W2q * h1W1l * h1W2l) + \
        (h2W1q * h2W2q * h1W1l * h1W2l) + \
        (h2W1q * h2W2q * h1W1l * h1W2l) + \

hh8l = h1W1l * h1W2l * h2W1l * h2W2l

ak.sum(hh8q) / len(evtDict["HH4V"]['weight'])
ak.sum(hh4q4l)
ak.sum(hh8l) / len(evtDict["HH4V"]['weight'])

tothh8q = ak.sum(hh8q)

hh8q


key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}


def get_vec(var, mass):
    dict1 = {}
    dict2 = {}
    for key, val in key_dict.items():
        dict1[key] = evtDict["HH4V"]['genHiggs1' + var + val]
        dict2[key] = evtDict["HH4V"]['genHiggs2' + var + val]
    dict1['mass'] = np.ones(len(evtDict["HH4V"]['genHiggs1' + var + val])) * mass
    dict2['mass'] = np.ones(len(evtDict["HH4V"]['genHiggs2' + var + val])) * mass

    return (ak.zip(dict1, with_name="PtEtaPhiMLorentzVector"), ak.zip(dict2, with_name="PtEtaPhiMLorentzVector"))


hvec = get_vec('', 125.1)
W1vec = get_vec('W1', 80.379)
W2vec = get_vec('W2', 80.379)
W1q1vec = get_vec('W1dau1', 0)
W1q2vec = get_vec('W1dau2', 0)
W2q1vec = get_vec('W2dau1', 0)
W2q2vec = get_vec('W2dau2', 0)


fatjet14q = (hvec[0].delta_r(W1q1vec[0]) < testr) * (hvec[0].delta_r(W1q2vec[0]) < testr) * (hvec[0].delta_r(W2q1vec[0]) < testr) * (hvec[0].delta_r(W2q2vec[0]) < testr)
fatjet24q = (hvec[1].delta_r(W1q1vec[1]) < testr) * (hvec[1].delta_r(W1q2vec[1]) < testr) * (hvec[1].delta_r(W2q1vec[1]) < testr) * (hvec[1].delta_r(W2q2vec[1]) < testr)
ak.sum(fatjet14q[hh8q]) / tothh8q
ak.sum(fatjet24q[hh8q]) / tothh8q

ak.sum(fatjet14q[hh8q] * fatjet24q[hh8q]) / tothh8q

evtDict["HH4V"]["fatJet14q"] = fatjet14q
evtDict["HH4V"]["fatJet24q"] = fatjet24q


ak.sum((hvec[0].delta_r(W1q1vec[0]) < testr) * (hvec[0].delta_r(W1q2vec[0]) < testr) * (hvec[0].delta_r(W2q1vec[0]) < testr) * (hvec[0].delta_r(W2q2vec[0]) < testr)) / totW


evtDict["HH4V"]["genHiggs2W2Decay"][hh8q]

plt.hist(ak.to_numpy(evtDict["HH4V"]["fatJet1DeepAK8MD_H4qvsQCD"][hh8q]), bins=np.linspace(0, 1, 101), histtype='step')
plt.hist(ak.to_numpy(evtDict["HH4V"]["fatJet2DeepAK8MD_H4qvsQCD"][hh8q]), bins=np.linspace(0, 1, 101), histtype='step')
evtDict["HH4V"]["fatJet2DeepAK8MD_H4qvsQCD"][hh8q]

ak.sum(h1W1l + h1W1q)

totW = ak.sum(h1W) + ak.sum(h2W)




















hists = {}

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


# fatJet1W = (evtDict['HHbbVV4q']['genHiggs1W1Pt'] > -99) * (evtDict['HHbbVV4q']['genHiggs1W2Pt'] > -99)
# fatJet2W = (evtDict['HHbbVV4q']['genHiggs2W1Pt'] > -99) * (evtDict['HHbbVV4q']['genHiggs2W2Pt'] > -99)
#
# fatJet1b = (evtDict['HHbbVV4q']['genHiggs1b1Pt'] > -99) * (evtDict['HHbbVV4q']['genHiggs1b2Pt'] > -99)
# fatJet2b = (evtDict['HHbbVV4q']['genHiggs2b1Pt'] > -99) * (evtDict['HHbbVV4q']['genHiggs2b2Pt'] > -99)
#
#
# evtDict['HHbbVV4q']['genHiggs1W1M'][fatJet1W]

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

# evtDict['HHbbVV4q'].fields

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
# H1W = evtDict['HHbbVV4q']['genHiggs1W1Decay'] == 0
# H2W = evtDict['HHbbVV4q']['genHiggs2W1Decay'] == 0
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


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True, pevtDict=None, useEDWeights=False, data=False, ak15=False):
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
                    if s == "HHbbVV4q" and fatJet and ak15: rnpf = "ak15" + npf
                    else: rnpf = npf
                    kwargs = {}
                    for j in range(nbins):
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
                if s == "HHbbVV4q" and fatJet and ak15: rnpf = "ak15" + npf
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


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True, hh4b=False, sepsig=True, stackall=False, log=False, lumilabel=None, data=False, rat_ylim=1, ak15=False):
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
        if not data: fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))
        else: fig, (axs, rax) = plt.subplots(2, nbins, figsize=(nbins * 9, 12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        if fatJet:
            if sepsig:
                rnpf = "ak15" + npf if ak15 else npf
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
                        axs[j].set_xlabel(None)
                        hist.plotratio(num=hists[name + var]['data'].project("sample", bpf + str(j + 1)).sum("sample"), denom=hists[name + var][regbg].project("sample", bpf + str(j + 1)).sum("sample"), ax=rax[j], error_opts=data_err_opts, unc='num')
                        rax[j].set_ylabel('data/MC')
                        rax[j].set_ylim(0, rat_ylim)
                        rax[j].grid()
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


init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=False)

tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


tot_events


scale_factor = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbVV4q': 1 / tot_events['HHbbVV4q'],
        'HH4b': 1 / tot_events['HH4b'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}



init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=True, hh4v=False)
plot_hists("Pt", "kin_tests", [40, 200, 2000], hh4v=False)


init_hists("MassSD", "Soft Drop Mass (GeV)", [50, 1, 400], name="data", scale=False)
fill_hists("MassSD", scale=False, data=True, name='data')
plot_hists("MassSD", "data_masssd_rat_pre_tf_qcd_only", [50, 1, 400], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)


init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], name="data", scale=False)
fill_hists("Pt", scale=False, data=True, name='data')
plot_hists("Pt", "data_pt_rat_post_tf_test", [40, 200, 2000], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=1.5, lumilabel=40)





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
bins = [[40, 200, 2000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True, ak15=True, hh4v=False)
plot_hists(vars, "jet_kin_ak15", bins, ak15=True, hh4v=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[40, 200, 2000], [50, 1, 400], [50, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True)
fill_hists(copy(vars), fatJet=True, scale=True, ak15=False, hh4v=False)
plot_hists(vars, "jet_kin", bins, ak15=False, hh4v=False)


hists

vars = {
    "HH4V": "PNetXbb_alt",
    "HHbbVV4q": "PNetXbb_alt",
    "HH4b": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
    "data": "PNetXbb",
}
disc_bin = [100, 0, 1]
init_hists("PNetXbb", "Particle Net Xbb", [100, 0, 1], fatJet=True)
fill_hists(vars, fatJet=True, scale=True, hh4v=False)
plot_hists(vars, "pnetxbb_new", [100, 0, 1], hh4v=False)



init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, hh4v=False)
plot_hists("DeepAK8MD_H4qvsQCD", "deepak8mdh4q", [100, 0, 1], hh4v=False)







# DeltaR study

h1W = (evtDict["HHbbVV4q"]['genHiggs1W1Decay'] == 1) * (evtDict["HHbbVV4q"]['genHiggs1W2Decay'] == 1)
h2W = (evtDict["HHbbVV4q"]['genHiggs2W1Decay'] == 1) * (evtDict["HHbbVV4q"]['genHiggs2W2Decay'] == 1)
totW = ak.sum(h1W) + ak.sum(h2W)
totW

len(evtDict["HHbbVV4q"]['genHiggs1W1Decay'])

weightsW = ak.flatten(ak.Array([weights["HHbbVV4q"][h1W], weights["HHbbVV4q"][h2W]]))


ak.sum(h1W * (evtDict["HHbbVV4q"]['genHiggs1W2Decay'] != 1))
ak.sum(~h1W * (evtDict["HHbbVV4q"]['genHiggs1W2Decay'] == 1))
ak.sum(h2W * (evtDict["HHbbVV4q"]['genHiggs2W2Decay'] != 1))
ak.sum(~h2W * (evtDict["HHbbVV4q"]['genHiggs2W2Decay'] == 1))


key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}


def get_vec(var, mass):
    dict = {}
    for key, val in key_dict.items():
        dict[key] = ak.flatten(ak.Array([evtDict["HHbbVV4q"]['genHiggs1' + var + val][h1W], evtDict["HHbbVV4q"]['genHiggs2' + var + val][h2W]]))
    dict['mass'] = np.ones(totW) * mass


    return ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


hvec = get_vec('', 125.1)
W1vec = get_vec('W1', 80.379)
W2vec = get_vec('W2', 80.379)
W1q1vec = get_vec('W1dau1', 0)
W1q2vec = get_vec('W1dau2', 0)
W2q1vec = get_vec('W2dau1', 0)
W2q2vec = get_vec('W2dau2', 0)

len(hvec)

testr = 0.8

ak.sum(hvec.delta_r(W1vec) < testr) / totW
ak.sum(hvec.delta_r(W2vec) < testr) / totW
ak.sum((hvec.delta_r(W1vec) < testr) * (hvec.delta_r(W2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr)) / totW
ak.sum((hvec.delta_r(W1q2vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q1vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr) * (hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)) / totW

ak8_4q = (hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr) * (hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)


testr = 1

ak.sum(hvec.delta_r(W1vec) < testr) / totW
ak.sum(hvec.delta_r(W2vec) < testr) / totW
ak.sum((hvec.delta_r(W1vec) < testr) * (hvec.delta_r(W2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr)) / totW
ak.sum((hvec.delta_r(W1q2vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q1vec) < testr)) / totW
ak.sum((hvec.delta_r(W2q2vec) < testr)) / totW

ak.sum((hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr) * (hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)) / totW
ak15_4q = (hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr) * (hvec.delta_r(W2q1vec) < testr) * (hvec.delta_r(W2q2vec) < testr)















events_bb_sorted = {}
fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD"]

for s, evts in evtDict.items():
    if s != "HH4V":
        pnet_key = "PNetXbb_alt" if s == "HHbbVV4q" else "PNetXbb"
        jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

        ak_dict =   {
                        "fatJet1PNetXbb": evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading,
                        "fatJet2PNetXbb": evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading,
                        "weight": weights[s],
                    }

        for var in fatjet_vars:
            ak_dict["fatJet1" + var] = evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading
            ak_dict["fatJet2" + var] = evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading
            if var != "DeepAK8MD_H4qvsQCD":
                ak_dict["ak15fatJet1" + var] = evts["ak15fatJet1" + var] * jet1_bb_leading + evts["ak15fatJet2" + var] * ~jet1_bb_leading
                ak_dict["ak15fatJet2" + var] = evts["ak15fatJet1" + var] * ~jet1_bb_leading + evts["ak15fatJet2" + var] * jet1_bb_leading

        events_bb_sorted[s] = ak.zip(ak_dict)


events_bb_sorted


# events_bb_sorted_ak15 = {}
# fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD"]
#
# for s, evts in evtDict.items():
#     if s != "HH4V":
#         pnet_key = "PNetXbb_alt" if s == "HHbbVV4q" else "PNetXbb"
#         jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]
#
#         ak_dict =   {
#                         "fatJet1PNetXbb": evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading,
#                         "fatJet2PNetXbb": evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading,
#                         "weight": weights[s],
#                     }
#
#         for var in fatjet_vars:
#             ak15 = "ak15" if s == 'HHbbVV4q' and var != "DeepAK8MD_H4qvsQCD" else ""
#             ak_dict["fatJet1" + var] = evts[ak15 + "fatJet1" + var] * jet1_bb_leading + evts[ak15 + "fatJet2" + var] * ~jet1_bb_leading
#             ak_dict["fatJet2" + var] = evts[ak15 + "fatJet1" + var] * ~jet1_bb_leading + evts[ak15 + "fatJet2" + var] * jet1_bb_leading
#
#         events_bb_sorted_ak15[s] = ak.zip(ak_dict)
#



var = "Pt"
varl = "$p_T$"
bins = [50, 250, 500]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted", pevtDict=events_bb_sorted, useEDWeights=True)
plot_hists(var, "jet_pt_bb_sorted_fine", bins, name="bb_sorted", sepsig=False)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins = [30, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted", pevtDict=events_bb_sorted, useEDWeights=True)
plot_hists(var, "jet_masssd_bb_sorted_fine", bins, name="bb_sorted", sepsig=False)


events_bb_sorted_ak15

var = "Pt"
varl = "$p_T$"
bins = [50, 250, 500]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted_ak15", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted_ak15", pevtDict=events_bb_sorted_ak15, useEDWeights=True)
plot_hists(var, "jet_pt_bb_sorted_fine_ak15", bins, name="bb_sorted_ak15", sepsig=False)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins = [30, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted_ak15", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted_ak15", pevtDict=events_bb_sorted_ak15, useEDWeights=True)
plot_hists(var, "jet_masssd_bb_sorted_fine_ak15", bins, name="bb_sorted_ak15", sepsig=False)


triggers17


events_bb_sorted_triggered = {}
fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD"]

for s, evts in evtDict.items():
    if s != "HH4V":
        triggered = evts[triggers17[0]]
        for i in range(1, len(triggers17)): triggered = triggered + evts[triggers17[i]]

        pnet_key = "PNetXbb_alt" if s == "HHbbVV4q" else "PNetXbb"
        jet1_bb_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

        ak_dict =   {
                        "fatJet1PNetXbb": (evts["fatJet1" + pnet_key] * jet1_bb_leading + evts["fatJet2" + pnet_key] * ~jet1_bb_leading)[triggered],
                        "fatJet2PNetXbb": (evts["fatJet1" + pnet_key] * ~jet1_bb_leading + evts["fatJet2" + pnet_key] * jet1_bb_leading)[triggered],
                        "weight": weights[s][triggered],
                    }

        for var in fatjet_vars:
            ak_dict["fatJet1" + var] = (evts["fatJet1" + var] * jet1_bb_leading + evts["fatJet2" + var] * ~jet1_bb_leading)[triggered]
            ak_dict["fatJet2" + var] = (evts["fatJet1" + var] * ~jet1_bb_leading + evts["fatJet2" + var] * jet1_bb_leading)[triggered]

        events_bb_sorted_triggered[s] = ak.zip(ak_dict)



events_bb_sorted_triggered



evtDict["HH4V"].fields

events_WW_sorted = {}
fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD", 'Tau2OverTau1', 'Tau4OverTau3', 'DeepAK8_H']

for s, evts in evtDict.items():
    print(s)
    if s != "HH4b" and s != "HHbbVV4q":
        if s == 'HH4V':
            pnet_key = "DeepAK8MD_H4qvsQCD"
            jet1_WW_leading = evts["fatJet1" + pnet_key][hh8q] > evts["fatJet2" + pnet_key][hh8q]

            ak_dict =   {
                            "weight": weights[s][hh8q],
                        }

            fatjet_vars.append("4q")

            for var in fatjet_vars:
                ak_dict["fatJet1" + var] = evts["fatJet1" + var][hh8q] * jet1_WW_leading + evts["fatJet2" + var][hh8q] * ~jet1_WW_leading
                ak_dict["fatJet2" + var] = evts["fatJet1" + var][hh8q] * ~jet1_WW_leading + evts["fatJet2" + var][hh8q] * jet1_WW_leading

            fatjet_vars = fatjet_vars[:-1]

            events_WW_sorted['HH4V8q'] = ak.zip(ak_dict)
        else:
            pnet_key = "DeepAK8MD_H4qvsQCD"
            jet1_WW_leading = evts["fatJet1" + pnet_key] > evts["fatJet2" + pnet_key]

            ak_dict =   {
                            "weight": weights[s],
                        }

            for var in fatjet_vars:
                ak_dict["fatJet1" + var] = evts["fatJet1" + var] * jet1_WW_leading + evts["fatJet2" + var] * ~jet1_WW_leading
                ak_dict["fatJet2" + var] = evts["fatJet1" + var] * ~jet1_WW_leading + evts["fatJet2" + var] * jet1_WW_leading

            events_WW_sorted[s] = ak.zip(ak_dict)


events_WW_sorted

ak.sum(events_WW_sorted["HH4V8q"]["fatJet14q"])
ak.sum(events_WW_sorted["HH4V8q"]["fatJet24q"])

# vars = ["Pt", "MassSD"]
# varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
# bins = [[50, 250, 600], [50, 30, 200]]
# init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
# fill_hists(copy(vars), fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
# plot_hists(vars, "jet_kin_bb_leading", bins, name="bb_sorted", hh4v=False, sepsig=False)
#
# init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD Score", [100, 0, 1], fatJet=True, name="bb_sorted", bbsorted=True)
# fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bb_sorted)
# plot_hists("DeepAK8MD_H4qvsQCD", "jet_deepak8mdh4q_bb_leading", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)
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


def cftable_func(cutflow, var_cuts, cut_labels=None, cut_idx=None):
    if cut_idx is None:
        i = 0
        cut_idx = []
        for brange in var_cuts.values():
            for j in brange:
                if j != 9999 and j != -9999:
                    cut_idx.append(i)
                i += 1

    return pandas.DataFrame(np.round(np.array(list(cutflow.values()))[:, cut_idx], 3), list(cutflow.keys()), cut_labels)






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















# kin cuts only

kin_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet2MassSD": [115, 135],
}

cutflow, events_bbs_kin_cuts = cutflow_func(kin_cuts)

var = {
    "HHbbVV4q": "PNetXbb_alt",
    "QCD": "PNetXbb",
    "tt": "PNetXbb",
}
var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [50, 0.9, 1]
init_hists("PNetXbb", varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_cut_pnet2", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=True, log=True)


init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="bb_sorted", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_cut_deepak8mdh4q_bb_leading2", [100, 0, 1], name="bb_sorted", hh4v=False, sepsig=False)


yields = {"sig": [], "bg": []}
for pnetcutoff in np.arange(0.99, 0.995, 0.001)[:-1]:
    sigy = []
    bgy = []
    for dakcutoff in np.arange(0.9, 1, 0.01):
        cuts = {}
        for s, evts in events_bbs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetXbb"] > pnetcutoff) * (evts["fatJet2DeepAK8MD_H4qvsQCD"] > dakcutoff)

        sig = np.sum(events_bbs_kin_cuts["HHbbVV4q"][cuts["HHbbVV4q"]].weight)
        bg = np.sum(events_bbs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_bbs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)



# kin + tagger cuts

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 300',
                'Jet2 pT > 300',
                'Jet1 MassSD > 75',
                'Jet1 MassSD < 150',
                'Jet1 MassSD > 50',
                'Jet1 MassSD < 150',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
            ]

cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
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


# 4b kin + bbVV tagger cuts

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 310',
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



# 4b kin + bbVV tagger cuts + mass2

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
    "fatJet2MassSD": [115, 135],
}

hhbbVV_cutflow = cutflow_func(var_cuts)[0]

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

cftable.to_csv('hhbbVV_cutflow.csv')



# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [105, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 145],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_cuts)[0]
msb2_cf




# AK15


# kin + bbVV tagger cuts - masssd

var_cuts = {
    "fatJet1Pt": [315, 9999],
    "fatJet2Pt": [320, 9999],
    # "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts, events_bb_sorted_ak15)

cut_labels = ['Jet1 pT > 315',
                'Jet2 pT > 320',
                # 'At least 1 jet pT > 350',
                'Jet1 MassSD > 110',
                'Jet1 MassSD < 140',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
            ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable



# 4b kin + bbVV tagger cuts + mass2

var_cuts = {
    "fatJet1Pt": [315, 9999],
    "fatJet2Pt": [320, 9999],
    # "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
    "fatJet2MassSD": [115, 135],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_bb_sorted_ak15)[0]

cut_labels = ['Jet1 pT > 315',
                'Jet2 pT > 320',
                # 'At least 1 jet pT > 350',
                'Jet1 MassSD > 110',
                'Jet1 MassSD < 140',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
                'Jet2 MassSD > 115',
                'Jet2 MassSD < 135',
            ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable

cftable.to_csv('hhbbVV_cutflow_ak15.csv')



# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [105, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 145],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_cuts)[0]
msb2_cf





# AK8


# kin + bbVV tagger cuts - masssd

var_cuts = {
    "fatJet1Pt": [315, 9999],
    "fatJet2Pt": [320, 9999],
    # "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts, events_bb_sorted)

cut_labels = ['Jet1 pT > 315',
                'Jet2 pT > 320',
                # 'At least 1 jet pT > 350',
                'Jet1 MassSD > 110',
                'Jet1 MassSD < 140',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
            ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable



# kin + bbVV tagger cuts + mass2

var_cuts = {
    "fatJet1Pt": [315, 9999],
    "fatJet2Pt": [320, 9999],
    # "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
    "fatJet2MassSD": [115, 135],
}

hhbbVV_cutflow = cutflow_func(var_cuts, events_bb_sorted)[0]

cut_labels = ['Jet1 pT > 315',
                'Jet2 pT > 320',
                # 'At least 1 jet pT > 350',
                'Jet1 MassSD > 110',
                'Jet1 MassSD < 140',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
                'Jet2 MassSD > 115',
                'Jet2 MassSD < 135',
            ]

del(hhbbVV_cutflow['HH4b'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable

cftable.to_csv('hhbbVV_cutflow.csv')



# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [105, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 145],
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
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
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
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
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





# bbVV kin + 4b tagger cuts

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
}

hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 300',
                'Jet2 pT > 300',
                'Jet1 MassSD > 75',
                'Jet1 MassSD < 150',
                'Jet1 MassSD > 50',
                'Jet1 MassSD < 150',
                'Jet1 PNetXbb > 0.99',
                'Jet2 PNetXbb > 0.94',
            ]

cftable = cftable_func(hhbbVV_cutflow, var_cuts, cut_labels)
cftable




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

init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("DeepAK8_H", "DeepAK8MD H", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("DeepAK8_H", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("DeepAK8_H", "jet_deepak8h_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau2OverTau1", "tau2/tau1", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau2OverTau1", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau2OverTau1", "jet_tau2otau1_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau4OverTau3", "tau4/tau3", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau4OverTau3", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau4OverTau3", "jet_tau4otau3_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)



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

init_hists("DeepAK8MD_H4qvsQCD", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="WW_sorted_cut", bbsorted=True)
fill_hists("DeepAK8MD_H4qvsQCD", fatJet=True, scale=True, name="WW_sorted_cut", pevtDict=events_WWs_kin_cuts, useEDWeights=True)
plot_hists("DeepAK8MD_H4qvsQCD", "jet_cut_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted_cut", hh4v=True, sepsig=False)

sig

hists['WW_sorted_cutDeepAK8MD_H4qvsQCD']['HH4V']


yields = {"sig": [], "bg": []}
for pnetcutoff1 in np.arange(0.8, 0.85, 0.01):
    sigy = []
    bgy = []
    for pnetcutoff2 in np.arange(0.8, 0.85, 0.01):
        cuts = {}
        for s, evts in events_WWs_kin_cuts.items():
            cuts[s] = (evts["fatJet1DeepAK8MD_H4qvsQCD"] > pnetcutoff1) * (evts["fatJet2DeepAK8MD_H4qvsQCD"] > pnetcutoff2)

        sig = np.sum(events_WWs_kin_cuts["HH4V"][cuts["HH4V"]].weight)
        bg = np.sum(events_WWs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_WWs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)


yields



#
# fatjet_vars = ["Pt", "MassSD"]
# var_names = ["pt", "msoftdrop"]
#
#
# pnetxbb = ak.pad_none(hh4b_samples[:]['FatJet_ParticleNetMD_probXbb'], 2, axis=1)[:, :2][~ak.is_none(hh4bfjpt[:, 1])]
# jet1_bb_leading = pnetxbb[:, 0] > pnetxbb[:, 1]
#
# ak_dict =   {
#                 "fatJet1PNetXbb": pnetxbb[:, 0] * jet1_bb_leading + pnetxbb[:, 1] * ~jet1_bb_leading,
#                 "fatJet2PNetXbb": pnetxbb[:, 0] * ~jet1_bb_leading + pnetxbb[:, 1] * jet1_bb_leading,
#                 "weight": np.ones(len(pnetxbb)) * RUN2LUMI * XSECHHBBBB / 372000,
#             }
#
# for var, name in zip(fatjet_vars, var_names):
#     evts = ak.pad_none(hh4b_samples[:]['FatJet_' + name], 2, axis=1)[:, :2][~ak.is_none(hh4bfjpt[:, 1])]
#     ak_dict["fatJet1" + var] = evts[:, 0] * jet1_bb_leading + evts[:, 1] * ~jet1_bb_leading
#     ak_dict["fatJet2" + var] = evts[:, 0] * ~jet1_bb_leading + evts[:, 1] * jet1_bb_leading
#
# ak_dict
# hh4b_events_bb_sorted = ak.zip(ak_dict)
# hh4b_events_bb_sorted
#


hh4b_kin_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    # "fatJet2MassSD": [50, 150],
    # "fatJet1PNetXbb": [0.985, 9999],
    # "fatJet1PNetXbb": [0.985, 9999],
    # "fatJet2DeepAK8MD_H4qvsQCD": [0.9, 9999],
}


hh4b_kincf, events_hh4b_kin_cuts = cutflow_func(hh4b_kin_var_cuts)
cftable = cftable_func(hh4b_kincf, hh4b_kin_var_cuts)
cftable


#
# s = 'hh4b'
# cuts = []
# for var, brange in hh4b_kin_var_cuts.items():
#     cuts.append(hh4b_events_bb_sorted[var] > brange[0])
#     cuts.append(hh4b_events_bb_sorted[var] < brange[1])
#
# cuts.append((hh4b_events_bb_sorted['fatJet1Pt'] > 350) + (hh4b_events_bb_sorted['fatJet2Pt'] > 350))
# cut = cuts[0]
# # hh4bcf[s] = []
# for i in np.arange(1, len(cuts)):
#     # hh4bcf[s].append(np.sum(evts[cut].weight))
#     cut = cut * cuts[i]
#
# # hh4bcf[s].append(np.sum(evts[cut].weight))
#
# events_hh4b_kin_cuts[s] = hh4b_events_bb_sorted[cut]
#
# hh4b_kincf






events_bbs_kin_cuts

var = "tagger2d_cut"
hists[var] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet1bb", r"Jet1 Tagger Score", 20, 0.9, 1),
                                hist.Bin("jet2WW", r"Jet2 Tagger Score", 20, 0.9, 1),
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


# s = 'hh4b'
# hists[var].fill(sample='HH4b',
#                  jet1bb = ak.to_numpy(events_hh4b_cuts[s]["fatJet1PNetXbb"]),
#                  jet2WW = ak.to_numpy(events_hh4b_cuts[s]["fatJet2PNetXbb"]),
#                  weight = ak.to_numpy(events_hh4b_cuts[s]["weight"]),
#                 )

patch_opts = {
    'cmap': 'jet',
    'vmin': 0,
    'vmax': 100,
}

# titles = ['HH4b', 'HHbbVV4q', 'QCD', 'tt']
titles = [['QCD bb1 vs WW2 tagger', 'QCD bb1 vs bb2 tagger'], ['HHbbVV4q', 'tt']]
fig, axs = plt.subplots(2, 2, figsize=(2 * 9, 2 * 9))
for i in range(2):
    if i == 1: del patch_opts['vmax']
    for j in range(2):
        hist.plot2d(hists[var][titles[i][j]].sum('sample'), 'jet1bb', ax=axs[i][j], patch_opts=patch_opts)
        axs[i][j].set_title(titles[i][j])

plt.tight_layout(pad=2)
plt.savefig("figs/tagger_2d_cut.pdf", bbox_inches='tight')
plt.show()



















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








# Tagger Plots/ROC Curves
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
