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
    "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWW4q": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
    "HH4b": "data/weighted/GluGluToHHTo4B_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8_1pb_weighted.root"
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")

evtDict["HH4V"] = uproot4.concatenate("data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root")
evtDict["HHbbWW4q"] = uproot4.concatenate("data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root")

data17 = uproot4.concatenate('data/weighted/JetHT*.root:tree')
evtDict["data"] = data17

evtDict["data"].fields

ak.sum((data17['fatJet1Pt'] > 250) * (data17['fatJet2Pt'] > 250) * (data17['fatJet1MassSD'] > 20) * (data17['fatJet2MassSD'] > 20))

full_hh4b_samples = uproot4.concatenate("data/hh4b/nano_*.root:Events")

evtDict["HHbbWW4q"].fields

data17.fields


evtDict["HH4V"].fields


evtDict["HH4V"]['genHiggs1W2Decay']


triggers = data17['HLT_PFJet500'] + \
data17['HLT_AK8PFJet500'] + \
data17['HLT_AK8PFJet360_TrimMass30'] + \
data17['HLT_AK8PFJet380_TrimMass30'] + \
data17['HLT_AK8PFJet400_TrimMass30'] + \
data17['HLT_AK8PFHT800_TrimMass50'] + \
data17['HLT_AK8PFJet330_PFAK8BTagCSV_p17']


ak.sum(triggers)

# DeltaR study

evtDict["HHbbWW4q"].fields

np.unique(ak.to_numpy(evtDict["HHbbWW4q"]['genHiggs1W1Decay']))

ak.sum(evtDict["HHbbWW4q"]['genHiggs1W1Decay'] > 1)
ak.sum(evtDict["HHbbWW4q"]['genHiggs2W1Decay'] > 1)
ak.sum(evtDict["HHbbWW4q"]['genHiggs2W1Decay'] == 5)
ak.sum(evtDict["HHbbWW4q"]['genHiggs1W1Decay'] == 5)


h1W = (evtDict["HHbbWW4q"]['genHiggs1W1Decay'] == 1) * (evtDict["HHbbWW4q"]['genHiggs1W2Decay'] == 1)
h2W = (evtDict["HHbbWW4q"]['genHiggs2W1Decay'] == 1) * (evtDict["HHbbWW4q"]['genHiggs2W2Decay'] == 1)
totW = ak.sum(h1W) + ak.sum(h2W)
totW

len(evtDict["HHbbWW4q"]['genHiggs1W1Decay'])

weightsW = ak.flatten(ak.Array([weights["HHbbWW4q"][h1W], weights["HHbbWW4q"][h2W]]))


ak.sum(h1W * (evtDict["HHbbWW4q"]['genHiggs1W2Decay'] != 1))
ak.sum(~h1W * (evtDict["HHbbWW4q"]['genHiggs1W2Decay'] == 1))
ak.sum(h2W * (evtDict["HHbbWW4q"]['genHiggs2W2Decay'] != 1))
ak.sum(~h2W * (evtDict["HHbbWW4q"]['genHiggs2W2Decay'] == 1))


key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}

def get_vec(var, mass):
    dict = {}
    for key, val in key_dict.items():
        dict[key] = ak.flatten(ak.Array([evtDict["HHbbWW4q"]['genHiggs1' + var + val][h1W], evtDict["HHbbWW4q"]['genHiggs2' + var + val][h2W]]))
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


testr = 1.5

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



data_err_opts_c = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}


hists['deltar'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total HWW4q Events'
str_8 = '0.8 $\Delta$R Captures All 4q'
str_15 = '1.5 $\Delta$R Captures All 4q'

hists['deltar'].fill(sample=str_all, h1=hvec.pt, weight=weightsW)
hists['deltar'].fill(sample=str_8, h1=hvec.pt[ak8_4q], weight=weightsW[ak8_4q])
hists['deltar'].fill(sample=str_15, h1=hvec.pt[ak15_4q], weight=weightsW[ak15_4q])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])



fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltar'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltar'][str_8].sum('sample'), hists['deltar'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltar'][str_15].sum('sample'), hists['deltar'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltar.pdf', bbox_inches='tight')
plt.show()









# delta r Hbb study

evtDict["HHbbWW4q"]['genHiggs1b1Pt']
evtDict["HHbbWW4q"]['genHiggs2b1Pt']

h1b = (evtDict["HHbbWW4q"]['genHiggs1b1Pt'] != -99) * (evtDict["HHbbWW4q"]['genHiggs1b2Pt'] != -99)
h2b = (evtDict["HHbbWW4q"]['genHiggs2b1Pt'] != -99) * (evtDict["HHbbWW4q"]['genHiggs2b2Pt'] != -99)
totb = ak.sum(h1b) + ak.sum(h2b)
totb
weightsb = ak.flatten(ak.Array([weights["HHbbWW4q"][h1b], weights["HHbbWW4q"][h2b]]))


key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}


def get_vec_bb(var, mass):
    dict = {}
    for key, val in key_dict.items():
        dict[key] = ak.flatten(ak.Array([evtDict["HHbbWW4q"]['genHiggs1' + var + val][h1b], evtDict["HHbbWW4q"]['genHiggs2' + var + val][h2b]]))
    dict['mass'] = np.ones(totb) * mass

    return ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


hbbvec = get_vec_bb('', 125.1)
b1vec = get_vec_bb('b1', 4.18)
b2vec = get_vec_bb('b2', 4.18)

testr = 0.8

ak.sum(hbbvec.delta_r(b1vec) < testr) / totb
ak.sum(hbbvec.delta_r(b2vec) < testr) / totb
ak.sum((hbbvec.delta_r(b1vec) < testr) * (hbbvec.delta_r(b2vec) < testr)) / totb

ak8_2b = (hbbvec.delta_r(b1vec) < testr) * (hbbvec.delta_r(b2vec) < testr)


testr = 1.5

ak.sum(hbbvec.delta_r(b1vec) < testr) / totb
ak.sum(hbbvec.delta_r(b2vec) < testr) / totb
ak.sum((hbbvec.delta_r(b1vec) < testr) * (hbbvec.delta_r(b2vec) < testr)) / totb

ak15_2b = (hbbvec.delta_r(b1vec) < testr) * (hbbvec.delta_r(b2vec) < testr)


hists['deltarbb'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total Hbb Events'
str_8 = '0.8 $\Delta$R Captures Both bb'
str_15 = '1.5 $\Delta$R Captures Both bb'

hists['deltarbb'].fill(sample=str_all, h1=hbbvec.pt, weight=weightsb)
hists['deltarbb'].fill(sample=str_8, h1=hbbvec.pt[ak8_2b], weight=weightsb[ak8_2b])
hists['deltarbb'].fill(sample=str_15, h1=hbbvec.pt[ak15_2b], weight=weightsb[ak15_2b])


fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltarbb'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})


data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltarbb'][str_8].sum('sample'), hists['deltarbb'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltarbb'][str_15].sum('sample'), hists['deltarbb'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltarbb.pdf', bbox_inches='tight')
plt.show()







# HH4V delta r

h1W1 = (evtDict["HH4V"]['genHiggs1W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W1Decay'] < 5)
h1W2 = (evtDict["HH4V"]['genHiggs1W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W2Decay'] < 5)
h2W1 = (evtDict["HH4V"]['genHiggs2W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W1Decay'] < 5)
h2W2 = (evtDict["HH4V"]['genHiggs2W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W2Decay'] < 5)

h1W = h1W1 * h1W2
h2W = h2W1 * h2W2

totW = ak.sum(h1W) + ak.sum(h2W)

weightsW = np.ones(totW) * weights["HH4V"][0]

key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}

def get_vec(var, mass):
    dict = {}
    for key, val in key_dict.items():
        dict[key] = ak.flatten(ak.Array([evtDict["HH4V"]['genHiggs1' + var + val][h1W], evtDict["HH4V"]['genHiggs2' + var + val][h2W]]))
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


testr = 1.5

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



data_err_opts_c = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}


hists['deltarl'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total HVV4l Events'
str_8 = '0.8 $\Delta$R Captures All 4l'
str_15 = '1.5 $\Delta$R Captures All 4l'

hists['deltarl'].fill(sample=str_all, h1=hvec.pt, weight=weightsW)
hists['deltarl'].fill(sample=str_8, h1=hvec.pt[ak8_4q], weight=weightsW[ak8_4q])
hists['deltarl'].fill(sample=str_15, h1=hvec.pt[ak15_4q], weight=weightsW[ak15_4q])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])



fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltarl'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltarl'][str_8].sum('sample'), hists['deltarl'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltarl'][str_15].sum('sample'), hists['deltarl'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltarl.pdf', bbox_inches='tight')
plt.show()







# HH4V delta r

h1W1l = (evtDict["HH4V"]['genHiggs1W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W1Decay'] < 5)
h1W2l = (evtDict["HH4V"]['genHiggs1W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W2Decay'] < 5)
h2W1l = (evtDict["HH4V"]['genHiggs2W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W1Decay'] < 5)
h2W2l = (evtDict["HH4V"]['genHiggs2W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W2Decay'] < 5)

h1W1q = (evtDict["HH4V"]['genHiggs1W1Decay'] == 1)
h1W2q = (evtDict["HH4V"]['genHiggs1W2Decay'] == 1)
h2W1q = (evtDict["HH4V"]['genHiggs2W1Decay'] == 1)
h2W2q = (evtDict["HH4V"]['genHiggs2W2Decay'] == 1)

h1W = h1W1l * h1W2q + h1W1q * h1W2l
h2W = h2W1l * h2W2q + h2W1q * h2W2l

totW = ak.sum(h1W) + ak.sum(h2W)
totW
weightsW = np.ones(totW) * weights["HH4V"][0]

key_dict = {'pt': 'Pt', 'eta': 'Eta', 'phi': 'Phi'}

def get_vec(var, mass):
    dict = {}
    for key, val in key_dict.items():
        dict[key] = ak.flatten(ak.Array([evtDict["HH4V"]['genHiggs1' + var + val][h1W], evtDict["HH4V"]['genHiggs2' + var + val][h2W]]))
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


testr = 1.5

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



data_err_opts_c = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}


hists['deltarql'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total HVV2ql$\\nu$ Events'
str_8 = '0.8 $\Delta$R Captures Both 2q'
str_15 = '1.5 $\Delta$R Captures Both 2q'

hists['deltarql'].fill(sample=str_all, h1=hvec.pt, weight=weightsW)
hists['deltarql'].fill(sample=str_8, h1=hvec.pt[ak8_4q], weight=weightsW[ak8_4q])
hists['deltarql'].fill(sample=str_15, h1=hvec.pt[ak15_4q], weight=weightsW[ak15_4q])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])



fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltarql'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltarql'][str_8].sum('sample'), hists['deltarql'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltarql'][str_15].sum('sample'), hists['deltarql'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltarql.pdf', bbox_inches='tight')
plt.show()










# acceptance calculation

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

jet_pt = ak.pad_none(full_hh4b_samples['FatJet_pt'], 2, axis=1)[:, :2]
jet_msd = ak.pad_none(full_hh4b_samples['FatJet_msoftdrop'], 2, axis=1)[:, :2]

ak.sum((jet_pt[:, 0] > 250) * (jet_pt[:, 1] > 250) * (jet_msd[:, 0] > 20) * (jet_msd[:, 1] > 20) )

fhiggs_pt
ntuple_cut_higgs_pt =





# Weights, all values in pb
RUN2LUMI = 137000
LUMI17 = 40000
XSECHHBBWWQQ = 31.05e-3 * 0.5807 * 0.2154 * 0.68**2 * 2
XSECHH4V = 31.05e-3 * (0.2154 + 0.2643)
XSECHHBBBB = 31.05e-3 * 0.5807**2
ACCEPTANCE = 0.08239891267414204
HHBBWW4Q_NEVENTS = 185417
HH4B_NEVENTS = 807198
DATA_MC_SF = 0.7709694064237897
QCD_MC_SF = 0.7677139264268172


weights = {}

weights["HH4V"] = evtDict["HH4V"]["weight"] * RUN2LUMI * XSECHH4V * ACCEPTANCE
weights["QCD"] = evtDict["QCD"]["totalWeight"] * RUN2LUMI * QCD_MC_SF
# weights["QCD"] = evtDict["QCD"]["totalWeight"] * LUMI17 # * QCD_MC_SF
weights["tt"] = evtDict["tt"]["totalWeight"] * RUN2LUMI # * DATA_MC_TF
weights["HHbbWW4q"] = evtDict["HHbbWW4q"]["weight"] * RUN2LUMI * XSECHHBBWWQQ * ACCEPTANCE
weights["HH4b"] = evtDict["HH4b"]["totalWeight"] * RUN2LUMI  # * ACCEPTANCE
weights["data"] = np.ones(len(data17))
# weights["fullHH4b"] = 0.0299 * np.ones(len(full_hh4b_samples['genWeight'])) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])
weights["fullHH4b"] = ((full_hh4b_samples['genWeight'] > 0) * 2 - 1) * RUN2LUMI * XSECHHBBBB / len(full_hh4b_samples['genWeight'])

ak.sum(weights["QCD"])
ak.sum(weights["HH4b"])
ak.sum(weights["HHbbWW4q"])
ak.sum(weights["tt"])







hists['genHiggsPt'] = hist.Hist("Events (Unit Normalized)",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Leading Higgs $p_T$ ", 51, 0, 1000),
                             )

hists['genHiggsPt'].fill(sample='HH4b', h1=evtDict['HH4b']['genHiggs1Pt'], weight=weights['HH4b'])
hists['genHiggsPt'].fill(sample='HHbbWW4q', h1=evtDict['HHbbWW4q']['genHiggs1Pt'], weight=weights['HHbbWW4q'])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])

hists['genHiggsPt'].scale(scale_factor, axis='sample')

hist.plot1d(hists['genHiggsPt'], line_opts=line_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.savefig('figs/hhbbWW_vs_hh4b_norm.pdf', bbox_inches='tight')
plt.show()


fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['genHiggsPt'], line_opts=line_opts, ax=axs, density=True)
# axs.set_yscale('log')
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
hist.plotratio(hists['genHiggsPt']['HHbbWW4q'].sum('sample'), hists['genHiggsPt']['HH4b'].sum('sample'), ax=rax, error_opts=data_err_opts, unc='num')
rax.set_ylabel("bbWW/4b")
rax.set_ylim(0, 2)
plt.savefig('figs/hhbbWW_vs_hh4b_norm.pdf', bbox_inches='tight')
plt.show()




hists['genHiggsMass'] = hist.Hist("Events (Unit Normalized)",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"HH Mass", 51, 200, 2500),
                             )


hists['genHiggsMass'].fill(sample='HH4b', h1=evtDict['HH4b']['genHH_mass'], weight=weights['HH4b'])
hists['genHiggsMass'].fill(sample='HHbbWW4q', h1=evtDict['HHbbWW4q']['genHH_mass'], weight=weights['HHbbWW4q'])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])
hists['genHiggsMass'].scale(scale_factor, axis='sample')

hist.plot1d(hists['genHiggsMass'], line_opts=line_opts)
plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.savefig('figs/hhbbWW_vs_hh4b_mass_norm.pdf', bbox_inches='tight')
plt.show()


fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['genHiggsMass'], line_opts=line_opts, ax=axs, density=True)
# axs.set_yscale('log')
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
hist.plotratio(hists['genHiggsMass']['HHbbWW4q'].sum('sample'), hists['genHiggsMass']['HH4b'].sum('sample'), ax=rax, error_opts=data_err_opts, unc='num')
rax.set_ylabel("bbWW/4b")
rax.set_ylim(0, 5)
plt.savefig('figs/hhbbWW_vs_hh4b_mass_norm.pdf', bbox_inches='tight')
plt.show()


tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


tot_events


scale_factor = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbWW4q': 1 / tot_events['HHbbWW4q'],
        'HH4b': 1 / tot_events['HH4b'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}





# fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
# hist.plot1d(hists['genHiggsPt'], line_opts=line_opts, ax=axs)
# # axs.set_yscale('log')
# axs.set_xlabel(None)
# axs.set_ylim(0, 5)
# axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
# hist.plotratio(hists['genHiggsPt']['HH4b gen-cut'].sum('sample'), hists['genHiggsPt']['HH4b full'].sum('sample'), ax=rax, error_opts=data_err_opts, unc='num')
# rax.set_ylabel("gen-cut/full")
# plt.savefig('figs/gen_cut_vs_full_no_acceptance.pdf', bbox_inches='tight')
# plt.show()





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

fatJet1 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet1Pt'], "eta": evtDict['HHbbWW4q']['fatJet1Eta'], "phi": evtDict['HHbbWW4q']['fatJet1Phi'], "mass": evtDict['HHbbWW4q']['fatJet1Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet2 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet2Pt'], "eta": evtDict['HHbbWW4q']['fatJet2Eta'], "phi": evtDict['HHbbWW4q']['fatJet2Phi'], "mass": evtDict['HHbbWW4q']['fatJet2Mass']}, with_name="PtEtaPhiMLorentzVector")
fatJet3 = ak.zip({"pt": evtDict['HHbbWW4q']['fatJet3Pt'], "eta": evtDict['HHbbWW4q']['fatJet3Eta'], "phi": evtDict['HHbbWW4q']['fatJet3Phi'], "mass": evtDict['HHbbWW4q']['fatJet3Mass']}, with_name="PtEtaPhiMLorentzVector")



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

fatJet1H1 = fatJet1.delta_r(genHiggs1) < testr
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


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True, pevtDict=None, useEDWeights=False, data=False):
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
                if (fatJet or "HH" in s) and (s != "HH4V" or hh4v) and (s != 'data' or data):
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

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}

regbg = re.compile('^(QCD|tt)$')
regsig = re.compile('^(HH4V|HHbbWW4q)$')
regnotdata = re.compile('^(?!data).*')
regmc = re.compile('^(QCD|tt|HHbbWW4q)$')


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True, hh4b=False, sepsig=True, stackall=False, log=False, lumilabel=None, data=False, rat_ylim=1):
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
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        if hh4v: hist.plot1d(hists[name + var]['hh4v'].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)

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

        var = vars["HHbbWW4q"]
        if not data: fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))
        else: fig, (axs, rax) = plt.subplots(2, nbins, figsize=(nbins * 9, 12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        if fatJet:
            if sepsig:
                _ = axs[0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[0].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(1, var)] * fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbWW4q - Hbb', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)
                _ = axs[1].hist(evtDict["HHbbWW4q"]["fatJet{}{}".format(2, var)] * fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbWW4q - HWW', weights=weights["HHbbWW4q"] * scale_factor["HHbbWW4q"], histtype='step', **line_opts)

                for j in range(nbins):
                    ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                    axs[j].set_prop_cycle(cycler(color=colors_cycle))
                    if hh4v: hist.plot1d(hists[name + var]['HH4V'].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)

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
                        hist.plot1d(hists[name + var][regnotdata].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbWW4q', 'tt', 'QCD'])
                    else:
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][sig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
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
        'HHbbWW4q': 1 / tot_events['HHbbWW4q'],
        'HH4b': 1 / tot_events['HH4b'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}




init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=True)
plot_hists("Pt", "kin_tests", [40, 200, 2000])


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
# vars = ["Pt", "Mass", "MassSD"]
# varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
# bins = [[50, 250, 500], [30, 50, 200], [30, 50, 200]]
# init_hists(vars, varsl, bins, fatJet=True, name="fine")
# fill_hists(copy(vars), fatJet=True, scale=True, name="fine")
# plot_hists(vars, "jet_kin_fine", bins, name="fine")
#

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





















events_bb_sorted = {}
fatjet_vars = ["Pt", "MassSD", "DeepAK8MD_H4qvsQCD"]

for s, evts in evtDict.items():
    if s != "HH4V":
        pnet_key = "PNetXbb_alt" if s == "HHbbWW4q" else "PNetXbb"
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


events_bb_sorted

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

    return pandas.DataFrame(np.round(np.array(list(cutflow.values()))[:, cut_idx], 2), list(cutflow.keys()), cut_labels)






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

cutflow

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

        sig = np.sum(events_bbs_kin_cuts["HHbbWW4q"][cuts["HHbbWW4q"]].weight)
        bg = np.sum(events_bbs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_bbs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)



yields['sig']
yields['']



# kin + tagger cuts

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbWW_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 300',
                'Jet2 pT > 300',
                'Jet1 MassSD > 75',
                'Jet1 MassSD < 150',
                'Jet1 MassSD > 50',
                'Jet1 MassSD < 150',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
            ]

cftable = cftable_func(hhbbWW_cutflow, var_cuts, cut_labels)
cftable


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [10, 50, 150]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_cut", pevtDict=events_bbs_cuts, useEDWeights=True)
plot_hists(var, "jet_cuts_masssd", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, lumilabel=True, stackall=True)


init_hists("MassSD", "Soft Drop Mass (GeV)", scale=False, [50, 1, 400], name="data", scale=False)
fill_hists("MassSD", scale=False, data=True, name='data')
plot_hists("MassSD", "data_masssd_rat_pre_tf", [50, 1, 400], data=True, name='data', stackall=True, sepsig=False, log=True, rat_ylim=0.5)


# 4b kin + bbWW tagger cuts

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
}

hhbbWW_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.94',
            ]

del(hhbbWW_cutflow['HH4b'])
cftable = cftable_func(hhbbWW_cutflow, var_cuts, cut_labels)
cftable

cftable.to_csv('hhbbWW_cutflow.csv')


# 4b kin + bbWW tagger cuts + mass2

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.94, 9999],
    "fatJet2MassSD": [115, 135],
}

hhbbWW_cutflow = cutflow_func(var_cuts)[0]

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

del(hhbbWW_cutflow['HH4b'])
cftable = cftable_func(hhbbWW_cutflow, var_cuts, cut_labels)
cftable

cftable.to_csv('hhbbWW_cutflow.csv')



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




# bbWW kin + 4b tagger cuts

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
}

hhbbWW_cutflow, events_bbs_cuts = cutflow_func(var_cuts)

cut_labels = ['Jet1 pT > 300',
                'Jet2 pT > 300',
                'Jet1 MassSD > 75',
                'Jet1 MassSD < 150',
                'Jet1 MassSD > 50',
                'Jet1 MassSD < 150',
                'Jet1 PNetXbb > 0.99',
                'Jet2 PNetXbb > 0.94',
            ]

cftable = cftable_func(hhbbWW_cutflow, var_cuts, cut_labels)
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

del(hh4bcf['HHbbWW4q'])

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

del(hh4bm2cf['HHbbWW4q'])

cftable = cftable_func(hh4bm2cf, hh4b_mass2_cuts, cut_labels)
cftable

cftable.to_csv('hh4b_cutflow_3.csv')


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





hh4b_cuts_evts['HH4b']









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

# titles = ['HH4b', 'HHbbWW4q', 'QCD', 'tt']
titles = [['QCD bb1 vs WW2 tagger', 'QCD bb1 vs bb2 tagger'], ['HHbbWW4q', 'tt']]
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
