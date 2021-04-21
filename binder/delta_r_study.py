import uproot4
import awkward1 as ak

# %matplotlib inline
import matplotlib.pyplot as plt
import mplhep as hep

import coffea.hist as hist

import numpy as np

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


hists = {}


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

data_err_opts_c = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}


# DeltaR study

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






# HH4V 2l2nu delta r

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






# HH4V 2q l\nu delta r

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

ak.sum((hvec.delta_r(W1q1vec) < testr) * (hvec.delta_r(W1q2vec) < testr) * (hvec.delta_r(W2q1vec) < testr)) / totW
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


hists['deltarlqboth'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total HVV2ql$\\nu$ Events'
str_8 = '0.8 $\Delta$R Captures All 4'
str_15 = '1.5 $\Delta$R Captures All 4'

hists['deltarlqboth'].fill(sample=str_all, h1=hvec.pt, weight=weightsW)
hists['deltarlqboth'].fill(sample=str_8, h1=hvec.pt[ak8_4q], weight=weightsW[ak8_4q])
hists['deltarlqboth'].fill(sample=str_15, h1=hvec.pt[ak15_4q], weight=weightsW[ak15_4q])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])



fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltarlqboth'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltarlqboth'][str_8].sum('sample'), hists['deltarlqboth'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltarlqboth'][str_15].sum('sample'), hists['deltarlqboth'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltarlqboth.pdf', bbox_inches='tight')
plt.show()








# HH4V 2q l\nu capturing 2q delta r

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

dict = {}
for key, val in key_dict.items():
    dict[key] = ak.flatten(ak.Array([
                                    (evtDict["HH4V"]['genHiggs1W1dau1' + val] * h1W1q + evtDict["HH4V"]['genHiggs1W2dau1' + val] * h1W2q)[h1W],
                                    (evtDict["HH4V"]['genHiggs2W1dau1' + val] * h2W1q + evtDict["HH4V"]['genHiggs2W2dau1' + val] * h2W2q)[h2W]]))
dict['mass'] = np.zeros(totW)

Wq1vec = ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


dict = {}
for key, val in key_dict.items():
    dict[key] = ak.flatten(ak.Array([
                                    (evtDict["HH4V"]['genHiggs1W1dau2' + val] * h1W1q + evtDict["HH4V"]['genHiggs1W2dau2' + val] * h1W2q)[h1W],
                                    (evtDict["HH4V"]['genHiggs2W1dau2' + val] * h2W1q + evtDict["HH4V"]['genHiggs2W2dau2' + val] * h2W2q)[h2W]]))
dict['mass'] = np.zeros(totW)

Wq2vec = ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


testr = 0.8

ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)) / totW
ak8_4q = (hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)


testr = 1.5

ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)) / totW
ak15_4q = (hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)


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







# HH4V  2q l\nu capturing 2q l delta r

h1W1l = (evtDict["HH4V"]['genHiggs1W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W1Decay'] < 4)
h1W2l = (evtDict["HH4V"]['genHiggs1W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs1W2Decay'] < 4)
h2W1l = (evtDict["HH4V"]['genHiggs2W1Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W1Decay'] < 4)
h2W2l = (evtDict["HH4V"]['genHiggs2W2Decay'] > 1) * (evtDict["HH4V"]['genHiggs2W2Decay'] < 4)

h1W1q = (evtDict["HH4V"]['genHiggs1W1Decay'] == 1)
h1W2q = (evtDict["HH4V"]['genHiggs1W2Decay'] == 1)
h2W1q = (evtDict["HH4V"]['genHiggs2W1Decay'] == 1)
h2W2q = (evtDict["HH4V"]['genHiggs2W2Decay'] == 1)

h1W = h1W1l * h1W2q + h1W1q * h1W2l
h2W = h2W1l * h2W2q + h2W1q * h2W2l

h1Wq1l = ((evtDict["HH4V"]['genHiggs1W1dau1Id'] == 11) + (evtDict["HH4V"]['genHiggs1W1dau1Id'] == 13)) * (h1W1l * h1W2q) + \
            ((evtDict["HH4V"]['genHiggs1W2dau1Id'] == 11) + (evtDict["HH4V"]['genHiggs1W2dau1Id'] == 13)) * (h1W1q * h1W2l)

h1Wq2l = ((evtDict["HH4V"]['genHiggs1W1dau2Id'] == 11) + (evtDict["HH4V"]['genHiggs1W1dau2Id'] == 13)) * (h1W1l * h1W2q) + \
            ((evtDict["HH4V"]['genHiggs1W2dau2Id'] == 11) + (evtDict["HH4V"]['genHiggs1W2dau2Id'] == 13)) * (h1W1q * h1W2l)

h2Wq1l = ((evtDict["HH4V"]['genHiggs2W1dau1Id'] == 11) + (evtDict["HH4V"]['genHiggs2W1dau1Id'] == 13)) * (h2W1l * h2W2q) + \
            ((evtDict["HH4V"]['genHiggs2W2dau1Id'] == 11) + (evtDict["HH4V"]['genHiggs2W2dau1Id'] == 13)) * (h2W1q * h2W2l)

h2Wq2l = ((evtDict["HH4V"]['genHiggs2W1dau2Id'] == 11) + (evtDict["HH4V"]['genHiggs2W1dau2Id'] == 13)) * (h2W1l * h2W2q) + \
            ((evtDict["HH4V"]['genHiggs2W2dau2Id'] == 11) + (evtDict["HH4V"]['genHiggs2W2dau2Id'] == 13)) * (h2W1q * h2W2l)


ak.sum(evtDict["HH4V"]['genHiggs1W1dau1Id'] == 11)
ak.sum(evtDict["HH4V"]['genHiggs1W1dau1Id'] == 13)

ak.sum(h1Wq1l)



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

dict = {}
for key, val in key_dict.items():
    dict[key] = ak.flatten(ak.Array([
                                    (evtDict["HH4V"]['genHiggs1W1dau1' + val] * h1W1q + evtDict["HH4V"]['genHiggs1W2dau1' + val] * h1W2q)[h1W],
                                    (evtDict["HH4V"]['genHiggs2W1dau1' + val] * h2W1q + evtDict["HH4V"]['genHiggs2W2dau1' + val] * h2W2q)[h2W]]))
dict['mass'] = np.zeros(totW)

Wq1vec = ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


dict = {}
for key, val in key_dict.items():
    dict[key] = ak.flatten(ak.Array([
                                    (evtDict["HH4V"]['genHiggs1W1dau2' + val] * h1W1q + evtDict["HH4V"]['genHiggs1W2dau2' + val] * h1W2q)[h1W],
                                    (evtDict["HH4V"]['genHiggs2W1dau2' + val] * h2W1q + evtDict["HH4V"]['genHiggs2W2dau2' + val] * h2W2q)[h2W]]))
dict['mass'] = np.zeros(totW)

Wq2vec = ak.zip(dict, with_name="PtEtaPhiMLorentzVector")


dict = {}
for key, val in key_dict.items():
    dict[key] = ak.flatten(ak.Array([
                                    ((evtDict["HH4V"]['genHiggs1W1dau1' + val] * h1W1l + evtDict["HH4V"]['genHiggs1W2dau1' + val] * h1W2l) * h1Wq1l + \
                                    (evtDict["HH4V"]['genHiggs1W1dau2' + val] * h1W1l + evtDict["HH4V"]['genHiggs1W2dau2' + val] * h1W2l) * h1Wq2l)[h1W],
                                    ((evtDict["HH4V"]['genHiggs2W1dau1' + val] * h2W1l + evtDict["HH4V"]['genHiggs2W2dau1' + val] * h2W2l) * h2Wq1l + \
                                    (evtDict["HH4V"]['genHiggs2W1dau2' + val] * h2W1l + evtDict["HH4V"]['genHiggs2W2dau2' + val] * h2W2l) * h2Wq2l)[h2W]
                                    ]))
dict['mass'] = np.zeros(totW)

Wlvec = ak.zip(dict, with_name="PtEtaPhiMLorentzVector")

len(hvec)

testr = 0.8

ak.sum((hvec.delta_r(Wlvec) < testr)) / totW
ak.sum((hvec.delta_r(Wq1vec) < testr)) / totW
ak.sum((hvec.delta_r(Wq2vec) < testr)) / totW

ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)) / totW
ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr) * (hvec.delta_r(Wlvec) < testr)) / totW
ak8_4q = (hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)


testr = 1.5

ak.sum((hvec.delta_r(Wlvec) < testr)) / totW
ak.sum((hvec.delta_r(Wq1vec) < testr)) / totW
ak.sum((hvec.delta_r(Wq2vec) < testr)) / totW

ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)) / totW
ak.sum((hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr) * (hvec.delta_r(Wlvec) < testr)) / totW
ak15_4q = (hvec.delta_r(Wq1vec) < testr) * (hvec.delta_r(Wq2vec) < testr)

hists['deltarql3'] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1", r"Higgs $p_T$ ", 51, 0, 1000),
                             )

str_all = 'Total HVV2ql$\\nu$ Events'
str_8 = '0.8 $\Delta$R Captures All 3 qql'
str_15 = '1.5 $\Delta$R Captures All 3 qql'

hists['deltarql3'].fill(sample=str_all, h1=hvec.pt, weight=weightsW)
hists['deltarql3'].fill(sample=str_8, h1=hvec.pt[ak8_4q], weight=weightsW[ak8_4q])
hists['deltarql3'].fill(sample=str_15, h1=hvec.pt[ak15_4q], weight=weightsW[ak15_4q])
# hists['genHiggsPt'].fill(sample='HH4b full', h1=higgs_pt[:, 0], weight=weights['fullHH4b'])



fig, (axs, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
hist.plot1d(hists['deltarql3'], line_opts=line_opts, ax=axs)
axs.set_xlabel(None)
# axs.set_ylim(0, 5)
axs.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

data_err_opts_c['color'] = 'blue'
hist.plotratio(hists['deltarql3'][str_8].sum('sample'), hists['deltarql3'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='0.8')
data_err_opts_c['color'] = 'orange'
hist.plotratio(hists['deltarql3'][str_15].sum('sample'), hists['deltarql3'][str_all].sum('sample'), ax=rax, error_opts=data_err_opts_c, unc='num', clear=False, label='1.5')
rax.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

rax.set_ylabel("captured/total")
rax.set_ylim(0, 1)
rax.grid(which='both')
plt.savefig('figs/deltarql3.pdf', bbox_inches='tight')
plt.show()


# HH4V

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
