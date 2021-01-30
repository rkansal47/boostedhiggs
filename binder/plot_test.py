import uproot4
import awkward1 as ak

# %matplotlib inline
import matplotlib.pyplot as plt
import mplhep as hep

import coffea.hist as hist

import numpy as np
from cycler import cycler

plt.style.use(hep.style.ROOT)


evts = uproot4.open("./data/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root")["tree"]
evts.show()


def getData(fnames=''):
    evts = uproot4.concatenate(fnames + ":tree")
    return evts


samples = {
    "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWWqq": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = getData(fname)

evts['fatJet1Pt']
evts['weight']

evtDict["QCD"]['fatJet1Pt']
evtDict["HHbbWWqq"]['fatJet1Pt']

evtDict["QCD"]['weight']
evtDict["HHbbWWqq"]['weight']

plt.hist(np.array(evtDict["HHbbWWqq"]['fatJet1Pt']))

hists = {}
hists["jet_pt"] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("jet1pt", r"Leading Jet $p_T$ ", 40, 200, 2000),
                             hist.Bin("jet2pt", r"Sub-Leading Jet $p_T$", 40, 200, 2000),
                             hist.Bin("jet3pt", r"3rd-Leading Jet $p_T$", 4, 200, 2000))

hists["genHiggs"] = hist.Hist("Events",
                             hist.Cat("sample", "Sample"),
                             hist.Bin("h1pt", r"Gen Higgs 1 $p_T$ ", 40, 200, 2000),
                             hist.Bin("h2pt", r"Gen Higgs 2 $p_T$", 40, 200, 2000))

for s, evts in evtDict.items():
    hists["jet_pt"].fill(sample=s,
                         jet1pt = evts["fatJet1Pt"],
                         jet2pt = evts["fatJet2Pt"],
                         jet3pt = evts["fatJet3Pt"],
                         # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                         )
    # for gen information we are only interested in looking at signal samples
    if("HH" in s):
        hists["genHiggs"].fill(sample=s,
                               h1pt = evts["genHiggs1Pt"],
                               h2pt = evts["genHiggs2Pt"])


# you can make comparison plots of the shape of histograms
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
hist.plot1d(hists["jet_pt"].sum("jet2pt", "jet3pt"), overlay='sample', ax=axs[0], density=True)
hist.plot1d(hists["jet_pt"].sum("jet1pt", "jet3pt"), overlay='sample', ax=axs[1], density=True)
fig.tight_layout()

fig, axs = plt.subplots(1, 2, figsize=(12, 6))
hist.plot1d(hists["genHiggs"].sum("h2pt"), overlay='sample', ax=axs[0], stack=True)
hist.plot1d(hists["genHiggs"].sum("h1pt"), overlay='sample', ax=axs[1], stack=True)
fig.tight_layout()

# or stack plots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
hist.plot1d(hists["jet_pt"].sum("jet2pt", "jet3pt"), overlay='sample', ax=axs[0], stack=True)
hist.plot1d(hists["jet_pt"].sum("jet1pt", "jet3pt"), overlay='sample', ax=axs[1], stack=True)
fig.tight_layout()



kin_vars = ['Pt', 'Eta', 'Phi', 'Mass']
kin_vars_labels = [r'$p_T$', r'$\eta$', r'$\varphi$', 'Mass']
bins = [[50, 0, 2000],
        [10, -3.5, 3.5],
        [10, -3.5, 3.5],
        [50, 0, 400]]



for i in range(len(kin_vars)):
    var = kin_vars[i]
    varl = kin_vars_labels[i]
    hists["jet" + var] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin("jet1", r"Leading Jet " + varl, bins[i][0], bins[i][1], bins[i][2]),
                                    hist.Bin("jet2", r"Sub-Leading Jet " + varl, bins[i][0], bins[i][1], bins[i][2]),
                                    hist.Bin("jet3", r"3rd-Leading Jet " + varl, bins[i][0], bins[i][1], bins[i][2])
                                    )


    # hists["genHiggs" + var] = hist.Hist("Events",
    #                                    hist.Cat("sample", "Sample"),
    #                                    hist.Bin("h1", r"Gen Higgs 1 " + varl, bins[i][0], bins[i][1], bins[i][2]),
    #                                    hist.Bin("h2", r"Gen Higgs 2 " + varl, bins[i][0], bins[i][1], bins[i][2])
    #                                    )


for s, evts in evtDict.items():
    for i in range(len(kin_vars)):
        var = kin_vars[i]
        hists["jet" + var].fill(sample=s,
                             jet1 = evts["fatJet1" + var],
                             jet2 = evts["fatJet2" + var],
                             jet3 = evts["fatJet3" + var],
                             weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
                             )
        # for gen information we are only interested in looking at signal samples
        # if("HH" in s):
        #     print(evts["fatJet1" + var])
        #     hists["genHiggs" + var].fill(sample=s,
        #                            h1 = evts["genHiggs1" + var],
        #                            h2 = evts["genHiggs2" + var]
        #                            )


# hists["jetEta"].rebin()

# hists["jetPt"].project("sample", "jet1")
# hists["jetPt"].sum("jet1", "jet2")
#
# hists['jetPt'].identifiers("sample")
#
# hists['jetPt'][:2].identifiers('sample')
#
# hists['jetPt'].scale(10)
#
#
#
# list(map(str, hists['jetPt'].identifiers("sample")))
#
hists['jetEta'][:2].values()

# colors = ['red', 'orange', 'blue', 'green']
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

fig, axs = plt.subplots(len(kin_vars), 3, figsize=(28, 36))
# axs.set_prop_cycle(cycler(color=colors))
for i in range(len(kin_vars)):
    var = kin_vars[i]
    for j in range(3):
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
        hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], density=True, clear=False, line_opts=line_opts)
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
        ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], density=True, stack=True, clear=False, fill_opts=fill_opts)
        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
plt.tight_layout(0.5)
plt.savefig("jet_kinematics.pdf", bbox_inches='tight')
plt.show()
