import uproot4
import awkward1 as ak

# %matplotlib inline
import matplotlib.pyplot as plt
import mplhep as hep

import coffea.hist as hist

import numpy as np
from cycler import cycler

from skhep.math.vectors import LorentzVector

from tqdm import tqdm

plt.style.use(hep.style.ROOT)


evts = uproot4.open("./data/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root")["tree"]
evts.show()

#
# def getData(fnames=''):
#     evts = uproot4.concatenate(fnames + ":tree")
#     return evts


samples = {
    "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWWqq": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot4.concatenate(fname + ":tree")


evts['fatJet1Pt']
evts['weight']

evtDict["QCD"]['fatJet1Pt']
evtDict["HHbbWWqq"]['fatJet1Pt']

evtDict["QCD"]['weight']
evtDict["HHbbWWqq"]['weight']


plt.hist(evtDict["HHbbWWqq"]['weight'])


plt.hist(evtDict["QCD"]['weight'][-100000:], bins=np.linspace(0, 0.001, 101))

evtDict["HHbbWWqq"]['fatJet1DeepAK8H']

xbb = np.array(evtDict["HHbbWWqq"]['fatJet1PNetXbb'])
xbb[xbb > 0.00001].shape

xbb[xbb > 0.001]



evts['fatJet1Mass']

np.nonzero(np.array(evtDict["HHbbWWqq"]['fatJet1PNetXbb']))

np.array(evtDict["HHbbWWqq"]['fatJet1PNetXbb']).shape
np.array(evtDict["HHbbWWqq"]['fatJet1PNetXbb']).nonzero()[0].shape


plt.hist(evtDict["HHbbWWqq"]['fatJet1PNetXbb'], bins=np.linspace(-1, 1, 100))
plt.title("HHbbWWqq")
plt.xlabel("Fat Jet 1 PNet Score")

plt.hist(evtDict["HHbbWWqq"]['fatJet1DeepAK8H'], bins=np.linspace(-1, 1, 100))
plt.title("HHbbWWqq")
plt.xlabel("Fat Jet 1 Deep AK8H Score")


plt.hist(np.array(evtDict["HHbbWWqq"]['fatJet1Pt']))
#
# hists = {}
# hists["jet_pt"] = hist.Hist("Events",
#                              hist.Cat("sample", "Sample"),
#                              hist.Bin("jet1pt", r"Leading Jet $p_T$ ", 40, 200, 2000),
#                              hist.Bin("jet2pt", r"Sub-Leading Jet $p_T$", 40, 200, 2000),
#                              hist.Bin("jet3pt", r"3rd-Leading Jet $p_T$", 4, 200, 2000))
#
# hists["genHiggs"] = hist.Hist("Events",
#                              hist.Cat("sample", "Sample"),
#                              hist.Bin("h1pt", r"Gen Higgs 1 $p_T$ ", 40, 200, 2000),
#                              hist.Bin("h2pt", r"Gen Higgs 2 $p_T$", 40, 200, 2000))
#
# for s, evts in evtDict.items():
#     hists["jet_pt"].fill(sample=s,
#                          jet1pt = evts["fatJet1Pt"],
#                          jet2pt = evts["fatJet2Pt"],
#                          jet3pt = evts["fatJet3Pt"],
#                          # weight = evts["totalWeight"],  # weight is a reserved keyword in Hist, and can be added to any fill() call
#                          )
#     # for gen information we are only interested in looking at signal samples
#     if("HH" in s):
#         hists["genHiggs"].fill(sample=s,
#                                h1pt = evts["genHiggs1Pt"],
#                                h2pt = evts["genHiggs2Pt"])
#
#
# # you can make comparison plots of the shape of histograms
# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# hist.plot1d(hists["jet_pt"].sum("jet2pt", "jet3pt"), overlay='sample', ax=axs[0], density=True)
# hist.plot1d(hists["jet_pt"].sum("jet1pt", "jet3pt"), overlay='sample', ax=axs[1], density=True)
# fig.tight_layout()
#
# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# hist.plot1d(hists["genHiggs"].sum("h2pt"), overlay='sample', ax=axs[0], stack=True)
# hist.plot1d(hists["genHiggs"].sum("h1pt"), overlay='sample', ax=axs[1], stack=True)
# fig.tight_layout()
#
# # or stack plots
# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# hist.plot1d(hists["jet_pt"].sum("jet2pt", "jet3pt"), overlay='sample', ax=axs[0], stack=True)
# hist.plot1d(hists["jet_pt"].sum("jet1pt", "jet3pt"), overlay='sample', ax=axs[1], stack=True)
# fig.tight_layout()
#

#
# kin_vars = ['Pt', 'Eta', 'Phi', 'Mass', 'MassSD']
# kin_vars_labels = [r'$p_T$', r'$\eta$', r'$\varphi$', 'Mass', 'Soft Drop Mass']
# bins = [[50, 0, 2000],
#         [10, -3.5, 3.5],
#         [10, -3.5, 3.5],
#         [50, 0, 400],
#         [50, 0, 400]]

hists = {}

kin_vars = ['Pt', 'Mass', 'MassSD']
kin_vars_labels = [r'$p_T$', 'Mass', 'Soft Drop Mass']
kin_bins = [[50, 0, 2000],
        [50, 0, 400],
        [50, 0, 400]]


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
        # print(evts["totalWeight"][:10])
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

scale = {'HH4V': 1 / tot_events['HH4V'],
        'HHbbWWqq': 1 / tot_events['HHbbWWqq'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}

for var in kin_vars:
    hists['jet' + var].scale(scale, axis='sample')


# list(map(str, hists['jetPt'].identifiers("sample")))

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

fig, axs = plt.subplots(len(kin_vars), 3, figsize=(28, len(kin_vars) * 9))
# axs.set_prop_cycle(cycler(color=colors))
for i in range(len(kin_vars)):
    var = kin_vars[i]
    for j in range(3):
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
        hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
        ylim1 = axs[i, j].get_ylim()[1]
        # print(ylim1)
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
        ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
        ylim2 = axs[i, j].get_ylim()[1]
        # print(ylim2)
        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
        # print(max(ylim1, ylim2))
        axs[i, j].set_ylim(0, 0.25)

plt.tight_layout(0.5)
plt.savefig("figs/jet_kinematics_weighted.pdf", bbox_inches='tight')
plt.show()


disc_var = 'DeepAK8MD_H4qvsQCD'
disc_var_label = 'DeepAK8MD_H4qvsQCD'
disc_var_bins = [100, -1, 1]





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




fig, axs = plt.subplots(len(disc_vars), 3, figsize=(28, len(disc_vars) * 9))
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')

for i in range(len(disc_vars)):
    var = disc_vars[i]
    for j in range(3):
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
        hist.plot1d(hists["jet" + var][:2].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
        ax = hist.plot1d(hists["jet" + var][2:].project("sample", "jet" + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)
plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
plt.savefig("figs/jet_disc_vars_weighted.pdf", bbox_inches='tight')
plt.show()


jet12mass = {}
jet123mass = {}

for s, evts in evtDict.items():
    print(s)
    jet12mass[s] = []
    jet123mass[s] = []
    for i in tqdm(range(len(evts["fatJet1Pt"]))):
        jet1 = LorentzVector()
        jet1.setptetaphim(evts["fatJet1Pt"][i], evts["fatJet1Eta"][i], evts["fatJet1Phi"][i], evts["fatJet1Mass"][i])
        jet2 = LorentzVector()
        jet2.setptetaphim(evts["fatJet2Pt"][i], evts["fatJet2Eta"][i], evts["fatJet2Phi"][i], evts["fatJet2Mass"][i])
        jet3 = LorentzVector()
        jet3.setptetaphim(evts["fatJet3Pt"][i], evts["fatJet3Eta"][i], evts["fatJet3Phi"][i], evts["fatJet3Mass"][i])

        jet12mass[s].append((jet1 + jet2).mass)
        jet123mass[s].append((jet1 + jet2 + jet3).mass)



hists["jets_mass"] = hist.Hist("Events",
                                hist.Cat("sample", "Sample"),
                                hist.Bin("jet12", r"Leading 2 Jets' Combined Mass", 50, 0, 800),
                                hist.Bin("jet123", r"Leading 3 Jets' Combined Mass", 50, 0, 800),
                                )


for s in evtDict.keys():
    hists["jets_mass"].fill(sample=s,
                         jet12 = ,
                         jet123 = evts["fatJet1Mass"] + evts["fatJet2Mass"] + evts["fatJet3Mass"]
                         )

tot_events = {}
vals = hists['jets_mass'].project("sample", "jet12").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


hists['jets_mass'].scale(scale, axis='sample')

fig, axs = plt.subplots(1, 2, figsize=(16, 9))

i = 0
axs[0].set_prop_cycle(cycler(color=colors_cycle))
hist.plot1d(hists["jets_mass"][:2].project("sample", "jet12"), overlay='sample', ax=axs[0], clear=False, line_opts=line_opts)
axs[0].set_prop_cycle(cycler(color=colors_cycle[2:]))
ax = hist.plot1d(hists["jets_mass"][2:].project("sample", "jet12"), overlay='sample', ax=axs[0], stack=True, clear=False, fill_opts=fill_opts)
axs[0].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

i = 1
axs[1].set_prop_cycle(cycler(color=colors_cycle))
hist.plot1d(hists["jets_mass"][:2].project("sample", "jet123"), overlay='sample', ax=axs[1], clear=False, line_opts=line_opts)
axs[1].set_prop_cycle(cycler(color=colors_cycle[2:]))
ax = hist.plot1d(hists["jets_mass"][2:].project("sample", "jet123"), overlay='sample', ax=axs[1], stack=True, clear=False, fill_opts=fill_opts)
axs[1].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

plt.tight_layout(0.5)

plt.savefig("figs/jets_mass.pdf", bbox_inches='tight')
plt.show()




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
