import uproot4
import awkward1 as ak

import matplotlib.pyplot as plt
import mplhep as hep

import coffea.hist as hist

import numpy as np
from cycler import cycler

from tqdm import tqdm
import pandas

from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
plt.style.use(hep.style.ROOT)


samples = {
    # "HH4V": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    # "HH4b": "data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWW4q": "data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "data/weighted/QCD_HT*.root",
    "tt": "data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    print("Loading {} samples".format(s))
    evtDict[s] = uproot4.concatenate(fname + ":tree")

print("Samples loaded")

# in fb
RUN2LUMI = 137
XSECHHBBWWQQ = 1.82
XSECHHBBBB = 31.05 * 0.58**2
ACCEPTANCE = 0.16928763440860214

weights = {}

# QCD, tt samples are already weighted with xsecs
# Acceptance calculated using BSM ggHH sample - need to check on SM one

# weights["HH4V"] = None
weights["QCD"] = evtDict["QCD"]["weight"] * RUN2LUMI
weights["tt"] = evtDict["tt"]["weight"] * RUN2LUMI
weights["HHbbWW4q"] = evtDict["HHbbWW4q"]["totalWeight"] * RUN2LUMI * XSECHHBBWWQQ * ACCEPTANCE
# weights["HH4b"] = np.ones(len(evtDict["HH4b"])) * RUN2LUMI * XSECHHBBBB * ACCEPTANCE / len(evtDict["HH4b"])



# fat jets by PNetXbb score
# output is a dict of all events with fatJet1 and 2 referring to the Hbb and HWW candidates resp., and storing the pT, soft drop mass, and tagger scores of each

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


# first value is minimum, second is maximum (9999 means no max)

var_cuts = {
    "fatJet1Pt": [300, 9999],
    "fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [75, 150],
    "fatJet2MassSD": [50, 150],
    "fatJet1PNetXbb": [0.99, 9999],
    "fatJet2DeepAK8MD_H4qvsQCD": [0.9, 9999],
}

# calculating the cuts and yields

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


print("Final event yields \nSignal: {:.2f}, BG: {:.2f}".format(cutflow['HHbbWW4q'][-1], cutflow['QCD'][-1] + cutflow['tt'][-1]))

# saving cut-flow table

cut_labels = ['Jet1 pT > 300',
                'Jet2 pT > 300',
                'Jet1 MassSD > 75',
                'Jet1 MassSD < 150',
                'Jet1 MassSD > 50',
                'Jet1 MassSD < 150',
                'Jet1 PNetXbb > 0.99',
                'Jet2 DeepAK8H4qvsQCD > 0.9',
            ]

cut_idx = [0, 2, 4, 5, 6, 7, 8, 10]

cftable = pandas.DataFrame(np.round(np.array(list(cutflow.values()))[:, cut_idx], 2), list(cutflow.keys()), cut_labels)
cftable.to_csv('cutflow.csv')

print("Cut-flow table saved")
