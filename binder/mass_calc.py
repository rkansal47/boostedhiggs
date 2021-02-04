import uproot
from skhep.math.vectors import LorentzVector
from tqdm import tqdm
import numpy as np

samples = {
    "HH4V": "/graphganvol/data/weighted/HHToVVVV_node_SM_Pt300_1pb_weighted.root",
    "HHbbWWqq": "/graphganvol/data/weighted/HHToBBVVToBBQQQQ_node_SM_1pb_weighted.root",
    "QCD": "/graphganvol/data/weighted/QCD_HT*.root",
    "tt": "/graphganvol/data/weighted/TTTo*.root",
}

evtDict = {}
for s, fname in samples.items():
    evtDict[s] = uproot.concatenate(fname + ":tree")

for s, evts in evtDict.items():
    print(s)
    jet12mass = []
    jet123mass = []
    for i in tqdm(range(len(evts["fatJet1Pt"]))):
        jet1 = LorentzVector()
        jet1.setptetaphim(evts["fatJet1Pt"][i], evts["fatJet1Eta"][i], evts["fatJet1Phi"][i], evts["fatJet1Mass"][i])
        jet2 = LorentzVector()
        jet2.setptetaphim(evts["fatJet2Pt"][i], evts["fatJet2Eta"][i], evts["fatJet2Phi"][i], evts["fatJet2Mass"][i])
        jet3 = LorentzVector()
        jet3.setptetaphim(evts["fatJet3Pt"][i], evts["fatJet3Eta"][i], evts["fatJet3Phi"][i], evts["fatJet3Mass"][i])

        jet12mass.append((jet1 + jet2).mass)
        print((jet1 + jet2).mass)

        jet123mass.append((jet1 + jet2 + jet3).mass)
        print((jet1 + jet2 + jet3).mass)

    np.save('/graphganvol/data/' + s + '_inv_mass.npy', np.array([jet12mass, jet123mass]))
