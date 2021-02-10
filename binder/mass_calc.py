import uproot
from skhep.math.vectors import LorentzVector
from tqdm import tqdm
import numpy as np
import sys
from os import listdir

j = int(sys.argv[1])
print(j)

j = 5

# dir = '/graphganvol/data/'
dir = 'data/'

fnames = listdir(dir + 'weighted/')

evts = uproot4.concatenate(dir + 'weighted/' + fnames[j] + ":tree")

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
    jet123mass.append((jet1 + jet2 + jet3).mass)

np.save(dir + fnames[j] + '_inv_mass.npy', np.array([jet12mass, jet123mass]))
