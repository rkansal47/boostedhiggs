import uproot4
import awkward1 as ak
import numpy as np

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
