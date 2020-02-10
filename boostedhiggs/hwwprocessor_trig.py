import numpy as np
from coffea import processor,hist,util
from uproot_methods import TLorentzVectorArray
from boostedhiggs.corrections import compiled
from .common import (getParticles)
import warnings
import argparse

class HwwProcessor_Trig(processor.ProcessorABC):
    def __init__(self, year):
        self._year = year
        self._corrections = compiled

        self._triggers = {
            '2017_had': [
                "AK8PFJet330_PFAK8BTagCSV_p17",
                "PFHT1050",
                "AK8PFJet400_TrimMass30",
                "AK8PFJet420_TrimMass30",
                "AK8PFHT800_TrimMass50",
                "PFJet500",
                "AK8PFJet500",
                ],
            '2017_muon': [
                "Mu50","Mu55",
                "Mu15_IsoVVVL_PFHT450_PFMET50","Mu15_IsoVVVL_PFHT600",
                ],
            '2017_electron': [
                "Ele27_WPTight_Gsf","Ele40_WPTight_Gsf","Ele20_WPLoose_Gsf",
                "Ele115_CaloIdVT_GsfTrkIdT",
                "Ele15_IsoVVVL_PFHT450_PFMET50","Ele15_IsoVVVL_PFHT600",
                ],
            '2017_met':[
                "PFMETNoMu110_PFMHTNoMu110_IDTight","PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","PFMETNoMu120_PFMHTNoMu120_IDTight",
                ],
            }

        self._triggers['2017_muonall'] = self._triggers['2017_muon']+self._triggers['2017_had']+self._triggers['2017_met']
        self._triggers['2017_electronall'] = self._triggers['2017_electron']+self._triggers['2017_had']+self._triggers['2017_met']

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        jetpt_axis = hist.Bin("jet_pt", r"Jet $p_T$", 30, 200, 800)
        genHpt_axis = hist.Bin("genH_pt", r"gen Higgs $p_T$",30, 200, 800)
        genElept_axis = hist.Bin("genEle_pt", r"gen Electron $p_T$",30, 20, 400)
        genMupt_axis = hist.Bin("genMu_pt", r"gen Muon $p_T$",30, 20, 400)
           
        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'
        hists['sumw'] = processor.defaultdict_accumulator(int)
        for sel in ['presel','noniso']:
            hists['Ele'+sel] = hist.Hist("Events / GeV",
                                          dataset_axis,
                                          jetpt_axis,
                                          genHpt_axis,
                                          genElept_axis,
                                          )
            hists['Mu'+sel] = hist.Hist("Events / GeV",
                                         dataset_axis,
                                         jetpt_axis,
                                         genHpt_axis,
                                         genMupt_axis,
                                         )
            for key in self._triggers.keys():
                key = sel+'_'+key.replace('_','')
                hists['Ele'+key] = hist.Hist("Events / GeV",
                                             dataset_axis,
                                             jetpt_axis,
                                             genHpt_axis,
                                             genElept_axis,
                                             )
                hists['Mu'+key] = hist.Hist("Events / GeV",
                                            dataset_axis,
                                            jetpt_axis,
                                            genHpt_axis,
                                            genMupt_axis,
                                            )
        self._accumulator = hists

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        dataset = df.metadata['dataset']
        output = self.accumulator.identity()
        selection = processor.PackedSelection()
        output = self.accumulator.identity()

        # pre-selection
        good = (
            (df.Muon.counts >= 1)
            & (df.FatJet.counts >= 1)
            )
        events = df[good]

        # trigger
        for key,item in self._triggers.items():
            trigger = np.ones(df.size, dtype='bool')
            for t in item:
                trigger = trigger | df.HLT[t]
            selection.add(key, trigger[good])

        # Muons
        goodMuon = (
            (events.Muon.pt > 30.)
            & (np.abs(events.Muon.eta) < 2.4)
            )
        nmuons = goodMuon.sum()

        # one Muon 
        selection.add('onemuon', (nmuons >= 1))

        # Electrons
        goodElectron = (
            (events.Electron.pt > 30)
            & (np.abs(events.Electron.eta) < 2.5)
            )
        nelectrons = goodElectron.sum()

        # one Electron
        selection.add('oneelectron', (nelectrons >= 1))

        # select FatJet with muon closest
        goodFatJet = (
            (events.FatJet.pt > 300.)
            & (np.abs(events.FatJet.eta) < 2.4)
            & (events.FatJet.msoftdrop > 10.)
            & (events.FatJet.jetId & 2)
            )
        leadingmuon = events.Muon[goodMuon][:, 0:1]
        leadingele = events.Electron[goodElectron][:, 0:1]
        leadingjet = events.FatJet[goodFatJet][:,0:1]

        selection.add('jetkin', (leadingjet.pt > 300).any())
        selection.add('nonisomuon', (leadingmuon.pt > 50).any())
        selection.add('nonisoelectron', (leadingmuon.pt > 115).any())

        # gen matching
        genH,genH_idx = getParticles(events,25,['fromHardProcess', 'isLastCopy'])
        genW,genW_idx = getParticles(events,24,['fromHardProcess', 'isLastCopy'])
        genE,genE_idx = getParticles(events,11,['fromHardProcess','isFirstCopy'],1)
        genM,genM_idx = getParticles(events,13,['fromHardProcess','isFirstCopy'],1)
        
        ishWW_qqelev = (genH.counts==1) & (genW.counts>=1) & (genE.counts==1)
        ishWW_qqmuv = (genH.counts==1) & (genW.counts>=1) & (genM.counts==1)
        
        selection.add('genEle',ishWW_qqelev)
        selection.add('genMu',ishWW_qqmuv)

        regions = {}
        for key in self._triggers.keys():
            regions[key.replace('_','')] = {key}
        regions['Elepresel'] = {'jetkin', 'oneelectron','genEle'}
        regions['Mupresel'] = {'jetkin', 'onemuon','genMu'}
        regions['Elenoniso'] = {'jetkin', 'oneelectron','genEle','nonisoelectron'}
        regions['Munoniso'] = {'jetkin', 'onemuon','genMu','nonisomuon'}

        for histname, h in output.items():
            if not isinstance(h, hist.Hist):
                continue
            if histname not in regions:
                continue
            region = []
            for r in regions.keys():
                if len(histname.split('_'))> 2:
                    histname2 = histname
                else:
                    histname2 = histname.split('_')
                if r in histname2:
                    region.extend(regions[r])
            cut = selection.all(*region)
            weight = cut

            def normalize(val):
                return val.pad(1, clip=True).fillna(0).flatten()

            if 'Mu' in histname:
                h.fill(jet_pt = normalize(leadingjet.pt),
                       genH_pt = normalize(genH.pt),
                       genMu_pt = normalize(genM.pt),
                       dataset=dataset,weight=weight)
            else:
                h.fill(jet_pt = normalize(leadingjet.pt),
                       genH_pt = normalize(genH.pt),
                       genEle_pt = normalize(genE.pt),
                       dataset=dataset,weight=weight)

        return output

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW processor for gen Trig studies')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    args = parser.parse_args()

    processor_instance = HwwProcessor_Trig(year=args.year)

    util.save(processor_instance, 'hwwprocessor_trig.coffea')
