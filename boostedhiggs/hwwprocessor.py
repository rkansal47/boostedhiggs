import numpy as np
from coffea import processor,hist,util
from uproot_methods import TLorentzVectorArray
from boostedhiggs.corrections import (
     compiled,
     corrected_msoftdrop,
     add_pileup_weight,
)
import warnings
import argparse

class HwwProcessor(processor.ProcessorABC):
    def __init__(self, year, trigger, channel, regions):
        self._year = year
        self._corrections = compiled
        self._btagWPs = {
            'med': {
                '2016': 0.6321,
                '2017': 0.4941,
                '2018': 0.4184,
            },
        }

        self._channel = channel
        self._trigger = trigger
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
                "Mu50",
                "Mu55",
                "OldMu100",
                "TkMu100",
                ],
            '2017_vvlmuon':[
                "Mu15_IsoVVVL_PFHT450_PFMET50",
                "Mu15_IsoVVVL_PFHT600",
                ],
            '2017_electron': [
                "Ele27_WPTight_Gsf",
                "Ele40_WPTight_Gsf",
                "Ele20_WPLoose_Gsf",
                "Ele115_CaloIdVT_GsfTrkIdT",
                ],
            '2017_vvlelectron': [
                "Ele15_IsoVVVL_PFHT450_PFMET50",
                "Ele15_IsoVVVL_PFHT600",
                ],
            '2017_met':[
                "PFMETNoMu110_PFMHTNoMu110_IDTight",
                "PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                "PFMETNoMu120_PFMHTNoMu120_IDTight",
                ],
            }
        self._triggers['2017_muonall'] = self._triggers['2017_muon']+self._triggers['2017_had']+self._triggers['2017_vvlmuon']
        self._triggers['2017_electronall'] = self._triggers['2017_electron']+self._triggers['2017_had']+self._triggers['2017_vvlelectron']

        self._fjetptMIN = 300.
        if trigger=='had':
            self._fjetptMIN = 450.

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        fjetpt_axis = hist.Bin("fjet_pt", r"Jet $p_T$", 20, 300, 1000)
        fjetlsf3_axis = hist.Bin("fjet_lsf3", r"Jet LSF$_3$", 20, 0, 1)
        fjetmsd_axis = hist.Bin("fjet_msd", r"Jet m$_{SD}$", 20, 20, 200)
        fjetmuondR_axis = hist.Bin("fjet_muondR", r"$\Delta R_{j,l}", 20, 0, 4.0)
        jetoppbtag_axis = hist.Bin("jet_oppbtag",r"Jet Opposite AK4 max b-tag", 20, 0, 1)
        jetawaybtag_axis = hist.Bin("jet_awaybtag",r"Jet Away AK4 max b-tag", 20, 0, 1)
        jetmmass_axis = hist.Bin("jet_mmass",r"Jet - Lep Mass", 20, 0, 120)
        flsjetpt_axis = hist.Bin("flsjet_pt", r"Lep Subtracted Jet $p_T$", 20, 300, 1000)
        flsjetmsd_axis = hist.Bin("flsjet_msd", r"Lep Subtracted Jet m$_{SD}$", 20, 0, 200)
        #flsjetbtag_axis = hist.Bin("flsjet_btag", r"Lep Subtracted Jet max b-tag", 20, 0, 1)
        muonpt_axis = hist.Bin("muon_pt", r"Muon $p_T$", 20, 0, 450)
        muonmiso_axis = hist.Bin("muon_miso", r"Muon mini PF ISO (total)", 20, 0, 1)
        metpt_axis = hist.Bin("met_pt", r"MET $p_T$", 20, 0, 350)
        metphi_axis = hist.Bin("met_phi", r"MET $\phi$", 20, -3.5, 3.5)
        
        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'
        hists['sumw'] = processor.defaultdict_accumulator(int)
        hists['cutflow'] = processor.defaultdict_accumulator(float)
        self._regions = regions.split(',')
        for region in self._regions:
            hists['%s_fjetprop'%region] = hist.Hist("Events / GeV",
                                                    dataset_axis,
                                                    fjetpt_axis,
                                                    fjetmsd_axis,
                                                    fjetlsf3_axis,
                                                    jetoppbtag_axis,
                                                )
            hists['%s_fmmjetprop'%region] = hist.Hist("Events / GeV",
                                                      dataset_axis,
                                                      jetmmass_axis)
            hists['%s_flsjetprop'%region] = hist.Hist("Events / GeV",
                                                      dataset_axis,
                                                      flsjetpt_axis,
                                                      flsjetmsd_axis,
                                                      fjetmuondR_axis,
                                                      jetawaybtag_axis,
                                                      #flsjetbtag_axis,
                                                  )
            hists['%s_muonprop'%region] = hist.Hist("Events / GeV",
                                                    dataset_axis,
                                                    muonpt_axis,
                                                    muonmiso_axis,
                                                )
            hists['%s_metprop'%region] = hist.Hist("Events / GeV",
                                                   dataset_axis,
                                                   metpt_axis,
                                                   metphi_axis)

        self._accumulator = hists

        print(self._fjetptMIN,self._trigger,self._channel)
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata['dataset']
        isRealData = 'genWeight' not in events.columns
        output = self.accumulator.identity()
        selection = processor.PackedSelection()
        output = self.accumulator.identity()

        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()

        # trigger
        trigger = np.zeros(events.size, dtype='bool')
        for t in self._triggers[self._year+'_'+self._trigger]:
            try:
                trigger = trigger | events.HLT[t]
            except:
                warnings.warn("Missing trigger %s" % t, RuntimeWarning)
        selection.add('trigger', trigger)

        # Muons
        goodMuon = (
            (events.Muon.pt > 27.)
            & (np.abs(events.Muon.eta) < 2.4)
            )
        nmuons = goodMuon.sum()

        # Electrons
        goodElectron = (
            (events.Electron.pt > 10)
            & (np.abs(events.Electron.eta) < 2.5)
            )
        nelectrons = goodElectron.sum()

        selection.add('onemuon', (nmuons >= 1)) 
        selection.add('oneelectron', (nelectrons >= 1))

        # FatJets
        fatjets = events.FatJet
        fatjets['msdcorr'] = corrected_msoftdrop(events.FatJet)
        goodFatJet = (
            (fatjets.pt > 300)
            & (np.abs(fatjets.eta) < 2.4)
            & (fatjets.msdcorr > 10.)
            & (fatjets.jetId & 2)
            )
        nfatjets = goodFatJet.sum()

        candidatejet = fatjets[goodFatJet][:,0:1]
        selection.add('jetkin', (candidatejet.pt > self._fjetptMIN).any())

        candidatemuon = events.Muon[goodMuon][:,0:1] 
        candidateelectron = events.Electron[goodElectron][:,0:1]

        # FatJets and Leptons
        selection.add('LSF3medium', (candidatejet.lsf3>0.7).any())

        ak8_muon_pair = candidatemuon.cross(events.FatJet)
        ak8_muon_dR = ak8_muon_pair.i0.delta_r(ak8_muon_pair.i1)

        print('nmuons>=1', (nmuons >= 1).any())
        print('nfatjets>=1', (nfatjets >= 1).any())

        good = (
            (events.Muon.counts >= 1)
            & (events.FatJet.counts >= 1)
            )
        eventsWmuon = events[good]
        print('eventsWmuon',len(events),len(eventsWmuon))
        leadingmuon = eventsWmuon.Muon[:,0:1]
        ak8_muon_pairWmuon = leadingmuon.cross(eventsWmuon.FatJet)
        ak8_muon_dRWmuon = ak8_muon_pairWmuon.i0.delta_r(ak8_muon_pairWmuon.i1)
        leadingjet = eventsWmuon.FatJet[ak8_muon_dRWmuon.argmin()]
        print('len jet muon',len(leadingjet.flatten()),len(leadingmuon.flatten()))

        print('len cand jet muon',len(candidatejet.flatten()),len(candidatemuon.flatten()))

        # FatJet substracted Lepton
        sj1_sj2_btagDeepB_pair = candidatejet.LSsubJet1btagDeepB.cross(candidatejet.LSsubJet2btagDeepB)
        #fls_btagDeepB_max = max(sj1_sj2_btagDeepB_pair.i0,sj1_sj2_btagDeepB_pair.i1)

        # Jets
        goodJet = (
            (events.Jet.pt > 30.)
            & (events.Jet.jetId & 2)
            )
        jets = events.Jet[goodJet]
        ak4_ak8_pair = jets.cross(candidatejet, nested=True)
        ak4_ak8_dphi = abs(ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1))
        ak4_opposite = jets[(ak4_ak8_dphi > np.pi / 2).all()]
        ak4_away = jets[(ak4_ak8_dphi > 0.8).all()]

        selection.add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < self._btagWPs['med'][self._year])
        selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() > self._btagWPs['med'][self._year])

        # MET
        MET = events.MET

        # MET eta with mass assumption
        met = eventsWmuon.MET
        mm = (leadingjet - leadingmuon).mass2
        jmass = (mm>0)*np.sqrt(np.maximum(0, mm)) + (mm<0)*leadingjet.mass
        joffshell = jmass < 62.5
        massassumption = 80.*joffshell + (125 - 80.)*~joffshell
        x = massassumption**2/(2*leadingmuon.pt*met.pt) + np.cos(leadingmuon.phi - met.phi)
        met_eta = (
            (x < 1)*np.arcsinh(x*np.sinh(leadingmuon.eta))
            + (x > 1)*(
                leadingmuon.eta - np.sign(leadingmuon.eta)*np.arccosh(leadingmuon.eta)
                )
            )

        met_p4 = TLorentzVectorArray.from_ptetaphim(np.array([0.]),np.array([0.]),np.array([0.]),np.array([0.]))
        if met.size > 0:
            met_p4 = TLorentzVectorArray.from_ptetaphim(met.pt, met_eta.fillna(0.), met.phi, np.zeros(met.size)) 

        # fill cutflow
        cutflow = ['trigger', 'jetkin', 'onemuon', 'LSF3medium','antiak4btagMediumOppHem','ak4btagMedium08']
        allcuts = set()
        output['cutflow']['none'] += len(events)
        for cut in cutflow:
            allcuts.add(cut)
            output['cutflow'][cut] += selection.all(*allcuts).sum()

        # weights
        weights = processor.Weights(len(events))
        weightsWmuon = processor.Weights(len(eventsWmuon))

        if not isRealData:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year)
            weightsWmuon.add('genweight', eventsWmuon.genWeight)
            add_pileup_weight(weightsWmuon, eventsWmuon.Pileup.nPU, self._year)

        regions = {}
        regions['presel'] = {'trigger','jetkin', 'onemuon'}
        regions['lsf07'] = {'trigger','jetkin', 'onemuon','LSF3medium'}
        regions['bopp'] =  {'trigger','jetkin', 'onemuon','LSF3medium','antiak4btagMediumOppHem'}

        for region in self._regions:
            selections = regions[region]
            cut = selection.all(*selections)
            weight = weights.weight()[cut] 
            
            def normalize(val):
                try:
                    return val[cut].pad(1, clip=True).fillna(0).flatten() 
                except:
                    return val[cut].flatten()

            output['%s_fjetprop'%region].fill(fjet_pt = normalize(candidatejet.pt),
                                              fjet_msd = normalize(candidatejet.msdcorr),
                                              fjet_lsf3 = normalize(candidatejet.lsf3),
                                              jet_oppbtag = normalize(ak4_opposite.btagDeepB.max()),
                                              dataset=dataset,
                                              weight=weight)
            output['%s_fmmjetprop'%region].fill(jet_mmass = jmass.flatten(),
                                                dataset=dataset,
                                                weight=weightsWmuon.weight())
            output['%s_flsjetprop'%region].fill(flsjet_pt = normalize(candidatejet.LSpt),
                                                flsjet_msd = normalize(candidatejet.LSmsoftdrop),
                                                jet_awaybtag = normalize(ak4_away.btagDeepB.max()),
                                                fjet_muondR = normalize(ak8_muon_dR.min()),
                                                dataset=dataset,
                                                weight=weight)
            output['%s_muonprop'%region].fill(muon_pt = normalize(candidatemuon.pt),
                                              muon_miso = normalize(candidatemuon.miniPFRelIso_all),
                                              dataset=dataset,
                                              weight=weight)
            output['%s_metprop'%region].fill(met_pt = normalize(MET.pt),
                                             met_phi = normalize(MET.phi),
                                             dataset=dataset,
                                             weight=weight)

        return output


    def postprocess(self, accumulator):

        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]                                                                                                                                                        
        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if dataset in self._corrections['xsections']:
                scale[dataset] = lumi*self._corrections['xsections'][dataset]/dataset_sumw
            else:
                warnings.warn("Missing cross section for dataset %s.  Normalizing to 1 pb" % dataset, RuntimeWarning)
                scale[dataset] = lumi / dataset_sumw

        for h in accumulator.values():
            if isinstance(h, hist.Hist):
                h.scale(scale, axis="dataset")

        return accumulator

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW processor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--trigger', choices=['muon','electron','had','muonall','electronall'], default='muon', help='trigger selection')
    parser.add_argument('--channel', choices=['muon','electron'], default='muon', help='channel')
    parser.add_argument('--regions', default='presel', help='regions')
    args = parser.parse_args()

    print('hww args',args)
    processor_instance = HwwProcessor(year=args.year,trigger=args.trigger,channel=args.channel,regions=args.regions)

    util.save(processor_instance, 'hwwprocessor.coffea')
