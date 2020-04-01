import numpy as np
import coffea
from coffea import processor,hist,util
from uproot_methods import TLorentzVectorArray
from boostedhiggs.common import (
    getParticles,
    matchedParticleFlavor,
    getFlavor,
)
from boostedhiggs.corrections import (
     compiled,
     corrected_msoftdrop,
     add_pileup_weight,
)
import warnings
import argparse

import logging
logger = logging.getLogger(__name__)

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

        genflavor_axis = hist.Bin('genflavor', 'Gen. jet flavor', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        fjetpt_axis = hist.Bin("fjet_pt", r"Cand Jet $p_T$", 20, 300, 1000)
        fjetlsf3_axis = hist.Bin("fjet_lsf3", r"Cand Jet LSF$_3$", 20, 0, 1)
        fjetmsd_axis = hist.Bin("fjet_msd", r"Cand Jet m$_{SD}$", 20, 20, 200)

        fjett41_axis = hist.Bin("fjet_t41", r"Cand Jet $\tau_{41}", 20, 0, 1.0)
        fjett42_axis = hist.Bin("fjet_t42", r"Cand Jet $\tau_{42}", 20, 0, 1.0)
        fjett31_axis = hist.Bin("fjet_t31", r"Cand Jet $\tau_{31}", 20, 0, 1.0)

        jetoppbtag_axis = hist.Bin("jet_oppbtag",r"Jet Opposite AK4 max b-tag", 20, 0, 1)

        fjetmmass_axis = hist.Bin("fjet_mmass",r"Jet - Lep Mass", 20, 0, 120)
        fjethmass_axis = hist.Bin("fjet_hmass",r"Jet + MET Mass", 25, 0, 500)

        flsjetpt_axis = hist.Bin("flsjet_pt", r"Lep Subtracted Jet $p_T$", 20, 300, 1000)
        flsjetmsd_axis = hist.Bin("flsjet_msd", r"Lep Subtracted Jet m$_{SD}$", 20, 0, 200)
        flsjetn2b1_axis = hist.Bin("flsjet_n2b1", r"Lep Subtracted Jet N$_{2}^{\beta = 1}$", 20, 0, 0.6)
        flsjetn3b1_axis = hist.Bin("flsjet_n3b1", r"Lep Subtracted Jet N$_{3}^{\beta = 1}$", 20, 0, 1.0)
        flsjett21_axis = hist.Bin("flsjet_t21", r"Lep Subtracted Jet $\tau_{21}", 20, 0, 1.0)
        flsjett32_axis = hist.Bin("flsjet_t32", r"Lep Subtracted Jet $\tau_{32}", 20, 0, 1.0)

        metpt_axis = hist.Bin("met_pt", r"MET $p_T$", 20, 0, 350)
        metphi_axis = hist.Bin("met_phi", r"MET $\phi$", 20, -3.5, 3.5)

        nmuons_axis = hist.Bin("nmuons", "N muons", 4, 0, 3)
        nelectrons_axis = hist.Bin("nelectrons", "N electrons", 4, 0, 3)

        #genweight_axis = hist.Bin("genweight", "gen weight", 20, 0, 500)
        #puweight_axis = hist.Bin("puweight", "pu weight", 20, 0, 500)

        leppt_axis = hist.Bin("lep_pt", r"Lepton $p_T$", 20, 0, 450)

        muonpt_axis = hist.Bin("muon_pt", r"Muon $p_T$", 20, 0, 450)
        muonmiso_axis = hist.Bin("muon_miso", r"Muon mini PF ISO (total)", 20, 0, 1)
        muonsip_axis = hist.Bin("muon_sip", r"Muon SIP", 20, 0, 22)
        muondz_axis = hist.Bin("muon_dz", r"Muon dz", 20, 0, 0.5)
        muondxy_axis = hist.Bin("muon_dxy", r"Muon dxy", 20, 0, 0.1)

        electronpt_axis = hist.Bin("electron_pt", r"Electron $p_T$", 20, 0, 450)
        electronmiso_axis = hist.Bin("electron_miso", r"Electron mini PF ISO (total)", 20, 0, 1)
        electronsip_axis = hist.Bin("electron_sip", r"Muon SIP", 20, 0,22)
        electrondz_axis = hist.Bin("electron_dz", r"Electron dz", 20, 0, 0.5)
        electrondxy_axis = hist.Bin("electron_dxy", r"Electron dxy", 20, 0, 0.1)

        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'
        hists['sumw'] = processor.defaultdict_accumulator(int)
        hists['cutflow'] = processor.defaultdict_accumulator(float)
        self._regions = regions.split(',')
        for region in self._regions:
            # hists['%s_fjetprop'%region] = hist.Hist("Events / GeV",
            #                                         dataset_axis,
            #                                         #fjetpt_axis,
            #                                         fjetmsd_axis,
            #                                         fjetlsf3_axis,
            #                                         #jetoppbtag_axis,
            #                                         #genflavor_axis,
            #                                     )
            # hists['%s_fjetextraprop'%region] =  hist.Hist("Events / GeV",
            #                                               dataset_axis,
            #                                               fjett41_axis,
            #                                               fjett42_axis,
            #                                               fjett31_axis,
            #                                           )
            # hists['%s_jetprop'%region] = hist.Hist("Events / GeV",
            #                                        dataset_axis,
            #                                        jetoppbtag_axis,
            #                                        #genflavor_axis,
            #                                   )
            hists['%s_fmmjetprop'%region] = hist.Hist("Events / GeV",
                                                      dataset_axis,
                                                      #fjetmmass_axis,
                                                      #fjethmass_axis,
                                                      fjetlsf3_axis,
                                                      fjetpt_axis,
                                                      leppt_axis,
                                                      genflavor_axis,
                                                  )
            hists['%s_fmmjetprop2'%region] = hist.Hist("Events / GeV",
                                                       dataset_axis,
                                                       fjetmmass_axis,                                                                                                           
                                                       fjetlsf3_axis,
                                                       genflavor_axis)
            # hists['%s_flsjetprop'%region] = hist.Hist("Events / GeV",
            #                                           dataset_axis,
            #                                           #flsjetpt_axis,
            #                                           flsjetmsd_axis,
            #                                           #flsjetn2b1_axis,
            #                                           #flsjetn3b1_axis,
            #                                           #flsjett21_axis,
            #                                           #flsjett32_axis
            #                                           genflavor_axis,
            #                                       )
            #hists['%s_metprop'%region] = hist.Hist("Events / GeV",
            #                                       dataset_axis,
            #                                       metpt_axis,
            #                                       metphi_axis)
            #hists['%s_weight'%region] = hist.Hist("Events / GeV",
            #                                      dataset_axis,
            #                                      genweight_axis,
            #                                      puweight_axis)
            # if self._channel=='muon':
            #     hists['%s_muonprop'%region] = hist.Hist("Events / GeV",
            #                                             dataset_axis,
            #                                             muonpt_axis,
            #                                             muonsip_axis,
            #                                             muonmiso_axis)
            #     hists['%s_muonextraprop'%region] = hist.Hist("Events / GeV",
            #                                                  dataset_axis,
            #                                                  nmuons_axis,
            #                                                  nelectrons_axis,
            #                                                  muondz_axis,
            #                                                  muondxy_axis)
            # else:
            #     hists['%s_electronprop'%region] = hist.Hist("Events / GeV",
            #                                                 dataset_axis,
            #                                                 electronpt_axis,
            #                                                 electronsip_axis,
            #                                                 electronmiso_axis)
            #     hists['%s_electronextraprop'%region] = hist.Hist("Events / GeV",
            #                                                      dataset_axis,
            #                                                      nmuons_axis,
            #                                                      nelectrons_axis,
            #                                                      electrondz_axis,
            #                                                      electrondxy_axis)


        self._accumulator = hists

        print(self._fjetptMIN,self._trigger,self._channel)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        dataset = df.metadata['dataset']
        isRealData = 'genWeight' not in df.columns
        output = self.accumulator.identity()
        selection = processor.PackedSelection()
        output = self.accumulator.identity()

        good = False
        goodMuon = (
            (df.Muon.pt > 27.)
            & (np.abs(df.Muon.eta) < 2.4)
            )
        nmuons = goodMuon.sum()

        goodElectron = (
            (df.Electron.pt > 30.)
            & (np.abs(df.Electron.eta) < 2.5)
            )
        nelectrons = goodElectron.sum()

        df.FatJet['msdcorr'] = corrected_msoftdrop(df.FatJet)

        goodFatJet = (
            (df.FatJet.pt > 300.)
            & (np.abs(df.FatJet.eta) < 2.4)
            & (df.FatJet.msdcorr > 10.)
            & (df.FatJet.isTight)
            )
        nfatjets = goodFatJet.sum()

        if self._channel == 'muon':
            good = (
                (nmuons >= 1)
                & (nfatjets >= 1)
            )
        else:
            good = (
                (nelectrons >= 1)
                & (nfatjets >= 1)
            )
        events = df[good]

        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()

        # trigger
        trigger = np.zeros(df.size, dtype='bool')
        for t in self._triggers[self._year+'_'+self._trigger]:
            try:
                trigger = trigger | df.HLT[t]
            except:
                warnings.warn("Missing trigger %s" % t, RuntimeWarning)
        selection.add('trigger', trigger[good])

        # Muons
        candidatemuon = events.Muon[:,0:1]
        nmuons = events.Muon.counts

        # Electrons
        candidateelectron = events.Electron[:,0:1]
        nelectrons = events.Electron.counts

        if self._channel == 'muon':
            candidatelep = candidatemuon
            selection.add('nootherlepton',(nelectrons == 0))
        else:
            candidatelep = candidateelectron
            selection.add('nootherlepton',(nmuons == 0))

        selection.add('iplepton', ((np.abs(candidatelep.dz) < 0.1)
                                   & (np.abs(candidatelep.dxy) < 0.05)).any())

        # FatJets
        ak8_lep_pair = candidatelep.cross(events.FatJet)
        ak8_lep_dR = ak8_lep_pair.i0.delta_r(ak8_lep_pair.i1)

        candidatejet = events.FatJet[ak8_lep_dR.argmin()]
        leadingjet = events.FatJet[:,0:1]
        
        ak8_lep_dR_closest = candidatelep.delta_r(candidatejet)

        selection.add('jetkin', (candidatejet.pt > self._fjetptMIN).any())
        selection.add('jetmsd', (candidatejet.msdcorr > 20).any())
        selection.add('LSF3medium', (candidatejet.lsf3 > 0.7).any())
        selection.add('LSF3tight', (candidatejet.lsf3 > 0.78).any())
        selection.add('lepnearjet', (ak8_lep_dR.min() < 1.5))
        selection.add('lepinjet', (ak8_lep_dR.min() < 0.8))

        # FatJet substracted Lepton
        # sj1_sj2_btagDeepB_pair = candidatejet.LSsubJet1btagDeepB.cross(candidatejet.LSsubJet2btagDeepB)
        # fls_btagDeepB_max = max(sj1_sj2_btagDeepB_pair.i0,sj1_sj2_btagDeepB_pair.i1)

        # Jets
        jets = events.Jet[
            (events.Jet.pt > 30.)
            & (abs(events.Jet.eta) < 2.5)
            & (events.Jet.isTight)
        ]
        ak4_ak8_pair = jets.cross(candidatejet, nested=True)
        ak4_ak8_dphi = abs(ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1))
        ak4_opposite = jets[(ak4_ak8_dphi > np.pi / 2).all()]
        ak4_away = jets[(ak4_ak8_dphi > 0.8).all()]

        selection.add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < self._btagWPs['med'][self._year])
        selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() < self._btagWPs['med'][self._year])

        # MET
        met = events.MET

        # MET eta with mass assumption
        mm = (candidatejet - candidatelep).mass2
        jmass = (mm>0)*np.sqrt(np.maximum(0, mm)) + (mm<0)*candidatejet.mass

        joffshell = jmass < 62.5
        massassumption = 80.*joffshell + (125 - 80.)*~joffshell
        x = massassumption**2/(2*candidatelep.pt*met.pt) + np.cos(candidatelep.phi - met.phi)
        met_eta = (
            (x < 1)*np.arcsinh(x*np.sinh(candidatelep.eta))
            + (x > 1)*(
                candidatelep.eta - np.sign(candidatelep.eta)*np.arccosh(candidatelep.eta)
                )
            )

        met_p4 = TLorentzVectorArray.from_ptetaphim(np.array([0.]),np.array([0.]),np.array([0.]),np.array([0.]))
        if met.size > 0:
            met_p4 = TLorentzVectorArray.from_ptetaphim(met.pt, met_eta.fillna(0.), met.phi, np.zeros(met.size)) 
            hmass = (candidatejet + met_p4).mass
        else:
            hmass = candidatejet.pt.zeros_like()

        # weights
        weights = processor.Weights(len(events),storeIndividual=True)
        if isRealData:
            genflavor = candidatejet.pt.zeros_like()
        else:
            try:
                weights.add('genweight', events.genWeight)
                add_pileup_weight(weights, events.Pileup.nPU, self._year)
                #print("Weight statistics: %r" % weights._weightStats)
            except:
                print('no gen weight')
            if 'TTTo' in dataset:
                genW,genW_idx = getParticles(events,24,['fromHardProcess', 'isLastCopy'])
                genb,genb_idx = getParticles(events,5,['fromHardProcess', 'isLastCopy'])
                genflavorW = matchedParticleFlavor(candidatelep, genW,'child', 0.4)
                genflavorb = matchedParticleFlavor(candidatelep, genb,'mom', 0.4)
                genflavor = getFlavor(genflavorW,genflavorb)
            elif( ('hww_2017' in dataset) or ('GluGluHToWW' in dataset)):
                genH,genH_idx = getParticles(events,25,['fromHardProcess', 'isLastCopy'])
                genW,genW_idx = getParticles(events,24,['fromHardProcess', 'isLastCopy'])
                genE,genE_idx = getParticles(events,11,['fromHardProcess','isFirstCopy'],1)
                genM,genM_idx = getParticles(events,13,['fromHardProcess','isFirstCopy'],1)
                genT,genT_idx = getParticles(events,15,['fromHardProcess','isFirstCopy'],1)
                genQ,genQ_idx = getParticles(events,[0,5],['fromHardProcess','isFirstCopy'])
                ishWW_qqelev = (genH.counts==1) & (genW.counts==2) & (genE.counts==1) & (genM.counts==0) & (genT.counts==0)
                ishWW_qqmuv = (genH.counts==1) & (genW.counts==2) & (genM.counts==1) & (genE.counts==0) & (genT.counts==0)
                ishWW_qqtauv = (genH.counts==1) & (genW.counts==2) & (genT.counts==1) & (genM.counts==0) & (genE.counts==0)
                ishWW_qqqq = (genH.counts==1) & (genW.counts==2) & (genQ.counts==4)  & (genM.counts==0) & (genE.counts==0)
                ishWW_muvelev = (genH.counts==1) & (genW.counts==2) & (genE.counts==1) & (genM.counts==1)
                ishWW_elevelev = (genH.counts==1) & (genW.counts==2) & (genE.counts==2) & (genM.counts==0)
                ishWW_tauvtauv = (genH.counts==1) & (genW.counts==2) & (genT.counts==2) & (genM.counts==0) & (genE.counts==0)
                ishWW_muvmuv = (genH.counts==1) & (genW.counts==2) & (genE.counts==0) & (genM.counts==2)
                genflavor = ((ishWW_qqelev) * 8 + (ishWW_qqmuv) * 9)
            else:
                genflavor = candidatejet.pt.zeros_like()

        # fill cutflow
        cutflow = ['trigger', 'jetkin', 'jetmsd', 'lepnearjet', 'lepinjet','antiak4btagMediumOppHem','nootherlepton','iplepton','LSF3medium', 'LSF3tight']
        allcuts = set()
        output['cutflow']['none'] += len(events)
        for cut in cutflow:
            allcuts.add(cut)
            output['cutflow'][cut] += selection.all(*allcuts).sum()

        regions = {}
        regions['presel'] = {'trigger','jetkin', 'jetmsd',  'lepinjet'}
        regions['antibtag'] = {'trigger','jetkin', 'jetmsd', 'antiak4btagMediumOppHem'}
        regions['noinjet'] = {'trigger','jetkin', 'jetmsd',  'lepnearjet', 'antiak4btagMediumOppHem'}
        regions['nolsf'] = {'trigger','jetkin', 'jetmsd',  'lepinjet', 'antiak4btagMediumOppHem'} #,'nootherlepton'}
        regions['lsf'] = {'trigger','jetkin', 'jetmsd',  'lepinjet', 'LSF3tight'}
        regions['bopp'] =  {'trigger','jetkin', 'jetmsd',  'lepinjet', 'LSF3tight','antiak4btagMediumOppHem'}
        regions['lep'] = {'trigger','jetkin', 'jetmsd',  'lepinjet', 'LSF3tight','antiak4btagMediumOppHem','nootherlepton','iplepton'}

        for region in self._regions:
            selections = regions[region]
            cut = selection.all(*selections)
            weight = weights.weight()[cut] 
            
            def normalize(val):
                try:
                    return val[cut].pad(1, clip=True).fillna(0).flatten() 
                except:
                    try:
                        return val[cut].flatten()
                    except:
                        return val[cut]

            # output['%s_fjetprop'%region].fill(#fjet_pt = normalize(candidatejet.pt),
            #                                   fjet_msd = normalize(candidatejet.msdcorr),
            #                                   fjet_lsf3 = normalize(candidatejet.lsf3),
            #                                   #jet_oppbtag = normalize(ak4_opposite.btagDeepB.max()),
            #                                   genflavor = normalize(genflavor),
            #                                   dataset=dataset,
            #                                   weight=weight
            # )
            # output['%s_fjetextraprop'%region].fill(fjet_t41 = normalize(candidatejet.tau4/candidatejet.tau1),
            #                                        fjet_t42 = normalize(candidatejet.tau4/candidatejet.tau2),
            #                                        fjet_t31 = normalize(candidatejet.tau3/candidatejet.tau1),
            #                                        dataset=dataset,
            #                                        weight=weight
            #                                    )
            # output['%s_jetprop'%region].fill(jet_oppbtag = normalize(ak4_opposite.btagDeepB.max()),
            #                                  genflavor = normalize(genflavor),
            #                                  dataset=dataset,
            #                                  weight=weight
            #                                 )
            output['%s_fmmjetprop'%region].fill(fjet_pt = normalize(candidatejet.pt),
                                                #fjet_mmass = normalize(jmass),
                                                #fjet_hmass = normalize(hmass),
                                                lep_pt = normalize(candidatelep.pt),
                                                fjet_lsf3 = normalize(candidatejet.lsf3),
                                                genflavor = normalize(genflavor),
                                                dataset=dataset,
                                                weight=weight)
            output['%s_fmmjetprop2'%region].fill(fjet_mmass = normalize(jmass),
                                                 fjet_lsf3 = normalize(candidatejet.lsf3),
                                                 genflavor = normalize(genflavor),
                                                 dataset=dataset,
                                                 weight=weight)
            # output['%s_flsjetprop'%region].fill(#flsjet_pt = normalize(candidatejet.LSpt),
            #                                     flsjet_msd = normalize(candidatejet.LSmsoftdrop),
            #                                     #flsjet_n2b1 = normalize(candidatejet.LSn2b1),
            #                                     #flsjet_n3b1 = normalize(candidatejet.LSn3b1),
            #                                     #flsjet_t21 = normalize(candidatejet.LStau2/candidatejet.LStau1),
            #                                     #flsjet_t32 = normalize(candidatejet.LStau3/candidatejet.LStau2),
            #                                     genflavor = normalize(genflavor),
            #                                     dataset=dataset,
            #                                     weight=weight)
            #output['%s_metprop'%region].fill(met_pt = normalize(met.pt),
            #                                 met_phi = normalize(met.phi),
            #                                 dataset=dataset,
            #                                 weight=weight)
            # output['%s_weight'%region].fill(puweight=weights.partial_weight(include=["pileup_weight"])[cut],
            #                                 genweight=weights.partial_weight(include=["genweight"])[cut],
            #                                 dataset=dataset,
            #                                 )
            # if self._channel=='muon':
            #     output['%s_muonprop'%region].fill(muon_pt = normalize(candidatemuon.pt),
            #                                       muon_miso = normalize(candidatemuon.miniPFRelIso_all),
            #                                       muon_sip = normalize(candidatemuon.sip3d),
            #                                       dataset=dataset,
            #                                       weight=weight)
            #     output['%s_muonextraprop'%region].fill(nmuons = normalize(nmuons),
            #                                            nelectrons = normalize(nelectrons),
            #                                            muon_dz = normalize(candidatemuon.dz),
            #                                            muon_dxy = normalize(candidatemuon.dxy),
            #                                            dataset=dataset,
            #                                            weight=weight)

            # else:
            #     output['%s_electronprop'%region].fill(electron_pt = normalize(candidateelectron.pt),
            #                                           electron_miso = normalize(candidateelectron.miniPFRelIso_all),
            #                                           electron_sip = normalize(candidateelectron.sip3d),
            #                                           dataset=dataset,
            #                                           weight=weight)
            #     output['%s_electronextraprop'%region].fill(nmuons = normalize(nmuons),
            #                                                nelectrons = normalize(nelectrons),
            #                                                electron_dz = normalize(candidateelectron.dz),
            #                                                electron_dxy = normalize(candidateelectron.dxy),
            #                                                dataset=dataset,
            #                                                weight=weight)

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
    parser.add_argument('--trigger', choices=['muon','vvlmuon','electron','vvlelectron','had','muonall','electronall'], default='muon', help='trigger selection')
    parser.add_argument('--channel', choices=['muon','electron'], default='muon', help='channel')
    parser.add_argument('--regions', default='presel', help='regions')
    args = parser.parse_args()

    print('hww args',args)
    processor_instance = HwwProcessor(year=args.year,trigger=args.trigger,channel=args.channel,regions=args.regions)

    util.save(processor_instance, 'hwwprocessor.coffea')
