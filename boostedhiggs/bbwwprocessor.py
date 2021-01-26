import numpy as np
from coffea import processor, hist
import awkward1 
from boostedhiggs.corrections import (
    add_pileup_weight
)
from boostedhiggs.btag import BTagEfficiency, BTagCorrector
from coffea.nanoaod.methods.leptons import LeptonCommon
from boostedhiggs.proc.util import *
import warnings

class HHbbWW(processor.ProcessorABC):
    #print(awkward1.__version__)
    def __init__(self, year):
        self._year = year
        self._trigger = {
	    2016: {
                "e": [
                    "Ele27_WPTight_Gsf",
                    "Ele45_WPLoose_Gsf",
                    "Ele25_eta2p1_WPTight_Gsf",
                    "Ele115_CaloIdVT_GsfTrkIdT",
                    "Ele15_IsoVVL_PFHT350",
                    "Ele15_IsoVVVL_PFHT400",
                    "Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50",
                    "Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
                    ],
                "mu": [
                    "IsoMu24",
                    "IsoTkMu24",
                    "Mu50",
                    "TkMu50",
                    "Mu15_IsoVVVL_PFHT400",
                    "Mu15_IsoVVVL_PFHT350",
	        ],
            }
        }
        self._trigger = self._trigger[int(self._year)]

        self._accumulator = processor.dict_accumulator({
            'sumw': processor.defaultdict_accumulator(float),
            'cutflow': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('channel', 'Channel'),
                hist.Bin('cut', 'Cut index', 9, 0, 9),
            ),
            })
        
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):

        # get meta infos
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        n_events = len(events)
        selection = processor.PackedSelection()
        weights = processor.Weights(n_events)
        output = self.accumulator.identity()

        # weights
        if not isRealData:
            output['sumw'][dataset] += awkward1.sum(events.genWeight)
        
        # trigger
        triggers = {}
        for channel in ["e","mu"]:
            trigger = np.zeros(len(events), dtype='bool')
            for t in self._trigger[channel]:
                try:
                    trigger = trigger | events.HLT[t]
                except:
                    warnings.warn("Missing trigger %s" % t, RuntimeWarning)
            triggers[channel] = trigger
            
        # met filter
        met_filters = ["goodVertices",
                       "globalSuperTightHalo2016Filter",
                       "HBHENoiseFilter",
                       "HBHENoiseIsoFilter",
                       "EcalDeadCellTriggerPrimitiveFilter",
                       "BadPFMuonFilter",
                       ]
        met_filters_mask = np.ones(len(events), dtype='bool')
        for t in met_filters:
            met_filters_mask = met_filters_mask & events.Flag[t]
        selection.add("met_filter", awkward1.to_numpy(met_filters_mask))
        
        # load objects
        muons = events.Muon
        electrons = events.Electron
        jets = events.Jet
        fatjets = events.FatJet
        subjets = events.SubJet
        fatjetsLS = events.FatJetLS
        met = events.MET
        
        # muons
        goodmuon = (
            (muons.mediumId)
            & (muons.miniPFRelIso_all <= 0.2)
            & (muons.pt >= 27)
            & (abs(muons.eta) <= 2.4)
            & (abs(muons.dz) < 0.1)
            & (abs(muons.dxy) < 0.05)
            & (muons.sip3d < 4)
        )
        good_muons = muons[goodmuon]
        ngood_muons = awkward1.sum(goodmuon, axis=1)

        # electrons
        goodelectron = (
            (electrons.mvaFall17V2noIso_WP90)
            & (electrons.pt >= 30)
            & (abs(electrons.eta) <= 1.479)
            & (abs(electrons.dz) < 0.1)
            & (abs(electrons.dxy) < 0.05)
            & (electrons.sip3d < 4)
        )
        good_electrons = electrons[goodelectron]
        ngood_electrons = awkward1.sum(goodelectron, axis=1)
        
        # good leptons
        good_leptons = awkward1.concatenate([good_muons, good_electrons], axis=1)
        good_leptons = good_leptons[awkward1.argsort(good_leptons.pt)]
        
        # lepton candidate
        candidatelep = awkward1.firsts(good_leptons)
        
        # lepton channel selection
        selection.add("ch_e", awkward1.to_numpy((triggers["e"]) & (ngood_electrons==1) & (ngood_muons==0))) # not sure if need to require 0 muons or 0 electrons in the next line
        selection.add("ch_mu", awkward1.to_numpy((triggers["mu"]) & (ngood_electrons==0) & (ngood_muons==1)))
        
        # jets
        ht = awkward1.sum(jets[jets.pt > 30].pt,axis=1)
        selection.add("ht_400", awkward1.to_numpy(ht>=400))
        goodjet = (
            (jets.isTight)
            & (jets.pt > 30)
            & (abs(jets.eta) <= 2.5)
            )
        good_jets = jets[goodjet]

        # fat jets
        jID = "isTight"
        # TODO: add mass correction

        # a way to get the first two subjets
        # cart = awkward1.cartesian([fatjets, subjets], nested=True)
        # idxes = awkward1.pad_none(awkward1.argsort(cart['0'].delta_r(cart['1'])), 2, axis=2)
        # sj1 = subjets[idxes[:,:,0]]
        # sj2 = subjets[idxes[:,:,1]]
        
        good_fatjet = (
            (getattr(fatjets, jID))
            & (abs(fatjets.eta) <= 2.4)
            & (fatjets.pt > 50)
            & (fatjets.msoftdrop > 30)
            & (fatjets.msoftdrop < 210)
            #& (fatjets.pt.copy(content=fatjets.subjets.content.counts) == 2) # TODO: require 2 subjets?
            # this can probably be done w FatJet_subJetIdx1 or FatJet_subJetIdx2
            & (awkward1.all(fatjets.subjets.pt >= 20))
            & (awkward1.all(abs(fatjets.subjets.eta) <= 2.4))
        )
        good_fatjets = fatjets[good_fatjet]

        # hbb candidate
        mask_hbb = (
            (good_fatjets.pt > 200)
            & (good_fatjets.delta_r(candidatelep) > 2.0)
            )
        candidateHbb = awkward1.firsts(good_fatjets[mask_hbb])

        # b-tag #& (good_fatjets.particleNetMD_Xbb > 0.9)
        selection.add('hbb_btag',awkward1.to_numpy(candidateHbb.deepTagMD_ZHbbvsQCD >= 0.8)) # score would be larger for tight category (0.97)  
        
        # No AK4 b-tagged jets away from bb jet
        jets_HbbV = jets[good_jets.delta_r(candidateHbb) >= 1.2]
        selection.add('hbb_vetobtagaway',  awkward1.to_numpy(awkward1.max(jets_HbbV.btagDeepB, axis=1, mask_identity=False) > BTagEfficiency.btagWPs[self._year]['medium']))
        
        # fat jets Lepton Subtracted
        # wjj candidate
        mask_wjj = (
            (fatjetsLS.pt > 50)
            & (fatjetsLS.delta_r(candidatelep) > 1.2)
            # need to add 2 subjets w pt > 20 & eta<2.4
            # need to add ID?
            )
        candidateWjj = awkward1.firsts(fatjetsLS[mask_wjj][awkward1.argmin(fatjetsLS[mask_wjj].delta_r(candidatelep),axis=1,keepdims=True)])
        # add t2/t1 <= 0.75 (0.45 HP)
        selection.add('hww_mass',  awkward1.to_numpy(candidateWjj.mass >= 10))

        print('met ',met)
        # wjjlnu info
        #HSolverLiInfo  hwwInfoLi;
        # qqSDmass = candidateWjj.msoftdrop
        # hwwLi   = hSolverLi->minimize(candidatelep.p4(), met.p4(), wjjcand.p4(), qqSDmass, hwwInfoLi)
        #neutrino = hwwInfoLi.neutrino;
        #wlnu     = hwwInfoLi.wlnu;
        #wqq      = hwwInfoLi.wqqjet;
        #hWW      = hwwInfoLi.hWW;
        #wwDM     = PhysicsUtilities::deltaR( wlnu,wqq) * hWW.pt()/2.0;
        # add dlvqq <= 11 (2.5 HP)
               
        # in the meantime let's add the mass
        '''
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
        
        # hh system
        candidateHH = candidateWjj + met_p4 + candidateHbb
        selection.add('hh_mass', candidateHH.mass >= 700)
        selection.add('hh_centrality', candidateHH.pt/candidateHH.mass >= 0.3)
        '''
        
        channels = {"e": ["met_filter","ch_e","ht_400","hbb_btag","hbb_vetobtagaway","hww_mass"], #,"hh_mass","hh_centrality"],
                    "mu": ["met_filter","ch_mu","ht_400","hbb_btag","hbb_vetobtagaway","hww_mass"] #,"hh_mass","hh_centrality"],
                    }

        # need to add gen info
        
        if not isRealData:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
            
        for channel, cuts in channels.items():
            allcuts = set()
            output['cutflow'].fill(dataset=dataset, channel=channel, cut=0, weight=weights.weight())
            for i, cut in enumerate(cuts):
                allcuts.add(cut)
                cut = selection.all(*allcuts)
                output['cutflow'].fill(dataset=dataset, channel=channel, cut=i + 1, weight=weights.weight()[cut])

        return output
        
    def postprocess(self, accumulator):
        return accumulator
