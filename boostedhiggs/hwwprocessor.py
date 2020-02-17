import numpy as np
from coffea import processor,hist,util
from uproot_methods import TLorentzVectorArray
from boostedhiggs.corrections import compiled
import warnings
import argparse

class HwwProcessor(processor.ProcessorABC):
    def __init__(self, year='2018', trigger='muon'):
        self._year = year
        self._corrections = compiled
        self._btagWPs = {
            'med': {
                '2016': 0.6321,
                '2017': 0.4941,
                '2018': 0.4184,
            },
        }

        self._trigger = trigger
        self._triggers = {
            '2016_had': [
                "PFHT800",
                "PFHT900",
                "AK8PFJet360_TrimMass30",
                "AK8PFHT700_TrimR0p1PT0p03Mass50",
                "PFHT650_WideJetMJJ950DEtaJJ1p5",
                "PFHT650_WideJetMJJ900DEtaJJ1p5",
                "AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",
                "PFJet450",
            ],
            '2017_had': [
                "AK8PFJet330_PFAK8BTagCSV_p17",
                "PFHT1050",
                "AK8PFJet400_TrimMass30",
                "AK8PFJet420_TrimMass30",
                "AK8PFHT800_TrimMass50",
                "PFJet500",
                "AK8PFJet500",
            ],
            '2018_had': [
                "AK8PFJet400_TrimMass30",
                "AK8PFJet420_TrimMass30",
                "AK8PFHT800_TrimMass50",
                "PFHT1050",
                "PFJet500",
                "AK8PFJet500",
                "AK8PFJet330_PFAK8BTagCSV_p17",
                "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",
            ],
            '2017_muon': [
                "Mu50","Mu55",
                "Mu15_IsoVVVL_PFHT450_PFMET50","Mu15_IsoVVVL_PFHT600",
            ],
            '2018_muon': [
                "Mu50","Mu55",
                "Mu15_IsoVVVL_PFHT450_PFMET50","Mu15_IsoVVVL_PFHT600",
            ],
            '2017_electron': [
                "Ele27_WPTight_Gsf","Ele40_WPTight_Gsf","Ele20_WPLoose_Gsf","Ele115_CaloIdVT_GsfTrkIdT"
                "Ele15_IsoVVVL_PFHT450_PFMET50","Ele15_IsoVVVL_PFHT600",
            ],
            '2018_electron': [
                "Ele27_WPTight_Gsf","Ele40_WPTight_Gsf","Ele20_WPLoose_Gsf","Ele115_CaloIdVT_GsfTrkIdT"
                "Ele15_IsoVVVL_PFHT450_PFMET50","Ele15_IsoVVVL_PFHT600",
            ],
            }

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        jetpt_axis = hist.Bin("jet_pt", r"Jet $p_T$", 20, 200, 1000)
        jetlsf3_axis = hist.Bin("jet_lsf3", r"Jet LSF$_3$", 20, 0, 1)
        jetmmass_axis = hist.Bin("jet_mmass",r"Jet - Lep Mass", 20, -100, 100)
        jetoppbtag_axis = hist.Bin("jet_oppbtag",r"Jet Opposite AK4 b-tag", 20, 0, 1)
        muonpt_axis = hist.Bin("muon_pt", r"Muon $p_T$", 20, 0, 400)
        muonmiso_axis = hist.Bin("muon_miso", r"Muon mini PF ISO (total)", 20, 0, 1)
        metpt_axis = hist.Bin("met_pt", r"MET $p_T$", 20, 0, 100)
        meteta_axis = hist.Bin("met_eta", r"MET $\eta$", 20, -4, 4)
        metphi_axis = hist.Bin("met_phi", r"MET $\phi$", 20, -4, 4)

        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'
        hists['sumw'] = processor.defaultdict_accumulator(int)
        hists['cutflow'] = processor.defaultdict_accumulator(float)        
        for region in ['trigger','presel','bopp','lepsel']:
            if region=='lepsel':
                hists['%s_fjetprop'%region] = hist.Hist("Events / GeV",
                                                        dataset_axis,
                                                        jetpt_axis,
                                                        jetlsf3_axis,
                                                        jetmmass_axis,
                                                        )
            hists['%s_trigprop'%region] = hist.Hist("Events / GeV",
                                                    dataset_axis,
                                                    jetpt_axis,
                                                    )
            hists['%s_jetprop'%region] = hist.Hist("Events / GeV",
                                                   dataset_axis,
                                                   jetoppbtag_axis)
            hists['%s_muonprop'%region] = hist.Hist("Events / GeV",
                                                    dataset_axis,
                                                    muonpt_axis,
                                                    muonmiso_axis)
            hists['%s_metprop'%region] = hist.Hist("Events / GeV",
                                                   dataset_axis,
                                                   metpt_axis,
                                                   meteta_axis,
                                                   metphi_axis,
                                                   )

        self._accumulator = hists

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        dataset = df.metadata['dataset']
        isRealData = 'genWeight' not in df.columns
        output = self.accumulator.identity()
        selection = processor.PackedSelection()
        output = self.accumulator.identity()

        # pre-selection
        good = (
            (df.Muon.counts >= 1)
            & (df.FatJet.counts >= 1)
            )
        events = df[good]

        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()

        # trigger
        trigger = np.zeros(df.size, dtype='bool')
        for t in self._triggers[self._year+'_'+self._trigger]:
            trigger = trigger | df.HLT[t]
        selection.add('trigger', trigger[good])

        # Muons
        goodMuon = (
            (events.Muon.pt > 27.)
            & (np.abs(events.Muon.eta) < 2.4)
            & (events.Muon.sip3d < 4)
            & (np.abs(events.Muon.dz) < 0.1)
            & (np.abs(events.Muon.dxy) < 0.05)
            #& (events.Muon.mvaId == 2)
            )
        nmuons = goodMuon.sum()

        # Electrons
        goodElectron = (
            (events.Electron.pt > 10)
            & (np.abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased & (1 << 2)).astype(bool)  # 2017V2 loose                                                                                                 
            )
        nelectrons = goodElectron.sum()

        # one Muon and zero Electrons   
        selection.add('onemuon', (nmuons == 1) & (nelectrons == 0))

        # select FatJet with muon closest
        goodFatJet = (
            (events.FatJet.pt > 300.)
            & (np.abs(events.FatJet.eta) < 2.4)
            & (events.FatJet.msoftdrop > 10.)
            & (events.FatJet.jetId & 2)
            )
        leadingmuon = events.Muon[goodMuon][:, 0:1]
        ak8_muon_pair = leadingmuon.cross(events.FatJet)
        ak8_muon_dR = ak8_muon_pair.i0.delta_r(ak8_muon_pair.i1)
        leadingjet = events.FatJet[goodFatJet][ak8_muon_dR.argmin()]
        #leadingjet = events.FatJet[goodFatJet][:,0:1]
        selection.add('jetkin', (leadingjet.pt > 300).any())

        # veto b-tag in opposite side of FatJet
        goodJet = (
            (events.Jet.pt > 30.)
            & (events.Jet.jetId & 2)
            )
        jets = events.Jet[goodJet]
        ak4_ak8_pair = jets.cross(leadingjet, nested=True)
        ak4_ak8_dphi = ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1)
        ak4_opposite = jets[(np.abs(ak4_ak8_dphi) > np.pi / 2).all()]
        selection.add('antibtag', ak4_opposite.btagDeepB.max() < self._btagWPs['med'][self._year])

        # lsf selection
        #selection.add('LSF3muinside', (leadingjet.muonIdx3SJ == 0).any())
        selection.add('LSF3medium', (leadingjet.lsf3>0.4).any())

        # building mass assumption
        mm = (leadingjet - leadingmuon).mass2
        jmass = (mm>0)*np.sqrt(np.maximum(0, mm)) + (mm<0)*leadingjet.mass
        
        #joffshell = leadingjet.msoftdrop < 125/2  # halfway point between offshell and onshell W
        met = events.MET
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
        cutflow = ['trigger', 'jetkin', 'onemuon', 'antibtag', 'LSF3medium']
        allcuts = set()
        output['cutflow']['none'] += len(events)
        for cut in cutflow:
            allcuts.add(cut)
            output['cutflow'][cut] += selection.all(*allcuts).sum()

        # weights
        weights = processor.Weights(len(events))
        if not isRealData:
            weights.add('genweight', events.genWeight)

        regions = {}
        regions['trigger'] = {'trigger'}
        regions['presel'] = {'trigger','jetkin', 'onemuon'}
        regions['bopp'] = {'trigger','jetkin', 'onemuon','antibtag'}
        regions['lepsel'] = {'trigger','jetkin', 'onemuon','antibtag','LSF3medium'}

        for histname, h in output.items():
            if not isinstance(h, hist.Hist):
                continue
            region = [r for r in regions.keys() if r in histname.split('_')]
            if len(region) == 1:
                region = region[0]
                selections = regions[region]
                cut = selection.all(*selections)
                #weight = weights.weight()[cut] 
                weight = cut

                def normalize(val):
                    #return val[cut].pad(1, clip=True).fillna(0).flatten()
                    return val.pad(1, clip=True).fillna(0).flatten()

                if '_fjetprop' in histname:
                    h.fill(jet_pt = normalize(leadingjet.pt),
                           jet_lsf3 = normalize(leadingjet.lsf3),
                           jet_mmass = normalize(jmass),
                           dataset=dataset,weight=weight)
                elif '_trigprop' in histname:
                    h.fill(jet_pt = normalize(leadingjet.pt),
                           dataset=dataset,weight=weight)
                elif '_muonprop' in histname:
                    h.fill(muon_pt = normalize(leadingmuon.pt),
                           muon_miso = normalize(leadingmuon.miniPFRelIso_all),
                           dataset=dataset,weight=weight)
                elif '_metprop' in histname:
                    if met.size > 0:
                        h.fill(met_pt = met.pt.flatten(),
                               met_eta = normalize(met_p4.eta),
                               met_phi = met.phi.flatten(),
                               dataset=dataset,weight=weight)
                elif '_jetprop' in histname:
                    h.fill(jet_oppbtag = ak4_opposite.btagDeepB.max().flatten(),
                           dataset=dataset,weight=weight)
            elif len(region) > 1:
                raise ValueError("Histogram '%s' has a name matching multiple region definitions: %r" % (histname, region))
            else:
                raise ValueError("Histogram '%s' does not fall into any region definitions." % (histname, ))

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
    parser.add_argument('--trigger', choices=['muon','electron','had'], default='muon', help='trigger selection')
    args = parser.parse_args()

    processor_instance = HwwProcessor(year=args.year,trigger=args.trigger)

    util.save(processor_instance, 'hwwprocessor.coffea')
