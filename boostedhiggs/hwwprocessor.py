import numpy as np
from coffea import processor,hist
from uproot_methods import TLorentzVectorArray

class HwwProcessor(processor.ProcessorABC):
    def __init__(self, year='2018', trigger='muon'):
        self._year = year

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
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
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
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                'AK8PFJet330_PFAK8BTagCSV_p17',
                "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",
            ],
            '2018_muon': [
                'Mu50','Mu55'
            ],
            '2018_electron': [
                'Ele27_WPTight_Gsf','Ele40_WPTight_Gsf','Ele20_WPLoose_Gsf','Ele115_CaloIdVT_GsfTrkIdT'
            ],

        }

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        jetpt_axis = hist.Bin("jet_pt", r"Jet $p_T$", 20, 200, 1000)
        jetlsf3_axis = hist.Bin("jet_lsf3", r"Jet LSF$_3$", 20, 0, 1)
        jetmmass_axis = hist.Bin("jet_mmass",r"Jet - Lep Mass", 20, -100, 100)
        jethmass_axis = hist.Bin("jet_hmass",r"Higgs Jet + MET Mass", 20, 0, 1000)
        jetoppbtag_axis = hist.Bin("jet_oppbtag",r"Jet Opposite AK4 b-tag", 20, 0, 1)
        muonpt_axis = hist.Bin("muon_pt", r"Muon $p_T$", 20, 0, 400)
        muonmiso_axis = hist.Bin("muon_miso", r"Muon mini PF ISO (total)", 20, 0, 1)
        metpt_axis = hist.Bin("met_pt", r"MET $p_T$", 20, 0, 100)
        meteta_axis = hist.Bin("met_eta", r"MET $\eta$", 20, -4, 4)
        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'
        hists['cutflow'] = processor.defaultdict_accumulator(float)
        for region in ['presel','muinjet']:
            hists['%s_jet'%region] = hist.Hist("Events / GeV",
                                               dataset_axis,
                                               jetpt_axis,
                                               jetlsf3_axis,
                                               jetmmass_axis,
                                               jethmass_axis,
                                               jetoppbtag_axis)
            hists['%s_muonevt'%region] = hist.Hist("Events / GeV",
                                                   muonpt_axis,
                                                   muonmiso_axis,
                                                   metpt_axis,
                                                   meteta_axis
                                                   )


        self._accumulator = hists

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        dataset = df.metadata['dataset']
        isRealData = 'genWeight' not in df.columns
        isSignal = 'hww' in dataset
        output = self.accumulator.identity()
        selection = processor.PackedSelection()

        # select at least one jet and one muon ( this is Pre-Selection! )                                                                                                       
        good = (
            (df.Muon.counts >= 1)
            & (df.CustomAK8Puppi.counts >= 1)
            )
        events = df[good]
        
        # trigger
        trigger = np.ones(df.size, dtype='bool')
        for t in self._triggers[self._year+'_'+self._trigger]:
            trigger = trigger & df.HLT[t]
        selection.add('trigger', trigger[good])

        # muon selection
        goodmuon = (
            (events.Muon.pt > 10)
            & (np.abs(events.Muon.eta) < 2.4)
            & (events.Muon.sip3d < 4)
            & (np.abs(events.Muon.dz) < 0.1)
            & (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.mvaId == 2)
        )
        nmuons = goodmuon.sum()
        leadingmuon = events.Muon[goodmuon].pad(1, clip=True).flatten()

        # fatjet closest to lepton 
        mujet_dR = leadingmuon.delta_r(events.CustomAK8Puppi)
        mu_in_cone = mujet_dR.min() < 0.8 # this I am not sure we have to put as a selection...
        mujet_bestidx = mujet_dR.argmin()
        leadingjet_mu = events.CustomAK8Puppi[mujet_bestidx].pad(1, clip=True).flatten()

        selection.add('jetkin', (
                (leadingjet_mu.pt > 300)
                & (leadingjet_mu.eta < 2.4)
                & (leadingjet_mu.msoftdrop > 10.)
                ).any())
        selection.add('jetid', (leadingjet_mu.jetId & 2).any())  # tight id 

        # lepton inside jet?
        selection.add('muinside', mu_in_cone.astype(bool))
        selection.add('LSF3muinside', (leadingjet_mu.electronIdx3SJ == 0).any())
        selection.add('LSF3medium', (leadingjet_mu.lsf3>0.78).any())

        # veto b-tag in opposite side
        jets = events.Jet[
            (events.Jet.pt > 30.)
            & (events.Jet.jetId & 2)  # tight id
            ]
        ak4_ak8_pair = jets.cross(leadingjet_mu, nested=True)
        dphi = ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1)
        ak4_opposite = jets[(np.abs(dphi) > np.pi / 2).all()]
        selection.add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < self._btagWPs['med'][self._year])

        # b-tag in same side
        #print(leadingjet_mu.columns)
        #subjets = events.CustomAK8PuppiSubJet[:, leadingjet_mu.subJetIdx1]

        # final lepton selection
        nelectrons = (
            (events.Electron.pt > 10)
            & (np.abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased & (1 << 2)).astype(bool)  # 2017V2 loose                                                                                                    
        ).sum()
        selection.add('onemuon', (nmuons == 1) & (nelectrons == 0)) # should we veto taus?                                                                                                          
        selection.add('muonkin', (
            (leadingmuon.pt > 27.)
            & (np.abs(leadingmuon.eta) < 2.4)
            ))

        # building variables
        leadingjet_mu_p4 = TLorentzVectorArray.from_ptetaphim(leadingjet_mu.pt, leadingjet_mu.eta, leadingjet_mu.phi, leadingjet_mu.msoftdrop)
        leadingmuon_p4 = TLorentzVectorArray.from_ptetaphim(leadingmuon.pt,leadingmuon.eta,leadingmuon.phi, leadingmuon.mass)
        print(leadingjet_mu_p4,leadingmuon_p4)

        mm = (leadingjet_mu_p4 - leadingmuon_p4).mass2 
        jmass = (mm>0)*np.sqrt(np.maximum(0, mm)) + (mm<0)*leadingjet_mu.mass # (jet - lep).M  

        met = events.MET
        joffshell = jmass < 125/2  # halfway point between offshell and onshell W
        massassumption = 80.*joffshell + (125 - 80.)*~joffshell
        x = massassumption**2/(2*leadingmuon.pt*met.pt) + np.cos(leadingmuon.phi - met.phi)
        met_eta = (
            (x < 1)*np.arcsinh(x*np.sinh(leadingmuon.eta))
            + (x > 1)*(
                leadingmuon.eta - np.sign(leadingmuon.eta)*np.arccosh(leadingmuon.eta)
                )
            )
        print(met_eta)
        met_p4 = TLorentzVectorArray.from_ptetaphim(met.pt, met_eta.fillna(0.), met.phi, np.zeros(met.size))

        # filling missing columns
        df['jet_pt'] = leadingjet_mu.pt
        df['jet_lsf3'] = leadingjet_mu.lsf3
        df['jet_mmass'] = jmass
        df['jet_hmass'] = (met_p4 + leadingjet_mu).mass
        df['jet_oppbtag'] = ak4_opposite.btagDeepB.max()
        df['muon_pt'] = leadingmuon.pt
        df['muon_miso'] = leadingmuon.miniPFRelIso_all
        df['met_pt'] = met.rho
        df['met_eta'] = met_eta

        # fill cutflow
        cutflow = ['trigger', 'jetkin', 'jetid', 'antiak4btagMediumOppHem', 'onemuon', 'muonkin', 'muinside', 'LSF3muinside','LSF3muinside']
        allcuts = set()
        output['cutflow']['none'] += len(events)
        for cut in cutflow:
            allcuts.add(cut)
            output['cutflow'][cut] += selection.all(*allcuts).sum()

        weights = processor.Weights(len(events))
        if not isRealData:
            weights.add('genweight', events.genWeight)

        regions = {}
        regions['presel'] = {'trigger', 'jetkin', 'jetid', 'antiak4btagMediumOppHem', 'onemuon', 'muonkin'}
        regions['muinjet'] = {'trigger', 'jetkin', 'jetid', 'antiak4btagMediumOppHem', 'onemuon', 'muonkin', 'muinside', 'LSF3muinside','LSF3muinside'}

        for histname, h in output.items():
            if not isinstance(h, hist.Hist):
                continue
            if not all(k in df or k == 'systematic' for k in h.fields):
                print("Missing fields %r from %r" % (set(h.fields) - set(df.keys()), h))
                continue
            fields = {k: df[k] for k in h.fields if k in df}
            region = [r for r in regions.keys() if r in histname.split('_')]
            if len(region) == 1:
                region = region[0]
                cut = selection.all(*regions[region])
                h.fill(**fields, weight=cut)
            elif len(region) > 1:
                raise ValueError("Histogram '%s' has a name matching multiple region definitions: %r" % (histname, region))
            else:
                raise ValueError("Histogram '%s' does not fall into any region definitions." % (histname, ))

        return output

    def postprocess(self, accumulator):
        return accumulator
