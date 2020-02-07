import os
import math
from array import array
import argparse
import sys
import json

#########################################
# To DO:
# 1. submit hadd to condor
#########################################

fdirLpcBacon = 'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/'
fdirLpcBaconUL = 'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/UL/'
fdataSamples = ['SingleMuon','SingleElectron','JetHT','MET']

def get2017files(idirLpcBacon,idirLpcBaconUL,sample='all'):

    tfiles = {
        'wqq': {'samples': ['WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
                            'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
                            'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8'],
                'dir': idirLpcBacon,
                'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                'tohadd': 20,
            },
        'wlnu-v1':{'samples': ['WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
                               'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
                               'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
                               'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
                               'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
                               'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
                           ],
                   'dir': idirLpcBacon,
                   'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                   'tohadd': 20,
               },
        'wlnu-v2': {'samples': ['WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
                            ],
                    'dir': idirLpcBacon,
                    'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v2',
                    'tohadd': 20,
                },
        'wlnu-v3': {'samples': ['WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                                'dir': idirLpcBacon,
                                'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v3',
                                'tohadd': 20,
                },
        'wlnu-inc': {'samples':['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
                            ],
                     'dir': idirLpcBaconUL,
                     'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
                     'tohadd': 20,
                 },
        'zll': {'samples':['DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
                           'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
                           'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
                           'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
                           'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',                                                              
                       ],
                'dir': idirLpcBacon,
                'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_v14-v2',
                'tohadd': 20,
            },
        'zll-ht200': {'samples':['DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                      'dir': idirLpcBacon,
                      'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v2',
                      'tohadd': 20,
                  },
        'zll-ht1200': {'samples': ['DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                       'dir': idirLpcBacon,
                       'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                       'tohadd': 20,
                  },
        'zll-inc': {'samples': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'],
                    'dir': idirLpcBaconUL,
                    'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
                    'tohadd': 20,
                },
        'zqq': {'samples': ['ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
                            'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
                            'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
                        ],
                'dir': idirLpcBacon,
                'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                'tohadd': 20,
            },
        'qcd': {'samples': ['QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            #'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
                            ],
                'dir': idirLpcBaconUL,
                'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
                'tohadd': 20,
            },
        'tt': {'samples': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
                           'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
                       ],
               'dir': idirLpcBaconUL,
               'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
               'tohadd': 20,
           },
        'tt-had':{'samples': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
                          ],
                  'dir': idirLpcBaconUL,
                  'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v4',
                  'tohadd': 20,
              },
        'st': {'samples': ['ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8',
                           'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8',
                           'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
                           'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
                           'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
                           'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
                       ],
               'dir': idirLpcBacon,
               'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
               'tohadd': 10,
           },
        'vv': {'samples': ['WW_TuneCP5_13TeV-pythia8',
                           'WZ_TuneCP5_13TeV-pythia8',
                           'ZZ_TuneCP5_13TeV-pythia8',
                       ],
               'dir': idirLpcBaconUL,
               'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
               'tohadd': 10,
           },
        'hwwlnu': {'samples': ['GluGluHToWWToLNuQQ_M125_TuneCP5_PSweight_13TeV-powheg2-jhugen727-pythia8'],
                   'dir': idirLpcBaconUL,
                   'path': 'pancakes-02_RunIISummer19UL17MiniAOD-106X_v6-v2',
                   'tohadd': 10,
               },
        'hh-ww':{'samples':['GluGluToHHTo2B2WToLNu2J_node_SM_TuneCP5_13TeV-madgraph_correctedcfg',
                        ],
                 'dir': idirLpcBacon,
                 'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                 'tohadd': 10,
             },
        'hh-vv':{'samples': ['BulkGravTohhTohVVhbb_narrow_M-600_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-800_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-900_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-1000_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-1200_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-1400_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-1600_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-1800_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-2000_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-2500_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-3000_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-3500_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-4000_TuneCP5_13TeV-madgraph-pythia8',
                             'BulkGravTohhTohVVhbb_narrow_M-4500_TuneCP5_13TeV-madgraph-pythia8',],
                 'dir': idirLpcBacon,
                 'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                 'tohadd': 2,
             },
        'hh-tautau':{'samples':['BulkGrav_hh_htatahbb_inclusive_narrow_M800_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M900_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M1000_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M1200_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M1400_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M1600_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M1800_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M2000_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M2500_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M3000_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M3500_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M4000_TuneCP5_13TeV-madgraph-pythia8',
                                'BulkGrav_hh_htatahbb_inclusive_narrow_M4500_TuneCP5_13TeV-madgraph-pythia8'
                            ],
                     'dir': idirLpcBacon,
                     'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                     'tohadd': 2,
                 },
        'zpqq-v1':{'samples':['VectorZPrimeGammaToQQGamma_flat_13TeV_madgraph_pythia8_TuneCP5',
                              'VectorZPrimeToQQ_M75_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M100_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M150_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M175_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M200_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M225_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M300_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M350_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M400_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M450_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M500_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              ],
                   'dir': idirLpcBacon,
                   'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1',
                   'tohadd': 1,
               },
        'zpqq-v2':{'samples':['VectorZPrimeToQQ_M50_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M125_pT300_TuneCP5_madgraph_pythia8_13TeV',
                              'VectorZPrimeToQQ_M250_pT300_TuneCP5_madgraph_pythia8_13TeV'],
                   'dir': idirLpcBacon,
                   'path': 'pancakes-02_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v2',
                   'tohadd': 1,
               },
        'SingleMuon': {'samples': ['SingleMuon/pancakes-02_Run2017B-09Aug2019_UL2017-v1',
                                   'SingleMuon/pancakes-02_Run2017C-09Aug2019_UL2017-v1',
                                   'SingleMuon/pancakes-02_Run2017D-09Aug2019_UL2017-v1',
                                   'SingleMuon/pancakes-02_Run2017E-09Aug2019_UL2017-v1',
                                   'SingleMuon/pancakes-02_Run2017F-09Aug2019_UL2017-v1',
                               ],
                       'dir':idirLpcBaconUL,
                       'path': '',
                       'tohadd': 20,
                   },
        'SingleElectron': {'samples':['SingleElectron/pancakes-02_Run2017B-09Aug2019_UL2017-v1',
                                      'SingleElectron/pancakes-02_Run2017C-09Aug2019_UL2017-v1',
                                      'SingleElectron/pancakes-02_Run2017D-09Aug2019_UL2017-v1',
                                      'SingleElectron/pancakes-02_Run2017E-09Aug2019_UL2017-v1',
                                      #'SingleElectron/pancakes-02_Run2017F-09Aug2019_UL2017-v1',
                                  ],
                           'dir':idirLpcBaconUL,
                           'path': '',
                           'tohadd': 20,
                       },
        'JetHT': {'samples':['JetHT/pancakes-02_Run2017B-09Aug2019_UL2017-v1',
                             'JetHT/pancakes-02_Run2017C-09Aug2019_UL2017-v1',
                             'JetHT/pancakes-02_Run2017D-09Aug2019_UL2017-v1',
                             'JetHT/pancakes-02_Run2017E-09Aug2019_UL2017-v1',
                             'JetHT/pancakes-02_Run2017F-09Aug2019_UL2017-v1',
                         ],
                  'dir':idirLpcBaconUL,
                  'path': '',
                  'tohadd': 20,
              },
        'MET': {'samples':['MET/pancakes-02_Run2017B-09Aug2019_UL2017-v1',
                           'MET/pancakes-02_Run2017C-09Aug2019_UL2017-v1',
                           'MET/pancakes-02_Run2017D-09Aug2019_UL2017-v1',
                           'MET/pancakes-02_Run2017E-09Aug2019_UL2017-v1',
                           'MET/pancakes-02_Run2017F-09Aug2019_UL2017-v1',
                       ],
                'dir':idirLpcBaconUL,
                'path': '',
                'tohadd': 20,
            },
    }

    if sample=='all':
        return tfiles
    else:
        newfiles = {}
        for s in sample.split(','):
            newfiles[s] = tfiles[s]
        return newfiles

def read(redirector,eosp):
    try:
        #print('xrdfs %s ls %s > tmp.txt'%(redirector,eosp))
        os.system('xrdfs %s ls %s > tmp.txt'%(redirector,eosp))
        with open("tmp.txt") as f: lineList = f.read().splitlines()
        return lineList
    except:
        return []

def readXroot(redirector,eosp):
    expandList = []
    lineList = read(redirector,eosp)
    if any(('.root' in it and 'NanoAOD' in it) for it in lineList):
        new_f = [it.replace('/store',redirector+'/store') for it in lineList]
        expandList.extend(new_f)
    elif any(('.root' in it and 'RunIIAutumn18MiniAOD-102X_v15-v1' in it) for it in lineList): # argh fix this next time
        new_f = [it.replace('/store',redirector+'/store') for it in lineList]                                              
        expandList.extend(new_f)
    else:
        for f in lineList:
            loop = True
            if any(('.root' in it and 'NanoAOD' in it) for it in f): 
                loop = False
            while loop:
                if type(f) == type([]):
                    newf = read(redirector,f[0])
                else:
                    newf = read(redirector,f)
                if len(newf)==0: 
                    loop = False
                else: 
                    f = newf
                    if any('.root' in it for it in f):
                        loop = False
            new_f = [it.replace('/store',redirector+'/store') for it in f]
            expandList.extend(new_f)
    return expandList

def expand(path,idir,midpath):
    expandedPaths = []
    redirector = 'root://cmseos.fnal.gov/'
    if 'cern' in idir: redirector = 'root://eoscms.cern.ch/'
    eosp = idir.replace(redirector,'')+'/'+path+'/'+midpath
    new_content = readXroot(redirector,eosp)
    expandedPaths.extend(new_content)
    return expandedPaths 

def expandPath(dicts,pathhadd=''):
    rdict = {}
    for sample,sampleDict in dicts.iteritems():
        d={} 
        for subSname in sampleDict['samples']:
            if pathhadd != '':
                expandedPath = expand(subSname+'-hadd',pathhadd,'')
            else:
                expandedPath = expand(subSname,sampleDict['dir'],sampleDict['path'])
            if len(expandedPath)==0:
                print "ERROR: %s has no files"%(subSname)
                print "Trying to expand path with %s"%sampleDict['path']
            d[subSname] = expandedPath 
        rdict[sample] =  d
    return rdict

def haddNano(output,listToHadd,outDir,sample):
    toHadd = ' '.join([str(i0) for i0 in listToHadd])
    outHadd = '%s/%s'%(outDir,output)
    os.system('mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/%s/'%sample)
    os.system('python haddnano.py /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/%s/%s %s'%(sample,output,toHadd))
    os.system('xrdcp /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/%s/%s %s'%(sample,output,outHadd))
    os.system('rm  /uscmst1b_scratch/lpc1/3DayLifetime/cmantill/%s/%s'%(sample,output))
    return outHadd

def slice_it(li, cols=2):
    start = 0
    for i in xrange(cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop

def main(args):
    year = args.year
    isample = args.sample
    idirLpcBacon = fdirLpcBacon.replace('/02/','/02/%s/'%year)
    idirLpcBaconUL = fdirLpcBaconUL.replace('/02/','/02/%s/'%year)

    lfiles = {}
    if year == '2017':
        lfiles = get2017files(idirLpcBacon,idirLpcBaconUL,isample)
        print(lfiles)

    finaljson = {}
    lnewjson = {}
    ltmpjson = {}

    if args.haddread:
        lnewjson = expandPath(lfiles,idirLpcBacon+'/hadd')
        lnewjson2 = expandPath(lfiles,idirLpcBaconUL+'/hadd')
        for key,item in lnewjson2.iteritems():
            lnewjson[key] = item
    else:
        ltmpjson = expandPath(lfiles)

    if args.hadd:
        for sample,sampleDict in lfiles.items():
            lnewjson[sample] = {}
            for ds in sampleDict['samples']:
                ds_FilesToHadd = ltmpjson[sample][ds]
                ds_base = 'nano_mc_%s'%year
                if sample in fdataSamples:
                    ds_base = ds_base.replace('mc','data')
                ds_outdir = sampleDict['dir'].replace('/pancakes-02_','_')+'/hadd/'+ds+'-hadd'
                os.system('mkdir -p %s'%ds_outdir.replace('root://cmseos.fnal.gov/','/eos/uscms/'))
                print(len(ds_FilesToHadd),sampleDict['tohadd'])
                if len(ds_FilesToHadd)>sampleDict['tohadd']:
                    splitFiles = slice_it(ds_FilesToHadd,sampleDict['tohadd'])
                else:
                    splitFiles = [ds_FilesToHadd]
                    print(splitFiles)
                lnewjson[sample][ds] = []
                for iL,iList in enumerate(splitFiles):
                    ds_output = ds_base+'_%i.root'%iL
                    ds_hadd = haddNano(ds_output,iList,ds_outdir,sample)
                    lnewjson[sample][ds].append(ds_hadd)
                #print(lnewjson[sample][ds])

    outname = 'hwwfiles_'+args.year
    if args.haddread or args.hadd:
        print('reading from hadd')
        finaljson[year] = lnewjson
        outname+='_hadd'
    else:
        print('reading from indiv dir')
        finaljson[year] = ltmpjson

    outf = open("metadata/%s.json"%outname,"w")
    outf.write((json.dumps(finaljson,indent=4)))
    for key,tfiles in sorted(finaljson.iteritems()):
        print "list of samples used by %s =  "%key, sorted(tfiles.keys())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='List files')
    parser.add_argument('--hadd', action='store_true', default=False, help='Hadd all files')
    parser.add_argument('--haddread', action='store_true', default=False, help='Read hadd files')
    parser.add_argument('--year', default='2017', help='year')
    parser.add_argument('--sample', default='all', help='sample')
    args = parser.parse_args()
    print(args)

    main(args)
