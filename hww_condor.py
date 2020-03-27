from coffea import util, processor
from boostedhiggs import HwwProcessor
import json
import sys
import argparse

def slice_it(li, cols=2):
    start = 0
    for i in range(0,cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop

def main(args):
    datasets = {}
    with open(args.fileset) as f:
        temp = json.load(f)
        for dsgroup,datasetlist in temp.items():
            if dsgroup != args.year: continue
            datasets = datasetlist

    ds = args.ds
    files = {}
    for process, processds in datasets.items():
        for ids, flist in processds.items():
            if ids != ds: continue
            if args.nsplit != 1:
                lSplit = slice_it(flist,args.nsplit)
                for iL,iList in enumerate(lSplit):
                    if iL == args.isplit:
                        files[ds] = {'files': iList, 'treename': 'Events'}
            else:
                files[ds] = {'files': flist, 'treename': 'Events'}
                
    p = HwwProcessor(year=args.year,trigger=args.trigger,channel=args.channel,regions=args.regions)

    exe_config = {
        'workers': 4,
        'savemetrics': True,
        'nano': True,
    }
    
    output, metrics = processor.run_uproot_job(files, 'Events', p, processor.iterative_executor, exe_config)
    util.save(output, 'output_%s_%iof%i_condor.coffea'%(ds.replace('/','-'),args.isplit,args.nsplit))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW processor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--trigger', choices=['muon','vvlmuon','electron','vvlelectron','had','muonall','electronall'], default='muon', help='trigger selection')
    parser.add_argument('--channel', choices=['muon','electron'], default='muon', help='channel')
    parser.add_argument('--regions', default='presel', help='regions')
    parser.add_argument('--fileset', default='../boostedhiggs/data/hwwfiles_2017.json', help='fileset')
    parser.add_argument('--ds', default='', help ='choose a dataset!')
    parser.add_argument("--isplit",type=int, default=0, help='split')
    parser.add_argument("--nsplit",type=int, default=1, help='number of jobs to split file')
  
    args = parser.parse_args()
    
    print(args)
    main(args)
