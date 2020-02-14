from coffea import util, processor
from boostedhiggs import HwwProcessor
import json
import sys
import argparse

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
        files[ds] = {'files': processds[ds], 'treename': 'Events'}
        
    p = HwwProcessor(year=args.year,trigger=args.trigger)

    exe_config = {
        'flatten': True,
        'workers': 4,
        'savemetrics': True,
        'compression': 0,
        'nano': True,
    }
    
    output, metrics = processor.run_uproot_job(files, 'Events', p, processor.iterative_executor, exe_config)
    util.save(output, 'output_%s_condor.coffea'%ds)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW processor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--trigger', choices=['muon','electron','had'], default='muon', help='trigger selection')
    parser.add_argument('--fileset', default='../boostedhiggs/data/hwwfiles_2017_hadd.json', help='fileset')
    parser.add_argument('--ds', default='', help ='choose a dataset!')
    args = parser.parse_args()
    
    main(args)
