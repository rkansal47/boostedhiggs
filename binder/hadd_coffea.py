from coffea import util, hist
import json
import os,sys
import argparse

def main(args):
    eosdir = args.eosdir
    idir = args.idir

    # read directory on eos
    lineList = []
    os.system('xrdfs root://cmseos.fnal.gov/ ls %s/%s > tmp_%s.txt'%(eosdir,idir,idir))
    with open('tmp_%s.txt'%idir) as f: lineList = f.read().splitlines()
    os.system('rm tmp_%s.txt'%idir) 

    flist = [ util.load('/eos/uscms/'+x) for x in lineList]
    for key in flist[0]:
      if isinstance(key, hist.Hist):
        for fi in range(1,len(flist)):
          flist[0][key].add(flist[fi][key])
      else:
        for fi in range(1,len(flist)):
          flist[0][key] = flist[0][key] + flist[fi][key]

    output = 'output_condor/%s_sum.coffea' % (idir)
    os.system('rm %s'%output)
    util.save(flist[0],output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW hadd coffea files in eos')
    parser.add_argument('--idir', default='', help='eos subdir')
    parser.add_argument('--eosdir', default='/store/user/cmantill/hww/', help='eos dir')
    args = parser.parse_args()

    main(args)

