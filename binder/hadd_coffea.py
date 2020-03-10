from coffea import util, hist
import json
import os,sys
import argparse

def addHist(flist,key):
    tmp = flist[0][key]
    if isinstance(key, hist.Hist):
        for fi in range(1,len(flist)):
            tmp.add(flist[fi][key])
    else:
        for fi in range(1,len(flist)):
            tmp = tmp + flist[fi][key]
    return tmp

def main(args):
    eosdir = args.eosdir
    idir = args.idir
    key = args.hist

    # read directory on eos
    lineList = []
    os.system('xrdfs root://cmseos.fnal.gov/ ls %s/%s > tmp_%s.txt'%(eosdir,idir,idir))
    with open('tmp_%s.txt'%idir) as f: lineList = f.read().splitlines()
    os.system('rm tmp_%s.txt'%idir) 

    flist = [ util.load('/eos/uscms/'+x) for x in lineList]
    os.system('mkdir -p output_condor/%s/'%idir)

    if key!='':
        print('adding ',key)
        out = addHist(flist,key)
        output = 'output_condor/%s/%s_sum.coffea' % (idir,key)
        os.system('rm %s'%output)
        util.save(out,output)
    else:
        out = {}
        for key in flist[0]:
            print('adding ',key)
            out[key] = addHist(flist,key)
            
        output = 'output_condor/%s_sum.coffea' % (idir)
        os.system('rm %s'%output)
        util.save(out,output)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted HWW hadd coffea files in eos')
    parser.add_argument('--idir', default='', help='eos subdir')
    parser.add_argument('--eosdir', default='/store/user/cmantill/hww/', help='eos dir')
    parser.add_argument('--hist', default='', help='add only this hist')
    args = parser.parse_args()

    main(args)

