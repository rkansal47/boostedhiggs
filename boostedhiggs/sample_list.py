import os
import math
from array import array
import argparse
import sys
import json
import shlex

from fileset import get2017files

fdirLpcBacon = 'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/'
fdirLpcBaconUL = 'root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/UL/'
fdataSamples = ['SingleMuon','SingleElectron','JetHT','MET']

cmssw = os.getenv('CMSSW_VERSION', 'CMSSW_10_2_6')
cmssw_base = os.getenv('CMSSW_BASE', 'CMSSW_10_2_6')

# submit nqueue jobs 
def write_condor(exe='runjob.sh', arguments = [], files = [],nqueue=1):
    job_name = exe.replace('.sh','.jdl')
    out = 'universe = vanilla\n'
    out += 'Executable = %s\n'%exe
    out += 'Should_Transfer_Files = YES\n'
    out += 'WhenToTransferOutput = ON_EXIT\n'
    out += 'Transfer_Input_Files = %s,%s\n'%(exe,','.join(files))
    out += 'Output = output/%s_$(Cluster)_$(Process).stdout\n'%job_name
    out += 'Error  = error/%s_$(Cluster)_$(Process).stderr\n'%job_name
    out += 'Log    = log/%s_$(Cluster)_$(Process).log\n'   %job_name
    out += 'Arguments = %s $(Process)\n'%(' '.join(arguments)) # last argument is number of splitting
    out += 'Queue %i\n'%nqueue
    with open(job_name, 'w') as f:
        f.write(out)
    os.system("condor_submit %s"%job_name)
    
# bash script: outdir must be in eos?
def write_bash(temp = 'runjob.sh', command = '', files = [],odir=''):
    out = '#!/bin/bash\n'
    # add arguments
    out += 'date\n'
    out += 'MAINDIR=`pwd`\n'
    out += 'ls\n'
    out += '#CMSSW from scratch (only need for root)\n'
    out += 'export CWD=${PWD}\n'
    out += 'export PATH=${PATH}:/cvmfs/cms.cern.ch/common\n'
    out += 'export CMS_PATH=/cvmfs/cms.cern.ch\n'
    out += 'export SCRAM_ARCH=slc7_amd64_gcc820\n'
    out += 'scramv1 project CMSSW %s\n'%cmssw
    out += 'cd %s/src\n'%cmssw
    out += 'eval `scramv1 runtime -sh` # cmsenv\n'
    out += 'cd - \n'
    out += command + '\n'
    out += 'echo "Inside $MAINDIR:"\n'
    out += 'ls -lrth \n'
    out += 'echo "DELETING..."\n'
    out += 'rm -rf %s\n'%cmssw
    out += 'rm -rf *.root\n'
    out += 'ls\n'
    out += 'date\n'
    with open(temp, 'w') as f:
        f.write(out)

# submit by list
def submitByList(iCommand,iJob,inSplit):
    lexe      = "runjob_%s.sh"%(iJob)
    files     = []
    arguments = []
    filesTransfer  = ['haddnano.py']
    write_bash(lexe,iCommand,files)
    os.system('mkdir -p output')
    os.system('mkdir -p error')
    os.system('mkdir -p log')
    write_condor(lexe,arguments,filesTransfer,inSplit)
    
# read directory on eos
def read(redirector,eosp,dryrun=False):
    try:
        #print('xrdfs %s ls %s > tmp.txt'%(redirector,eosp))
        os.system('xrdfs %s ls %s > tmp.txt'%(redirector,eosp))
        with open("tmp.txt") as f: lineList = f.read().splitlines()
        os.system('rm tmp.txt')
        return lineList
    except:
        return []

# find last directory and call read
def readXroot(redirector,eosp,dryrun=False):
    expandList = []
    lineList = read(redirector,eosp,dryrun)
    if any(('.root' in it and 'nano' in it) for it in lineList):
        new_f = [it.replace('/store',redirector+'/store') for it in lineList]
        expandList.extend(new_f)
    else:
        for f in lineList:
            loop = True
            if any(('.root' in it and 'nano' in it) for it in f): 
                loop = False
            while loop:
                if type(f) == type([]):
                    newf = read(redirector,f[0],dryrun)
                else:
                    newf = read(redirector,f,dryrun)
                if len(newf)==0: 
                    loop = False
                else: 
                    f = newf
                    if any('.root' in it for it in f):
                        loop = False
            new_f = [it.replace('/store',redirector+'/store') for it in f]
            expandList.extend(new_f)
    return expandList

# expand path into files
def expand(path,idir,midpath,dryrun=False):
    expandedPaths = []
    redirector = 'root://cmseos.fnal.gov/'
    if 'cern' in idir: redirector = 'root://eoscms.cern.ch/'
    eosp = idir.replace(redirector,'')+'/'+path+'/'+midpath
    new_content = readXroot(redirector,eosp,dryrun)
    expandedPaths.extend(new_content)
    return expandedPaths 

# expand path and sub dirs (if not hadd)
def expandPath(dicts,pathhadd='',dryrun=False):
    rdict = {}
    for sample,sampleDict in dicts.items():
        d={} 
        for subSname in sampleDict['samples']:
            if pathhadd != '':
                expandedPath = expand(subSname+'-hadd',pathhadd,'',dryrun)
            else:
                expandedPath = expand(subSname,sampleDict['dir'],sampleDict['path'],dryrun)
            if len(expandedPath)==0:
                print("ERROR: %s has no files"%subSname)
                print("Trying to expand path with %s"%sampleDict['path'])
            d[subSname] = expandedPath 
        rdict[sample] =  d
    return rdict

# build hadd command
def haddNano(output,listToHadd,outDir,sample,idir=''):
    toHadd = ' '.join([str(i0) for i0 in listToHadd])
    outHadd = '%s/%s'%(outDir,output)
    tmpdir = '%s%s/'%(idir,sample)
    command = 'mkdir -p %s \n'%(tmpdir)
    command += 'python haddnano.py %s/%s %s \n'%(tmpdir,output,toHadd)
    command += 'xrdcp %s/%s %s \n'%(tmpdir,output,outHadd)
    command += 'rm  %s/%s \n'%(tmpdir,output)
    return command,outHadd

# my fav utility
def slice_it(li, cols=2):
    start = 0
    for i in range(0,cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop

def main(args):
    year = args.year
    isample = args.sample
    idirLpcBacon = fdirLpcBacon.replace('/02/','/02/%s/'%year)
    idirLpcBaconUL = fdirLpcBaconUL.replace('/02/','/02/%s/'%year)
    dryrun = args.year

    lfiles = {}
    if year == '2017':
        lfiles = get2017files(idirLpcBacon,idirLpcBaconUL,isample)

    finaljson = {}
    lnewjson = {}
    ltmpjson = {}

    # get json of files
    if args.haddread:
        lnewjson = expandPath(lfiles,idirLpcBacon+'/hadd')
        lnewjson1 = expandPath(lfiles,idirLpcBaconUL+'/hadd')
        for key,value in lnewjson1.items():
            for key1,value1 in value.items():
                if( len(lnewjson[key][key1]) == 0):
                    lnewjson[key][key1] = value1
    else:
        ltmpjson = expandPath(lfiles,'',dryrun)

    # hadd list from files
    if args.hadd:
        pwd = os.getcwd()
        if args.condor:
            os.system('mkdir -p condor_hadd/')
            os.system('cp haddnano.py condor_hadd/')
        for sample,sampleDict in lfiles.items():
            lnewjson[sample] = {}
            for ds in sampleDict['samples']:
                lcommands = []

                ds_FilesToHadd = ltmpjson[sample][ds]
                ds_base = 'nano_mc_%s'%year
                if sample in fdataSamples: ds_base = ds_base.replace('mc','data')
                ds_outdir = sampleDict['dir'].replace('/pancakes-02_','_')+'/hadd/'+ds+'-hadd'
                os.system('mkdir -p %s'%ds_outdir.replace('root://cmseos.fnal.gov/','/eos/uscms/'))
                if len(ds_FilesToHadd)>sampleDict['tohadd']:
                    splitFiles = slice_it(ds_FilesToHadd,sampleDict['tohadd'])
                else:
                    splitFiles = [ds_FilesToHadd]

                lnewjson[sample][ds] = []
                for iL,iList in enumerate(splitFiles):
                    ds_output = ds_base+'_%i.root'%iL
                    idircondor = '/uscmst1b_scratch/lpc1/3DayLifetime/cmantill/'
                    if args.condor: idircondor = ''
                    command_hadd,ds_hadd = haddNano(ds_output,iList,ds_outdir,sample,idircondor)
                    lnewjson[sample][ds].append(ds_hadd)
                    lcommands.append(command_hadd)

                # hadd for real
                if args.condor:
                    os.chdir('condor_hadd/')
                    for i0,l in enumerate(lcommands):
                        submitByList(l,ds+'_%i'%i0,1)
                    os.chdir(pwd)
                else:
                    for l in lcommands:
                        for x in l.split('\n'):
                            os.system(x)
    
    # create json
    outname = 'hwwfiles_'+args.year
    if args.haddread or args.hadd:
        print('reading from hadd')
        finaljson[year] = lnewjson
        outname+='_hadd'
    else:
        print('reading from indiv dir')
        finaljson[year] = ltmpjson

    outf = open("data/%s.json"%outname,"w")
    outf.write((json.dumps(finaljson,indent=4)))
    for key,tfiles in sorted(finaljson.items()):
        print("list of samples used by %s =  "%key, sorted(tfiles.keys()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='List files')
    parser.add_argument('--year', default='2017', help='year')
    parser.add_argument('--sample', default='all', help='sample')
    parser.add_argument('--hadd', action='store_true', default=False, help='Hadd all files')
    parser.add_argument('--haddread', action='store_true', default=False, help='Read hadd files')
    parser.add_argument('--condor', action='store_true', default=False, help='condor')
    parser.add_argument('--dryrun', action='store_true', default=False, help='print commands')

    args = parser.parse_args()
    print(args)

    main(args)
