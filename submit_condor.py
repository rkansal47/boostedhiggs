import json
import os,sys
import argparse

# submit nqueue jobs 
def write_condor(exe='runjob.sh', arguments = [], files = [],nqueue=1):
    job_name = exe.replace('.sh','.jdl')
    out = 'universe = vanilla\n'
    out += 'Executable = %s\n'%exe
    out += 'Should_Transfer_Files = YES\n'
    out += 'WhenToTransferOutput = ON_EXIT\n'
    out += 'Transfer_Input_Files = %s,%s\n'%(exe,','.join(files))
    out += 'request_memory = 3.3GB\n'
    out += 'Output = output/%s_$(Cluster)_$(Process).stdout\n'%job_name
    out += 'Error  = error/%s_$(Cluster)_$(Process).stderr\n'%job_name
    out += 'Log    = log/%s_$(Cluster)_$(Process).log\n'   %job_name
    out += 'Arguments = %s $(Process)\n'%(' '.join(arguments)) # last argument is number of splitting
    out += 'Queue %i\n'%nqueue
    with open(job_name, 'w') as f:
        f.write(out)
    os.system("condor_submit %s"%job_name)

# bash script
def write_bash(temp = 'runjob.sh', command = '', files = [],odir=''):
    out = '#!/bin/bash\n'
    # add arguments
    out += 'date\n'
    out += 'ls\n'
    for f in files:
        if '.tgz' in f or '.gz'  in f:
            out += 'tar -zxf %s \n'%f
    out += 'source coffeaenv/bin/activate \n'
    out += "python -c 'import coffea' \n"

    out += command + '\n'
    out += 'ls -lrth \n'
    out += 'xrdcp output*.coffea root://cmseos.fnal.gov//store/user/cmantill/hww/%s/ \n'%odir
    out += 'rm output*.coffea \n'
    out += 'date\n'
    with open(temp, 'w') as f:
        f.write(out)

def main(args):
    datasets = {}
    with open(args.fileset) as f:
        temp = json.load(f)
        for dsgroup,datasetlist in temp.items():
            if dsgroup != args.year: continue
            datasets = datasetlist

    mc_hww = {'qcd': 30,
              'tt': 60,
              'vv': 10,
              'st': 20,
              'zll': 10,
              'zll-ht200': 10,
              'zll-ht1200': 10,
              'wlnu-v1': 10,
              'wlnu-v2': 10,
              'wlnu-v3': 10,
              'wqq': 10,
              'zqq': 10,
              'hwwlnu': 2, 
              #'gghwwlnu': 1,
              'tt-had': 100,
              #'JetHT': 40,
              #'SingleMuon': 40,
              #'SingleElectron': 40,
              'hww_private': 30,
          }

    ds_torun = []
    for process, processds in datasets.items():
        if process not in mc_hww.keys(): continue
        if args.ds != 'all':
            if process not in args.ds.split(','): continue
        if args.channel=='muon': 
            if process=='SingleElectron': continue
            if process=='JetHT' and args.trigger=='muon': continue
            if process=='SingleMuon' and args.trigger=='had': continue
        if args.channel=='electron':
            if process=='SingleMuon': continue
            if process=='JetHT' and args.trigger=='electron': continue
            if process=='SingleElectron' and args.trigger=='had': continue
        for ds, flist in processds.items():
            ds_torun.append([ds,mc_hww[process]])

    print('ds_torun ',ds_torun)
    odir = 'hww%s_trigger%s_channel%s_region%s_%s'%(args.year,args.trigger,args.channel,args.regions,args.tag) 
    os.system('mkdir -p /eos/uscms/store/user/cmantill/hww/%s/'%odir)                                                                                                     
    os.system('mkdir -p condor/')                                                                                                                                               
    os.system('mkdir -p condor/error')                                                                                                                                           
    os.system('mkdir -p condor/output')                                                                                                                                           
    os.system('mkdir -p condor/log')
    os.system('tar -zcf boostedhiggs.tgz boostedhiggs/')
    cwd = os.getcwd()
    os.chdir('condor/')
    os.system('cp ../coffeaenv.tar.gz .')
    os.system('cp ../boostedhiggs.tgz .')
    exePython = 'hww_condor.py'
    os.system('cp ../%s .'%exePython)
    files = ['coffeaenv.tar.gz','boostedhiggs.tgz',exePython]

    for ds,nsplit in ds_torun:
        command='python %s --year %s --trigger %s --channel %s --regions %s --fileset %s --ds %s --nsplit %i --isplit ${1}'%(exePython,args.year,args.trigger,args.channel,args.regions,args.fileset,ds,nsplit)
        write_bash('runjob%s_%s.sh'%(args.jobtag,ds.replace('/','')), command=command,files=files,odir=odir)
        write_condor('runjob%s_%s.sh'%(args.jobtag,ds.replace('/','')),arguments = [], files =files,nqueue=nsplit)
        
    os.chdir(cwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Submit to condor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--trigger', choices=['muon','vvlmuon','electron','vvlelectron','had','muonall','electronall'], default='muon', help='trigger selection')
    parser.add_argument('--channel', choices=['muon','electron'], default='muon', help='channel')
    parser.add_argument('--regions', default='presel', help='regions')
    parser.add_argument('--fileset', default='boostedhiggs/data/hwwfiles_2017.json', help='fileset')
    parser.add_argument('--ds', default='all', help ='choose a dataset or run on all')
    parser.add_argument('--tag', default='', help ='output tag')
    parser.add_argument('--jobtag', default='', help ='job tag')
    args = parser.parse_args()
    
    main(args)
