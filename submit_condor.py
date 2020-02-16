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
    out += 'source coffeaenv/bin/activate.csh \n'
    out += "python -c 'import coffea' \n"
    out += 'tcsh coffeaenv/bin/activate.csh \n'
    out += "python -c 'import coffea' \n"
    out += 'bash coffeaenv/bin/activate \n'
    out += "python -c 'import coffea' \n"
    out += 'source coffeaenv/bin/activate \n'
    out += "python -c 'import coffea' \n"

    out += command + '\n'
    out += 'ls -lrth \n'
    out += 'xrdcp output*.coffea root://cmseos.fnal.gov//store/user/cmantill/hww/%s/ \n'%odir
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

    mc_hww = [
        'qcd','tt','tt-had','vv','st',
        'zll','zll-ht200','zll-ht1200','wlnu-v1','wlnu-v2','wlnu-v3','wqq','zqq',
        'hwwlnu', #'gghwwlnu'
    ]

    ds_torun = []
    if args.ds == 'all':
        for process, processds in datasets.items():
            if process not in mc_hww: continue
            for ds, flist in processds.items():
                ds_torun.append(ds)
    else:
        ds_torun = args.ds.split(',')


    odir = 'hww_%s_%s'%(args.year,args.trigger)                                                                                                                                 
    os.system('mkdir -p /eos/uscms/store/user/cmantill/hww/%s/'%odir)                                                                                                     
    os.system('mkdir -p condor/')                                                                                                                                               
    os.system('mkdir -p condor/error')                                                                                                                                           
    os.system('mkdir -p condor/output')                                                                                                                                           
    os.system('mkdir -p condor/log')
    os.system('tar -zcf boostedhiggs.tgz boostedhiggs/')
    cwd = os.getcwd()
    os.chdir('condor/')
    print(os.getcwd())
    os.system('cp ../coffeaenv.tar.gz .')
    os.system('cp ../boostedhiggs.tgz .')
    exePython = 'hww_condor.py'
    os.system('cp ../binder/%s .'%exePython)
    files = ['coffeaenv.tar.gz','boostedhiggs.tgz',exePython]
    odir = 'hww_%s_%s'%(args.year,args.trigger)

    for ds in ds_torun:
        command='python %s --year %s --trigger %s --fileset %s --ds %s'%(exePython,args.year,args.trigger,args.fileset,ds)
        write_bash('runjob_%s.sh'%ds, command=command,files=files,odir=odir)
        write_condor('runjob_%s.sh'%ds,arguments = [], files =files,nqueue=1)
        
    os.chdir(cwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Submit to condor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--trigger', choices=['muon','electron','had'], default='muon', help='trigger selection')
    parser.add_argument('--fileset', default='boostedhiggs/data/hwwfiles_2017_hadd.json', help='fileset')
    parser.add_argument('--ds', default='all', help ='choose a dataset or run on all')
    args = parser.parse_args()
    
    main(args)
