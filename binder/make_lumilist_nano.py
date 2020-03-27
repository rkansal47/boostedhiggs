#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import subprocess

import uproot
import numpy as np
from coffea import lumi_tools

def slice_it(li, cols=2):
    start = 0
    for i in range(0,cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop

year = '2017'
fileset = '../boostedhiggs/data/hwwfiles_%s.json'%year
alldatasets = {}
with open(fileset) as f:
     temp = json.load(f)
     for dsgroup,datasetlist in temp.items():
          if dsgroup != year: continue
          alldatasets = datasetlist

ds = alldatasets['JetHT']
samples = {}
for ids, flist in ds.items():
    samples[ids] = {'files': flist, 'treename': 'Events'}

slices=100


lumivalues = lumi_tools.LumiData("../boostedhiggs/data/lumi%s.csv.gz"%year)

def get_lumilist(dataset, filename, treename):
     file = uproot.open(filename)
     if treename not in file:
          print("Bad file:", filename)
          return dataset, lumi_tools.LumiList()
     tree = file[treename]
     if tree.numentries == 0:
          return dataset, lumi_tools.LumiList()
     run, lumi = tree["run"].array(), tree["luminosityBlock"].array()
     if run.size==0: return dataset, lumi_tools.LumiList()
     lumilist = lumi_tools.LumiList(run, lumi)
     return dataset, lumilist


dataset_lumi = {}
nworkers = 12

with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
     futures = set()
     print(samples.keys())
     for dataset, files in samples.items():
         splitFiles = slice_it(files['files'],slices)
         for iL,iList in enumerate(splitFiles):
             futures.update(executor.submit(get_lumilist, dataset, file, files['treename']) for file in iList)
     try:
         total = len(futures)
         processed = 0
         while len(futures) > 0:
             finished = set(job for job in futures if job.done())
             for job in finished:
                 dataset, accumulator = job.result()
                 if dataset in dataset_lumi:
                     dataset_lumi[dataset] += accumulator
                 else:
                     dataset_lumi[dataset] = accumulator
                 processed += 1
                 if processed % 10 == 0:
                     print("Processing: done with % 4d / % 4d files" % (processed, total))
             futures -= finished
             del finished
     except KeyboardInterrupt:
          print("Ok quitter")
          for job in futures: job.cancel()
     except:
          for job in futures: job.cancel()
          raise


print("dataset, lumi [/pb], lumisections, unique lumisections")
for ds, ll in dataset_lumi.items():
     lumi = lumivalues.get_lumi(ll.array)
     nunique = np.unique(ll.array, axis=0).shape[0]
     ntot = ll.array.shape[0]
     print("%50s %6.0f %6d %6d" % (ds, lumi, ntot, nunique))
     
with gzip.open("lumilist_pancakes_2017.pkl.gz", "wb") as fout:
     pickle.dump(dataset_lumi, fout)
