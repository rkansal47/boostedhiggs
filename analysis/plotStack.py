from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from coffea import hist
from coffea.util import load

import pickle
import gzip
import math

import argparse
from boostedhiggs import processmap
from hists_map import *

plt.style.use(hep.style.CMS)
plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': False,
        })

fill_opts = {
    'edgecolor': (0,0,0,0.3),
    'alpha': 0.8
    }
err_opts = {
    'label':'Stat. Unc.',
    'hatch':'///',
    'facecolor':'none',
    'edgecolor':(0,0,0,.5),
    'linewidth': 0
    }
data_err_opts = {
    'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
}

def drawStack(h,var_name,var_label,lumifb,sel,vars_cut,save=False,drawData=False):
    # make cuts
    exceptions = ['process', var_name]
    for var,val in vars_cut.items():
        exceptions.append(var)
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions])
    for var,val in vars_cut.items():
        if var!=var_name:
            x = x.integrate(var,slice(val[0],val[1]))
    if var_name in vars_cut.keys():
        x = x[:, vars_cut[var_name][0]:vars_cut[var_name][1]]

    # numbers
    #x[privsig].scale(0.3799) # xsec?                                                                                                                                            
    #x[privsig].scale(0.29) # higgs pT?                                       

    scalesig = 5000.
    scalehpt = 0.29
    hsig = x[privsig]
    hsig.scale(scalesig)
    hsig.scale(scalehpt)

    all_bkg = 0
    for key,val in x[bkg].values().items():
        all_bkg+=val.sum()
    
    all_sig = 0
    for key,val in x[privsig].values().items():
        all_sig +=val.sum()
    
    if all_bkg>0:
        print('allbkg %.2f, allsig %.2f, s/sqrt(b) %.4f'%(all_bkg,all_sig,all_sig/np.math.sqrt(all_bkg)))

    # save
    if save:
        out = "templates_%s.root"%var_name
        if os.path.exists(out):
            os.remove(out)
            
        fout = uproot.create(out)
        for ih,proc in enumerate(x.identifiers('process')):
            pass_template = (x.integrate('process', proc))
            name = "%s_pass" % (proc.name)
            if 'hww_private' in proc.name:
                print('scaling by hpt',proc)
                pass_template.scale(scalehpt)
            fout[name] = hist.export1d(pass_template)
        fout.close()

    # identifiers
    xaxis = var_name
    x.axis(xaxis).label = var_label
    for ih,hkey in enumerate(x.identifiers('process')):
        name = process_latex[hkey.name]
        if 'hww_private' in hkey.name:
            name += ' x %i'%scalesig
        x.identifiers('process')[ih].label = name

    x.axis('process').sorting = 'integral'

    fig,ax = plt.subplots(1,1, figsize=(8,8))
    hist.plot1d(x[bkg],
                overlay='process',ax=ax,
                clear=False,
                stack=True,
                fill_opts=fill_opts,
                error_opts=err_opts,
                )
    if drawData:
        hist.plot1d(x['data_obs_jetht'],overlay="process",ax=ax,clear=False,error_opts=data_err_opts)

    hist.plot1d(x[offsig],ax=ax,overlay='process',clear=False,error_opts={'color': 'aquamarine','linewidth':2})
    hist.plot1d(hsig,ax=ax,overlay='process',clear=False,error_opts={'color': 'greenyellow','linewidth':2})

    ax.autoscale(axis='x', tight=True)
    maxval = 0
    for newkey,val in x[bkg].values().items():
        maxval = max(maxval,np.amax(val))

    ax.set_ylim(0, None)
    ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    if 'triggermuon' in sel:
        leg = ax.legend(handles=old_handles,labels=old_labels,title='Muon trigger')
    elif 'triggerhad' in sel:
        leg = ax.legend(handles=old_handles,labels=old_labels,title='Hadronic trigger')
    elif 'triggermuonall' in sel:
        leg = ax.legend(handles=old_handles,labels=old_labels,title='All trigger (muon)')
    elif 'triggerelectronall' in sel:
        leg = ax.legend(handles=old_handles,labels=old_labels,title='All trigger (electron)')
    else:
        leg = ax.legend(handles=old_handles,labels=old_labels,title='All trigger')

    lumi = plt.text(1., 1., r"%i fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("stack_%s_lumi%i.pdf"%(var_name,lumifb))
    ax.set_ylim(0.01, 100*maxval)
    ax.set_yscale('log')   
    fig.savefig("stack_%s_lumi%i_log.pdf"%(var_name,lumifb))

def getPlots(args):
    lumifb = args.lumi
    hist_name = args.sel
    tag = args.tag+'_'+hist_name

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s'%args.hists)

    # properties                                                                                                                                                                                                               
    var_name = args.var
    var_label = var_properties[args.var]['var_label']
    vars_cut = {}
    if args.cuts is not None:
        for cut in args.cuts.split(","):
            label = cut.split(":")[0]
            vals_s = cut.split(":")[1].split("-")
            vars_cut[label] = [float(i) for i in vals_s]
    try:
        #print(hists_unmapped)
        h = hists_unmapped[hist_name]
    except:
        h = hists_unmapped

    # map to hists
    if isinstance(h, hist.Hist):
        h_mapped = processmap.apply(h)

    # normalize to lumi (except for data) 
    h_mapped[nodata].scale({p: lumifb for p in h_mapped.identifiers('process')}, axis="process")

    os.chdir(odir)
    drawStack(h_mapped,var_name,var_label,lumifb,hist_name,vars_cut,args.save)
    os.chdir(pwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",  help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",       help="tag")
    parser.add_argument('--lumi',       dest='lumi',     default=41.1,     help="lumi")
    parser.add_argument('--sel',        dest='sel',      default='',       help='hist name')
    parser.add_argument('--var',        dest='var',      default='',       help='var ')
    parser.add_argument('--cuts',       dest='cuts',     default=None,     help='cuts defined by [,]')
    parser.add_argument('--save',       dest='save',     action='store_true', default=False)
    args = parser.parse_args()

    getPlots(args)

