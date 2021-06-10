import uproot4
import awkward1 as ak

import pickle

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import mplhep as hep
from matplotlib.patches import Rectangle

import coffea.hist as hist
import coffea
from coffea.lookup_tools import extractor
from coffea.nanoevents.methods.vector import PtEtaPhiELorentzVector

import numpy as np
from cycler import cycler

from tqdm import tqdm
import pandas

import re

import math

from sklearn import metrics

from os import listdir
from copy import copy, deepcopy
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
plt.style.use(hep.style.ROOT)
%matplotlib inline

plt.rcParams.update({'font.size': 16})
plt.style.use(hep.style.CMS)

evtDict = {'HHbbVV4q': [], 'QCD': [], 'tt': [], 'data': []}
hists = {}


def init_hists(vars, labels, binrs, fatJet=True, name=None, scale=True, bbsorted=False, wwsorted=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    lpf = ["Leading Jet", "Sub-Leading Jet", "3rd Leading Jet"] if fatJet else ["Gen Higgs 1", "Gen Higgs 2"]
    if bbsorted: lpf = ["Hbb Candidate Jet", "HWW Candidate Jet"]
    if wwsorted: lpf = ["Leading H4q Candidate Jet", "Sub-Leading H4q Candidate Jet"]
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    ename = "Events (Unit Normalized)" if scale else "Events"

    if type(vars) is list:
        if type(binrs[0]) is not list:
            temp = binrs
            binrs = []
            for i in range(len(vars)):
                binrs.append(temp)

        for i in range(len(vars)):
            bins = []
            for j in range(nbins):
                bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels[i]), *binrs[i]))

            hists[name + vars[i]] = hist.Hist(ename, hist.Cat("sample", "Sample"), *bins)
    else:
        bins = []
        for j in range(nbins):
            if type(binrs) == list:
                bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels), *binrs))
            else:
                bins.append(hist.Bin("{}{}".format(bpf, j + 1), r"{} {}".format(lpf[j], labels), binrs))

        hists[name + vars] = hist.Hist(ename, hist.Cat("sample", "Sample"), *bins)


def fill_hists(vars, fatJet=True, name=None, hh4v=True, scale=True, pevtDict=None, useEDWeights=True, data=False, ak15=False, jet2ak15=False, blinding=None, data_blinded_only=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf
    if pevtDict is None: pevtDict = evtDict
    if useEDWeights and pevtDict is not None:
        pweights = {}
        for s, evts in pevtDict.items():
            pweights[s] = evts.weight
    else: pweights = weights

    if type(vars) is list:
        for i in range(len(vars)):
            if not type(vars[i]) is dict:
                temp = vars[i]
                vars[i] = {}
                for s in pevtDict.keys():
                    vars[i][s] = temp

            for s, evts in pevtDict.items():
                if (fatJet or "HH" in s) and (s != "HH4V" or hh4v) and (s != 'data' or data):
                    kwargs = {}
                    for j in range(nbins):
                        if fatJet and (ak15 or (j == 1 and jet2ak15)): rnpf = "ak15" + npf
                        else: rnpf = npf
                        kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(rnpf, j + 1, vars[i][s])]
                    hists[name + vars[i]["HHbbVV4q"]].fill(sample=s, weight=pweights[s], **kwargs)

            if scale: hists[name + vars[i]["HHbbVV4q"]].scale(scale_factor, axis='sample')
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in pevtDict.keys():
                vars[s] = temp

        for s, evts in pevtDict.items():
            if (fatJet or "HH" in s) and (s != "HH4V" or hh4v):
                if fatJet and ak15: rnpf = "ak15" + npf
                else: rnpf = npf
                kwargs = {}
                for j in range(nbins):
                    kwargs["{}{}".format(bpf, j + 1)] = evts["{}{}{}".format(rnpf, j + 1, vars[s])]
                hists[name + vars["QCD"]].fill(sample=s, weight=pweights[s], **kwargs)

        if scale: hists[name + vars["QCD"]].scale(scale_factor, axis='sample')
        if blinding is not None:
            hists[name + vars["QCD"]].values()[('data', )][:, blinding[0]:blinding[1]] = 0
            if not data_blinded_only: hists[name + vars["QCD"]].values()[('QCD', )][:, blinding[0]:blinding[1]] = 0


colors_cycle = ['#e31a1c', '#33a02c', '#1f78b4', '#a6cee3', '#b2df8a', '#fb9a99', '#f1005f', '#d09d09']

fill_opts = {
    'edgecolor': (0, 0, 0, 0.3),
    'linewidth': 2,
    'alpha': 0.8
}

line_opts = {
    'linestyle': '-',
    'linewidth': 3,
    'alpha': 0.8
}

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}

regbg = re.compile('^(QCD|tt)$')
regsig = re.compile('^(HH4V.*|HHbbVV4q)$')
regnotdata = re.compile('^(?!data).*')
regmc = re.compile('^(QCD|tt|HHbbVV4q)$')
reghh4v = re.compile('HH4V.*')


def plot_hists(vars, fname, binrs, fatJet=True, name=None, hh4v=True, hh4b=False, sepsig=True, stackall=False, log=False, log_limits=None, lumilabel=None, data=False, ratio=True, rat_ylim=1, ak15=False, blinding=None, legend_loc=None, incl_signal=False):
    npf = "fatJet" if fatJet else "genHiggs"
    bpf = "jet" if fatJet else "h"
    nbins = 2 if fatJet else 2
    if name is None: name = npf

    if type(vars) is list:
        fig, axs = plt.subplots(len(vars), nbins, figsize=(nbins * 9, len(vars) * 9))
        for i in range(len(vars)):
            if not type(vars[i]) is dict:
                temp = vars[i]
                vars[i] = {}
                for s in evtDict.keys():
                    vars[i][s] = temp

            var = vars[i]["HHbbVV4q"]
            if type(binrs[0]) is not list:
                temp = binrs
                binrs = []
                for k in range(len(vars)):
                    binrs.append(temp)

            if fatJet:
                if sepsig:
                    rnpf = "ak15" + npf if ak15 else npf
                    if ak15:
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    else:
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2b, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                        _ = axs[i, 1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2W, bins=np.linspace(binrs[i][1], binrs[i][2], binrs[i][0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)


                    for j in range(nbins):
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        if hh4v: hist.plot1d(hists[name + var][reghh4v].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)

                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                        axs[i, j].set_ylim(0, ylim)
                else:
                    for j in range(nbins):
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], stack=True, clear=False, fill_opts=fill_opts)
                        axs[i, j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][regsig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[i, j], clear=False, line_opts=line_opts)
                        axs[i, j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                        axs[i, j].set_ylim(0, ylim)

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()
    else:
        if not type(vars) is dict:
            temp = vars
            vars = {}
            for s in evtDict.keys():
                vars[s] = temp

        var = vars["HHbbVV4q"]
        if data and ratio: fig, (axs, rax) = plt.subplots(2, nbins, figsize=(nbins * 9, 12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        else: fig, axs = plt.subplots(1, nbins, figsize=(nbins * 9, 9))

        if fatJet:
            if sepsig:
                rnpf = "ak15" + npf if ak15 else npf
                if ak15:
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * ak15fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * ak15fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                else:
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[0].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 1, var)] * fatJet1W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2b, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[1], label='HHbbVV4q - Hbb', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)
                    _ = axs[1].hist(evtDict["HHbbVV4q"]["{}{}{}".format(rnpf, 2, var)] * fatJet2W, bins=np.linspace(binrs[1], binrs[2], binrs[0] + 1), color=colors_cycle[5], label='HHbbVV4q - HVV', weights=weights["HHbbVV4q"] * scale_factor["HHbbVV4q"], histtype='step', **line_opts)

                for j in range(nbins):
                    ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                    axs[j].set_prop_cycle(cycler(color=colors_cycle))
                    if hh4v: hist.plot1d(hists[name + var][reghh4v].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)

                    axs[j].set_prop_cycle(cycler(color=colors_cycle[2:]))
                    hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts)
                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})
                    axs[j].set_ylim(0, ylim)
            else:
                for j in range(nbins):
                    if data:
                        # hist.plot1d(hists[name + var][regbg].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['tt', 'QCD'])
                        if stackall: hist.plot1d(hists[name + var][regnotdata].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbVV4q', 'tt', 'QCD'])
                        elif incl_signal:
                            ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                            axs[j].set_prop_cycle(cycler(color=colors_cycle[3::-1]))
                            hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts, order=['tt', 'QCD'])
                            axs[j].set_prop_cycle(cycler(color=colors_cycle))
                            hist.plot1d(hists[name + var][regsig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
                            axs[j].set_ylim(0, ylim)
                        else:
                            axs[j].set_prop_cycle(cycler(color=colors_cycle[3:]))
                            hist.plot1d(hists[name + var]['QCD'].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['QCD'])


                        hist.plot1d(hists[name + var]['data'].project("sample", bpf + str(j + 1)), ax=axs[j], clear=False, error_opts=data_err_opts)
                        ylim = np.max(list(hists[name + var].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_ylim(0, ylim)
                        # axs[j].ticklabel_format(useMathText=True)

                        if ratio:
                            axs[j].set_xlabel(None)
                            hist.plotratio(num=hists[name + var]['data'].project("sample", bpf + str(j + 1)).sum("sample"), denom=hists[name + var][regbg].project("sample", bpf + str(j + 1)).sum("sample"), ax=rax[j], error_opts=data_err_opts, unc='num')
                            rax[j].set_ylabel('data/MC bg')
                            rax[j].set_ylim(0, rat_ylim)
                            rax[j].grid()

                    elif stackall:
                        hist.plot1d(hists[name + var][regnotdata].project('sample', bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, stack=True, fill_opts=fill_opts, order=['HHbbVV4q', 'tt', 'QCD'])
                    else:
                        ylim = np.max(list(hists[name + var][regnotdata].project("sample", bpf + str(j + 1)).values().values())) * 1.1
                        axs[j].set_prop_cycle(cycler(color=colors_cycle[3::-1]))
                        hist.plot1d(hists[name + var][regbg].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], stack=True, clear=False, fill_opts=fill_opts, order=['tt', 'QCD'])
                        axs[j].set_prop_cycle(cycler(color=colors_cycle))
                        hist.plot1d(hists[name + var][regsig].project("sample", bpf + str(j + 1)), overlay='sample', ax=axs[j], clear=False, line_opts=line_opts)
                        axs[j].set_ylim(0, ylim)

                    axs[j].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16}, loc=legend_loc)

                    if log:
                        if log_limits is None: axs[j].set_ylim(1)
                        else: axs[j].set_ylim(log_limits)
                        axs[j].set_yscale('log')

                    if lumilabel is not None: hep.label.lumitext(str(lumilabel) + "fb$^{-1}$", ax=axs[j])

        plt.tight_layout(0.5)
        plt.savefig("figs/{}.pdf".format(fname), bbox_inches='tight')
        plt.show()


def ak815_comp_plot(vars, varsl, bins, ak8events, ak15events, name, qcd=True, lumilabel=40):
    for i in range(len(vars)):
        histname = "ak8_ak15_{}_comp".format(vars[i])
        hists[histname] = hist.Hist("Events",
                                        hist.Cat("sample", "Sample"),
                                        hist.Bin("jet1", r"Jet {}".format(varsl[i]), *bins[i]),
                                        hist.Bin("jet2", r"Jet {}".format(varsl[i]), *bins[i]),
                                        )

        for s in ["HHbbVV4q", "QCD"]:
            hists[histname].fill(sample=s + " AK8",
                                 jet1=ak8events[s]["fatJet1" + vars[i]],
                                 jet2=ak8events[s]["fatJet2" + vars[i]],
                                 weight=ak8events[s]["weight"]
                                 )
            hists[histname].fill(sample=s + " AK15",
                                 jet1=ak15events[s]["ak15fatJet1" + vars[i]],
                                 jet2=ak15events[s]["ak15fatJet2" + vars[i]],
                                 weight=ak15events[s]["weight"]
                                 )


    fig, axs = plt.subplots(len(vars), 2, figsize=(len(vars) * 9, 18))

    xlabs = ['Hbb', 'HWW']

    loc = [7, 2] if qcd else [1, 1]

    for i in range(2):
        for j in range(len(vars)):
            histname = "ak8_ak15_{}_comp".format(vars[j])
            # print(str(i) + ", " + str(j))
            hist.plot1d(hists[histname]["HHbbVV4q AK8"].project("sample", "jet" + str(i + 1)), clear=False, line_opts=line_opts, ax=axs[j, i])
            hist.plot1d(hists[histname]["HHbbVV4q AK15"].project("sample", "jet" + str(i + 1)), clear=False, line_opts=line_opts, ax=axs[j, i])
            # axs[j, i].hist(events_ak15bb_sorted[sample]["ak15fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='red', weights=events_ak15bb_sorted[sample]['weight'], label='AK15', linewidth=2)
            if qcd:
                axqcd = axs[j, i].twinx()
                axqcd.set_prop_cycle(cycler(color=colors_cycle[0:]))
                hist.plot1d(hists[histname]["QCD AK15"].project("sample", "jet" + str(i + 1)), clear=False, line_opts={"linestyle": "dotted", "linewidth": 5}, ax=axqcd)
                hist.plot1d(hists[histname]["QCD AK8"].project("sample", "jet" + str(i + 1)), clear=False, line_opts={"linestyle": "dotted", "linewidth": 5}, ax=axqcd)
                axqcd.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16}, loc=1)
                axqcd.set_ylabel("QCD Events")
            # axqcd.hist(events_ak8bb_sorted['QCD']["fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='green', weights=events_ak8bb_sorted['QCD']['weight'], label='AK8', linewidth=3, linestyle='dashed')
            # axqcd.hist(events_ak15bb_sorted['QCD']["ak15fatJet{}{}".format(i + 1, vars[j])], bins=np.linspace(bins[j][1], bins[j][2], bins[j][0] + 1), histtype='step', color='red', weights=events_ak15bb_sorted['QCD']['weight'], label='AK15', linewidth=3, linestyle='dashed')
            axs[j, i].legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16}, loc=loc[j])
            axs[j, i].set_ylabel("HHbbVV4q Events")
            axs[j, i].set_xlabel("{} Candidate Jet {}".format(xlabs[i], varsl[j]))
            axs[j, i].set_xlim([bins[j][1], bins[j][2]])
            # hep.label.lumitext(, ax=axs[j, i])

            axs[j, i].text(
                x=0.15,
                y=1.005,
                s=str(lumilabel) + "fb$^{-1}$",
                transform=axs[j, i].transAxes,
                ha="right",
                va="bottom",
                fontsize=plt.rcParams["font.size"] * 0.95,
                fontweight="normal",
                # fontname=fontname,
            )

    plt.tight_layout(pad=2.0)
    plt.savefig("figs/{}.pdf".format(name))
    plt.show()


def tagger2d_plots(histname, ak8events, ak15events, vmaxhh=0.01, vmaxqcd=250, name=None, lumilabel=40):
    if name is None: name = histname
    titles = ['HHbbVV4q bb1 vs VV2 tagger', 'QCD bb1 vs VV2 tagger', 'QCD bb1 vs bb2 tagger']

    hists[histname] = hist.Hist("Events",
                                    hist.Cat("sample", "Sample"),
                                    hist.Bin("jet1bb", r"Jet1 Tagger Score", 20, 0.9, 1),
                                    hist.Bin("jet2WW", r"Jet2 Tagger Score", 20, 0.9, 1),
                                    )


    s = "HHbbVV4q"
    hists[histname].fill(sample="AK8 " + titles[0],
                         jet1bb = ak8events[s]["fatJet1PNetXbb"],
                         jet2WW = ak8events[s]["fatJet2PNetHqqqq"],
                         weight = ak8events[s]["weight"],
                        )

    s = 'QCD'
    hists[histname].fill(sample="AK8 " + titles[1],
                                 jet1bb = ak8events[s]["fatJet1PNetXbb"],
                                 jet2WW = ak8events[s]["fatJet2PNetHqqqq"],
                                 weight = ak8events[s]["weight"],
                                )

    hists[histname].fill(sample="AK8 " + titles[2],
                                 jet1bb = ak8events[s]["fatJet1PNetXbb"],
                                 jet2WW = ak8events[s]["fatJet2PNetXbb"],
                                 weight = ak8events[s]["weight"],
                                )


    s = "HHbbVV4q"
    hists[histname].fill(sample="AK15 " + titles[0],
                         jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                         jet2WW = ak15events[s]["ak15fatJet2PNetHqqqq"],
                         weight = ak15events[s]["weight"],
                        )

    s = 'QCD'
    hists[histname].fill(sample="AK15 " + titles[1],
                                 jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                                 jet2WW = ak15events[s]["ak15fatJet2PNetHqqqq"],
                                 weight = ak15events[s]["weight"],
                                )

    hists[histname].fill(sample="AK15 " + titles[2],
                                 jet1bb = ak15events[s]["ak15fatJet1PNetMDXbb"],
                                 jet2WW = ak15events[s]["ak15fatJet2PNetMDXbb"],
                                 weight = ak15events[s]["weight"],
                                )

    patch_opts = {
        'cmap': 'jet',
        'vmin': 0,
        'vmax': vmaxhh,
    }

    samps = ["AK8 ", "AK15 "]

    fig, axs = plt.subplots(2, 3, figsize=(3*9, 2*8))
    for j in range(3):
        if j == 1:
            patch_opts['vmax'] = vmaxqcd
            # patch_opts['vmin'] = 1
            # patch_opts['norm'] = mpl.colors.LogNorm()
        for i in range(2):
            hist.plot2d(hists[histname][samps[i] + titles[j]].sum('sample'), 'jet1bb', ax=axs[i][j], patch_opts=patch_opts)
            axs[i][j].set_title(samps[i] + titles[j], size=24)

            axs[i][j].text(
                x=1.2,
                y=1.01,
                s=str(lumilabel) + "fb$^{-1}$",
                transform=axs[i][j].transAxes,
                ha="right",
                va="bottom",
                fontsize=plt.rcParams["font.size"] * 0.95,
                fontweight="normal",
                # fontname=fontname,
            )

    plt.tight_layout(pad=1)
    plt.savefig("figs/{}.pdf".format(name), bbox_inches='tight')
    plt.show()



def single_roc_curves(tagger, tagger_name, name, jet, events, refill=True, ak15=True):
    jk = 'jet' + str(jet)

    init_hists(ak8tagger, tagger_name, disc_bin, name=name + "ak8", fatJet=True)
    fill_hists(ak8tagger, pevtDict=ak8events, name=name + "ak8", fatJet=True, scale=False, hh4v=False, ak15=False)

    disc_vals_ak8 = {'sig': [], 'bg': []}
    for s in ak8events.keys():
        if s == 'QCD' or s == 'tt':
            disc_vals_ak8['bg'].append(hists[name + "ak8" + ak8tagger].project("sample", jk).values()[(s, )])
        if s == 'HHbbVV4q':
            disc_vals_ak8['sig'].append(hists[name + "ak8" + ak8tagger].project("sample", jk).values()[(s, )])


    disc_vals_ak8['sig'] = np.sum(np.array(disc_vals_ak8['sig']), axis=0)
    disc_vals_ak8['bg'] = np.sum(np.array(disc_vals_ak8['bg']), axis=0)
    disc_vals_ak8['tpr'] = np.cumsum(disc_vals_ak8['sig'][::-1])[::-1] / np.sum(np.array(disc_vals_ak8['sig']))
    disc_vals_ak8['fpr'] = np.cumsum(disc_vals_ak8['bg'][::-1])[::-1] / np.sum(np.array(disc_vals_ak8['bg']))
    disc_vals_ak8['tpr'] = np.append(disc_vals_ak8['tpr'], 0)
    disc_vals_ak8['fpr'] = np.append(disc_vals_ak8['fpr'], 0)

    tpr_aves = (disc_vals_ak8['tpr'][:-1] + disc_vals_ak8['tpr'][1:]) / 2
    fpr_diffs = (disc_vals_ak8['fpr'][:-1] - disc_vals_ak8['fpr'][1:])
    aucak8 = np.sum(tpr_aves * fpr_diffs)


    init_hists(ak15tagger, tagger_name, disc_bin, name=name + "ak15", fatJet=True)
    fill_hists(ak15tagger, pevtDict=ak15events, name=name + "ak15", fatJet=True, scale=False, hh4v=False, ak15=True)

    disc_vals_ak15 = {'sig': [], 'bg': []}
    for s in ak15events.keys():
        if s == 'QCD' or s == 'tt':
            disc_vals_ak15['bg'].append(hists[name + "ak15" + ak15tagger].project("sample", jk).values()[(s, )])
        if s == 'HHbbVV4q':
            disc_vals_ak15['sig'].append(hists[name + "ak15" + ak15tagger].project("sample", jk).values()[(s, )])


    disc_vals_ak15['sig'] = np.sum(np.array(disc_vals_ak15['sig']), axis=0)
    disc_vals_ak15['bg'] = np.sum(np.array(disc_vals_ak15['bg']), axis=0)
    disc_vals_ak15['tpr'] = np.cumsum(disc_vals_ak15['sig'][::-1])[::-1] / np.sum(np.array(disc_vals_ak15['sig']))
    disc_vals_ak15['fpr'] = np.cumsum(disc_vals_ak15['bg'][::-1])[::-1] / np.sum(np.array(disc_vals_ak15['bg']))
    disc_vals_ak15['tpr'] = np.append(disc_vals_ak15['tpr'], 0)
    disc_vals_ak15['fpr'] = np.append(disc_vals_ak15['fpr'], 0)

    tpr_aves = (disc_vals_ak15['tpr'][:-1] + disc_vals_ak15['tpr'][1:]) / 2
    fpr_diffs = (disc_vals_ak15['fpr'][:-1] - disc_vals_ak15['fpr'][1:])
    aucak15 = np.sum(tpr_aves * fpr_diffs)


    if ak8: plt.plot(disc_vals_ak8['fpr'], disc_vals_ak8['tpr'], label='AK8 AUC = {:.2f}'.format(aucak8))
    if ak15: plt.plot(disc_vals_ak15['fpr'], disc_vals_ak15['tpr'], label='AK15 AUC = {:.2f}'.format(aucak15))

    plt.title(f'{tagger_name} ROC Curve')
    plt.xlabel('Background Efficiency')
    plt.ylabel('Signal Efficiency')
    plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

    plt.tight_layout(0.5)
    plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
    plt.savefig(f"figs/{name}.pdf", bbox_inches='tight')
    plt.show()


    if ak8: plt.semilogx(disc_vals_ak8['fpr'], disc_vals_ak8['tpr'], label='AK8 AUC = {:.2f}'.format(aucak8))
    if ak15: plt.semilogx(disc_vals_ak15['fpr'], disc_vals_ak15['tpr'], label='AK15 AUC = {:.2f}'.format(aucak15))

    plt.title(f'{tagger_name} Semilog ROC Curve')
    plt.xlabel('Background Efficiency')
    plt.ylabel('Signal Efficiency')
    plt.xlim(0.0001, 1)
    plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

    plt.tight_layout(0.5)
    plt.savefig(f"figs/{name}_semilog.pdf", bbox_inches='tight')
    plt.show()




def roc_curves(disc_bin, ak8tagger, ak15tagger, tagger_name, name, jet, ak8events, ak15events, refill=True, ak8=True, ak15=True):
    jk = 'jet' + str(jet)

    init_hists(ak8tagger, tagger_name, disc_bin, name=name + "ak8", fatJet=True)
    fill_hists(ak8tagger, pevtDict=ak8events, name=name + "ak8", fatJet=True, scale=False, hh4v=False, ak15=False)

    disc_vals_ak8 = {'sig': [], 'bg': []}
    for s in ak8events.keys():
        if s == 'QCD' or s == 'tt':
            disc_vals_ak8['bg'].append(hists[name + "ak8" + ak8tagger].project("sample", jk).values()[(s, )])
        if s == 'HHbbVV4q':
            disc_vals_ak8['sig'].append(hists[name + "ak8" + ak8tagger].project("sample", jk).values()[(s, )])


    disc_vals_ak8['sig'] = np.sum(np.array(disc_vals_ak8['sig']), axis=0)
    disc_vals_ak8['bg'] = np.sum(np.array(disc_vals_ak8['bg']), axis=0)
    disc_vals_ak8['tpr'] = np.cumsum(disc_vals_ak8['sig'][::-1])[::-1] / np.sum(np.array(disc_vals_ak8['sig']))
    disc_vals_ak8['fpr'] = np.cumsum(disc_vals_ak8['bg'][::-1])[::-1] / np.sum(np.array(disc_vals_ak8['bg']))
    disc_vals_ak8['tpr'] = np.append(disc_vals_ak8['tpr'], 0)
    disc_vals_ak8['fpr'] = np.append(disc_vals_ak8['fpr'], 0)

    tpr_aves = (disc_vals_ak8['tpr'][:-1] + disc_vals_ak8['tpr'][1:]) / 2
    fpr_diffs = (disc_vals_ak8['fpr'][:-1] - disc_vals_ak8['fpr'][1:])
    aucak8 = np.sum(tpr_aves * fpr_diffs)


    init_hists(ak15tagger, tagger_name, disc_bin, name=name + "ak15", fatJet=True)
    fill_hists(ak15tagger, pevtDict=ak15events, name=name + "ak15", fatJet=True, scale=False, hh4v=False, ak15=True)

    disc_vals_ak15 = {'sig': [], 'bg': []}
    for s in ak15events.keys():
        if s == 'QCD' or s == 'tt':
            disc_vals_ak15['bg'].append(hists[name + "ak15" + ak15tagger].project("sample", jk).values()[(s, )])
        if s == 'HHbbVV4q':
            disc_vals_ak15['sig'].append(hists[name + "ak15" + ak15tagger].project("sample", jk).values()[(s, )])


    disc_vals_ak15['sig'] = np.sum(np.array(disc_vals_ak15['sig']), axis=0)
    disc_vals_ak15['bg'] = np.sum(np.array(disc_vals_ak15['bg']), axis=0)
    disc_vals_ak15['tpr'] = np.cumsum(disc_vals_ak15['sig'][::-1])[::-1] / np.sum(np.array(disc_vals_ak15['sig']))
    disc_vals_ak15['fpr'] = np.cumsum(disc_vals_ak15['bg'][::-1])[::-1] / np.sum(np.array(disc_vals_ak15['bg']))
    disc_vals_ak15['tpr'] = np.append(disc_vals_ak15['tpr'], 0)
    disc_vals_ak15['fpr'] = np.append(disc_vals_ak15['fpr'], 0)

    tpr_aves = (disc_vals_ak15['tpr'][:-1] + disc_vals_ak15['tpr'][1:]) / 2
    fpr_diffs = (disc_vals_ak15['fpr'][:-1] - disc_vals_ak15['fpr'][1:])
    aucak15 = np.sum(tpr_aves * fpr_diffs)


    if ak8: plt.plot(disc_vals_ak8['fpr'], disc_vals_ak8['tpr'], label='AK8 AUC = {:.2f}'.format(aucak8))
    if ak15: plt.plot(disc_vals_ak15['fpr'], disc_vals_ak15['tpr'], label='AK15 AUC = {:.2f}'.format(aucak15))

    plt.title(f'{tagger_name} ROC Curve')
    plt.xlabel('Background Efficiency')
    plt.ylabel('Signal Efficiency')
    plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

    plt.tight_layout(0.5)
    plt.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True, style='sci')
    plt.savefig(f"figs/{name}.pdf", bbox_inches='tight')
    plt.show()


    if ak8: plt.semilogx(disc_vals_ak8['fpr'], disc_vals_ak8['tpr'], label='AK8 AUC = {:.2f}'.format(aucak8))
    if ak15: plt.semilogx(disc_vals_ak15['fpr'], disc_vals_ak15['tpr'], label='AK15 AUC = {:.2f}'.format(aucak15))

    plt.title(f'{tagger_name} Semilog ROC Curve')
    plt.xlabel('Background Efficiency')
    plt.ylabel('Signal Efficiency')
    plt.xlim(0.0001, 1)
    plt.legend(fancybox=True, shadow=True, frameon=True, prop={'size': 16})

    plt.tight_layout(0.5)
    plt.savefig(f"figs/{name}_semilog.pdf", bbox_inches='tight')
    plt.show()


def ddttemplate(evts, ak15=False, ak815=False, eff=0.99, name="", ddt=False, lo_threshold=0):
    key = 'PNetH4qvsQCD' if ak15 else 'PNetHqqqq'
    key_name = 'pneth4qvsQCD'

    ak15str = "ak15" if (ak15 or ak815) else ""
    hists[ak15str + 'ddthist'] = hist.Hist("Events",
                        hist.Cat("sample", "Sample"),
                        hist.Bin("pt", r"$p_T$ (GeV)", 100, 250, 1000),
                        hist.Bin("rho", r"$\rho = \ln (m_{SD}^2 / p_T^2)$", 100, -8, -0.5),
                        hist.Bin(key_name, key_name, 1000, lo_threshold, 1),
                        )

    hists[ak15str + 'ddthist'].fill(sample='QCD',
                         pt = evts['QCD'][ak15str + 'fatJet2Pt'],
                         rho = np.log(evts['QCD'][ak15str + 'fatJet2MassSD'] ** 2 / evts['QCD'][ak15str + 'fatJet2Pt'] ** 2),
                         pneth4qvsQCD = evts['QCD'][ak15str + f'fatJet2{key}'],
                         weight = evts['QCD']["weight"]
                         )

    val_QCD = hists[ak15str + 'ddthist']['QCD'].values(overflow='allnan')[('QCD',)]
    qcd_maxval_temp = np.cumsum(val_QCD, axis=2)
    qcd_maxval = qcd_maxval_temp[:, :, -1]
    norma = qcd_maxval_temp / np.maximum(1e-10, qcd_maxval[:, :, np.newaxis])

    hist_y_QCD = deepcopy(hists[ak15str + 'ddthist'])
    template = hist_y_QCD.sum(key_name, )  # pt, rho base
    hist_y_QCD.clear()
    hist_y_QCD._sumw = {():norma}

    res = np.apply_along_axis(lambda norma: norma.searchsorted(eff), axis = 2, arr = norma)
    res[res > 1000] = 0

    # evaluation of GRU cut from quantile (trivial if GRU has 100 bins)
    def bineval(a):
        return hist_y_QCD.identifiers(key_name,overflow='allnan')[a].lo

    binfunc = np.vectorize(bineval)
    qmap = binfunc(res)

    qmap[qmap == -math.inf] = lo_threshold
    qmap

    template.clear()
    template._sumw = {():qmap}
    template.label = 'Cut for {}% B effciency'.format(np.round((1 - eff) * 100, 2))

    if ak815: ak15str = "ak815"

    hist.plot2d(template.sum('sample'), xaxis = "rho", patch_opts={'vmax': 1, 'cmap': 'jet'})
    plt.savefig("figs/{}{}{}_ddtmap.pdf".format(name, ak15str, int(eff * 10000)))
    plt.show()

    coffea.util.save(template, "ddtmaps/{}{}ddtmap_{}_cut.coffea".format(name, ak15str, int(eff * 10000)))


def ddttagger(evts, ak15=False, ak815=False, eff=0.99, name=""):
    key = 'PNetH4qvsQCD' if ak15 else 'PNetHqqqq'

    if ak15: ak15str = "ak15"
    elif ak815: ak15str = "ak815"
    else: ak15str = ""

    ddtmap = uproot4.open("ddtmaps/{}{}ddtmap_{}_cut.root".format(name, ak15str, int(eff * 10000)))
    print(ddtmap['h1'].values())

    ext = extractor()
    ext.add_weight_sets(["ddtmap h1 ddtmaps/{}{}ddtmap_{}_cut.root".format(name, ak15str, int(eff * 10000))])
    ext.finalize()
    evaluator = ext.make_evaluator()

    for s in evts.keys():
        rho = np.log(evts[s][ak15str + 'fatJet2MassSD'] ** 2 / evts[s][ak15str + 'fatJet2Pt'] ** 2)
        pnet_cut = evaluator["ddtmap"](evts[s][ak15str + "fatJet2Pt"], rho)
        pnetddt = evts[s][ak15str + f'fatJet2{key}'] - pnet_cut
        evts[s] = ak.zip(dict(zip(evts[s].fields + [ak15str + 'fatJet2{}{}{}Eff'.format(key, name, int(eff * 10000)), ak15str + 'fatJet2{}{}DDT{}'.format(key, name, int(eff * 10000))], ak.unzip(evts[s]) + (ak.Array(pnet_cut), ak.Array(pnetddt)))))


def triggered_events(events_in):
    events_triggered = {}
    for s, evts in events_in.items():
        ak_dict = {}
        for field in evts.fields:
            ak_dict[field] = evts[field][evts['triggered']]
        events_triggered[s] = ak.zip(ak_dict)

    return events_triggered


def cutflow_func(var_cuts, events, ret_sb=False):
    cutflow = {}
    events_cuts = {}
    for s, evts in events.items():
        cuts = []
        for var, brange in var_cuts.items():
            if '+' in var:
                vars = var.split('+')
                cut1 = evts[vars[0]] > brange[0]
                cut2 = evts[vars[0]] < brange[1]
                for tvars in vars[1:]:
                    cut1 = cut1 + (evts[tvars] > brange[0])
                    cut2 = cut2 + (evts[tvars] < brange[1])

                cuts.append(cut1)
                cuts.append(cut2)
            else:
                cuts.append(evts[var] > brange[0])
                cuts.append(evts[var] < brange[1])

        cut = cuts[0]
        cutflow[s] = []
        for i in np.arange(1, len(cuts)):
            cutflow[s].append(np.sum(evts[cut].weight))
            cut = cut * cuts[i]

        cutflow[s].append(np.sum(evts[cut].weight))
        events_cuts[s] = evts[cut]

    if ret_sb:
        return [cutflow["HHbbVV4q"][-1], cutflow["QCD"][-1], cutflow["tt"][-1]]

    return cutflow, events_cuts


def cftable_func(cutflow, var_cuts, cut_labels=None, cut_idx=None):
    if cut_labels is None:
        cut_labels = []

        for var, cuts in var_cuts.items():
            varname = var.split('fat')[-1]
            if cuts[0] > 0 or "DDT" in varname: cut_labels.append("{} > {}".format(varname, cuts[0]))
            if cuts[1] < 9999: cut_labels.append("{} < {}".format(varname, cuts[1]))

    if cut_idx is None:
        i = 0
        cut_idx = []
        for brange in var_cuts.values():
            for j in brange:
                if j != 9999 and j != -9999:
                    cut_idx.append(i)
                i += 1

    return pandas.DataFrame(np.round(np.array(list(cutflow.values()))[:, cut_idx], 3), list(cutflow.keys()), cut_labels)


def plot_single_roc(evts, keys, sig_label, tagger_name, name, sig='HHbbVV4q'):
    if type(keys) is str: keys = [keys]

    y_score = np.concatenate((evts[sig][keys[0]], evts['QCD'][keys[0]]))
    y_true = np.concatenate((np.ones(len(evts[sig][keys[0]])), np.zeros(len(evts['QCD'][keys[0]]))))
    y_weight = np.concatenate((evts[sig]['weight'], evts['QCD']['weight']))

    for key in keys[1:]:
        y_score = np.concatenate((y_score, evts[sig][key], evts['QCD'][key]))
        y_true = np.concatenate((y_true, np.ones(len(evts[sig][key])), np.zeros(len(evts['QCD'][key]))))
        y_weight = np.concatenate((y_weight, evts[sig]['weight'], evts['QCD']['weight']))

    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score, sample_weight=y_weight)
    tpr_aves = (tpr[:-1] + tpr[1:]) / 2
    fpr_diffs = (fpr[1:] - fpr[:-1])
    auc = np.sum(tpr_aves * fpr_diffs)

    plt.semilogx(fpr, tpr, label=f'AUC = {auc:.2f}')
    plt.legend(loc='upper left', fancybox=True)
    plt.xlabel('QCD Background Efficiency')
    plt.ylabel(sig_label + ' Signal Efficiency')
    plt.xlim(1e-4, 1)
    plt.title(f"ParticleNet {tagger_name} ROC Curve")
    plt.savefig(f"figs/{name}.pdf")
    plt.show()
    return fpr, tpr, auc


# bVfpr, bVtpr, bVauc = plot_single_roc(hh4b_kin_cuts_evts, "fatJet1PNetXbb", 'Hbb', 'Txbb', 'pnetxbb_kin_cuts_roc')
# bbfpr, bbtpr, bbauc = plot_single_roc(hh4b_kin_cuts_evts, "fatJet1PNetXbb", 'Hbb', 'Txbb', 'pnetxbb_hh4b_kin_cuts_roc', sig='HH4b')




filehandler = open('events_ak15bb_sorted_triggered_effs.obj', 'rb')
events_ak15bb_sorted = pickle.load(filehandler)
filehandler.close()

filehandler = open('events_bb_sorted_triggered_effs.obj', 'rb')
events_ak8bb_sorted = pickle.load(filehandler)
filehandler.close()

events_ak15bb_sorted



scan = {
    "fatJet1Pt": [250, 275],
    "fatJet2Pt": [250, 275],
    "fatJet1MassSD": [100, 110, 115],
    "fatJet1MassSD_UPPER": [30, 40],
    "fatJet2MassSD": [100, 110, 115],
    "fatJet2MassSD_UPPER": [30, 40],
    "fatJet1PNetXbb": [0.9, 0.95],
}

all_cuts = []


def recursive_scanner(scan, keys, cut={}):
    key = keys[0]

    for val in scan[key]:
        tcut = deepcopy(cut)
        tcut[key] = [val]

        if key + "_UPPER" in scan:
            for val2 in scan[key + "_UPPER"]:
                tcut2 = deepcopy(tcut)
                tcut2[key].append(val + val2)

                if len(keys) == 2: all_cuts.append(tcut2)
                else: recursive_scanner(scan, keys[2:], tcut2)
        else:
            tcut[key].append(9999)
            if len(keys) == 1: all_cuts.append(tcut)
            else: recursive_scanner(scan, keys[1:], tcut)


recursive_scanner(scan, list(scan.keys()), cut={"fatJet2PNetH4qvsQCDNtupleCutsDDT9900": [0, 9999]})

len(all_cuts)

all_cuts

sbs = np.array(sbs)

significance = sbs[:, 0] / np.sqrt(np.sum(sbs, axis=1))

np.argmax(significance)

np.sort(significance)[-20:]

sbs[259]
sbs[283]
sbs[137]
all_cuts[259]


scan = {
    "ak15fatJet1Pt": [250, 275],
    "ak15fatJet2Pt": [250, 275],
    "ak15fatJet1MassSD": [100, 110],
    "ak15fatJet1MassSD_UPPER": [40, 50],
    "ak15fatJet2MassSD": [100, 110],
    "ak15fatJet2MassSD_UPPER": [40, 50],
    "ak15fatJet1PNetMDXbb": [0.9, 0.925, 0.95],
}



ak15all_cuts = []

def ak15recursive_scanner(scan, keys, cut={}):
    key = keys[0]
    for val in scan[key]:
        tcut = deepcopy(cut)
        tcut[key] = [val]

        if key + "_UPPER" in scan:
            for val2 in scan[key + "_UPPER"]:
                tcut2 = deepcopy(tcut)
                tcut2[key].append(val + val2)

                if len(keys) == 2: ak15all_cuts.append(tcut2)
                else: ak15recursive_scanner(scan, keys[2:], tcut2)
        else:
            tcut[key].append(9999)
            if len(keys) == 1: ak15all_cuts.append(tcut)
            else: ak15recursive_scanner(scan, keys[1:], tcut)


ak15recursive_scanner(scan, list(scan.keys()), cut={"ak15fatJet2PNetH4qvsQCDNtupleCutsDDT9900": [0, 9999]})


len(ak15all_cuts)

ak15sbs = []


for i, cuts in tqdm(enumerate(ak15all_cuts), total=len(ak15all_cuts)):
    ak15sbs.append(cutflow_func(cuts, events_bbs_ak15_cuts, ret_sb=True))

ak15sbs




ak15sbs = np.array(ak15sbs)


ak15sbs[:, 0].max()


ak15significance = ak15sbs[:, 0] / np.sqrt(np.sum(ak15sbs, axis=1))

np.argmax(ak15significance)
np.sort(ak15significance)[-20:]

np.sort(ak15significance)[-20:]

ak15sbs[np.argsort(ak15significance)[-10:]]

np.argsort(ak15significance)[-10:]


ak15sbs[261]

ak15all_cuts[261]












# filehandler = open('events_bb_sorted.obj', 'rb')
# events_ak8bb_sorted_no_trigger = pickle.load(filehandler)
# filehandler.close()


var = "PNetMDXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists(var, varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(var, "jet_pnet_ak15", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False, ak15=True)


var = "PNetMDXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists(var, varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak8bb_sorted, useEDWeights=True, ak15=True)
plot_hists(var, "jet_pnet_ak15_ak8_sorted", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False, ak15=True)


init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak15bb_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_ak15bb_leading", [100, 0, 1], name="bb_sorted_cut", hh4v=False, sepsig=False)


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak8bb_sorted, useEDWeights=True)
plot_hists(vars, "jet_kin_bb_sorted", bins, name="bb_sorted", sepsig=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_bb_sorted_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak8bb_sorted, useEDWeights=True)
plot_hists(vars, "jet_kin_bb_sorted_fine", bins, name="bb_sorted", sepsig=False)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_ak15bb_sorted, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_bb_sorted_fine_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
ak815_comp_plot(vars, varsl, bins, events_ak8bb_sorted, events_ak15bb_sorted, "ak815_or_cut_triggered")


var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists("PNetXbb", varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_ak8bb_sorted, useEDWeights=True)
plot_hists(var, "jet_cut_pnet2", bins, name="bb_sorted_cut", hh4v=False, sepsig=False)



vars = ["Pt", "Mass", "MassSD"]
bins = [[50, 200, 600], [50, 1, 250], [50, 1, 250]]

fig, axs = plt.subplots(2, len(vars), figsize=(len(vars) * 10, 18))

for i in range(len(vars)):
    for j in range(2):
        if j == 0:
            ak8hists = events_ak8bb_sorted["HHbbVV4q"]["fatJet1{}".format(vars[i])] * ak8fj1match
            ak15hists = events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet1{}".format(vars[i])] * ak15fj1ak8fj1 + events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet2{}".format(vars[i])] * ak15fj2ak8fj1
        else:
            ak8hists = events_ak8bb_sorted["HHbbVV4q"]["fatJet2{}".format(vars[i])] * ak8fj2match
            ak15hists = events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet1{}".format(vars[i])] * ak15fj1ak8fj2 + events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet2{}".format(vars[i])] * ak15fj2ak8fj2

        h2d = axs[j, i].hist2d(ak.to_numpy(ak8hists), ak.to_numpy(ak15hists), bins=bins[i][0], range=[[bins[i][1], bins[i][2]], [bins[i][1], bins[i][2]]])
        ax_divider = make_axes_locatable(axs[j, i])
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        fig.colorbar(h2d[3], cax=cax)
        axs[j, i].set_xlabel("FatJet{} {} (GeV)".format(j + 1, vars[i]))
        axs[j, i].set_ylabel("AK15FatJet {} {} (GeV)".format(j + 1, vars[i]))

fig.tight_layout(pad=2)
plt.savefig("figs/ak15ak82dplots.pdf")
plt.show()


plt.hist2d(np.concatenate([ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["fatJet1PNetXbb"] * ak8fj1match), ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["fatJet2PNetXbb"] * ak8fj2match)]),
            np.concatenate([ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet1PNetMDXbb"] * ak15fj1ak8fj1 + events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet2PNetMDXbb"] * ak15fj2ak8fj1),
                            ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet1PNetMDXbb"] * ak15fj1ak8fj2 + events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet2PNetMDXbb"] * ak15fj2ak8fj2)]),
            bins=25, range=[[0, 1], [0, 1]], norm=mpl.colors.LogNorm())
plt.colorbar()
plt.xlabel("AK8 PNet Score")
plt.ylabel("AK15 PNet Score")
plt.savefig("figs/pnet_2d.pdf")


vars = ["PNet"]
bins = [[50, 200, 600], [50, 1, 250], [50, 1, 250]]

fig, axs = plt.subplots(1, 2, figsize=(len(vars) * 10, 18))

for j in range(2):
    h2d = axs[j, i].hist2d(ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["fatJet{}{}".format(j + 1, vars[i])]), ak.to_numpy(events_ak8bb_sorted["HHbbVV4q"]["ak15fatJet{}{}".format(j + 1, vars[i])]), bins=bins[i][0], range=[[bins[i][1], bins[i][2]], [bins[i][1], bins[i][2]]])
    ax_divider = make_axes_locatable(axs[j, i])
    cax = ax_divider.append_axes("right", size="7%", pad="2%")
    fig.colorbar(h2d[3], cax=cax)
    axs[j, i].set_xlabel("FatJet{} {} (GeV)".format(j + 1, vars[i]))
    axs[j, i].set_ylabel("AK15FatJet {} {} (GeV)".format(j + 1, vars[i]))

fig.tight_layout(pad=2)
plt.savefig("figs/ak15ak82dplots.pdf")
plt.show()






# ak8 kin ntuple cuts

var_cuts = {
    "fatJet1Pt": [250, 9999],
    "fatJet2Pt": [250, 9999],
    "fatJet1MassSD": [20, 9999],
    "fatJet2MassSD": [20, 9999],
}

hhbbVV_cutflow, events_bbs_ak8_cuts = cutflow_func(var_cuts, events_ak8bb_sorted)

ak.sum(events_bbs_ak8_cuts['HHbbVV4q']['weight'])
ak.sum(events_bbs_ak8_cuts['data']['weight'])

del(events_bbs_ak8_cuts['HH4b'])

# for i in [0.99, 0.997, 0.999, ]
ddttemplate(events_bbs_ak8_cuts, ak15=False, name="NtupleCuts", eff=0.99)
ddttagger(events_bbs_ak8_cuts, ak15=False, name="NtupleCuts", eff=0.99)
events_bbs_ak8_cuts['HHbbVV4q'].fields


cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable


vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak8_cuts, useEDWeights=True)
plot_hists(vars, "jet_kin_ak8_cuts_bb_sorted", bins, name="bb_sorted", sepsig=False)

vars = ["Pt", "Mass", "MassSD"]
varsl = ["$p_T$ (GeV)", "Mass (GeV)", "Soft Drop Mass (GeV)"]
bins = [[80, 200, 1000], [80, 1, 400], [80, 1, 400]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak8_cuts, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_ak8_cuts_bb_sorted_ak15", bins, name="bb_sorted", sepsig=False, ak15=True)



# ak15 kin ntuple cuts

var_cuts = {
    "ak15fatJet1Pt": [250, 9999],
    "ak15fatJet2Pt": [250, 9999],
    "ak15fatJet1MassSD": [20, 9999],
    "ak15fatJet2MassSD": [20, 9999],
}

hhbbVV_cutflow, events_bbs_ak15_cuts = cutflow_func(var_cuts, events_ak15bb_sorted)

ak.sum(events_bbs_ak15_cuts['HHbbVV4q']['weight'])
ak.sum(events_bbs_ak15_cuts['data']['weight'])

ddttemplate(events_bbs_ak15_cuts, ak15=True, name="NtupleCuts", lo_threshold=0, eff=0.9900)
ddttagger(events_bbs_ak15_cuts, ak15=True, name="NtupleCuts", eff=0.99)


events_bbs_ak15_cuts['QCD'].fields

cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
ak815_comp_plot(vars, varsl, bins, events_bbs_ak8_cuts, events_bbs_ak15_cuts, "ak815_resp_cut")


tagger2d_plots('tagger2d_ntuple_cuts', events_bbs_ak8_cuts, events_bbs_ak15_cuts, vmaxhh=0.01, vmaxqcd=250)

ak.sum(events_bbs_ak8_cuts['data']['weight'])
ak.sum(events_bbs_ak15_cuts['data']['weight'])

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_ak8_cuts, useEDWeights=True)
plot_hists(var, "jetak8_ntuple_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_bbs_ak15_cuts, useEDWeights=True)
plot_hists(var, "jetak15_ntuple_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])

del(events_bbs_ak8_cuts['HH4b'])
del(events_bbs_ak15_cuts['HH4b'])

events_bbs_ak8_cuts["HHbb"]

roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetXbb", ak15tagger="PNetMDXbb", tagger_name="ParticleNet Xbb", name="PNetXbb_no_cuts_no_trigger_roc", jet=1, ak8events=events_bbs_ak8_cuts, ak15events=events_bbs_ak15_cuts, ak15=False)
roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetHqqqq", ak15tagger="PNetH4qvsQCD", tagger_name="ParticleNet Th4q", name="PNetH4qvsQCD_kin_cuts_roc", jet=2, ak8events=events_bbs_ak8_cuts, ak15events=events_bbs_ak15_cuts, ak8=False)
roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetHqqqq", ak15tagger="PNetMDXbb", tagger_name="ParticleNet Txbb", name="ak15PNetTxbb_kin_cuts_roc", jet=1, ak8events=events_bbs_ak8_cuts, ak15events=events_bbs_ak15_cuts, ak8=False)
roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetXbb", ak15tagger="PNetMDXbb", tagger_name="ParticleNet Txbb", name="ak8PNetTxbb_kin_cuts_roc", jet=1, ak8events=events_bbs_ak8_cuts, ak15events=events_bbs_ak15_cuts, ak15=False)



roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetHqqqq", ak15tagger="PNetHqqqq", tagger_name="ParticleNet Hqqqq", name="PNetHqqqq_kin_cuts_roc", jet=2, ak8events=events_bbs_ak8_cuts, ak15events=events_bbs_ak15_cuts, ak8=False)

name = "test"
init_hists("PNetHqqqq", "PNetHqqqq", [1000, 0, 1], name=name + "ak8", fatJet=True)
fill_hists("PNetHqqqq", pevtDict=events_bbs_ak8_cuts, name=name + "ak8", fatJet=True, scale=False, hh4v=False, ak15=False)




init_hists("Pt", "$p_T$ (GeV)", [40, 200, 2000], fatJet=True)
fill_hists("Pt", fatJet=True, scale=False, pevtDict=events_bbs_ak15_cuts, useEDWeights=True, ak15=True)

tot_events = {}
vals = hists['fatJetPt'].project("sample", "jet1").values()
for s in evtDict.keys():
    tot_events[s] = np.sum(vals[(s,)])


tot_events


scale_factor = {
        # 'HH4V': 1 / tot_events['HH4V'],
        'HHbbVV4q': 1 / tot_events['HHbbVV4q'],
        # 'HH4b': 1 / tot_events['HH4b'],
        'QCD': 1 / (tot_events['QCD'] + tot_events['tt']),
        'tt': 1 / (tot_events['QCD'] + tot_events['tt'])}





vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_kin_ak15_ntuple_cuts", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = "PNetMDXbb"
varsl = "MD PNet $T_{Xbb}$"
bins = [100, 0, 1]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_pnetxbb_ak15_ntuple_cuts", bins, name="bb_sorted", sepsig=False, ak15=True)


vars = "PNetH4qvsQCD"
varsl = "Non-MD PNet $T_{H4q}$ Score"
bins = [100, 0, 1]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_pneth4qvsqcd_ak15_ntuple_cuts", bins, name="bb_sorted", sepsig=False, ak15=True)




vars = "Pt"
varsl = "$p_T$ (GeV)"
bins = [50, 250, 500]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", scale=False, bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, scale=False, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_pt_ak15_ntuple_cuts_datamc", bins, name="bb_sorted", sepsig=False, ak15=True, data=True, ratio=True, rat_ylim=2, lumilabel=40)


vars = "MassSD"
varsl = "Soft Drop Mass (GeV)"
bins = [50, 50, 200]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", scale=False, bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, scale=False, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_mass_ak15_ntuple_cuts_datamc", bins, name="bb_sorted", sepsig=False, ak15=True, data=True, ratio=True, rat_ylim=2, lumilabel=40)


vars = "PNetMDXbb"
varsl = "MD PNet $T_{Xbb}$"
bins = [100, 0, 1]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", scale=False, bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, scale=False, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_pnetxbb_ak15_ntuple_cuts_datamc", bins, name="bb_sorted", sepsig=False, ak15=True, data=True, ratio=True, rat_ylim=2, lumilabel=40)


vars = "PNetH4qvsQCD"
varsl = "Non-MD PNet $T_{h4q}$ Score"
bins = [100, 0, 1]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted", scale=False, bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted", pevtDict=events_bbs_ak15_cuts, scale=False, useEDWeights=True, ak15=True)
plot_hists(vars, "jet_pneth4qvsqcd_ak15_ntuple_cuts_datamc", bins, name="bb_sorted", sepsig=False, ak15=True, data=True, ratio=True, rat_ylim=2, lumilabel=40)





events_bbs_ak8_cuts
events_bbs_ak8_cuts['data'].fields

# ak815 kin ntuple cuts

var_cuts = {
    "fatJet1Pt": [250, 9999],
    "ak15fatJet2Pt": [250, 9999],
    "fatJet1MassSD": [20, 9999],
    "ak15fatJet2MassSD": [20, 9999],
}

hhbbVV_cutflow, events_bbs_ak815_cuts = cutflow_func(var_cuts, events_8bb_15VV_sorted_triggered)

ddttemplate(events_bbs_ak815_cuts, ak815=True, name="NtupleCuts", lo_threshold=0, eff=0.9900)
ddttagger(events_bbs_ak15_cuts, ak15=True, name="NtupleCuts", eff=0.99)




# tagger cuts only ak15

tagger_cuts = {
    "ak15fatJet1PNetMDXbb": [0, 9999],
    "ak15fatJet2PNetHqqqq": [0.94, 9999],
}

cutflow, events_ak15bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak15_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_tagger_cuts, useEDWeights=True, ak15=True)
plot_hists(var, "jetak15_loose_tagger_cuts_{}_{}_masssd_rat".format(tagger_cuts["ak15fatJet1PNetMDXbb"][0], tagger_cuts["ak15fatJet2PNetHqqqq"][0]), bins, ak15=True, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, rat_ylim=2.5, blinding=[115, 145])

# tagger cuts only ak8

tagger_cuts = {
    "fatJet1PNetXbb": [0, 9999],
    "fatJet2PNetHqqqq": [0.94, 9999],
}

cutflow, events_ak8bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak8_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_tagger_cut", pevtDict=events_ak8bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak8_loose_tagger_cuts_{}_{}_masssd".format(tagger_cuts["fatJet1PNetXbb"][0], tagger_cuts["fatJet2PNetHqqqq"][0]), bins, name="bb_sorted_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, rat_ylim=2.5, blinding=[115, 145])





# tagger cuts only ak15

tagger_cuts = {
    "ak15fatJet1PNetMDXbb": [0, 9999],
    "ak15fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}

cutflow, events_ak15bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak15_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

events_bbs_ak15_cuts['QCD'].fields

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [9, 70, 205]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_tagger_cuts, useEDWeights=True, ak15=True)
plot_hists(var, "jetak15_loose_taggerddt_cuts_{}_9900_masssd_rat".format(tagger_cuts["ak15fatJet1PNetMDXbb"][0]), bins, ak15=True, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, rat_ylim=2.5, blinding=[115, 145])



# tagger cuts only ak8

tagger_cuts = {
    "fatJet1PNetXbb": [0, 9999],
    "fatJet2PNetHqqqqNtupleCutsDDT990": [0, 9999],
    # "fatJet2PNetHqqqqNtupleCutsDDT990": [0, 9999],
}

cutflow, events_ak8bbs_tagger_cuts = cutflow_func(tagger_cuts, events_bbs_ak8_cuts)
# cftable = cftable_func(cutflow, tagger_cuts)
# cftable

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_tagger_cut", pevtDict=events_ak8bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak8_loose_taggerddt_cuts_{}_masssd".format(tagger_cuts["fatJet1PNetXbb"][0]), bins, name="bb_sorted_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, rat_ylim=2.5, blinding=[115, 145])







# pt cuts only

kin_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    # "fatJet1PNetXbb": [0.99, 9999],
    # "fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}

cutflow, events_bbs_ak8_pt_cuts = cutflow_func(kin_cuts, events_bbs_ak8_cuts)

# cftable = cftable_func(cutflow, kin_cuts)
# cftable
#
# cftable.to_csv('hhbbVV_cutflow_pt_cuts.csv')
#
# ddttemplate(events_bbs_pt_cuts, ak15=False, name="PtCuts", lo_threshold=0, eff=0.99)
# ddttagger(events_bbs_pt_cuts, ak15=False, name="PtCuts", eff=0.99)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [10, 55, 205]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_pt_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_pt_tagger_cut", pevtDict=events_bbs_ak8_pt_cuts, useEDWeights=True, blinding=[4, 6])
plot_hists(var, "jetak8_pt_cuts__masssd", bins, name="bb_sorted_pt_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False)



# pt cuts only

kin_cuts = {
    "ak15fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "ak15fatJet1MassSD": [110, 140],
    # "ak15fatJet1PNetMDXbb": [0.99, 9999],
    # "ak15fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}

cutflow, events_bbs_ak15_pt_cuts = cutflow_func(kin_cuts, events_bbs_ak15_cuts)

# cftable = cftable_func(cutflow, kin_cuts)
# cftable
#
# cftable.to_csv('hhbbVV_cutflow_ak15_pt_cuts.csv')
#
# ddttemplate(events_bbs_ak15_pt_cuts, ak15=True, name="PtCuts", eff=0.9997)
# ddttagger(events_bbs_ak15_pt_cuts, ak15=True, name="PtCuts", eff=0.9997)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [10, 55, 205]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_pt_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_pt_tagger_cut", pevtDict=events_bbs_ak15_pt_cuts, useEDWeights=True, blinding=[4, 6], ak15=True)
plot_hists(var, "jetak15_pt_cuts_masssd", bins, name="bb_sorted_pt_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, ak15=True)


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[50, 250, 500], [50, 50, 200]]
ak815_comp_plot(vars, varsl, bins, events_bbs_ak8_pt_cuts, events_bbs_ak15_pt_cuts, "ak815_pt_cut")



roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetXbb", ak15tagger="PNetMDXbb", tagger_name="ParticleNet Xbb", name="PNetXbb_pt_cuts_roc", jet=1, ak8events=events_bbs_ak8_pt_cuts, ak15events=events_bbs_ak15_pt_cuts)
roc_curves(disc_bin=[1000, 0, 1], ak8tagger="PNetHqqqq", ak15tagger="PNetHqqqq", tagger_name="ParticleNet Hqqqq", name="PNetHqqqq_pt_cuts_roc", jet=2, ak8events=events_bbs_ak8_pt_cuts, ak15events=events_bbs_ak15_pt_cuts)




# kin cuts only ak15

kin_cuts = {
    "ak15fatJet1Pt": [275, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "ak15fatJet1MassSD": [110, 140],
    "ak15fatJet2MassSD": [115, 145],
}

cutflow, events_bbs_ak15_kin_cuts = cutflow_func(kin_cuts, events_bbs_ak15_cuts)

ddttemplate(events_ak15bbs_kin_cuts, ak15=True)


# kin cuts only ak8

kin_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet2MassSD": [115, 145],
}

cutflow, events_bbs_ak8_kin_cuts = cutflow_func(kin_cuts, events_bbs_ak8_cuts)

ddttemplate(events_bbs_kin_cuts, ak15=False)
ddttagger(events_bbs_kin_cuts, ak15=False)

tagger2d_plots('tagger2d_kin_cuts', events_bbs_kin_cuts, events_ak15bbs_kin_cuts, vmaxhh=0.005, vmaxqcd=10)

roc_curves(disc_bin=[10000, 0, 1], ak8tagger="PNetXbb", ak15tagger="PNetMDXbb", tagger_name="ParticleNet Xbb", name="PNetXbb_kin_cuts_roc", jet=1, ak8events=events_bbs_ak8_kin_cuts, ak15events=events_bbs_ak15_kin_cuts)
roc_curves(disc_bin=[10000, 0, 1], ak8tagger="PNetHqqqq", ak15tagger="PNetHqqqq", tagger_name="ParticleNet Hqqqq", name="PNetHqqqq_kin_cuts_roc", jet=2, ak8events=events_bbs_ak8_kin_cuts, ak15events=events_bbs_ak15_kin_cuts)



var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jetak8_kin_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_ak15bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jetak15_kin_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_ak8_cuts_masssd", bins, name="bb_sorted", sepsig=False)

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [50, 50, 200]
init_hists(var, varl, bins, fatJet=True, name="bb_sorted_ak15", bbsorted=True)
fill_hists(var, fatJet=True, name="bb_sorted_ak15", pevtDict=events_bbs_kin_cuts, useEDWeights=True, ak15=True)
plot_hists(var, "jet_ak8_cuts_masssd_ak15", bins, name="bb_sorted_ak15", sepsig=False, ak15=True)


var = "PNetXbb"
varl = "PNetXbb Score"
bins =  [100, 0, 1]
init_hists("PNetXbb", varl, bins, scale=True, fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists(var, "jet_cut_pnet2", bins, name="bb_sorted_cut", hh4v=False, sepsig=False, stackall=False)

init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="bb_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="bb_sorted_cut", pevtDict=events_bbs_kin_cuts, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_bb_leading2", [100, 0, 1], name="bb_sorted_cut", hh4v=False, sepsig=False)


yields = {"sig": [], "bg": []}
for pnetcutoff in np.arange(0.99, 0.995, 0.001)[:-1]:
    sigy = []
    bgy = []
    for dakcutoff in np.arange(0.9, 1, 0.01):
        cuts = {}
        for s, evts in events_bbs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetXbb"] > pnetcutoff) * (evts["fatJet2PNetHqqqq"] > dakcutoff)

        sig = np.sum(events_bbs_kin_cuts["HHbbVV4q"][cuts["HHbbVV4q"]].weight)
        bg = np.sum(events_bbs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_bbs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)







# tagger cuts only ak8 + ak15

tagger_cuts = {
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqq": [0.935, 9999],
}

cutflow, events_ak815bbs_tagger_cuts = cutflow_func(tagger_cuts, events_8bb_15VV_sorted)
cftable = cftable_func(cutflow, tagger_cuts)
cftable


vars = ["Pt", "MassSD"]
varsl = ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"]
bins = [[20, 250, 500], [25, 50, 200]]
init_hists(vars, varsl, bins, fatJet=True, name="bb_sorted815", bbsorted=True)
fill_hists(vars, fatJet=True, name="bb_sorted815", pevtDict=events_ak815bbs_tagger_cuts, useEDWeights=True, jet2ak15=True)
plot_hists(vars, "jet_kin_tagger_cuts_ak815", bins, name="bb_sorted815", sepsig=False, ak15=True)


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_tagger_cut", pevtDict=events_ak8bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak8_tagger_cuts_masssd", bins, name="bb_sorted_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_tagger_cuts, useEDWeights=True)
plot_hists(var, "jetak15_tagger_cuts_masssd", bins, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])






# AK15

# ak15 all cuts

var_cuts = {
    'ak15fatJet1Pt': [250, 9999],
    'ak15fatJet2Pt': [250, 9999],
    'ak15fatJet1MassSD': [100, 150],
    'ak15fatJet2MassSD': [100, 150],
    'ak15fatJet1PNetMDXbb': [0.9, 9999],
    'ak15fatJet2PNetH4qvsQCDNtupleCutsDDT9900': [0, 9999],
}

hhbbVV_cutflow, events_bbs_ak15_all_cuts = cutflow_func(var_cuts, events_bbs_ak15_cuts)

# del(hhbbVV_cutflow['data'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable


1.035/1.37
3209 / 6.7e5

1.8/6e6

(1.8/6e6) / 6.2e-7

0.68/8.3e5

2e-6 / (0.68/8.3e5)

2.44/2

0.532 / 16000

3.2e-4 / (0.532 / 16000)

0.5/0.683

# sig eff of pnet xbb at 0.5% bg eff
0.356 / 0.683

# comparing to 4b xbb
(1.035 / 1.37) / (0.356 / 0.683)


220 / 16000

# bbVV second cut
0.261 / 0.532
400 / 13600

# 4b second cut
0.33 / 1
14 / 3200

(0.33 / 0.0044) / (0.5 / 0.03)


# branching ratio * mass resolution
2.6 * 1.22
9.62 / (2.6 * 1.22)  # bb tagger factor




# branching ratio * mass resolution * bb tagger * lower bb cut * VV tagger * lower VV cut
2.3 * 1.22 * 1.45 * 2 * 2 * 2.25



# effect of increasing signal yield to 0.3 (mass resolution * worse bb tagger * lower bb cut)
1.22 * 1.45 * 2

4 / 3.172 # bb tagger factor

1.22 * 1.45 * 2 * (1.22 ** 2)

# equal taggers + mass regression
2.3 * 2 * 2.25 / (1.22 ** 2)


494/14


0.3 / np.sqrt(0.3 + 494 / 5)
1 / np.sqrt(1 + 494 * 3/ 5)

1 / np.sqrt(1 + 14 * 3)

cftable.to_csv('hhbbVV_cutflow_ak15.csv')


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_pt_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_pt_tagger_cut", pevtDict=events_bbs_ak15_pt_tagger_cuts, useEDWeights=True, blinding=[3, 5], ak15=True)
plot_hists(var, "jetak15_pt_tagger_cuts_masssd", bins, name="bb_sorted_pt_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, ak15=True)



# bg estimatation using mass2 sidebands

var_cuts = {
    'ak15fatJet1Pt': [275, 9999],
    'ak15fatJet2Pt': [275, 9999],
    'ak15fatJet1MassSD': [110, 140],
    # 'ak15fatJet2MassSD': [110, 150],
    'ak15fatJet1PNetMDXbb': [0.95, 9999],
    'ak15fatJet2PNetH4qvsQCDNtupleCutsDDT9900': [0, 9999],
}


cutflow, events_ak15bbs_nomass2_cuts = cutflow_func(var_cuts, events_bbs_ak15_cuts)

del(cutflow['data'])

cftable = cftable_func(cutflow, var_cuts)
cftable


var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [6, 70, 190]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_tagger_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_tagger_cut", pevtDict=events_ak15bbs_nomass2_cuts, useEDWeights=True, ak15=True, blinding=[2, 4], data_blinded_only=True)
plot_hists(var, "jetak15_post_all_cuts_ddt_0bb_cuts", bins, ak15=True, name="bb_sortedak15_tagger_cut", hh4v=False, sepsig=False, lumilabel=40, data=True, ratio=True, rat_ylim=4, log=True, log_limits=[0.01, 1e3], legend_loc='center right', incl_signal=True)




# mass sidebands

mass_sb1 = {
    "ak15fatJet2MassSD": [90, 110],
}

msb1_cf = cutflow_func(mass_sb1, events_ak15bbs_nomass2_cuts)[0]
msb1_cf


mass_sb2 = {
    "ak15fatJet2MassSD": [150, 170],
}

msb2_cf = cutflow_func(mass_sb2, events_ak15bbs_nomass2_cuts)[0]
msb2_cf






# AK8

# ak8 all cuts

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    "fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.985, 9999],
    # "fatJet2PNetHqqqq": [0.94, 9999],
    "fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}

# del(events_bbs_ak8_cuts['HH4b'])

hhbbVV_cutflow, events_nopneth4q = cutflow_func(var_cuts, events_bbs_ak8_cuts)

del(hhbbVV_cutflow['data'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

# non md tagger, all kin cuts, 0.94 cut:
0.3/0.438
11940/441592

# ddt, all kin cuts, 0.01 bg eff
0.23/0.438
4630/441592

cftable.to_csv('hhbbVV_cutflow.csv')


# bg estimatation using mass2 sidebands

var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [110, 140],
    # "fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    # "fatJet2PNetHqqqq": [0.94, 9999],
    "fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}

cutflow, events_bbs_nomass2_cuts = cutflow_func(var_cuts, events_bbs_ak8_cuts)

del(cutflow['data'])

cftable = cftable_func(cutflow, var_cuts)
cftable

20000/3600000
0.47/1.1




# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [100, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_nomass2_cuts)[0]
msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [145, 160],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_nomass2_cuts)[0]
msb2_cf



ak815_comp_plot(["Pt", "MassSD"], ["$p_T$ (GeV)", "Soft Drop Mass (GeV)"], [[20, 250, 500], [25, 50, 200]], events_bbs_cuts, events_ak15bbs_cuts, "ak815_nomassd_cut", qcd=False)

var = "MassSD"
varl = "Soft Drop Mass (GeV)"
bins =  [7, 70, 175]
init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sorted_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sorted_all_cut", pevtDict=events_bbs_nomass2_cuts, useEDWeights=True)
plot_hists(var, "jetak8_all-mass2_cuts_masssd", bins, name="bb_sorted_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])


init_hists(var, varl, bins, scale=False, fatJet=True, name="bb_sortedak15_all_cut", bbsorted=True)
fill_hists(var, fatJet=True, scale=False, name="bb_sortedak15_all_cut", pevtDict=events_ak15bbs_nomass2_cuts, useEDWeights=True)
plot_hists(var, "jetak15_all-mass2_cuts_masssd", bins, name="bb_sortedak15_all_cut", hh4v=False, sepsig=False, lumilabel=40, stackall=True, data=True, ratio=False, blinding=[115, 145])




# AK8 for bb, AK15 for WW


# + mass2

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [110, 140],
    "ak15fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}



hhbbVV_cutflow = cutflow_func(var_cuts, events_bbs_ak815_cuts)[0]

del(hhbbVV_cutflow['data'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

cftable.to_csv('hhbbVV_ak8bb15VV_cutflow.csv')


# kin + bbVV tagger cuts - masssd

var_cuts = {
    "fatJet1Pt": [340, 9999],
    "ak15fatJet2Pt": [300, 9999],
    "fatJet1MassSD": [110, 140],
    # "ak15fatJet2MassSD": [115, 145],
    "fatJet1PNetXbb": [0.99, 9999],
    "ak15fatJet2PNetHqqqqNtupleCutsDDT9900": [0, 9999],
}


hhbbVV_cutflow, events_bbs_cuts = cutflow_func(var_cuts, events_bbs_ak815_cuts)

del(hhbbVV_cutflow['data'])
cftable = cftable_func(hhbbVV_cutflow, var_cuts)
cftable

# mass sidebands

mass_sb1 = {
    "ak15fatJet2MassSD": [100, 115],
}

msb1_cf = cutflow_func(mass_sb1, events_bbs_cuts)[0]
msb1_cf


mass_sb2 = {
    "ak15fatJet2MassSD": [145, 160],
}

msb2_cf = cutflow_func(mass_sb2, events_bbs_cuts)[0]
msb2_cf












# 4b cuts

hh4b_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet2MassSD": [105, 135],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
}

hh4bcf, hh4b_cuts_evts = cutflow_func(hh4b_var_cuts, events_bbs_ak8_cuts)

ak.sum(events_bbs_ak8_cuts["HH4b"]['weight'])

cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet2 MassSD > 105',
                'Jet2 MassSD < 135',
                'Jet1 PNetXbb > 0.985',
                'Jet2 PNetXbb > 0.985',
            ]

del(hh4bcf['HHbbVV4q'])
del(hh4bcf['data'])

cftable = cftable_func(hh4bcf, hh4b_var_cuts, cut_labels)
cftable

cftable.to_csv('hh4b_cutflow.csv')


# 4b cut - jet2 MassSD

hh4b_mass2_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet1PNetXbb": [0.985, 9999],
    "fatJet2PNetXbb": [0.985, 9999],
}

hh4bm2cf, hh4b_no_mass2_cuts_evts = cutflow_func(hh4b_mass2_cuts, events_bbs_ak8_cuts)


# mass sidebands

mass_sb1 = {
    "fatJet2MassSD": [90, 105],
}

hh4b_msb1_cf = cutflow_func(mass_sb1, hh4b_no_mass2_cuts_evts)[0]
hh4b_msb1_cf


mass_sb2 = {
    "fatJet2MassSD": [135, 150],
}

hh4b_msb2_cf = cutflow_func(mass_sb2, hh4b_no_mass2_cuts_evts)[0]
hh4b_msb2_cf




evtDict["HH4V"].fields


init_hists("PNetHqqqq", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("DeepAK8_H", "DeepAK8MD H", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("DeepAK8_H", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("DeepAK8_H", "jet_deepak8h_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau2OverTau1", "tau2/tau1", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau2OverTau1", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau2OverTau1", "jet_tau2otau1_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("Tau4OverTau3", "tau4/tau3", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Tau4OverTau3", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Tau4OverTau3", "jet_tau4otau3_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)


init_hists("PNetHqqqq", "DeepAK8MD H4q", [100, 0, 1], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted", hh4v=True, sepsig=False)



init_hists("Pt", "$p_T$", [100, 250, 1200], fatJet=True, name="WW_sorted", wwsorted=True)
fill_hists("Pt", fatJet=True, scale=True, name="WW_sorted", pevtDict=events_WW_sorted, useEDWeights=True)
plot_hists("Pt", "jet_pt_WW_leading", [100, 250, 1200], name="WW_sorted", hh4v=True, sepsig=False)


# WW kin cuts only

kin_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [75, 140],
    "fatJet2MassSD": [75, 140],
}

cutflow, events_WWs_kin_cuts = cutflow_func(kin_cuts, events_WW_sorted)

cutflow

init_hists("PNetHqqqq", "DeepAK8MD H4q vs QCD", [100, 0, 1], fatJet=True, name="WW_sorted_cut", bbsorted=True)
fill_hists("PNetHqqqq", fatJet=True, scale=True, name="WW_sorted_cut", pevtDict=events_WWs_kin_cuts, useEDWeights=True)
plot_hists("PNetHqqqq", "jet_cut_deepak8mdh4q_WW_leading", [100, 0, 1], name="WW_sorted_cut", hh4v=True, sepsig=False)

sig

hists['WW_sorted_cutPNetHqqqq']['HH4V']


yields = {"sig": [], "bg": []}
for pnetcutoff1 in np.arange(0.8, 0.85, 0.01):
    sigy = []
    bgy = []
    for pnetcutoff2 in np.arange(0.8, 0.85, 0.01):
        cuts = {}
        for s, evts in events_WWs_kin_cuts.items():
            cuts[s] = (evts["fatJet1PNetHqqqq"] > pnetcutoff1) * (evts["fatJet2PNetHqqqq"] > pnetcutoff2)

        sig = np.sum(events_WWs_kin_cuts["HH4V"][cuts["HH4V"]].weight)
        bg = np.sum(events_WWs_kin_cuts["QCD"][cuts["QCD"]].weight) + np.sum(events_WWs_kin_cuts["tt"][cuts["tt"]].weight)
        sigy.append(sig)
        bgy.append(bg)

    yields["sig"].append(sigy)
    yields["bg"].append(bgy)


yields







# 4b cuts only kin

hh4b_var_cuts = {
    "fatJet1Pt": [310, 9999],
    "fatJet2Pt": [310, 9999],
    "fatJet1Pt+fatJet2Pt": [350, 9999],
    "fatJet1MassSD": [105, 135],
    "fatJet2MassSD": [105, 135],
}

hh4bcf, hh4b_kin_cuts_evts = cutflow_func(hh4b_var_cuts, events_bbs_ak8_cuts)

ak.sum(events_bbs_ak8_cuts["HH4b"]['weight'])

cut_labels = ['Jet1 pT > 310',
                'Jet2 pT > 310',
                'At least 1 jet pT > 350',
                'Jet1 MassSD > 105',
                'Jet1 MassSD < 135',
                'Jet2 MassSD > 105',
                'Jet2 MassSD < 135',
            ]

# del(hh4bcf['HHbbVV4q'])
del(hh4bcf['data'])

cftable = cftable_func(hh4bcf, hh4b_var_cuts, cut_labels)
cftable


bVfpr, bVtpr, bVauc = plot_single_roc(hh4b_kin_cuts_evts, "fatJet1PNetXbb", 'Hbb', 'Txbb', 'pnetxbb_kin_cuts_roc')
bbfpr, bbtpr, bbauc = plot_single_roc(hh4b_kin_cuts_evts, ["fatJet1PNetXbb", "fatJet2PNetXbb"], 'Hbb', 'Txbb', 'pnetxbb_hh4b_kin_cuts_roc', sig='HH4b')


plt.semilogx(bbfpr, bbtpr, label=f'4b AUC = {bbauc:.2f}')
plt.semilogx(bVfpr, bVtpr, label=f'bbVV AUC = {bVauc:.2f}')
plt.legend(loc='upper left', fancybox=True)
plt.xlabel('QCD Background Efficiency')
plt.ylabel('Hbb Signal Efficiency')
plt.xlim(1e-4, 1)
plt.title(f"ParticleNet PNetXbb ROC Curve")
plt.savefig(f"figs/pnetxbb_comp_bbVV_4b_roc.pdf")
plt.show()
