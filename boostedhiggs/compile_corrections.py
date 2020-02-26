#!/usr/bin/env python
import json
import gzip
import uproot
import numexpr
import numpy as np
from coffea import hist, lookup_tools
from coffea.util import load, save
from coffea.hist import plot

corrections = {}

def read_xsections(filename):
    out = {}
    with open(filename) as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            dataset, xsexpr, *_ = line.split()
            try:
                xs = float(numexpr.evaluate(xsexpr))
            except:
                print("numexpr evaluation failed for line: %s" % line)
                raise
            if xs <= 0:
                warnings.warn("Cross section is <= 0 in line: %s" % line, RuntimeWarning)
            out[dataset] = xs
    return out

corrections['xsections'] = read_xsections("data/xSections.dat")

pileup_corr = load('data/pileup_mc.coffea')
del pileup_corr['data_obs_jet']
del pileup_corr['data_obs_mu']
with uproot.open("data/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root") as fin_pileup:
    norm = lambda x: x / x.sum()
    data_pu = norm(fin_pileup["pileup"].values)
    data_pu_puUp = norm(fin_pileup["pileup_plus"].values)
    data_pu_puDown = norm(fin_pileup["pileup_minus"].values)

    # https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py
    mc_pu = np.array([
            1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
            0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
            0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
            0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
            0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
            0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
            0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
            0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
            0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
            0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
            0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
            0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
            0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
            3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
            1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 
            ])

    # https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
    '''
    mc_pu = np.array([
              3.39597497605e-05,
              6.63688402133e-06,
              1.39533611284e-05,
              3.64963078209e-05,
              6.00872171664e-05,
              9.33932578027e-05,
              0.000120591524486,
              0.000128694546198,
              0.000361697233219,
              0.000361796847553,
              0.000702474896113,
              0.00133766053707,
              0.00237817050805,
              0.00389825605651,
              0.00594546732588,
              0.00856825906255,
              0.0116627396044,
              0.0148793350787,
              0.0179897368379,
              0.0208723871946,
              0.0232564170641,
              0.0249826433945,
              0.0262245860346,
              0.0272704617569,
              0.0283301107549,
              0.0294006137386,
              0.0303026836965,
              0.0309692426278,
              0.0308818046328,
              0.0310566806228,
              0.0309692426278,
              0.0310566806228,
              0.0310566806228,
              0.0310566806228,
              0.0307696426944,
              0.0300103336052,
              0.0288355370103,
              0.0273233309106,
              0.0264343533951,
              0.0255453758796,
              0.0235877272306,
              0.0215627588047,
              0.0195825559393,
              0.0177296309658,
              0.0160560731931,
              0.0146022004183,
              0.0134080690078,
              0.0129586991411,
              0.0125093292745,
              0.0124360740539,
              0.0123547104433,
              0.0123953922486,
              0.0124360740539,
              0.0124360740539,
              0.0123547104433,
              0.0124360740539,
              0.0123387597772,
              0.0122414455005,
              0.011705203844,
              0.0108187105305,
              0.00963985508986,
              0.00827210065136,
              0.00683770076341,
              0.00545237697118,
              0.00420456901556,
              0.00367513566191,
              0.00314570230825,
              0.0022917978982,
              0.00163221454973,
              0.00114065309494,
              0.000784838366118,
              0.000533204105387,
              0.000358474034915,
              0.000238881117601,
              0.0001984254989,
              0.000157969880198,
              0.00010375646169,
              6.77366175538e-05,
              4.39850477645e-05,
              2.84298066026e-05,
              1.83041729561e-05,
              1.17473542058e-05,
              7.51982735129e-06,
              6.16160108867e-06,
              4.80337482605e-06,
              3.06235473369e-06,
              1.94863396999e-06,
              1.23726800704e-06,
              7.83538083774e-07,
              4.94602064224e-07,
              3.10989480331e-07,
              1.94628487765e-07,
              1.57888581037e-07,
              1.2114867431e-07,
              7.49518929908e-08,
              4.6060444984e-08,
              2.81008884326e-08,
              1.70121486128e-08,
              1.02159894812e-08])
              mc_pu = np.r_[mc_pu, np.zeros(1)]
              '''
    mc_pu = np.r_[mc_pu, np.zeros(26)] 
    mask = mc_pu > 0.
    corr = data_pu.copy()
    corr_puUp = data_pu_puUp.copy()
    corr_puDown = data_pu_puDown.copy()
    corr[mask] /= mc_pu[mask]
    corr_puUp[mask] /= mc_pu[mask]
    corr_puDown[mask] /= mc_pu[mask]
    pileup_corr = lookup_tools.dense_lookup.dense_lookup(corr, fin_pileup["pileup"].edges)
    pileup_corr_puUp = lookup_tools.dense_lookup.dense_lookup(corr_puUp, fin_pileup["pileup"].edges)
    pileup_corr_puDown = lookup_tools.dense_lookup.dense_lookup(corr_puDown, fin_pileup["pileup"].edges)

corrections['2017_pileupweight'] = pileup_corr
corrections['2017_pileupweight_puUp'] = pileup_corr_puUp
corrections['2017_pileupweight_puDown'] = pileup_corr_puDown

save(corrections, 'data/corrections.coffea')
