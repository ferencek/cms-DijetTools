#!/usr/bin/env python

import os
from ROOT import *

path_prefix = 'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets'
#path_prefix = 'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets'

histo_name_postfix = "bbbar"
#histo_name_postfix = "ccbar"
#histo_name_postfix = "qqbarlight"
#histo_name_postfix = "gg"

histo_names = [
    'myAnalyzer/x_dist_500',
    'myAnalyzer/x_dist_700',
    'myAnalyzer/x_dist_1200',
    'myAnalyzer/x_dist_2000',
    'myAnalyzer/x_dist_3500'
]

files = [
    'RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root'
]

array_names = [
    'double y500[50]',
    'double y700[50]',
    'double y1200[50]',
    'double y2000[50]',
    'double y3500[50]'
]

for i, fl in enumerate(files):
  file = TFile(os.path.join(path_prefix,fl))
  histo = file.Get(histo_names[i] + '_' + histo_name_postfix)

  array_string = array_names[i] + ' = {'
  
  for b in range(1,histo.GetNbinsX()+1):
    array_string += str(histo.GetBinContent(b))
    if b<histo.GetNbinsX(): array_string += ', '

  array_string += '};'
  print array_string
