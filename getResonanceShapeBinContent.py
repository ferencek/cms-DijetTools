#!/usr/bin/env python

import os
from ROOT import *

path_prefix = 'RSG_resonance_shape_files'
histo_name = 'myAnalyzer/cutHisto_allPreviousCuts________x_bbbar'

files = [
    'RSG_M_500__histograms.root',
    'RSG_M_700__histograms.root',
    'RSG_M_1200__histograms.root',
    'RSG_M_2000__histograms.root',
    'RSG_M_3500__histograms.root'
]

array_names = [
    'y500[50]',
    'y700[50]',
    'y1200[50]',
    'y2000[50]',
    'y3500[50]'
]

for i, fl in enumerate(files):
  file = TFile(os.path.join(path_prefix,fl))
  histo = file.Get(histo_name)

  array_string = array_names[i] + ' = {'
  
  for b in range(1,histo.GetNbinsX()+1):
    array_string += str(histo.GetBinContent(b))
    if b<histo.GetNbinsX(): array_string += ', '

  array_string += '}'
  print array_string
