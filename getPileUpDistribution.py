#!/usr/bin/env python

import sys, os, string, re, copy
from optparse import OptionParser
from ROOT import *


def main():
  # usage description
  usage = "Usage: getPileUpDistribution.py [options] \nExample: ./getPileUpDistribution.py -i PileUp_dist.root"

  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-i", "--input_file", dest="input_file",
                    help="Input pile-up distribution file",
                    metavar="INPUT_FILE")

  parser.add_option("-o", "--output_file", dest="output_file",
                    help="Output pile-up distribution file",
                    metavar="OUTPUT_FILE")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not options.input_file:
    print usage
    sys.exit()

  file = TFile(options.input_file)
  histo = file.Get('pileup')
  bin_content = []

  for i in range(1,histo.GetNbinsX()):
    bin_content.append(histo.GetBinContent(i))

  print bin_content


if __name__ == "__main__":
  main()
