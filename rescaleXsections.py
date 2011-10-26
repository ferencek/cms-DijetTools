#!/usr/bin/env python

import sys, os, string, re
from optparse import OptionParser
from math import *


def main():
  # usage description
  usage = "Usage: rescaleXsections.py [options] \nExample: ./rescaleXsections.py -i xsections.txt -o xsections_rescaled.txt -r 1.255 -d QCD_Pt.+TuneZ2_7TeV_pythia6"

  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-i", "--input_xsection_file", dest="input_xsection_file",
                    help="Input cross section file",
                    metavar="INPUT_XSECTION_FILE")

  parser.add_option("-o", "--output_xsection_file", dest="output_xsection_file",
                    help="Output cross section file",
                    metavar="OUTPUT_XSECTION_FILE")

  parser.add_option("-r", "--rescale_factor", dest="rescale_factor",
                    help="Rescale factor",
                    metavar="RESCALE_FACTOR")

  parser.add_option("-d", "--dataset_regex", dest="dataset_regex",
                    help="Regular expression matching dataset(s) to be rescaled",
                    metavar="DATASET_REGEX")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.input_xsection_file and options.output_xsection_file and options.rescale_factor and options.dataset_regex):
    print usage
    sys.exit()

  # open and read the input_xsection_file
  input_xsection_file = open(options.input_xsection_file,'r')
  input_xsection_lines = input_xsection_file.readlines()

  # create new cross section file
  output_xsection_file = open(options.output_xsection_file,'w')
  
  for line in input_xsection_lines:

    if( re.search(options.dataset_regex, line) ):
      line_elements = line.split()
      new_xsection = float(line_elements[1])*float(options.rescale_factor)
      new_line = line.replace(line_elements[1],str(new_xsection))
      output_xsection_file.write(new_line) 
    else:
      output_xsection_file.write(line)

  input_xsection_file.close()
  output_xsection_file.close()

  print ''
  print "New cross section file: " + options.output_xsection_file
  print ''


if __name__ == "__main__":
  main()

