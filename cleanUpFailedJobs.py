#!/usr/bin/env python

import sys, os, subprocess, shutil, string, re
from optparse import OptionParser


def main():
  # usage description
  usage = "Usage: cleanUpFailedJobs.py [options] \nExample: ./cleanUpFailedJobs.py -w crab_workdir -j 1,5-20"

  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-w", "--workdir", dest="workdir",
                    help="CRAB working directory",
                    metavar="WORKDIR")

  parser.add_option("-j", "--jobs", dest="jobs",
                    help="List of jobs or job ranges",
                    metavar="JOBS")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.workdir and options.jobs):
    print usage
    sys.exit()

  workdir = options.workdir

  # redefine workdir as an absolute path (if not defined in such form already)
  if not re.search("^/", workdir):
    workdir = os.path.join(os.getcwd(),workdir)

  if not os.path.isdir(workdir):
    print '%s not found'%workdir
    sys.exit(1)

  # path to the res/ subdirectory
  res_path = os.path.join(workdir,'res')

  # change directory to the res/ subdirectory
  os.chdir(res_path)

  input_jobs = options.jobs
  final_jobs = []

  for j in input_jobs.split(','):
    job_range = j.split('-')
    if len(job_range) > 2:
      print '%s is not a valid job range'%j
      sys.exit(1)
    for i in range(0,len(job_range)):
      if not job_range[i].isdigit():
        print '%s is not a valid job number'%job_range[i]
        sys.exit(1)
    for k in range(int(job_range[0]),int(job_range[-1])+1):
      final_jobs.append(k)
    
  #print final_jobs
  
  print 'Starting cleanup...'

  files = 'CMSSW_JOBNUMBER.stdout CMSSW_JOBNUMBER.stderr crab_fjr_JOBNUMBER.xml'

  for job in final_jobs:
    cmd = ('rm -f ' + files).replace('JOBNUMBER',str(job))
    proc = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
    output = proc.communicate()[0]
    if proc.returncode != 0:
      print output
      sys.exit(1)

  print 'Done'


if __name__ == "__main__":
  main()

