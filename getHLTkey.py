#!/usr/bin/env python 
import os
from sys import stderr, exit
from RecoLuminosity.LumiDB import selectionParser
import commands
import fileinput
import string

#
# HLT key and run query: E.P., 27 July 2010
# Specific trigger information added: bdahmes, 23 Aug 2010 
# query to the Run Registry taken from a script by Giovanni Petrucianni
#


from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options]")
parser.add_option("--firstRun",  dest="firstRun",  help="first run", type="int", metavar="RUN", default="1")
parser.add_option("--lastRun",   dest="lastRun",   help="last run",  type="int", metavar="RUN", default="9999999")
parser.add_option("--groupName", dest="groupName", help="select runs of name like NAME", metavar="NAME", default="Collisions%")
parser.add_option("--rrurl",     dest="rrurl",     help="run registry xmlrpc url", metavar="URL", default="http://pccmsdqm04.cern.ch/runregistry/xmlrpc")
parser.add_option("--HLTkey",    dest="HLTkey",    help="name of the HLTkey e.g. /cdaq/physics/Run2010/v3.1/HLT_1.6E30/V1",metavar="HLT")
parser.add_option("--perKey",    action="store_true",default=False,dest="perKey",help="list the runs per HLT key",metavar="perKey")
parser.add_option("--json",      dest="json",      help="JSON file for good runs", type="string", default="None")
parser.add_option("--HLTpaths",  dest="HLTpaths",  help="Comma-separated list of HLT paths (a,b,c)", default="None")
parser.add_option("--latex",     action="store_true", default=False, dest="latex",help="table output in LaTeX format",metavar="latex")
parser.add_option("--runsPerRow",dest="runsPerRow",default=-1, type="int",help="Runs printed before new table row",metavar="NUMBER")
(options, args) = parser.parse_args()

def max(i,j):
    if i > j:
        return i
    else:
        return j

def queryRR():
    stderr.write("Querying run registry for range [%d, %d], group name like %s ...\n" % (options.firstRun, options.lastRun, options.groupName))
    import xmlrpclib
    import xml.dom.minidom
    server = xmlrpclib.ServerProxy(options.rrurl)
    run_data = server.DataExporter.export('RUN', 'GLOBAL', 'xml_datasets', "{runNumber} >= %d AND {runNumber} <= %d AND {groupName} like '%s' AND {datasetName} = '/Global/Online/ALL'"  % (options.firstRun, options.lastRun, options.groupName))
    ret = {}
    xml_data = xml.dom.minidom.parseString(run_data)
    xml_runs = xml_data.documentElement.getElementsByTagName("RUN_DATASET")
    for xml_run in xml_runs:
        ret[xml_run.getElementsByTagName("RUN_NUMBER")[0].firstChild.nodeValue] = xml_run.getElementsByTagName("RUN_HLTKEY")[0].firstChild.nodeValue
    return ret

tempHLTfile = 'tempHLT.log'
def pathInfo(path):
    searchString = " " + path + " "
    result = [line for line in
              open(tempHLTfile) if searchString in line]
    if len(result) > 0:
        return result[0]
    else:
        return ""

runKeys = queryRR()
prescaleTable = {}
runs = runKeys.keys(); runs.sort()
paths = []
runsPerRow = int(options.runsPerRow)

if options.HLTpaths != "None":
    pathString = options.HLTpaths
    beginString = 0
    ctr = 0 
    while beginString < len(pathString):
        ctr += 1
        nextString = pathString[beginString:-1].find(',')
        endString = beginString + nextString
        if endString < beginString:
            endString = len(pathString)
        newPath = pathString[beginString:endString]
        beginString = endString + 1
        paths.append(newPath)

if options.json != "None":
    jfile  = open(options.json,'r')
    filteredRunList = []
    goodLS = ''
    parsingResult = '' 
    goodLS = jfile.read()
    parsingResult = selectionParser.selectionParser(goodLS)
    if not parsingResult:
        print 'Failed to parse the input JSON file',ifilename
        raise 
    goodRuns = parsingResult.runs()
    for run in runs:
        key = runKeys[run]
        for anotherRun in goodRuns:
            if int(run) == anotherRun:
                filteredRunList.append(run)
    runs = filteredRunList

if options.perKey:
    runsPerKey={}
    for run in runs:
        key = runKeys[run]
        if not key in runsPerKey.keys():
            tmpruns=[]
            tmpruns.append(run)
            runsPerKey[key] = tmpruns
        else:
            runsPerKey[key].append(run)

	theKeys = runsPerKey.keys()
    # Resort the menus based on run number
    resortedKeys = []
    totalKeys = len(theKeys)
    while len(resortedKeys) != totalKeys:
        earliestRun = 999999
        earliestKey = ''
        for key in theKeys:
            first = int(runsPerKey[key][0])
            if first < earliestRun:
                earliestKey = key
                earliestRun = first
        resortedKeys.append(earliestKey)
        theKeys.remove(earliestKey)
    theKeys = resortedKeys

    # Look for special features based on interesting paths
    # NOTE: DQM is ignored (avoids seeing path appear more than once)
    featuresPerKey = {}
    if options.HLTpaths != "None":
        tempHLTfile = "tempHLT.log"
        for key in theKeys:
            mySetup = "touch " + tempHLTfile + " ; rm " + tempHLTfile
            os.system(mySetup)
            myGet = "edmConfigFromDB --orcoff --configName " + key + " --noes --services PrescaleService | hltDumpStream > " + tempHLTfile
            os.system(myGet)
            infoString = ''
            for path in paths:
                info = pathInfo(path).split()
                i = 1
                foundPrescale = False
                while i < len(info):
                    if info[i].isdigit():
                        if foundPrescale:
                            info.pop(i)
                        else:
                            foundPrescale = True
                            i += 1
                    else:
                        i += 1
                if len(info) > 4:
                    lastColumn = info[len(info)-1].split()
                    info[3] = lastColumn[len(lastColumn)-1]
                    while len(info) > 4:
                        dummy = info.pop()
                for item in info:
                    infoString += item.strip() + ':'
                infoString = infoString[:len(infoString)-1] + ';'
            # if options.latex:
            #     infoString = string.replace(infoString,'_','\\_')
            # print infoString
            featuresPerKey[key] = infoString
            myCleanup = "rm " + tempHLTfile
            os.system(myCleanup)
    
    pathFeatures = {}
    pathPrescale = {}
    pathL1seed = {}
    updatedFeatures = {}
    updatedPrescale = {}
    updatedL1seed = {}
    for path in paths:
        pathFeatures[path]    = ''
        pathPrescale[path]    = ''
        pathL1seed[path]      = ''
        updatedFeatures[path] = 'X'
        updatedPrescale[path] = ''
        updatedL1seed[path]   = ''

    beginString = ''
    separatorString = ''
    endString = ''
    if options.latex:
        separatorString = '&'
        endString = '\\\\'

    maxLen = 0 
    for key in theKeys:
        theLen = len(key)
        if theLen > maxLen:
            maxLen = theLen

    for key in theKeys:
        runList = []
        theruns = runsPerKey[key]
        topr = ""
        for r in theruns:
            topr += r + ","
            if runsPerRow > 0 and topr.count(',') == runsPerRow and theruns.index(r) != (len(theruns) - 1):
                runList.append(topr)
                topr = ""
        topr = topr[:len(topr)-1]
        runList.append(topr)

        # Always assume a path has been disabled, and prepare to be wrong
        for path in paths:
            updatedFeatures[path] = 'X'

        if options.HLTpaths != "None":
            parseMenuInfo = featuresPerKey[key].split(';',len(paths))
            for pathBlock in parseMenuInfo:
                parsePathInfo = pathBlock.split(':')
                if len(parsePathInfo) > 1:
                    thePath = parsePathInfo[len(parsePathInfo)-3]
                    # thePath = string.replace(thePath,'\\_','_')
                    updatedFeatures[thePath] = 'U'
                    if parsePathInfo[0] == '~':
                        updatedFeatures[thePath] = 'G'
                    updatedPrescale[thePath] = parsePathInfo[len(parsePathInfo)-2]
                    updatedL1seed[thePath] = parsePathInfo[len(parsePathInfo)-1]

        xtraSpace = ''
        for i in range(0,maxLen-len(key)):
            xtraSpace += ' '
        outputKey = key
        if options.latex:
            outputKey = string.replace(key,'_','\\_')
        print beginString,outputKey,xtraSpace,separatorString,
        featureString = []
        for path in paths:
            stringPrefix = ''
            if updatedFeatures[path] != pathFeatures[path]:
                pathFeatures[path] = updatedFeatures[path]
                stringPrefix += '(' + pathFeatures[path] + ')' 
            L1andPrescale = ''
            if updatedPrescale[path] != pathPrescale[path] or updatedL1seed[path] != pathL1seed[path]:
                if updatedPrescale[path] == pathPrescale[path]:
                    L1andPrescale = '(' + updatedL1seed[path] + ')'
                elif updatedL1seed[path] == pathL1seed[path]:
                    L1andPrescale = '(' + updatedPrescale[path] + ')'
                else:
                    L1andPrescale = '(' + updatedL1seed[path] + ',' + updatedPrescale[path] + ')'
                pathPrescale[path] = updatedPrescale[path]
                pathL1seed[path]   = updatedL1seed[path]
            pathString = ''
            if len(stringPrefix) > 0 or len(L1andPrescale) > 0:
                pathString = path
                if len(stringPrefix) > 0:
                    pathString = stringPrefix + pathString
                if len(L1andPrescale) > 0:
                    pathString += L1andPrescale
                if options.latex:
                    pathString = string.replace(pathString,'_','\\_')
            featureString.append(pathString)
        oldLength = len(featureString)
        for i in range(0,oldLength):
            idx = oldLength - 1 - i
            if len(featureString[idx]) == 0:
                featureString.pop(idx)
        for i in range(0,max(len(featureString),len(runList))):
            if i > 0:
                print beginString,' ', 
                for j in range(0,maxLen):
                    print '',
                print separatorString,
            if i >= len(runList):
                print featureString[i],separatorString,' ',
            elif i >= len(featureString):
                print ' ', separatorString,runList[i],
            else:
                print featureString[i],separatorString,runList[i],
            print endString
            
    exit(1)
			
if options.HLTkey:
	HLTkey = options.HLTkey
	print "List of runs taken with HLT key = ",HLTkey 
for run in runs:
    key = runKeys[run]

    if not options.HLTkey:
       print run,key	
    else:
	if key == options.HLTkey:
	   print run
