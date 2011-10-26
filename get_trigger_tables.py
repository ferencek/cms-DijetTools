#!/usr/bin/env python

import os, sys, subprocess, string, re

trigger_menus = {
 # To get the trigger menu history in a given run range, run
 # ./getHLTkey.py --firstRun=170249 --lastRun=173692 --perKey
 # For more info, see https://twiki.cern.ch/twiki/bin/view/CMS/HLTriggerTools and
 #                    http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/CommonToolsV2/test/TriggerInfo/README
 
 ## /Jet/Run2011A-May10ReReco-v1/AOD
 '/cdaq/physics/Run2011/5e32/v4.2/HLT/V6':     [160431,160432,160433,160442,160443,160444,160445,160446,160447,160449,160450,160454,160455,160456,160462,160463,160466,160467],
 '/cdaq/physics/Run2011/5e32/v4.2/HLT/V7':     [160497,160498,160499,160577,160578,160579],
 '/cdaq/physics/Run2011/5e32/v4.2/HLT/V8':     [160808,160815],
 '/cdaq/physics/Run2011/5e32/v4.3/HLT/V3':     [160819],
 '/cdaq/physics/Run2011/5e32/v4.3/HLT/V4':     [160827,160835],
 '/cdaq/physics/Run2011/5e32/v5.1/HLT/V3':     [160853,160871,160872,160873,160874,160875,160876,160877],
 '/cdaq/physics/Run2011/5e32/v5.2/HLT/V2':     [160888],
 '/cdaq/physics/Run2011/5e32/v5.2/HLT/V5':     [160890],
 '/cdaq/physics/Run2011/5e32/v5.2/HLT/V6':     [160894,160898],
 '/cdaq/physics/Run2011/5e32/v5.2/HLT/V7':     [160907,160911,160913,160914,160915,160916,160935,160936,160937,160938,160939,160940,160942,160943],
 '/cdaq/physics/Run2011/5e32/v5.3/HLT/V1':     [160954,160955,160956,160957,160994,160998,161008,161016,161020],
 '/cdaq/physics/Run2011/5e32/v5.3/HLT/V2':     [161103,161106,161107,161113,161116,161117,161119,161156,161165,161176],
 '/cdaq/physics/Run2011/5e32/v6.1/HLT/V1':     [161216],
 '/cdaq/physics/Run2011/5e32/v6.1/HLT/V3':     [161217],
 '/cdaq/physics/Run2011/5e32/v6.1/HLT/V5':     [161222,161223,161224,161233],
 '/cdaq/physics/Run2011/5e32/v6.1/HLT/V6':     [161310,161311,161312],
 '/cdaq/physics/Run2011/5e32/v6.2/HLT/V2':     [162718,162733,162739,162742,162760,162762,162765],
 '/cdaq/physics/Run2011/5e32/v6.2/HLT/V3':     [162803,162808,162810,162811,162822,162825,162826,162827,162828,162909,162917,162924,162925,162926],
 '/cdaq/physics/Run2011/5e32/v6.2/HLT/V4':     [162929,163045,163046,163069,163071,163072,163078,163232,163233,163234,163235,163237,163238,163252,163255,163261],
 '/cdaq/physics/Run2011/5e32/v8.1/HLT/V5':     [163269,163270],
 '/cdaq/physics/Run2011/5e32/v8.1/HLT/V6':     [163286,163289],
 '/cdaq/physics/Run2011/5e32/v8.1/HLT/V8':     [163296,163297,163300,163301,163302,163308,163332,163333,163334],
 '/cdaq/physics/Run2011/5e32/v8.2/HLT/V3':     [163337,163338,163339,163340,163358,163369,163370,163371,163372,163374,163375,163376,163378,163385,163387,163402,163475,163476,163478,163479,163480,163481,163482,163483,163581,163582,163583,163584,163585,163586,163587,163588,163589,163596,163630,163655,163657,163658,163659,163660,163661,163662,163663,163664,163665,163668],
 '/cdaq/physics/Run2011/5e32/v8.3/HLT/V2':     [163738,163753],
 '/cdaq/physics/Run2011/5e32/v8.3/HLT/V3':     [163754],
 '/cdaq/physics/Run2011/5e32/v8.3/HLT/V4':     [163757,163758,163759,163760,163761,163763,163765,163795,163796,163817,163869],

 ## /Jet/Run2011A-PromptReco-v4/AOD
 '/cdaq/physics/Run2011/1e33/v1.3/HLT/V2':     [165088,165098,165099,165102,165103,165120,165121],
 '/cdaq/physics/Run2011/1e33/v1.3/HLT/V6':     [165205],
 '/cdaq/physics/Run2011/1e33/v1.3/HLT/V7':     [165208],
 '/cdaq/physics/Run2011/1e33/v1.3/HLT/V12':    [165364],
 '/cdaq/physics/Run2011/1e33/v1.3/HLT/V13':    [165400,165402,165415,165467,165472,165486,165487,165506,165514,165523,165525,165529,165537,165542,165548,165558,165567,165570,165617,165619,165620,165633],
 '/cdaq/physics/Run2011/1e33/v2.3/HLT/V1':     [165970,165993],
 '/cdaq/physics/Run2011/1e33/v2.3/HLT/V3':     [166010,166011,166033,166034,166049,166051,166052,166149,166150],
 '/cdaq/physics/Run2011/1e33/v2.4/HLT/V2':     [166161,166163,166164],
 '/cdaq/physics/Run2011/1e33/v2.5/HLT/V1':     [166346],
 '/cdaq/physics/Run2011/1e33/v2.4/HLT/V4':     [166374,166377,166379,166380,166394,166408,166429,166438,166462,166486,166502,166512],
 '/cdaq/physics/Run2011/1e33/v2.4/HLT/V5':     [166514,166530,166554,166563,166565,166699,166701,166763,166781,166782],
 '/cdaq/physics/Run2011/1e33/v2.4/HLT/V6':     [166784,166787],
 '/cdaq/physics/Run2011/1e33/v2.4/HLT/V8':     [166839,166841,166842,166859,166860,166861,166862,166863,166864,166888,166889,166890,166893,166894,166895,166911,166921,166922,166923,166946,166950,166960,166966,166967],
 '/cdaq/physics/Run2011/1.4e33/v1.1/HLT/V1':   [167039,167041,167043],
 '/cdaq/physics/Run2011/1.4e33/v1.2/HLT/V1':   [167078,167098,167102,167103,167151,167281,167282,167284],
 '/cdaq/physics/Run2011/1.4e33/v1.2/HLT/V3':   [167551,167673,167674,167675,167676,167740,167746,167754,167784,167786,167807,167830,167898,167913],

 ## /Jet/Run2011A-05Aug2011-v1/AOD and /Jet/Run2011A-PromptReco-v6/AOD
 '/cdaq/physics/Run2011/2e33/v1.1/HLT/V1':     [170249,170255,170286,170292,170298,170303,170304,170307,170348,170354,170376,170378,170380,170382,170397,170406],
 '/cdaq/physics/Run2011/2e33/v1.1/HLT/V2':     [170452,170527,170722,170759],
 '/cdaq/physics/Run2011/2e33/v1.2/HLT/V1':     [170826,170842,170854,170876,170896,170899,170901,171050,171091,171098,171102,171106,171116,171117,171156,171178],
 '/cdaq/physics/Run2011/2e33/v1.2/HLT/V4':     [171219,171274,171282,171315,171369,171446,171484,171578,171812,171875,171876,171879,171880,171890,171895,171897,171921,171926,172014,172024,172033],
 '/cdaq/physics/Run2011/2e33/v1.2/HLT/V5':     [172163,172208,172252,172254,172255,172268,172286,172389,172399,172400,172401,172411],
 '/cdaq/physics/Run2011/2e33/v1.2/HLT/V7':     [172478,172485,172488,172495,172497,172507,172516,172619,172620,172630,172635,172778,172791,172798,172799,172801,172802,172819,172822,172824,172847,172865,172868,172949,172951,172952,172992,172998,172999,173198],
 '/cdaq/physics/Run2011/3e33/v1.1/HLT/V1':     [173236,173240,173241],
 '/cdaq/physics/Run2011/3e33/v1.1/HLT/V3':     [173243,173244],
 '/cdaq/physics/Run2011/3e33/v1.1/HLT/V4':     [173380,173381,173389,173406,173430,173431,173438,173439],
 '/cdaq/physics/Run2011/3e33/v1.2/HLT/V1':     [173657,173658,173659,173660,173661,173662,173663,173664,173692],

 ## /Jet/Run2011B-PromptReco-v1/AOD
 '/cdaq/physics/Run2011/3e33/v2.0/HLT/V4':     [175832,175833,175834,175835,175837],
 '/cdaq/physics/Run2011/3e33/v2.0/HLT/V7':     [175857,175858,175860,175863,175865,175866,175877,175881,175886,175887,175888,175906,175910,175921],
 '/cdaq/physics/Run2011/3e33/v2.1/HLT/V1':     [175971,175973,175974,175975,175976,175990,176023],
 '/cdaq/physics/Run2011/3e33/v2.1/HLT/V2':     [176161,176163,176165,176166,176167,176169,176201,176202,176206,176207],
 '/cdaq/physics/Run2011/3e33/v2.2/HLT/V3':     [176286,176289,176304,176308,176309],
 '/cdaq/physics/Run2011/3e33/v2.3/HLT/V2':     [176461,176463,176464,176465,176466,176467,176468,176469,176470],
 '/cdaq/physics/Run2011/3e33/v3.0/HLT/V2':     [176545,176547,176548],
 '/cdaq/physics/Run2011/3e33/v3.1/HLT/V1':     [176697,176701,176702,176765,176771,176795,176796,176797,176799,176801,176805,176807,176841,176842,176844,176848,176850,176860,176868,176885,176886,176889,176928,176929,176933,176959,176982,177053],
 '/cdaq/physics/Run2011/3e33/v4.0/HLT/V2':     [177074,177088,177095,177096,177131,177138,177139,177140,177141,177183,177184],
 '/cdaq/physics/Run2011/3e33/v4.0/HLT/V3':     [177201],
 '/cdaq/physics/Run2011/3e33/v4.0/HLT/V5':     [177222,177313,177317,177318,177319,177449,177452,177507,177509,177515,177730,177776,177782,177783,177785,177786,177788,177789,177790,177791,177875,177878],
 '/cdaq/physics/Run2011/3e33/v4.0/HLT/V6':     [177718,177719],
 '/cdaq/physics/Run2011/3e33/v5.0/HLT/V1':     [178078,178098,178099,178100,178101,178102,178110,178116,178151,178160,178162],
 '/cdaq/physics/Run2011/5e33/v1.4/HLT/V3':     [178420,178421,178424,178479],
 '/cdaq/physics/Run2011/5e33/v1.4/HLT/V4':     [178667,178675,178677,178703,178708],
 '/cdaq/physics/Run2011/5e33/v1.4/HLT/V5':     [178712,178724,178731,178738,178786,178803,178840,178854,178866,178871,178920,178970,178985],
}

#dataset = 'Jet'
dataset = 'HT'

for key in trigger_menus:
  print key
  idx = trigger_menus.keys().index(key)
  fileName = '%i-%i_%s%s.txt'%(trigger_menus[key][0],trigger_menus[key][-1],dataset,string.replace(key,'/','__'))
  outputFile = open(fileName,'w')

  cmd = 'edmConfigFromDB --orcoff --configName %s | hltDumpStream'%key

  proc = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
  output = proc.communicate()[0]
  if proc.returncode != 0:
    print output
    sys.exit(1)
  #print output

  doPrintOut = False
  firstLine = True

  for line in output.splitlines():
    if ( doPrintOut and (string.find(line,'dataset') != -1 or string.find(line,'stream') != -1)): doPrintOut = False
    if ( doPrintOut and firstLine ):
      trigger = line.split()[0]
      if trigger == '~':
        trigger = line.split()[1]
      cmd_trg = 'edmConfigFromDB --orcoff --configName %s --format summary.ascii --paths %s'%(key,trigger)
      proc_trg = subprocess.Popen( cmd_trg, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
      output_trg = proc_trg.communicate()[0]
      if proc_trg.returncode != 0:
        print output_trg
        sys.exit(1)
      outputFile.write(output_trg)
      firstLine = False
    if doPrintOut: outputFile.write(line + '\n')
    if re.search(('dataset %s$'%(dataset)),line): doPrintOut = True

  outputFile.close()
