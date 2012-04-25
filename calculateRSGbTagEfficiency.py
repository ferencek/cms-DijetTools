#!/usr/bin/env python

import os, sys
from ROOT import *
from array import array
import copy


gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadTopMargin(0.07)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.13)
gStyle.SetPadRightMargin(0.07)

graph_eff0_array = []
graph_eff0_syst_array = []
graph_eff1_array = []
graph_eff1_syst_array = []
graph_eff2_array = []
graph_eff2_syst_array = []

def eff(inputDir, final_state, title, outputFilename, outputFormat = 'png'):

  files = [
    'RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root'
  ]

  masses = [500., 700., 1200., 2000., 3500.]

  mass_array = array('d')
  mass_err_array = array('d')
  eff0_array = array('d')
  eff0_err_array = array('d')
  eff1_array = array('d')
  eff1_err_array = array('d')
  eff2_array = array('d')
  eff2_err_array = array('d')

  for i, fl in enumerate(files):

    file = TFile(os.path.join(inputDir,fl))

    histo_N0 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_0tag')
    histo_N1 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_1tag')
    histo_N2 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_2tag')

    binMin = histo_N0.GetXaxis().FindBin(0.3*masses[i])
    binMax = histo_N0.GetXaxis().FindBin(1.3*masses[i])

    N0 = histo_N0.Integral(binMin,binMax)
    N1 = histo_N1.Integral(binMin,binMax)
    N2 = histo_N2.Integral(binMin,binMax)

    N = N0 + N1 + N2

    eff0 = N0/N
    eff1 = N1/N
    eff2 = N2/N
    
    mass_array.append(masses[i])
    mass_err_array.append(0.)

    eff0_array.append( eff0 )
    eff0_err_array.append( sqrt( (eff0*(1.-eff0))/N ) )
    eff1_array.append( eff1 )
    eff1_err_array.append( sqrt( (eff1*(1.-eff1))/N ) )
    eff2_array.append( eff2 )
    eff2_err_array.append( sqrt( (eff2*(1.-eff2))/N ) )

  print "masses:"
  print mass_array
  print "eff 0-tag:"
  print eff0_array
  #print "eff 1-tag:"
  #print eff1_array
  print "eff 2-tag:"
  print eff2_array
    
  graph_eff0 = TGraphErrors(len(mass_array),mass_array,eff0_array,mass_err_array,eff0_err_array)
  graph_eff0.SetMarkerStyle(24)
  graph_eff0.SetLineWidth(2)
  graph_eff0.SetLineStyle(1)
  graph_eff0.SetLineColor(1)
  graph_eff0.SetTitle(title)
  graph_eff0.GetXaxis().SetTitle("Resonance Mass [GeV]")
  graph_eff0.GetYaxis().SetTitle("Efficiency")
  graph_eff0.GetYaxis().SetRangeUser(0.,1.)

  graph_eff1 = TGraphErrors(len(mass_array),mass_array,eff1_array,mass_err_array,eff1_err_array)
  graph_eff1.SetMarkerStyle(25)
  graph_eff1.SetLineWidth(2)
  graph_eff1.SetLineStyle(5)
  graph_eff1.SetLineColor(2)

  graph_eff2 = TGraphErrors(len(mass_array),mass_array,eff2_array,mass_err_array,eff2_err_array)
  graph_eff2.SetMarkerStyle(26)
  graph_eff2.SetLineWidth(2)
  graph_eff2.SetLineStyle(7)
  graph_eff2.SetLineColor(4)

  c = TCanvas("c", "",800,800)
  c.cd()

  graph_eff0.Draw("ALP")
  graph_eff1.Draw("LP")
  graph_eff2.Draw("LP")

  legend = TLegend(.65,.7,.85,.85)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.AddEntry(graph_eff0,"0Tag","lp")
  legend.AddEntry(graph_eff1,"1Tag","lp")
  legend.AddEntry(graph_eff2,"2Tag","lp")
  legend.Draw()

  c.SaveAs(outputFilename + '.' + outputFormat)


def eff_syst(inputDirs, final_state, title, outputFilename, outputFormat = 'png', pos = 0.80):

  files = [
    'RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root',
    'RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root'
  ]

  masses = [500., 700., 1200., 2000., 3500.]

  eff0_aa = []
  eff1_aa = []
  eff2_aa = []
    
  for inDir in inputDirs:
  
    eff0 = []
    eff1 = []
    eff2 = []

    for i, fl in enumerate(files):

      file = TFile(os.path.join(inDir,fl))

      histo_N0 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_0tag')
      histo_N1 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_1tag')
      histo_N2 = file.Get('myAnalyzer/h1_DijetMass_' + final_state + '_2tag')

      binMin = histo_N0.GetXaxis().FindBin(0.3*masses[i])
      binMax = histo_N0.GetXaxis().FindBin(1.3*masses[i])

      N0 = histo_N0.Integral(binMin,binMax)
      N1 = histo_N1.Integral(binMin,binMax)
      N2 = histo_N2.Integral(binMin,binMax)

      N = N0 + N1 + N2

      eff0.append( N0/N )
      eff1.append( N1/N )
      eff2.append( N2/N )

    eff0_aa.append(eff0)
    eff1_aa.append(eff1)
    eff2_aa.append(eff2)

  #print eff2_aa[0]
  #print eff2_aa[1]
  #print eff2_aa[2]
  #print eff2_aa[3]
  #print eff2_aa[4]
    
  eff0_syst_up = []
  eff0_syst_down = []
  eff1_syst_up = []
  eff1_syst_down = []
  eff2_syst_up = []
  eff2_syst_down = []
    
  for i in range(0,len(masses)):
    eff0_syst_up_SFb = max(eff0_aa[1][i],eff0_aa[2][i]) - eff0_aa[0][i]
    if( eff0_syst_up_SFb < 0 ): eff0_syst_up_SFb = 0.
    eff0_syst_down_SFb = min(eff0_aa[1][i],eff0_aa[2][i]) - eff0_aa[0][i]
    if( eff0_syst_down_SFb > 0 ): eff0_syst_down_SFb = 0.
    eff0_syst_up_SFl = max(eff0_aa[3][i],eff0_aa[4][i]) - eff0_aa[0][i]
    if( eff0_syst_up_SFl < 0 ): eff0_syst_up_SFl = 0.
    eff0_syst_down_SFl = min(eff0_aa[3][i],eff0_aa[4][i]) - eff0_aa[0][i]
    if( eff0_syst_down_SFl > 0 ): eff0_syst_down_SFl = 0.

    eff0_syst_up.append( sqrt( eff0_syst_up_SFb**2 + eff0_syst_up_SFl**2 ) )
    eff0_syst_down.append( sqrt( eff0_syst_down_SFb**2 + eff0_syst_down_SFl**2 ) )

    eff1_syst_up_SFb = max(eff1_aa[1][i],eff1_aa[2][i]) - eff1_aa[0][i]
    if( eff1_syst_up_SFb < 0 ): eff1_syst_up_SFb = 0.
    eff1_syst_down_SFb = min(eff1_aa[1][i],eff1_aa[2][i]) - eff1_aa[0][i]
    if( eff1_syst_down_SFb > 0 ): eff1_syst_down_SFb = 0.
    eff1_syst_up_SFl = max(eff1_aa[3][i],eff1_aa[4][i]) - eff1_aa[0][i]
    if( eff1_syst_up_SFl < 0 ): eff1_syst_up_SFl = 0.
    eff1_syst_down_SFl = min(eff1_aa[3][i],eff1_aa[4][i]) - eff1_aa[0][i]
    if( eff1_syst_down_SFl > 0 ): eff1_syst_down_SFl = 0.

    eff1_syst_up.append( sqrt( eff1_syst_up_SFb**2 + eff1_syst_up_SFl**2 ) )
    eff1_syst_down.append( sqrt( eff1_syst_down_SFb**2 + eff1_syst_down_SFl**2 ) )

    eff2_syst_up_SFb = max(eff2_aa[1][i],eff2_aa[2][i]) - eff2_aa[0][i]
    if( eff2_syst_up_SFb < 0 ): eff2_syst_up_SFb = 0.
    eff2_syst_down_SFb = min(eff2_aa[1][i],eff2_aa[2][i]) - eff2_aa[0][i]
    if( eff2_syst_down_SFb > 0 ): eff2_syst_down_SFb = 0.
    eff2_syst_up_SFl = max(eff2_aa[3][i],eff2_aa[4][i]) - eff2_aa[0][i]
    if( eff2_syst_up_SFl < 0 ): eff2_syst_up_SFl = 0.
    eff2_syst_down_SFl = min(eff2_aa[3][i],eff2_aa[4][i]) - eff2_aa[0][i]
    if( eff2_syst_down_SFl > 0 ): eff2_syst_down_SFl = 0.

    eff2_syst_up.append( sqrt( eff2_syst_up_SFb**2 + eff2_syst_up_SFl**2 ) )
    eff2_syst_down.append( sqrt( eff2_syst_down_SFb**2 + eff2_syst_down_SFl**2 ) )

  print "eff0_sys_up:"
  print eff0_syst_up
  print "eff0_sys_down:"
  print eff0_syst_down

  print "eff2_sys_up:"
  print eff2_syst_up
  print "eff2_sys_down:"
  print eff2_syst_down
    
  mass_array = array('d')
  mass_array_syst = array('d')
  eff0_array = array('d')
  eff0_array_syst = array('d')
  eff1_array = array('d')
  eff1_array_syst = array('d')
  eff2_array = array('d')
  eff2_array_syst = array('d')
    
  for i in range(0,len(masses)):
    mass_array.append(masses[i])
    mass_array_syst.append(masses[i])
    eff0_array.append( eff0_aa[0][i] )
    eff0_array_syst.append( eff0_aa[0][i] - eff0_syst_down[i] )
    eff1_array.append( eff1_aa[0][i] )
    eff1_array_syst.append( eff1_aa[0][i] - eff1_syst_down[i] )
    eff2_array.append( eff2_aa[0][i] )
    eff2_array_syst.append( eff2_aa[0][i] - eff2_syst_down[i] )

  for i in range(0,len(masses)):
    mass_array_syst.append(masses[len(masses)-i-1])
    eff0_array_syst.append( eff0_aa[0][len(masses)-i-1] + eff0_syst_up[len(masses)-i-1] )
    eff1_array_syst.append( eff1_aa[0][len(masses)-i-1] + eff1_syst_up[len(masses)-i-1] )
    eff2_array_syst.append( eff2_aa[0][len(masses)-i-1] + eff2_syst_up[len(masses)-i-1] )
    
  graph_eff0 = TGraph(len(mass_array),mass_array,eff0_array)
  graph_eff0.SetMarkerStyle(24)
  graph_eff0.SetMarkerColor(1)
  graph_eff0.SetLineWidth(2)
  graph_eff0.SetLineStyle(1)
  graph_eff0.SetLineColor(1)
  graph_eff0.SetFillColor(1)
  graph_eff0.SetFillStyle(3013)
  graph_eff0.SetTitle(title)
  graph_eff0.GetXaxis().SetTitle("Resonance Mass [GeV]")
  graph_eff0.GetYaxis().SetTitle("Efficiency")
  graph_eff0.GetYaxis().SetRangeUser(0.,1.)

  graph_eff0_syst = TGraph(len(mass_array_syst),mass_array_syst,eff0_array_syst)
  graph_eff0_syst.SetFillColor(1)
  graph_eff0_syst.SetFillStyle(3013)

  graph_eff1 = TGraph(len(mass_array),mass_array,eff1_array)
  graph_eff1.SetMarkerStyle(25)
  graph_eff1.SetMarkerColor(2)
  graph_eff1.SetLineWidth(2)
  graph_eff1.SetLineStyle(5)
  graph_eff1.SetLineColor(2)
  graph_eff1.SetFillColor(2)
  graph_eff1.SetFillStyle(3004)

  graph_eff1_syst = TGraph(len(mass_array_syst),mass_array_syst,eff1_array_syst)
  graph_eff1_syst.SetFillColor(2)
  graph_eff1_syst.SetFillStyle(3004)
  
  graph_eff2 = TGraph(len(mass_array),mass_array,eff2_array)
  graph_eff2.SetMarkerStyle(26)
  graph_eff2.SetMarkerColor(4)
  graph_eff2.SetLineWidth(2)
  graph_eff2.SetLineStyle(7)
  graph_eff2.SetLineColor(4)
  graph_eff2.SetFillColor(4)
  graph_eff2.SetFillStyle(3005)

  graph_eff2_syst = TGraph(len(mass_array_syst),mass_array_syst,eff2_array_syst)
  graph_eff2_syst.SetFillColor(4)
  graph_eff2_syst.SetFillStyle(3005)

  c = TCanvas("c", "",800,800)
  c.cd()

  graph_eff0.Draw("APL")
  graph_eff0_syst.Draw("F")
  graph_eff1.Draw("PL")
  graph_eff1_syst.Draw("F")
  graph_eff2.Draw("PL")
  graph_eff2_syst.Draw("F")

  graph_eff0_array.append( copy.deepcopy(graph_eff0) )
  graph_eff0_syst_array.append( copy.deepcopy(graph_eff0_syst) )
  graph_eff1_array.append( copy.deepcopy(graph_eff1) )
  graph_eff1_syst_array.append( copy.deepcopy(graph_eff1_syst) )
  graph_eff2_array.append( copy.deepcopy(graph_eff2) )
  graph_eff2_syst_array.append( copy.deepcopy(graph_eff2_syst) )

  legend = TLegend(.63,.72,.83,.87)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  legend.AddEntry(graph_eff0,"CSVL 0-tag","lfp")
  legend.AddEntry(graph_eff1,"CSVL 1-tag","lfp")
  legend.AddEntry(graph_eff2,"CSVL 2-tag","lfp")
  legend.Draw()

  l1 = TLatex()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetNDC()
  l1.SetTextSize(0.05)
  l1.DrawLatex(0.17,0.87, title)
  l1.SetTextSize(0.04)
  l1.DrawLatex(0.17,pos, "CMS Simulation")
  l1.DrawLatex(0.18,pos-0.05, "#sqrt{s} = 7 TeV")
  l1.DrawLatex(0.17,pos-0.10, "|#eta| < 2.5, |#Delta#eta| < 1.3")
  l1.DrawLatex(0.17,pos-0.15, "Wide Jets")
  
  c.SaveAs(outputFilename + '.' + outputFormat)  

  
if __name__ == "__main__":
 
  eff('CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets', 'bbbar', 'G#rightarrowb#bar{b}, CSVL', 'RSGToBBbar_bTagEfficiency_CSVL_PUSFReweighted')
  #eff('CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets', 'ccbar', 'G#rightarrowc#bar{c}, CSVL', 'RSGToCCbar_bTagEfficiency_CSVL_PUSFReweighted')
  #eff('CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets', 'qqbarlight', 'G#rightarrowq#bar{q} (q=u,d,s), CSVL', 'RSGToQQbarLight_bTagEfficiency_CSVL_PUSFReweighted')
  #eff('CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets', 'gg', 'G#rightarrowgg, CSVL', 'RSGToGG_bTagEfficiency_CSVL_PUSFReweighted')

  #eff('DijetBBTag_RSGraviton/CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted', 'bbbar', 'G#rightarrowb#bar{b}, CSVM', 'RSGToBBbar_bTagEfficiency_CSVM_PUSFReweighted')
  #eff('DijetBBTag_RSGraviton/CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted', 'ccbar', 'G#rightarrowc#bar{c}, CSVM', 'RSGToCCbar_bTagEfficiency_CSVM_PUSFReweighted')
  #eff('DijetBBTag_RSGraviton/CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted', 'qqbarlight', 'G#rightarrowq#bar{q} (q=u,d,s), CSVM', 'RSGToQQbarLight_bTagEfficiency_CSVM_PUSFReweighted')
  #eff('DijetBBTag_RSGraviton/CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted', 'gg', 'G#rightarrowgg, CSVM', 'RSGToGG_bTagEfficiency_CSVM_PUSFReweighted')
  
  # efficiency plots with systematic uncertainty bands
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbUp',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlUp'],
           #'bbbar', 'G#rightarrowb#bar{b}', 'RSGToBBbar_bTagEfficiency_CSVL_PUSFReweighted_syst','png')
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbUp',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlUp'],
           #'ccbar', 'G#rightarrowc#bar{c}', 'RSGToCCbar_bTagEfficiency_CSVL_PUSFReweighted_syst','png')
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbUp',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlUp'],
           #'qqbarlight', 'G#rightarrowq#bar{q} (q=u,d,s)', 'RSGToQQbarLight_bTagEfficiency_CSVL_PUSFReweighted_syst','png',0.6)
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFbUp',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlDown',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets_SFlUp'],
           #'gg', 'G#rightarrowgg', 'RSGToGG_bTagEfficiency_CSVL_PUSFReweighted_syst','png',0.6)

  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbP1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlP1'],
           #'bbbar', 'G#rightarrowb#bar{b}, CSVM', 'RSGToBBbar_bTagEfficiency_CSVM_PUSFReweighted_syst')
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbP1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlP1'],
           #'ccbar', 'G#rightarrowc#bar{c}, CSVM', 'RSGToCCbar_bTagEfficiency_CSVM_PUSFReweighted_syst')
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbP1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlP1'],
           #'qqbarlight', 'G#rightarrowq#bar{q} (q=u,d,s), CSVM', 'RSGToQQbarLight_bTagEfficiency_CSVM_PUSFReweighted_syst')
  #eff_syst(['CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFbP1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlM1',
            #'CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVM_PUSFReweighted_SFlP1'],
           #'gg', 'G#rightarrowgg, CSVM', 'RSGToGG_bTagEfficiency_CSVM_PUSFReweighted_syst')

  #c = TCanvas("c", "",800,800)
  #c.cd()

  ##graph_eff0_array[2].SetTitle('RSG, CSVL')
  #graph_eff0_array[2].Draw("APL")
  #graph_eff0_syst_array[2].Draw("F")
  #graph_eff1_array[2].Draw("PL")
  #graph_eff1_syst_array[2].Draw("F")
  #graph_eff2_array[2].Draw("PL")
  #graph_eff2_syst_array[2].Draw("F")

  #graph_eff0_array[3].SetLineColor(6)
  #graph_eff0_array[3].SetFillColor(6)
  #graph_eff0_array[3].SetMarkerColor(6)
  #graph_eff0_syst_array[3].SetFillColor(6)
  #graph_eff1_array[3].SetLineColor(kGreen+2)
  #graph_eff1_array[3].SetFillColor(kGreen+2)
  #graph_eff1_array[3].SetMarkerColor(kGreen+2)
  #graph_eff1_syst_array[3].SetFillColor(kGreen+2)
  #graph_eff1_syst_array[3].SetFillStyle(3005)
  #graph_eff2_array[3].SetLineColor(kOrange+7)
  #graph_eff2_array[3].SetFillColor(kOrange+7)
  #graph_eff2_array[3].SetMarkerColor(kOrange+7)
  #graph_eff2_syst_array[3].SetFillColor(kOrange+7)
  #graph_eff2_syst_array[3].SetFillStyle(3004)
  
  #graph_eff0_array[3].Draw("PL")
  #graph_eff0_syst_array[3].Draw("F")
  #graph_eff1_array[3].Draw("PL")
  #graph_eff1_syst_array[3].Draw("F")
  #graph_eff2_array[3].Draw("PL")
  #graph_eff2_syst_array[3].Draw("F")

  #legend1 = TLegend(.18,.41,.48,.59)
  #legend1.SetBorderSize(0)
  #legend1.SetFillColor(0)
  #legend1.SetFillStyle(0)
  #legend1.SetTextFont(42)
  #legend1.SetTextSize(0.03)
  #legend1.SetHeader("G#rightarrowq#bar{q} (q=u,d,s)")
  #legend1.AddEntry(graph_eff0_array[2],"CSVL 0-tag","lfp")
  #legend1.AddEntry(graph_eff1_array[2],"CSVL 1-tag","lfp")
  #legend1.AddEntry(graph_eff2_array[2],"CSVL 2-tag","lfp")
  #legend1.Draw()

  #legend2 = TLegend(.5,.41,.8,.59)
  #legend2.SetBorderSize(0)
  #legend2.SetFillColor(0)
  #legend2.SetFillStyle(0)
  #legend2.SetTextFont(42)
  #legend2.SetTextSize(0.03)
  #legend2.SetHeader("G#rightarrowgg")
  #legend2.AddEntry(graph_eff0_array[3],"CSVL 0-tag","lfp")
  #legend2.AddEntry(graph_eff1_array[3],"CSVL 1-tag","lfp")
  #legend2.AddEntry(graph_eff2_array[3],"CSVL 2-tag","lfp")
  #legend2.Draw()

  #l1 = TLatex()
  #l1.SetTextAlign(12)
  #l1.SetTextFont(42)
  #l1.SetNDC()
  #l1.SetTextSize(0.04)
  #l1.DrawLatex(0.60,0.88, "CMS Simulation")
  #l1.DrawLatex(0.61,0.83, "#sqrt{s} = 7 TeV")
  #l1.DrawLatex(0.60,0.78, "|#eta| < 2.5, |#Delta#eta| < 1.3")
  #l1.DrawLatex(0.60,0.73, "Wide Jets")
  
  #c.SaveAs('RSGToGGandQQbarLight_bTagEfficiency_CSVL_PUSFReweighted_syst.eps')

  ##graph_eff0_array[6].SetTitle('RSG, CSVM')
  #graph_eff0_array[6].Draw("APL")
  #graph_eff0_syst_array[6].Draw("F")
  #graph_eff1_array[6].Draw("PL")
  #graph_eff1_syst_array[6].Draw("F")
  #graph_eff2_array[6].Draw("PL")
  #graph_eff2_syst_array[6].Draw("F")

  #graph_eff0_array[7].SetLineColor(7)
  #graph_eff0_array[7].SetFillColor(7)
  #graph_eff0_array[7].SetMarkerColor(7)
  #graph_eff0_syst_array[7].SetFillColor(7)
  #graph_eff1_array[7].SetLineColor(8)
  #graph_eff1_array[7].SetFillColor(8)
  #graph_eff1_array[7].SetMarkerColor(8)
  #graph_eff1_syst_array[7].SetFillColor(8)
  #graph_eff1_syst_array[7].SetFillStyle(3005)
  #graph_eff2_array[7].SetLineColor(kOrange-3)
  #graph_eff2_array[7].SetFillColor(kOrange-3)
  #graph_eff2_array[7].SetMarkerColor(kOrange-3)
  #graph_eff2_syst_array[7].SetFillColor(kOrange-3)
  #graph_eff2_syst_array[7].SetFillStyle(3004)

  #graph_eff0_array[7].Draw("PL")
  #graph_eff0_syst_array[7].Draw("F")
  #graph_eff1_array[7].Draw("PL")
  #graph_eff1_syst_array[7].Draw("F")
  #graph_eff2_array[7].Draw("PL")
  #graph_eff2_syst_array[7].Draw("F")

  #legend3 = TLegend(.15,.38,.45,.55)
  #legend3.SetBorderSize(0)
  #legend3.SetFillColor(0)
  #legend3.SetFillStyle(0)
  #legend3.SetTextSize(0.03)
  #legend3.AddEntry(graph_eff0_array[6],"0Tag, q#bar{q} (q=u,d,s)","lfp")
  #legend3.AddEntry(graph_eff1_array[6],"1Tag, q#bar{q} (q=u,d,s)","lfp")
  #legend3.AddEntry(graph_eff2_array[6],"2Tag, q#bar{q} (q=u,d,s)","lfp")
  #legend3.Draw()

  #legend4 = TLegend(.5,.38,.8,.55)
  #legend4.SetBorderSize(0)
  #legend4.SetFillColor(0)
  #legend4.SetFillStyle(0)
  #legend4.SetTextSize(0.03)
  #legend4.AddEntry(graph_eff0_array[7],"0Tag, gg","lfp")
  #legend4.AddEntry(graph_eff1_array[7],"1Tag, gg","lfp")
  #legend4.AddEntry(graph_eff2_array[7],"2Tag, gg","lfp")
  #legend4.Draw()

  #c.SaveAs('RSGToGGandQQbarLight_bTagEfficiency_CSVM_PUSFReweighted_syst.png')
