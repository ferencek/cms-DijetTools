#!/usr/bin/env python

import os, sys
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
#gStyle.SetPadLeftMargin(0.13);
#gStyle.SetPadRightMargin(0.07);


def main(inputDir, tagger, BR, outputFilename):

  if BR > 1.:
    print 'bbbar branching ratio has to be less than 1'
    sys.exit(1)

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

    histo_N0_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_0tag')
    histo_N1_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_1tag')
    histo_N2_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_2tag')

    histo_N0_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_0tag')
    histo_N1_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_1tag')
    histo_N2_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_2tag')

    binMin = histo_N0_bbbar.GetXaxis().FindBin(0.3*masses[i]);
    binMax = histo_N0_bbbar.GetXaxis().FindBin(1.3*masses[i]);

    N0_bbbar = histo_N0_bbbar.Integral(binMin,binMax)
    N1_bbbar = histo_N1_bbbar.Integral(binMin,binMax)
    N2_bbbar = histo_N2_bbbar.Integral(binMin,binMax)

    N0_nonbbbar = histo_N0_nonbbbar.Integral(binMin,binMax)
    N1_nonbbbar = histo_N1_nonbbbar.Integral(binMin,binMax)
    N2_nonbbbar = histo_N2_nonbbbar.Integral(binMin,binMax)

    N_bbbar = N0_bbbar + N1_bbbar + N2_bbbar
    N_nonbbbar = N0_nonbbbar + N1_nonbbbar + N2_nonbbbar

    eff0_bbbar = N0_bbbar/N_bbbar
    eff1_bbbar = N1_bbbar/N_bbbar
    eff2_bbbar = N2_bbbar/N_bbbar
    
    eff0_nonbbbar = N0_nonbbbar/N_nonbbbar
    eff1_nonbbbar = N1_nonbbbar/N_nonbbbar
    eff2_nonbbbar = N2_nonbbbar/N_nonbbbar
   
    r_nonbbbar = (N_bbbar*(1./BR - 1.))/(N_nonbbbar)
    
    mass_array.append(masses[i])
    mass_err_array.append(0.)

    eff0_array.append( (N0_bbbar + r_nonbbbar*N0_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )
    eff0_err_array.append( sqrt(((1./(r_nonbbbar*N_nonbbbar + N_bbbar))**2)*(eff0_bbbar*(1.-eff0_bbbar)*N_bbbar + (r_nonbbbar**2)*eff0_nonbbbar*(1.-eff0_nonbbbar)*N_nonbbbar)) )
    eff1_array.append( (N1_bbbar + r_nonbbbar*N1_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )
    eff1_err_array.append( sqrt(((1./(r_nonbbbar*N_nonbbbar + N_bbbar))**2)*(eff1_bbbar*(1.-eff1_bbbar)*N_bbbar + (r_nonbbbar**2)*eff1_nonbbbar*(1.-eff1_nonbbbar)*N_nonbbbar)) )
    eff2_array.append( (N2_bbbar + r_nonbbbar*N2_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )
    eff2_err_array.append( sqrt(((1./(r_nonbbbar*N_nonbbbar + N_bbbar))**2)*(eff2_bbbar*(1.-eff2_bbbar)*N_bbbar + (r_nonbbbar**2)*eff2_nonbbbar*(1.-eff2_nonbbbar)*N_nonbbbar)) )
    
  graph_eff0 = TGraphErrors(len(mass_array),mass_array,eff0_array,mass_err_array,eff0_err_array)
  graph_eff0.SetMarkerStyle(24)
  graph_eff0.SetLineWidth(2)
  graph_eff0.SetLineStyle(1)
  graph_eff0.SetLineColor(1)
  graph_eff0.SetTitle(tagger + ', BR(RSG#rightarrowb#bar{b})=' + str(BR))
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

  c = TCanvas("c", "",1000,1000)
  c.cd()

  graph_eff0.Draw("APC")
  graph_eff1.Draw("PC")
  graph_eff2.Draw("PC")

  legend = TLegend(.65,.7,.85,.85)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.AddEntry(graph_eff0,"0Tag","lp")
  legend.AddEntry(graph_eff1,"1Tag","lp")
  legend.AddEntry(graph_eff2,"2Tag","lp")
  legend.Draw()

  c.SaveAs(outputFilename.replace('.png','_BR' + str(BR) + '.png'))


def main_syst(inputDirs, tagger, BR, outputFilename):

  if BR > 1.:
    print 'bbbar branching ratio has to be less than 1'
    sys.exit(1)

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

      histo_N0_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_0tag')
      histo_N1_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_1tag')
      histo_N2_bbbar = file.Get('myAnalyzer/h1_DijetMass_bbbar_2tag')

      histo_N0_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_0tag')
      histo_N1_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_1tag')
      histo_N2_nonbbbar = file.Get('myAnalyzer/h1_DijetMass_nonbbbar_2tag')

      binMin = histo_N0_bbbar.GetXaxis().FindBin(0.3*masses[i]);
      binMax = histo_N0_bbbar.GetXaxis().FindBin(1.3*masses[i]);

      N0_bbbar = histo_N0_bbbar.Integral(binMin,binMax)
      N1_bbbar = histo_N1_bbbar.Integral(binMin,binMax)
      N2_bbbar = histo_N2_bbbar.Integral(binMin,binMax)

      N0_nonbbbar = histo_N0_nonbbbar.Integral(binMin,binMax)
      N1_nonbbbar = histo_N1_nonbbbar.Integral(binMin,binMax)
      N2_nonbbbar = histo_N2_nonbbbar.Integral(binMin,binMax)

      N_bbbar = N0_bbbar + N1_bbbar + N2_bbbar
      N_nonbbbar = N0_nonbbbar + N1_nonbbbar + N2_nonbbbar

      r_nonbbbar = (N_bbbar*(1./BR - 1.))/(N_nonbbbar)

      eff0.append( (N0_bbbar + r_nonbbbar*N0_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )
      eff1.append( (N1_bbbar + r_nonbbbar*N1_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )
      eff2.append( (N2_bbbar + r_nonbbbar*N2_nonbbbar)/(r_nonbbbar*N_nonbbbar + N_bbbar) )

    eff0_aa.append(eff0)
    eff1_aa.append(eff1)
    eff2_aa.append(eff2)

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
    eff0_array_syst.append( min(eff0_aa[1][i],eff0_aa[2][i]) )
    eff1_array.append( eff1_aa[0][i] )
    eff1_array_syst.append( min(eff1_aa[1][i],eff1_aa[2][i]) )
    eff2_array.append( eff2_aa[0][i] )
    eff2_array_syst.append( min(eff2_aa[1][i],eff2_aa[2][i]) )

  for i in range(0,len(masses)):
    mass_array_syst.append(masses[len(masses)-i-1])
    eff0_array_syst.append( max(eff0_aa[1][len(masses)-i-1],eff0_aa[2][len(masses)-i-1]) )
    eff1_array_syst.append( max(eff1_aa[1][len(masses)-i-1],eff1_aa[2][len(masses)-i-1]) )
    eff2_array_syst.append( max(eff2_aa[1][len(masses)-i-1],eff2_aa[2][len(masses)-i-1]) )
    
  graph_eff0 = TGraph(len(mass_array),mass_array,eff0_array)
  graph_eff0.SetMarkerStyle(24)
  graph_eff0.SetMarkerColor(1)
  graph_eff0.SetLineWidth(2)
  graph_eff0.SetLineStyle(1)
  graph_eff0.SetLineColor(1)
  graph_eff0.SetFillColor(1)
  graph_eff0.SetFillStyle(3013)
  graph_eff0.SetTitle(tagger + ', BR(RSG#rightarrowb#bar{b})=' + str(BR))
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

  c = TCanvas("c", "",1000,1000)
  c.cd()

  graph_eff0.Draw("APL")
  graph_eff0_syst.Draw("F")
  graph_eff1.Draw("PL")
  graph_eff1_syst.Draw("F")
  graph_eff2.Draw("PL")
  graph_eff2_syst.Draw("F")

  legend = TLegend(.65,.7,.85,.85)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.AddEntry(graph_eff0,"0Tag","lfp")
  legend.AddEntry(graph_eff1,"1Tag","lfp")
  legend.AddEntry(graph_eff2,"2Tag","lfp")
  legend.Draw()

  c.SaveAs(outputFilename.replace('.png','_BR' + str(BR) + '_syst.png'))  

  
if __name__ == "__main__":
 
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted', 'TCHEL', 1., 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted', 'TCHEL', 0.75, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted', 'TCHEL', 0.5, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted', 'TCHEL', 0.1, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')

  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp', 'TCHEL', 1., 'RSG_bTagEfficiency_PUSFReweighted_SFUp_TCHEL.png')
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp', 'TCHEL', 0.5, 'RSG_bTagEfficiency_PUSFReweighted_SFUp_TCHEL.png')

  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown', 'TCHEL', 1., 'RSG_bTagEfficiency_PUSFReweighted_SFDown_TCHEL.png')
  #main('CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown', 'TCHEL', 0.5, 'RSG_bTagEfficiency_PUSFReweighted_SFDown_TCHEL.png')

  
  # efficiency plots with systematic uncertainty bands
  main_syst(['CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown'],
             'TCHEL', 1., 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')

  main_syst(['CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown'],
             'TCHEL', 0.75, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')

  main_syst(['CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown'],
             'TCHEL', 0.5, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')

  main_syst(['CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown'],
             'TCHEL', 0.1, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')

  main_syst(['CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFUp',
             'CRAB_Jobs_RSGraviton_bTagEfficiency_TCHEL_PUSFReweighted_SFDown'],
             'TCHEL', 0.01, 'RSG_bTagEfficiency_PUSFReweighted_TCHEL.png')
            
      