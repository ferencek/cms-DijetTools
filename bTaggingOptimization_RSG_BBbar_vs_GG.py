#!/usr/bin/env python

import os, sys
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetPaintTextFormat("1.2g")
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)

def plot_JP_CSV(xMax,yMax,fileName,linesAlgo=""):

  dijetMassBins = [
    "DijetMass3500to4000GeV",
    "DijetMass3000to3500GeV",
    "DijetMass2500to3000GeV",
    "DijetMass2000to2500GeV",
    "DijetMass1500to2000GeV",
    "DijetMass1000to1500GeV",
    "DijetMass500to1000GeV"
  ]

  file1 = TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoBBbar_240212/Final__histograms.root")
  file2 = TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoJJ_010312/Final__histograms.root")

  algoBins = {"JPL": 10, "JPM": 11, "JPT": 12, "CSVL": 13, "CSVM": 14, "CSVT": 15}

  algoArrays = {}

  for algo in algoBins.keys():
    algoArrays[algo] = [array('d'), array('d')]

  for dmb in dijetMassBins:

    h2_RSGravitonToBBbar_eff_denom = file1.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_denom")
    h2_RSGravitonToBBbar_eff_num = file1.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_num")

    h2_RSGravitonToGG_eff_denom = file2.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_denom")
    h2_RSGravitonToGG_eff_num = file2.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_num")

    h2_RSGravitonToBBbar_eff = TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 15, 0.5, 15.5, 15, 0.5, 15.5)
    h2_RSGravitonToBBbar_eff.Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom)

    h2_RSGravitonToGG_eff = TH2D("h2_RSGravitonToGG_eff", "h2_RSGravitonToGG_eff", 15, 0.5, 15.5, 15, 0.5, 15.5)
    h2_RSGravitonToGG_eff.Divide(h2_RSGravitonToGG_eff_num,h2_RSGravitonToGG_eff_denom)

    for algo in algoBins.keys():

      b = algoBins[algo]
      algoArrays[algo][0].append(h2_RSGravitonToBBbar_eff.GetBinContent(b,b))
      algoArrays[algo][1].append(h2_RSGravitonToGG_eff.GetBinContent(b,b))

    del h2_RSGravitonToBBbar_eff
    del h2_RSGravitonToGG_eff

    
  c = TCanvas("c","",1000,1000)
  c.cd()

  bkg = TH2D("","",100,0,xMax,100,0,yMax)
  bkg.GetXaxis().SetTitle("#epsilon_{RSG#rightarrowb#bar{b}}")
  bkg.GetYaxis().SetTitle("#epsilon_{RSG#rightarrowgg}")
  bkg.SetTitleOffset(1.5,"Y")
  bkg.Draw()

  markerColor = {"JPL": [24, 28],  "CSVL": [20, 4], "JPM": [25, 8], "CSVM": [21, 2], "JPT": [26, 6], "CSVT": [22, 1]}

  legend = TLegend(.6,.6,.85,.85)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)

  algoGraphs = {}
  algoFuncs = {}
  
  for algo in ["JPL", "CSVL", "JPM", "CSVM", "JPT", "CSVT"]:
    
    algoGraphs[algo] = TGraph(7,algoArrays[algo][0],algoArrays[algo][1])
    algoGraphs[algo].SetMarkerStyle(markerColor[algo][0])
    algoGraphs[algo].SetMarkerColor(markerColor[algo][1])
    algoGraphs[algo].SetLineWidth(2)
    algoGraphs[algo].SetLineColor(markerColor[algo][1])
    algoGraphs[algo].Draw("PLsame")
    legend.AddEntry(algoGraphs[algo],algo,"lp")

    if algo==linesAlgo:
      for i in range(0,len(dijetMassBins)):
        algoFuncs[algo+str(i)] = TF1(algo+str(i),"("+str(algoArrays[algo][1][i])+"/"+str(algoArrays[algo][0][i])+")*x",0,1)
        algoFuncs[algo+str(i)].SetLineColor(markerColor[algo][1])
        algoFuncs[algo+str(i)].SetLineStyle(3)
        algoFuncs[algo+str(i)].SetLineWidth(1)
        algoFuncs[algo+str(i)].Draw("same")
      
  legend.Draw()
  c.SaveAs(fileName)


def plot_TCHEL_CSVL(xMax,yMax,fileName,linesAlgo=""):

  dijetMassBins = [
    "DijetMass3500to4000GeV",
    "DijetMass3000to3500GeV",
    "DijetMass2500to3000GeV",
    "DijetMass2000to2500GeV",
    "DijetMass1500to2000GeV",
    "DijetMass1000to1500GeV",
    "DijetMass500to1000GeV"
  ]

  file1 = TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoBBbar_240212/Final__histograms.root")
  file2 = TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoJJ_010312/Final__histograms.root")

  algoBins = {"TCHEL": 1, "CSVL": 13}

  algoArrays = {}

  for algo in algoBins.keys():
    algoArrays[algo] = [array('d'), array('d')]

  for dmb in dijetMassBins:

    h2_RSGravitonToBBbar_eff_denom = file1.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_denom")
    h2_RSGravitonToBBbar_eff_num = file1.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_num")

    h2_RSGravitonToGG_eff_denom = file2.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_denom")
    h2_RSGravitonToGG_eff_num = file2.Get("RSGravitonJJ__h2_BTaggers_" + dmb + "_eff_num")

    h2_RSGravitonToBBbar_eff = TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 15, 0.5, 15.5, 15, 0.5, 15.5)
    h2_RSGravitonToBBbar_eff.Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom)

    h2_RSGravitonToGG_eff = TH2D("h2_RSGravitonToGG_eff", "h2_RSGravitonToGG_eff", 15, 0.5, 15.5, 15, 0.5, 15.5)
    h2_RSGravitonToGG_eff.Divide(h2_RSGravitonToGG_eff_num,h2_RSGravitonToGG_eff_denom)

    for algo in algoBins.keys():

      b = algoBins[algo]
      algoArrays[algo][0].append(h2_RSGravitonToBBbar_eff.GetBinContent(b,b))
      algoArrays[algo][1].append(h2_RSGravitonToGG_eff.GetBinContent(b,b))

    del h2_RSGravitonToBBbar_eff
    del h2_RSGravitonToGG_eff


  c = TCanvas("c","",1000,1000)
  c.cd()

  bkg = TH2D("","",100,0,xMax,100,0,yMax)
  bkg.GetXaxis().SetTitle("#epsilon_{G#rightarrowb#bar{b}}")
  bkg.GetYaxis().SetTitle("#epsilon_{G#rightarrowgg}")
  bkg.SetTitleOffset(1.1,"Y")
  bkg.Draw()

  markerColor = {"TCHEL": [24, kRed],  "CSVL": [26, kGreen+2]}

  legend = TLegend(.2,.4,.45,.55)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)

  l1 = TLatex()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetNDC()
  l1.SetTextSize(0.04)
  l1.DrawLatex(0.19,0.90, "CMS Simulation")
  l1.DrawLatex(0.20,0.85, "#sqrt{s} = 7 TeV")
  l1.DrawLatex(0.19,0.80, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.19,0.75, "Anti-k_{T} R = 0.7 PF Jets");
  l1.SetTextSize(0.03)
  l1.SetTextColor(kGreen+2)
  l1.DrawLatex(0.7,0.2, "0.5<m_{jj}<1 TeV")
  l1.DrawLatex(0.19,0.27, "3.5<m_{jj}<4 TeV")
  l1.SetTextColor(kRed)
  l1.DrawLatex(0.76,0.42, "0.5<m_{jj}<1 TeV")
  l1.DrawLatex(0.61,0.87, "3.5<m_{jj}<4 TeV")

  algoGraphs = {}
  algoFuncs = {}

  for algo in ["TCHEL", "CSVL"]:

    algoGraphs[algo] = TGraph(7,algoArrays[algo][0],algoArrays[algo][1])
    algoGraphs[algo].SetMarkerStyle(markerColor[algo][0])
    algoGraphs[algo].SetMarkerColor(markerColor[algo][1])
    algoGraphs[algo].SetLineWidth(2)
    algoGraphs[algo].SetLineColor(markerColor[algo][1])
    algoGraphs[algo].Draw("PLsame")
    legend.AddEntry(algoGraphs[algo],algo + ", 2 b-tags","lp")

    if algo==linesAlgo:
      for i in range(0,len(dijetMassBins)):
        algoFuncs[algo+str(i)] = TF1(algo+str(i),"("+str(algoArrays[algo][1][i])+"/"+str(algoArrays[algo][0][i])+")*x",0,1)
        algoFuncs[algo+str(i)].SetLineColor(markerColor[algo][1])
        algoFuncs[algo+str(i)].SetLineStyle(3)
        algoFuncs[algo+str(i)].SetLineWidth(1)
        algoFuncs[algo+str(i)].Draw("same")

  legend.Draw()
  c.SaveAs(fileName)

  
if __name__ == "__main__":

  #plot_JP_CSV(0.6,0.23,"DoubleTag_eff_RSG_GG_vs_BBbar.png","CSVM")
  #plot_JP_CSV(0.36,0.045,"DoubleTag_eff_RSG_GG_vs_BBbar_zoomed.png")

  plot_TCHEL_CSVL(0.75,0.85,"DoubleTag_eff_TCHEL_CSVL_RSG_GG_vs_BBbar.eps")
