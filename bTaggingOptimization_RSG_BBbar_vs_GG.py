#!/usr/bin/env python

import os, sys
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0)
gStyle.SetLabelSize(0.035, "XYZ");
gStyle.SetPaintTextFormat("1.2g");
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
gStyle.SetPadLeftMargin(0.13);
gStyle.SetPadRightMargin(0.07);


def main(xMax,yMax,fileName,linesAlgo=""):

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
  bkg.SetTitleOffset(1.5,"Y");
  bkg.Draw()

  markerColor = {"JPL": [24, 28],  "CSVL": [20, 4], "JPM": [25, 8], "CSVM": [21, 2], "JPT": [26, 6], "CSVT": [22, 1]}

  legend = TLegend(.6,.6,.85,.85)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

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

  
if __name__ == "__main__":

  main(0.6,0.23,"DoubleTag_eff_RSG_BBbar_vs_GG.png","CSVM")
  main(0.36,0.045,"DoubleTag_eff_RSG_BBbar_vs_GG_zoomed.png")
  