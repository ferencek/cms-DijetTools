#!/usr/bin/env python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
from math import *


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

# stat only
#xs_limits_0p2 = array('d', [1.26494, 0.51744500000000004, 0.24494099999999999, 0.15290899999999999, 0.11716699999999999, 0.127417, 0.0980293, 0.074709399999999995, 0.052263499999999997, 0.047051599999999999, 0.077230499999999994, 0.074166700000000002, 0.047022000000000001, 0.033457899999999999, 0.034281499999999999, 0.038926099999999998, 0.030869500000000001, 0.0180441, 0.0115011, 0.0081844799999999992, 0.00617397, 0.0043310700000000002, 0.0035545199999999998, 0.0033220699999999999, 0.0034729100000000001, 0.0033162199999999999, 0.0031368400000000001, 0.0028479400000000002, 0.0025428299999999998, 0.0021667800000000001, 0.0018278700000000001])
#xs_limits_1p0 = array('d', [0.56668799999999997, 0.58086000000000004, 0.17114699999999999, 0.089566499999999993, 0.059730699999999998, 0.056915300000000002, 0.064923800000000004, 0.091563800000000001, 0.066524299999999995, 0.051226899999999999, 0.056009000000000003, 0.053176099999999997, 0.041119599999999999, 0.033465799999999997, 0.034853299999999997, 0.038552099999999999, 0.032492, 0.0215675, 0.0135388, 0.010174300000000001, 0.0076462500000000003, 0.0060970299999999998, 0.0046936, 0.0043732399999999996, 0.0043788999999999998, 0.0040874099999999997, 0.0036084099999999998, 0.0032902000000000001, 0.0031122300000000001, 0.0026434700000000002, 0.0022806300000000001])

# syst+stat
xs_limits_0p2 = array('d', [1.3648100000000001, 0.58711899999999995, 0.27155200000000002, 0.183139, 0.13209099999999999, 0.12962000000000001, 0.108418, 0.0842802, 0.060856199999999999, 0.0575132, 0.077734600000000001, 0.078043399999999999, 0.060779100000000003, 0.040989100000000001, 0.0360913, 0.038837299999999998, 0.035022699999999997, 0.026342399999999998, 0.0167986, 0.0113588, 0.0075952099999999998, 0.0055416399999999996, 0.0040954900000000002, 0.0035792699999999998, 0.00344086, 0.0033120699999999999, 0.0031159500000000001, 0.0029412100000000001, 0.0026468500000000001, 0.0023606199999999999, 0.00205522])
xs_limits_1p0 = array('d', [0.61651699999999998, 0.64604499999999998, 0.21223, 0.106056, 0.0684279, 0.063028699999999993, 0.069857000000000002, 0.093030399999999999, 0.073984900000000006, 0.055434499999999998, 0.058917700000000003, 0.056365800000000001, 0.046353699999999998, 0.037593399999999999, 0.036586, 0.038987599999999997, 0.036139200000000003, 0.029617399999999999, 0.018850599999999999, 0.013448099999999999, 0.0097049700000000003, 0.0070975300000000003, 0.0057316299999999997, 0.0046952699999999997, 0.0043628499999999997, 0.0040803899999999997, 0.0037947599999999999, 0.003447, 0.0032466299999999999, 0.0028686499999999999, 0.0025190299999999998])

masses = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0])

m_zprime = array('d', [700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800.])
m_s8 = array('d', [700., 1000., 1200., 1400., 1700., 2100., 2700.])
zprime = array('d', [0.6246E+01, 0.3427E+01, 0.1969E+01, 0.1172E+01, 0.7171E+00, 0.4486E+00, 0.2857E+00, 0.1845E+00, 0.1206E+00, 0.7961E-01, 0.5295E-01, 0.3545E-01, 0.2386E-01, 0.1611E-01, 0.1092E-01, 0.7413E-02, 0.5039E-02, 0.3426E-02, 0.2329E-02, 0.1580E-02, 0.1070E-02, 0.7231E-03, 0.4867E-03, 0.3261E-03, 0.2174E-03, 0.1440E-03, 0.9477E-04, 0.6190E-04, 0.4007E-04])
s8 = array('d', [6.52288, 0.995265, 0.319478, 0.108527, 0.0235598, 0.00331955, 0.000188328])

for i in range(0,len(masses)):
 xs_limits_0p2[i] = log(xs_limits_0p2[i])
 xs_limits_1p0[i] = log(xs_limits_1p0[i])

for i in range(0,len(m_zprime)):
 zprime[i] = log(zprime[i])

for i in range(0,len(m_s8)):
 s8[i] = log(s8[i])

graph_limit_0p2 = TGraph(len(masses),masses,xs_limits_0p2)
graph_limit_1p0 = TGraph(len(masses),masses,xs_limits_1p0)

graph_zprime = TGraph(len(m_zprime),m_zprime,zprime)
graph_s8 = TGraph(len(m_s8),m_s8,s8)


m_final = array('d')
zprime_ratio = array('d')
s8_ratio = array('d')

for i in range(0,18):

 mass = 1000. + 100.*i

 m_final.append(mass)

 zprime_ratio.append( exp(graph_zprime.Eval(mass)) / exp(graph_limit_0p2.Eval(mass)) )
 s8_ratio.append( exp(graph_s8.Eval(mass)) / exp(graph_limit_1p0.Eval(mass)) )

#print m_final
#print zprime_ratio
#print s8_ratio
 
graph_zprime_ratio = TGraph(len(m_final),m_final,zprime_ratio)
graph_zprime_ratio.GetXaxis().SetTitle("Resonance Mass [GeV]")
graph_zprime_ratio.GetYaxis().SetTitle("(Theory #sigma) / (95% CL Limit on #sigma)")
graph_zprime_ratio.GetYaxis().SetRangeUser(0,2)
graph_zprime_ratio.GetXaxis().SetNdivisions(505)
graph_zprime_ratio.SetLineWidth(2)
graph_zprime_ratio.SetLineStyle(8)
graph_zprime_ratio.SetLineColor(kTeal+3)

graph_s8_ratio = TGraph(len(m_final),m_final,s8_ratio)
graph_s8_ratio.SetLineWidth(2)
graph_s8_ratio.SetLineStyle(6)
graph_s8_ratio.SetLineColor(kBlue)


c = TCanvas("c", "",800,800)
c.cd()

graph_zprime_ratio.Draw("AL")
graph_s8_ratio.Draw("L")

line = TLine(1000,1,2700,1);
line.SetLineWidth(2)
line.SetLineColor(kBlack)
line.SetLineStyle(1)
line.Draw("same")

legend = TLegend(.50,.63,.80,.76)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
#legend.SetHeader('95% CL Upper Limits (stat. only)')
legend.AddEntry(graph_zprime_ratio,"Z' xs / Obs. upper limit","lp")
legend.AddEntry(graph_s8_ratio,"S8 xs / Obs. upper limit","lp")
legend.Draw()

#l1 = TLatex()
#l1.SetTextAlign(12)
#l1.SetTextFont(42)
#l1.SetNDC()
#l1.SetTextSize(0.04)
#l1.DrawLatex(0.18,0.89, "qq/bb, M=1 TeV, f_{b#bar{b}}=0.5")
#l1.SetTextSize(0.04)
#l1.DrawLatex(0.18,0.43, "CMS Preliminary")
#l1.DrawLatex(0.18,0.35, "#intLdt = 5 fb^{-1}")
#l1.DrawLatex(0.19,0.30, "#sqrt{s} = 7 TeV")
#l1.DrawLatex(0.18,0.25, "|#eta| < 2.5, |#Delta#eta| < 1.3")
#l1.DrawLatex(0.18,0.20, "Wide Jets")
#l1.SetTextSize(0.035)
#l1.DrawLatex(0.70,0.52, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}")
#l1.SetTextSize(0.055)
#l1.DrawLatex(0.50,0.80, "0, 1 and 2 b-tags")

#gPad.RedrawAxis();

c.SetGridx()
c.SetGridy()
c.SaveAs('theory_xs_obs_limit_ratio.eps')

