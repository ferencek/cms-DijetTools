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


# syst+stat
masses = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0])
xs_obs_limits = array('d', [0.61651699999999998, 0.64604499999999998, 0.21223, 0.106056, 0.0684279, 0.063028699999999993, 0.069857000000000002, 0.093030399999999999, 0.073984900000000006, 0.055434499999999998, 0.058917700000000003, 0.056365800000000001, 0.046353699999999998, 0.037593399999999999, 0.036586, 0.038987599999999997, 0.036139200000000003, 0.029617399999999999, 0.018850599999999999, 0.013448099999999999, 0.0097049700000000003, 0.0070975300000000003, 0.0057316299999999997, 0.0046952699999999997, 0.0043628499999999997, 0.0040803899999999997, 0.0037947599999999999, 0.003447, 0.0032466299999999999, 0.0028686499999999999, 0.0025190299999999998])
xs_exp_limits = array('d', [0.24473600000000001, 0.19852, 0.15547900000000001, 0.12864600000000001, 0.10428900000000001, 0.086690299999999998, 0.078473000000000001, 0.065750100000000006, 0.056944399999999999, 0.049186100000000003, 0.043140299999999999, 0.036887700000000002, 0.032925299999999998, 0.0278639, 0.024965299999999999, 0.020008499999999999, 0.018001699999999999, 0.0163096, 0.014044299999999999, 0.0121666, 0.0107533, 0.0096080000000000002, 0.0082403600000000004, 0.0073010599999999998, 0.00645394, 0.0056334699999999998, 0.0048711500000000003, 0.0043512400000000001, 0.0037534999999999999, 0.0032833900000000002, 0.0028720999999999998])
xs_exp_limits_1sigma = array('d', [0.230077, 0.172379, 0.138352, 0.108768, 0.087911500000000004, 0.071576600000000004, 0.062381499999999999, 0.054009300000000003, 0.049493000000000002, 0.0406264, 0.034007500000000003, 0.029669899999999999, 0.026186899999999999, 0.0229819, 0.019173699999999998, 0.016801699999999999, 0.014201800000000001, 0.0127403, 0.0110461, 0.00954348, 0.0082169900000000004, 0.0073561700000000004, 0.00662862, 0.0057504100000000001, 0.0048199799999999998, 0.0041103499999999996, 0.0037769000000000001, 0.00329934, 0.00288594, 0.0024359300000000002, 0.0020966299999999999, 0.00386171, 0.0043015400000000004, 0.00505012, 0.0058275899999999997, 0.0066554600000000002, 0.0080018599999999995, 0.0089053799999999992, 0.0100944, 0.0114018, 0.01311, 0.014928800000000001, 0.0164914, 0.0202256, 0.0226237, 0.025392399999999999, 0.028238300000000001, 0.0366719, 0.037945199999999998, 0.045357300000000003, 0.052951900000000003, 0.063769999999999993, 0.073031200000000004, 0.085376499999999994, 0.095834100000000005, 0.12121700000000001, 0.143344, 0.166016, 0.20974300000000001, 0.22234899999999999, 0.30954100000000001, 0.48728700000000003])
xs_exp_limits_1sigma_up = array('d', [0.48728700000000003, 0.30954100000000001, 0.22234899999999999, 0.20974300000000001, 0.166016, 0.143344, 0.12121700000000001, 0.095834100000000005, 0.085376499999999994, 0.073031200000000004, 0.063769999999999993, 0.052951900000000003, 0.045357300000000003, 0.037945199999999998, 0.0366719, 0.028238300000000001, 0.025392399999999999, 0.0226237, 0.0202256, 0.0164914, 0.014928800000000001, 0.01311, 0.0114018, 0.0100944, 0.0089053799999999992, 0.0080018599999999995, 0.0066554600000000002, 0.0058275899999999997, 0.00505012, 0.0043015400000000004, 0.00386171])

masses_exp = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0, 1100.0, 1000.0])


m_signal = array('d', [700., 1000., 1200., 1400., 1700., 2100., 2700.])
signal = array('d', [6.52288, 0.995265, 0.319478, 0.108527, 0.0235598, 0.00331955, 0.000188328])


for i in range(0,len(masses)):
 xs_obs_limits[i] = log(xs_obs_limits[i])
 xs_exp_limits[i] = log(xs_exp_limits[i])
 xs_exp_limits_1sigma[i] = log(xs_exp_limits_1sigma[i])
 xs_exp_limits_1sigma_up[i] = log(xs_exp_limits_1sigma_up[i])

for i in range(0,len(m_signal)):
 signal[i] = log(signal[i])

graph_obs_limit = TGraph(len(masses),masses,xs_obs_limits)
graph_exp_limit = TGraph(len(masses),masses,xs_exp_limits)
graph_exp_1sigma_limit = TGraph(len(masses),masses,xs_exp_limits_1sigma)
graph_exp_1sigma_up_limit = TGraph(len(masses),masses,xs_exp_limits_1sigma_up)

graph_signal = TGraph(len(m_signal),m_signal,signal)


m_final = array('d')
m_exp_final = array('d')
signal_obs_ratio = array('d')
signal_exp_ratio = array('d')
signal_exp_1sigma_ratio = array('d')
signal_exp_1sigma_up_ratio = array('d')


for i in range(0,18):

 mass = 1000. + 100.*i

 m_final.append(mass)
 m_exp_final.append(mass)

 signal_obs_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_obs_limit.Eval(mass)) )
 signal_exp_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_limit.Eval(mass)) )
 signal_exp_1sigma_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_1sigma_limit.Eval(mass)) )
 signal_exp_1sigma_up_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_1sigma_up_limit.Eval(mass)) )

for i in range(0,18):
  m_exp_final.append( masses[18-i-1] )
  signal_exp_1sigma_ratio.append( signal_exp_1sigma_up_ratio[18-i-1] )

#print m_final
#print signal_obs_ratio

graph_signal_obs_ratio = TGraph(len(m_final),m_final,signal_obs_ratio)
graph_signal_obs_ratio.SetLineWidth(2)
graph_signal_obs_ratio.SetLineStyle(1)
graph_signal_obs_ratio.SetLineColor(kBlack)

graph_signal_exp_ratio = TGraph(len(m_final),m_final,signal_exp_ratio)
graph_signal_exp_ratio.SetLineWidth(2)
graph_signal_exp_ratio.SetLineStyle(2)
graph_signal_exp_ratio.SetLineColor(kBlue)
graph_signal_exp_ratio.SetFillColor(kGreen+1)

graph_signal_exp_1sigma_ratio = TGraph(len(m_exp_final),m_exp_final,signal_exp_1sigma_ratio)
graph_signal_exp_1sigma_ratio.GetXaxis().SetTitle("Resonance Mass [GeV]")
graph_signal_exp_1sigma_ratio.GetYaxis().SetTitle("(Theory #sigma) / (95% CL Limit on #sigma)")
graph_signal_exp_1sigma_ratio.GetYaxis().SetRangeUser(0,2.2)
graph_signal_exp_1sigma_ratio.GetXaxis().SetNdivisions(505)
graph_signal_exp_1sigma_ratio.SetFillColor(kGreen+1)

c = TCanvas("c", "",800,800)
c.cd()

graph_signal_exp_1sigma_ratio.Draw("AF")
graph_signal_obs_ratio.Draw("L")
graph_signal_exp_ratio.Draw("L")


line = TLine(1000,1,2700,1);
line.SetLineWidth(2)
line.SetLineColor(kRed)
line.SetLineStyle(1)
line.Draw("same")

legend = TLegend(.50,.63,.90,.76)
#legend.SetBorderSize(0)
legend.SetFillColor(0)
#legend.SetFillStyle(0)
legend.SetShadowColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
#legend.SetHeader('95% CL Upper Limits (stat. only)')
legend.AddEntry(graph_signal_obs_ratio,"Observed","l")
legend.AddEntry(graph_signal_exp_ratio,"Expected with 1#sigma band","lf")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(12)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)
l1.DrawLatex(0.70,0.89, "S8 (f_{b#bar{b}} #approx 1.0)")
#l1.SetTextSize(0.04)
#l1.DrawLatex(0.18,0.43, "CMS Preliminary")
#l1.DrawLatex(0.18,0.35, "#intLdt = 5 fb^{-1}")
#l1.DrawLatex(0.19,0.30, "#sqrt{s} = 7 TeV")
#l1.DrawLatex(0.18,0.25, "|#eta| < 2.5, |#Delta#eta| < 1.3")
#l1.DrawLatex(0.18,0.20, "Wide Jets")
l1.SetTextSize(0.035)
l1.DrawLatex(0.70,0.40, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}")
l1.SetTextSize(0.05)
l1.DrawLatex(0.18,0.20, "0, 1 and 2 b-tags")

gPad.RedrawAxis();

c.SetGridx()
c.SetGridy()
c.SaveAs('theory_xs_obsexp_limit_ratio_S8.eps')

