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
xs_obs_limits = array('d', [1.3648100000000001, 0.58711899999999995, 0.27155200000000002, 0.183139, 0.13209099999999999, 0.12962000000000001, 0.108418, 0.0842802, 0.060856199999999999, 0.0575132, 0.077734600000000001, 0.078043399999999999, 0.060779100000000003, 0.040989100000000001, 0.0360913, 0.038837299999999998, 0.035022699999999997, 0.026342399999999998, 0.0167986, 0.0113588, 0.0075952099999999998, 0.0055416399999999996, 0.0040954900000000002, 0.0035792699999999998, 0.00344086, 0.0033120699999999999, 0.0031159500000000001, 0.0029412100000000001, 0.0026468500000000001, 0.0023606199999999999, 0.00205522])
xs_exp_limits = array('d', [0.485738, 0.38718799999999998, 0.29228599999999999, 0.21132999999999999, 0.18193000000000001, 0.13775399999999999, 0.11717, 0.093375100000000003, 0.075917499999999999, 0.061835599999999998, 0.053033700000000003, 0.041665800000000003, 0.037294800000000003, 0.029347999999999999, 0.024883499999999999, 0.0215517, 0.0182154, 0.014678699999999999, 0.0129204, 0.011520000000000001, 0.0097334699999999993, 0.0082745100000000005, 0.0072972200000000001, 0.0062461000000000001, 0.0052646500000000001, 0.0044042899999999999, 0.00389595, 0.0034369499999999998, 0.00294407, 0.0025829, 0.00217556])
xs_exp_limits_1sigma = array('d', [0.43630400000000003, 0.32172899999999999, 0.24041899999999999, 0.18268999999999999, 0.138153, 0.110932, 0.0981123, 0.071031399999999995, 0.058592499999999999, 0.050565400000000003, 0.040761800000000001, 0.032824699999999998, 0.027845200000000001, 0.022961499999999999, 0.018748500000000001, 0.017023, 0.014075799999999999, 0.011798599999999999, 0.0098583799999999999, 0.0087295500000000009, 0.0076932600000000004, 0.0065515199999999999, 0.0057051799999999998, 0.0046342500000000003, 0.0039584399999999997, 0.0033864899999999998, 0.0029631900000000001, 0.0024761399999999999, 0.0021318399999999999, 0.00190876, 0.00169425, 0.0030458600000000001, 0.0036301200000000001, 0.0037631100000000001, 0.0043824800000000002, 0.0051964400000000001, 0.0063956400000000002, 0.0075242299999999998, 0.0087268399999999996, 0.0100743, 0.010973200000000001, 0.013247, 0.015426499999999999, 0.017418099999999999, 0.021967199999999999, 0.025088900000000001, 0.0297872, 0.034244400000000001, 0.039367600000000003, 0.053964499999999999, 0.060651400000000001, 0.075624499999999997, 0.086407499999999998, 0.108255, 0.12942400000000001, 0.170566, 0.19799600000000001, 0.26902300000000001, 0.31684400000000001, 0.46920099999999998, 0.56133900000000003, 1.0550200000000001])
xs_exp_limits_1sigma_up = array('d', [1.0550200000000001, 0.56133900000000003, 0.46920099999999998, 0.31684400000000001, 0.26902300000000001, 0.19799600000000001, 0.170566, 0.12942400000000001, 0.108255, 0.086407499999999998, 0.075624499999999997, 0.060651400000000001, 0.053964499999999999, 0.039367600000000003, 0.034244400000000001, 0.0297872, 0.025088900000000001, 0.021967199999999999, 0.017418099999999999, 0.015426499999999999, 0.013247, 0.010973200000000001, 0.0100743, 0.0087268399999999996, 0.0075242299999999998, 0.0063956400000000002, 0.0051964400000000001, 0.0043824800000000002, 0.0037631100000000001, 0.0036301200000000001, 0.0030458600000000001])
xs_exp_limits_2sigma = array('d', [0.41592699999999999, 0.28797299999999998, 0.215001, 0.15284800000000001, 0.116769, 0.097762299999999996, 0.082964099999999999, 0.0566831, 0.050276099999999997, 0.042489300000000001, 0.0308177, 0.026813699999999999, 0.022752999999999999, 0.0167624, 0.014708199999999999, 0.0133565, 0.012002799999999999, 0.0099016999999999994, 0.0083029499999999999, 0.0069091100000000004, 0.0055788699999999997, 0.0050905899999999999, 0.0041207400000000003, 0.0037357200000000001, 0.0032398800000000001, 0.0025992300000000001, 0.0021666599999999999, 0.00189145, 0.00175124, 0.0016241599999999999, 0.0013480499999999999, 0.0040998700000000003, 0.0046907700000000004, 0.0049422900000000002, 0.0058606800000000001, 0.0067569800000000001, 0.0078149699999999992, 0.0091602799999999998, 0.010872100000000001, 0.0131544, 0.014439, 0.016490899999999999, 0.0191985, 0.023836099999999999, 0.0274966, 0.0316607, 0.037979400000000003, 0.043410400000000002, 0.0506413, 0.067983299999999997, 0.0782254, 0.090525099999999997, 0.11908299999999999, 0.142204, 0.18047299999999999, 0.23049500000000001, 0.28344999999999998, 0.37257899999999999, 0.46100400000000002, 0.63282799999999995, 0.89756899999999995, 1.33544])
xs_exp_limits_2sigma_up = array('d', [1.33544, 0.89756899999999995, 0.63282799999999995, 0.46100400000000002, 0.37257899999999999, 0.28344999999999998, 0.23049500000000001, 0.18047299999999999, 0.142204, 0.11908299999999999, 0.090525099999999997, 0.0782254, 0.067983299999999997, 0.0506413, 0.043410400000000002, 0.037979400000000003, 0.0316607, 0.0274966, 0.023836099999999999, 0.0191985, 0.016490899999999999, 0.014439, 0.0131544, 0.010872100000000001, 0.0091602799999999998, 0.0078149699999999992, 0.0067569800000000001, 0.0058606800000000001, 0.0049422900000000002, 0.0046907700000000004, 0.0040998700000000003])

masses_exp = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0, 1100.0, 1000.0])

m_signal = array('d', [700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800.])
signal = array('d', [0.6246E+01, 0.3427E+01, 0.1969E+01, 0.1172E+01, 0.7171E+00, 0.4486E+00, 0.2857E+00, 0.1845E+00, 0.1206E+00, 0.7961E-01, 0.5295E-01, 0.3545E-01, 0.2386E-01, 0.1611E-01, 0.1092E-01, 0.7413E-02, 0.5039E-02, 0.3426E-02, 0.2329E-02, 0.1580E-02, 0.1070E-02, 0.7231E-03, 0.4867E-03, 0.3261E-03, 0.2174E-03, 0.1440E-03, 0.9477E-04, 0.6190E-04, 0.4007E-04])


for i in range(0,len(masses)):
 xs_obs_limits[i] = log(xs_obs_limits[i])
 xs_exp_limits[i] = log(xs_exp_limits[i])
 xs_exp_limits_1sigma[i] = log(xs_exp_limits_1sigma[i])
 xs_exp_limits_1sigma_up[i] = log(xs_exp_limits_1sigma_up[i])
 xs_exp_limits_2sigma[i] = log(xs_exp_limits_2sigma[i])
 xs_exp_limits_2sigma_up[i] = log(xs_exp_limits_2sigma_up[i])

for i in range(0,len(m_signal)):
 signal[i] = log(signal[i])

graph_obs_limit = TGraph(len(masses),masses,xs_obs_limits)
graph_exp_limit = TGraph(len(masses),masses,xs_exp_limits)
graph_exp_1sigma_limit = TGraph(len(masses),masses,xs_exp_limits_1sigma)
graph_exp_1sigma_up_limit = TGraph(len(masses),masses,xs_exp_limits_1sigma_up)
graph_exp_2sigma_limit = TGraph(len(masses),masses,xs_exp_limits_2sigma)
graph_exp_2sigma_up_limit = TGraph(len(masses),masses,xs_exp_limits_2sigma_up)

graph_signal = TGraph(len(m_signal),m_signal,signal)


m_final = array('d')
m_exp_final = array('d')
signal_obs_ratio = array('d')
signal_exp_ratio = array('d')
signal_exp_1sigma_ratio = array('d')
signal_exp_1sigma_up_ratio = array('d')
signal_exp_2sigma_ratio = array('d')
signal_exp_2sigma_up_ratio = array('d')


for i in range(0,18):

 mass = 1000. + 100.*i

 m_final.append(mass)
 m_exp_final.append(mass)

 signal_obs_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_obs_limit.Eval(mass)) )
 signal_exp_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_limit.Eval(mass)) )
 signal_exp_1sigma_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_1sigma_limit.Eval(mass)) )
 signal_exp_1sigma_up_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_1sigma_up_limit.Eval(mass)) )
 signal_exp_2sigma_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_2sigma_limit.Eval(mass)) )
 signal_exp_2sigma_up_ratio.append( exp(graph_signal.Eval(mass)) / exp(graph_exp_2sigma_up_limit.Eval(mass)) )

for i in range(0,18):
  m_exp_final.append( masses[18-i-1] )
  signal_exp_1sigma_ratio.append( signal_exp_1sigma_up_ratio[18-i-1] )
  signal_exp_2sigma_ratio.append( signal_exp_2sigma_up_ratio[18-i-1] )

#print m_final
#print signal_obs_ratio
#print signal_exp_ratio

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
graph_signal_exp_1sigma_ratio.SetFillColor(kGreen+1)

graph_signal_exp_2sigma_ratio = TGraph(len(m_exp_final),m_exp_final,signal_exp_2sigma_ratio)
graph_signal_exp_2sigma_ratio.GetXaxis().SetTitle("Resonance Mass [GeV]")
graph_signal_exp_2sigma_ratio.GetYaxis().SetTitle("(Theory #sigma) / (95% CL Limit on #sigma)")
graph_signal_exp_2sigma_ratio.GetYaxis().SetRangeUser(0,2.3)
graph_signal_exp_2sigma_ratio.GetXaxis().SetNdivisions(505)
graph_signal_exp_2sigma_ratio.SetFillColor(kYellow)

c = TCanvas("c", "",800,800)
c.cd()

graph_signal_exp_2sigma_ratio.Draw("AF")
graph_signal_exp_1sigma_ratio.Draw("F")
graph_signal_obs_ratio.Draw("L")
graph_signal_exp_ratio.Draw("L")


line = TLine(1000,1,2700,1);
line.SetLineWidth(2)
line.SetLineColor(kRed)
line.SetLineStyle(1)
line.Draw("same")

legend = TLegend(.36,.79,.56,.92)
#legend.SetBorderSize(0)
legend.SetFillColor(0)
#legend.SetFillStyle(0)
legend.SetShadowColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
#legend.SetHeader('95% CL Upper Limits (stat. only)')
legend.AddEntry(graph_signal_obs_ratio,"Observed","l")
legend.AddEntry(graph_signal_exp_ratio,"Expected","l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(12)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)
l1.DrawLatex(0.70,0.89, "Z' (f_{b#bar{b}} #approx 0.2)")
l1.SetTextSize(0.04)
l1.DrawLatex(0.62,0.76, "CMS Preliminary")
l1.DrawLatex(0.62,0.68, "#intLdt = 5 fb^{-1}")
l1.DrawLatex(0.63,0.63, "#sqrt{s} = 7 TeV")
l1.DrawLatex(0.62,0.58, "|#eta| < 2.5, |#Delta#eta| < 1.3")
l1.DrawLatex(0.62,0.53, "Wide Jets")
l1.SetTextSize(0.035)
l1.DrawLatex(0.70,0.40, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}")
l1.SetTextSize(0.05)
l1.DrawLatex(0.18,0.20, "0, 1 and 2 b-tags")

gPad.RedrawAxis();

c.SetGridx()
c.SetGridy()
c.SaveAs('theory_xs_obsexp_limit_ratio_Zprime.eps')

