#!/usr/bin/env python

import string, re
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE)
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


PtBins_b = array('d',[30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670])

# CSVL
CSVL_SFb_0to2p4 = TF1("CSVL_SFb_0to2p4","1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x)))", 30.,670.)

CSVL_SFb_errors = TH1D("CSVL_SFb_errors", "CSVL_SFb_errors", 14, PtBins_b)
CSVL_SFb_errors.SetBinContent( 0,0.12)
CSVL_SFb_errors.SetBinContent( 1,0.0188743)
CSVL_SFb_errors.SetBinContent( 2,0.0161816)
CSVL_SFb_errors.SetBinContent( 3,0.0139824)
CSVL_SFb_errors.SetBinContent( 4,0.0152644)
CSVL_SFb_errors.SetBinContent( 5,0.0161226)
CSVL_SFb_errors.SetBinContent( 6,0.0157396)
CSVL_SFb_errors.SetBinContent( 7,0.0161619)
CSVL_SFb_errors.SetBinContent( 8,0.0168747)
CSVL_SFb_errors.SetBinContent( 9,0.0257175)
CSVL_SFb_errors.SetBinContent(10,0.026424)
CSVL_SFb_errors.SetBinContent(11,0.0264928)
CSVL_SFb_errors.SetBinContent(12,0.0315127)
CSVL_SFb_errors.SetBinContent(13,0.030734)
CSVL_SFb_errors.SetBinContent(14,0.0438259)
CSVL_SFb_errors.SetBinContent(15,(2*0.0438259))

CSVL_SFl_0to2p4 =   TF1("CSVL_SFl_0to2p4","((1.0344+(0.000962994*x))+(-3.65392e-06*(x*x)))+(3.23525e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_0to0p5 =   TF1("CSVL_SFl_0to0p5","((1.07536+(0.000175506*x))+(-8.63317e-07*(x*x)))+(3.27516e-10*(x*(x*x)))", 20.,670.)
CSVL_SFl_0p5to1p0 = TF1("CSVL_SFl_0p5to1p0","((1.07846+(0.00032458*x))+(-1.30258e-06*(x*x)))+(8.50608e-10*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p0to1p5 = TF1("CSVL_SFl_1p0to1p5","((1.08294+(0.000474818*x))+(-1.43857e-06*(x*x)))+(1.13308e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p5to2p4 = TF1("CSVL_SFl_1p5to2p4","((1.0617+(0.000173654*x))+(-5.29009e-07*(x*x)))+(5.55931e-10*(x*(x*x)))", 20.,670.)

CSVL_SFl_0to2p4_min =   TF1("CSVL_SFl_0to2p4_min","((0.956023+(0.000825106*x))+(-3.18828e-06*(x*x)))+(2.81787e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_0to0p5_min =   TF1("CSVL_SFl_0to0p5_min","((0.994425+(-8.66392e-05*x))+(-3.03813e-08*(x*x)))+(-3.52151e-10*(x*(x*x)))", 20.,670.)
CSVL_SFl_0p5to1p0_min = TF1("CSVL_SFl_0p5to1p0_min","((0.998088+(6.94916e-05*x))+(-4.82731e-07*(x*x)))+(1.63506e-10*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p0to1p5_min = TF1("CSVL_SFl_1p0to1p5_min","((1.00294+(0.000289844*x))+(-7.9845e-07*(x*x)))+(5.38525e-10*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p5to2p4_min = TF1("CSVL_SFl_1p5to2p4_min","((0.979816+(0.000138797*x))+(-3.14503e-07*(x*x)))+(2.38124e-10*(x*(x*x)))", 20.,670.)

CSVL_SFl_0to2p4_max =   TF1("CSVL_SFl_0to2p4_max","((1.11272+(0.00110104*x))+(-4.11956e-06*(x*x)))+(3.65263e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_0to0p5_max =   TF1("CSVL_SFl_0to0p5_max","((1.15628+(0.000437668*x))+(-1.69625e-06*(x*x)))+(1.00718e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_0p5to1p0_max = TF1("CSVL_SFl_0p5to1p0_max","((1.15882+(0.000579711*x))+(-2.12243e-06*(x*x)))+(1.53771e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p0to1p5_max = TF1("CSVL_SFl_1p0to1p5_max","((1.16292+(0.000659848*x))+(-2.07868e-06*(x*x)))+(1.72763e-09*(x*(x*x)))", 20.,670.)
CSVL_SFl_1p5to2p4_max = TF1("CSVL_SFl_1p5to2p4_max","((1.14357+(0.00020854*x))+(-7.43519e-07*(x*x)))+(8.73742e-10*(x*(x*x)))", 20.,670.)

# CSVM
CSVM_SFb_0to2p4 = TF1("CSVM_SFb_0to2p4","0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)))", 30.,670.)

CSVM_SFb_errors = TH1D("CSVM_SFb_errors", "CSVM_SFb_errors", 14, PtBins_b)
CSVM_SFb_errors.SetBinContent( 0,0.12)
CSVM_SFb_errors.SetBinContent( 1,0.0295675)
CSVM_SFb_errors.SetBinContent( 2,0.0295095)
CSVM_SFb_errors.SetBinContent( 3,0.0210867)
CSVM_SFb_errors.SetBinContent( 4,0.0219349)
CSVM_SFb_errors.SetBinContent( 5,0.0227033)
CSVM_SFb_errors.SetBinContent( 6,0.0204062)
CSVM_SFb_errors.SetBinContent( 7,0.0185857)
CSVM_SFb_errors.SetBinContent( 8,0.0256242)
CSVM_SFb_errors.SetBinContent( 9,0.0383341)
CSVM_SFb_errors.SetBinContent(10,0.0409675)
CSVM_SFb_errors.SetBinContent(11,0.0420284)
CSVM_SFb_errors.SetBinContent(12,0.0541299)
CSVM_SFb_errors.SetBinContent(13,0.0578761)
CSVM_SFb_errors.SetBinContent(14,0.0655432)
CSVM_SFb_errors.SetBinContent(15,(2*0.0655432))

CSVM_SFl_0to2p4 =   TF1("CSVM_SFl_0to2p4","((1.04318+(0.000848162*x))+(-2.5795e-06*(x*x)))+(1.64156e-09*(x*(x*x)))", 20.,670.)
CSVM_SFl_0to0p8 =   TF1("CSVM_SFl_0to0p8","((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)))", 20.,670.)
CSVM_SFl_0p8to1p6 = TF1("CSVM_SFl_0p8to1p6","((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)))", 20.,670.)
CSVM_SFl_1p6to2p4 = TF1("CSVM_SFl_1p6to2p4","((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)))", 20.,670.)

CSVM_SFl_0to2p4_min =   TF1("CSVM_SFl_0to2p4_min","((0.962627+(0.000448344*x))+(-1.25579e-06*(x*x)))+(4.82283e-10*(x*(x*x)))", 20.,670.)
CSVM_SFl_0to0p8_min =   TF1("CSVM_SFl_0to0p8_min","((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)))", 20.,670.)
CSVM_SFl_0p8to1p6_min = TF1("CSVM_SFl_0p8to1p6_min","((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)))", 20.,670.)
CSVM_SFl_1p6to2p4_min = TF1("CSVM_SFl_1p6to2p4_min","((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)))", 20.,670.)

CSVM_SFl_0to2p4_max =   TF1("CSVM_SFl_0to2p4_max","((1.12368+(0.00124806*x))+(-3.9032e-06*(x*x)))+(2.80083e-09*(x*(x*x)))", 20.,670.)
CSVM_SFl_0to0p8_max =   TF1("CSVM_SFl_0to0p8_max","((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)))", 20.,670.)
CSVM_SFl_0p8to1p6_max = TF1("CSVM_SFl_0p8to1p6_max","((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)))", 20.,670.)
CSVM_SFl_1p6to2p4_max = TF1("CSVM_SFl_1p6to2p4_max","((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)))", 20.,670.)


def plot_SFl(mass_min, mass_max, n_steps, SFl_mean, SFl_min, SFl_max, xAxisTitle, yAxisTitle, label, outputFilename):

  masses = array('d')
  masses_unc = array('d')

  step = (mass_max-mass_min)/float(n_steps)

  for i in range(0,n_steps+1):
    mass = mass_min+float(i)*step
    masses.append(mass)
    masses_unc.append(mass)

  for i in range(0,n_steps+1):
    masses_unc.append(masses[n_steps-i])

  a_SFl = array('d')
  a_SFl_unc = array('d')

  for i in range(0,len(masses)):
    a_SFl.append(SFl_mean.Eval(masses[i]))

  for i in range(0,len(masses_unc)):
    if(i <= n_steps):
      a_SFl_unc.append(SFl_min.Eval(masses_unc[i]))
    else:
      a_SFl_unc.append(SFl_max.Eval(masses_unc[i]))

  g_SFl = TGraph(len(masses),masses,a_SFl)
  g_SFl.SetLineWidth(2)
  g_SFl.SetLineColor(kRed)
  g_SFl.GetXaxis().SetTitle(xAxisTitle)
  g_SFl.GetYaxis().SetTitle(yAxisTitle)
  g_SFl.GetYaxis().SetRangeUser(0.,1.3)

  g_SFl_unc = TGraph(len(masses_unc),masses_unc,a_SFl_unc)
  g_SFl_unc.SetLineColor(0)
  g_SFl_unc.SetFillColor(kRed)
  g_SFl_unc.SetFillStyle(3013)

  c = TCanvas("c", "",800,600)
  c.cd()

  g_SFl.Draw("AL")
  g_SFl_unc.Draw("F")

  legend = TLegend(.55,.55,.85,.60)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  legend.AddEntry(g_SFl_unc, "Uncertainty (stat \oplus syst)","f")
  legend.Draw()

  l1 = TLatex()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetNDC()
  l1.SetTextSize(0.04)
  l1.DrawLatex(0.18,0.34, "CMS Preliminary")
  l1.DrawLatex(0.18,0.26, "#intLdt = 5 fb^{-1}")
  l1.DrawLatex(0.19,0.21, "#sqrt{s} = 7 TeV")
  l1.SetTextSize(0.05)
  l1.DrawLatex(0.57,0.40, label)

  c.SaveAs(outputFilename)

def plot_SFb(mass_min, mass_max, n_steps, SFb_mean, SFb_error, xAxisTitle, yAxisTitle, label, outputFilename):

  masses = array('d')
  masses_unc = array('d')

  step = (mass_max-mass_min)/float(n_steps)

  for i in range(0,n_steps+1):
    mass = mass_min+float(i)*step
    masses.append(mass)
    masses_unc.append(mass)

  for i in range(0,n_steps+1):
    masses_unc.append(masses[n_steps-i])

  a_SFb = array('d')
  a_SFb_unc = array('d')

  for i in range(0,len(masses)):
    a_SFb.append(SFb_mean.Eval(masses[i]))

  for i in range(0,len(masses_unc)):
    if(i <= n_steps):
      a_SFb_unc.append(SFb_mean.Eval(masses_unc[i]) - SFb_error.GetBinContent(SFb_error.GetXaxis().FindBin(masses_unc[i])))
    else:
      a_SFb_unc.append(SFb_mean.Eval(masses_unc[i]) + SFb_error.GetBinContent(SFb_error.GetXaxis().FindBin(masses_unc[i])))

  g_SFb = TGraph(len(masses),masses,a_SFb)
  g_SFb.SetLineWidth(2)
  g_SFb.SetLineColor(kRed)
  g_SFb.GetXaxis().SetTitle(xAxisTitle)
  g_SFb.GetYaxis().SetTitle(yAxisTitle)
  g_SFb.GetYaxis().SetRangeUser(0.,1.3)

  g_SFb_unc = TGraph(len(masses_unc),masses_unc,a_SFb_unc)
  g_SFb_unc.SetLineColor(0)
  g_SFb_unc.SetFillColor(kRed)
  g_SFb_unc.SetFillStyle(3013)

  c = TCanvas("c", "",800,600)
  c.cd()

  g_SFb.Draw("AL")
  g_SFb_unc.Draw("F")

  legend = TLegend(.55,.55,.85,.60)
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  legend.SetTextSize(0.04)
  legend.AddEntry(g_SFb_unc, "Uncertainty (stat \oplus syst)","f")
  legend.Draw()

  l1 = TLatex()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetNDC()
  l1.SetTextSize(0.04)
  l1.DrawLatex(0.18,0.34, "CMS Preliminary")
  l1.DrawLatex(0.18,0.26, "#intLdt = 5 fb^{-1}")
  l1.DrawLatex(0.19,0.21, "#sqrt{s} = 7 TeV")
  l1.SetTextSize(0.05)
  l1.DrawLatex(0.57,0.40, label)

  c.SaveAs(outputFilename)
  
  
if __name__ == "__main__":
  # CSVL
  plot_SFl(30., 669.99, 100, CSVL_SFl_0to2p4, CSVL_SFl_0to2p4_min, CSVL_SFl_0to2p4_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVL}{|#eta| < 2.4}", "SFl_eta0to2p4.eps")
  plot_SFl(30., 669.99, 100, CSVL_SFl_0to0p5, CSVL_SFl_0to0p5_min, CSVL_SFl_0to0p5_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVL}{|#eta| < 0.5}", "SFl_eta0to0p5.eps")
  plot_SFl(30., 669.99, 100, CSVL_SFl_0p5to1p0, CSVL_SFl_0p5to1p0_min, CSVL_SFl_0p5to1p0_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVL}{0.5 #leq |#eta| < 1.0}", "SFl_eta0p5to1p0.eps")
  plot_SFl(30., 669.99, 100, CSVL_SFl_1p0to1p5, CSVL_SFl_1p0to1p5_min, CSVL_SFl_1p0to1p5_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVL}{1.0 #leq |#eta| < 1.5}", "SFl_eta1p0to1p5.eps")
  plot_SFl(30., 669.99, 100, CSVL_SFl_1p5to2p4, CSVL_SFl_1p5to2p4_min, CSVL_SFl_1p5to2p4_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVL}{1.5 #leq |#eta| < 2.4}", "SFl_eta1p5to2p4.eps")

  plot_SFb(30., 669.99, 100, CSVL_SFb_0to2p4, CSVL_SFb_errors, "p_{T} [GeV]", "SF_{b}", "#splitline{CSVL}{|#eta| < 2.4}", "SFb_eta0to2p4.eps")

  # CSVM
  plot_SFl(30., 669.99, 100, CSVM_SFl_0to2p4, CSVM_SFl_0to2p4_min, CSVM_SFl_0to2p4_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVM}{|#eta| < 2.4}", "SFl_CSVM_eta0to2p4.eps")
  plot_SFl(30., 669.99, 100, CSVM_SFl_0to0p8, CSVM_SFl_0to0p8_min, CSVM_SFl_0to0p8_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVM}{|#eta| < 0.8}", "SFl_CSVM_eta0to0p8.eps")
  plot_SFl(30., 669.99, 100, CSVM_SFl_0p8to1p6, CSVM_SFl_0p8to1p6_min, CSVM_SFl_0p8to1p6_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVM}{0.8 #leq |#eta| < 1.6}", "SFl_CSVM_eta0p8to1p6.eps")
  plot_SFl(30., 669.99, 100, CSVM_SFl_1p6to2p4, CSVM_SFl_1p6to2p4_min, CSVM_SFl_1p6to2p4_max, "p_{T} [GeV]", "SF_{l}", "#splitline{CSVM}{1.6 #leq |#eta| < 2.4}", "SFl_CSVM_eta1p6to2p4.eps")

  plot_SFb(30., 669.99, 100, CSVM_SFb_0to2p4, CSVM_SFb_errors, "p_{T} [GeV]", "SF_{b}", "#splitline{CSVM}{|#eta| < 2.4}", "SFb_CSVM_eta0to2p4.eps")
  