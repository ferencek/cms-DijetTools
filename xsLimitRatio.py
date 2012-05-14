#!/usr/bin/env python

import string, re
from ROOT import *
from array import array


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

comb = True
#comb = False

# 2-tag, f_bbbar=0.5
xs_stat = [0.50567200000000001, 0.48394500000000001, 0.28312199999999998, 0.178814, 0.13411200000000001, 0.119211, 0.13411000000000001, 0.165992, 0.15706500000000001, 0.14152300000000001, 0.124318, 0.10662199999999999, 0.090940300000000002, 0.081103900000000007, 0.079797099999999996, 0.074583700000000003, 0.0618661, 0.050253100000000002, 0.0391156, 0.0298023, 0.022351699999999999, 0.018626500000000001, 0.018626500000000001, 0.018626500000000001, 0.0204891, 0.0213805, 0.019093300000000001, 0.018626500000000001, 0.018626500000000001, 0.016763799999999999, 0.016763799999999999]
xs_sys_all = [0.55037499999999995, 0.51374799999999998, 0.357628, 0.208616, 0.14901300000000001, 0.14156299999999999, 0.14901200000000001, 0.17344300000000001, 0.17196600000000001, 0.16014999999999999, 0.14666999999999999, 0.12152300000000001, 0.105841, 0.096005099999999996, 0.087247599999999995, 0.082034300000000004, 0.074904600000000002, 0.063291600000000003, 0.050291500000000003, 0.037252899999999999, 0.0298023, 0.022351699999999999, 0.0204891, 0.0204891, 0.022351699999999999, 0.0213805, 0.0209559, 0.020489199999999999, 0.020489199999999999, 0.018626500000000001, 0.018626500000000001]
xs_sys_lumi = [0.50567200000000001, 0.48394500000000001, 0.28312199999999998, 0.178814, 0.13411200000000001, 0.119211, 0.13411000000000001, 0.165992, 0.15706500000000001, 0.14152300000000001, 0.124318, 0.10662199999999999, 0.090940300000000002, 0.081103900000000007, 0.079797099999999996, 0.076446299999999995, 0.0618661, 0.050253100000000002, 0.037253000000000001, 0.0298023, 0.022351699999999999, 0.018626500000000001, 0.018626500000000001, 0.018626500000000001, 0.0204891, 0.0213805, 0.019093300000000001, 0.018626500000000001, 0.018626500000000001, 0.016763799999999999, 0.016763799999999999]
xs_sys_JES = [0.52057299999999995, 0.48394500000000001, 0.342727, 0.208616, 0.13411200000000001, 0.126661, 0.14156099999999999, 0.15854199999999999, 0.16078999999999999, 0.14524799999999999, 0.131769, 0.110348, 0.094665600000000003, 0.084829199999999993, 0.079797099999999996, 0.076446299999999995, 0.069316600000000006, 0.057703699999999997, 0.048428899999999997, 0.037252899999999999, 0.026076999999999999, 0.022351699999999999, 0.0204891, 0.018626500000000001, 0.0204891, 0.019517799999999998, 0.019093300000000001, 0.018626500000000001, 0.018626500000000001, 0.018626500000000001, 0.016763799999999999]
xs_sys_JER = [0.50567200000000001, 0.46904400000000002, 0.28312199999999998, 0.178814, 0.13411200000000001, 0.119211, 0.14156099999999999, 0.165992, 0.15706500000000001, 0.14524799999999999, 0.12804399999999999, 0.10662199999999999, 0.094665600000000003, 0.081103900000000007, 0.081659700000000002, 0.074583700000000003, 0.0618661, 0.050253100000000002, 0.037253000000000001, 0.0298023, 0.022351699999999999, 0.0204891, 0.018626500000000001, 0.018626500000000001, 0.0204891, 0.0213805, 0.019093300000000001, 0.018626500000000001, 0.018626500000000001, 0.016763799999999999, 0.016763799999999999]
xs_sys_btag = [0.50567200000000001, 0.48394500000000001, 0.28312199999999998, 0.178814, 0.13411200000000001, 0.119211, 0.14156099999999999, 0.165992, 0.16451499999999999, 0.14524799999999999, 0.131769, 0.110348, 0.094665600000000003, 0.084829199999999993, 0.083522399999999997, 0.080171599999999996, 0.065591300000000005, 0.050253100000000002, 0.040978300000000002, 0.0298023, 0.022351699999999999, 0.0204891, 0.018626500000000001, 0.0204891, 0.022351699999999999, 0.0213805, 0.0209559, 0.020489199999999999, 0.018626500000000001, 0.018626500000000001, 0.018626500000000001]
xs_sys_bkg = [0.52057299999999995, 0.49884600000000001, 0.29802299999999998, 0.178814, 0.13411200000000001, 0.119211, 0.14156099999999999, 0.16971800000000001, 0.16078999999999999, 0.14524799999999999, 0.12804399999999999, 0.10662199999999999, 0.094665600000000003, 0.081103900000000007, 0.081659700000000002, 0.076446299999999995, 0.063728699999999999, 0.050253100000000002, 0.040978300000000002, 0.0298023, 0.022351699999999999, 0.018626500000000001, 0.018626500000000001, 0.018626500000000001, 0.0204891, 0.0213805, 0.019093300000000001, 0.018626500000000001, 0.018626500000000001, 0.016763799999999999, 0.016763799999999999]

if comb:
 # Combined, f_bbbar=0.5
 xs_stat = [0.43779699999999999, 0.38328499999999999, 0.24687100000000001, 0.14901200000000001, 0.096857499999999999, 0.096857499999999999, 0.104308, 0.104522, 0.081644700000000001, 0.063330300000000006, 0.065327399999999994, 0.065512200000000007, 0.053470700000000003, 0.042758999999999998, 0.038424899999999998, 0.038156099999999998, 0.033347300000000003, 0.0249803, 0.017419, 0.0121072, 0.0093132400000000004, 0.00745058, 0.0055879399999999996, 0.0051222699999999999, 0.0051222699999999999, 0.0046566100000000003, 0.0041909499999999997, 0.0041909499999999997, 0.0037253, 0.0032596299999999999, 0.0027939699999999998]
 xs_sys_all = [0.46759899999999999, 0.42053800000000002, 0.291574, 0.178814, 0.111759, 0.104308, 0.111759, 0.108247, 0.092820600000000003, 0.070780899999999994, 0.069052699999999995, 0.067374799999999999, 0.059058699999999999, 0.049278299999999997, 0.041218900000000003, 0.038156099999999998, 0.035209900000000002, 0.0305682, 0.023938299999999999, 0.016763799999999999, 0.011175900000000001, 0.0083818999999999994, 0.0065192599999999998, 0.0055879399999999996, 0.0051222699999999999, 0.0046566100000000003, 0.0046566100000000003, 0.0041909499999999997, 0.0037253, 0.0034924600000000002, 0.0032596299999999999]
 xs_sys_lumi = [0.43779699999999999, 0.38328499999999999, 0.24687100000000001, 0.14901200000000001, 0.096857499999999999, 0.096857499999999999, 0.104308, 0.104522, 0.081644700000000001, 0.063330300000000006, 0.065327399999999994, 0.065512200000000007, 0.053470700000000003, 0.041827700000000002, 0.038424899999999998, 0.038156099999999998, 0.033347300000000003, 0.0249803, 0.017419, 0.0121072, 0.0093132400000000004, 0.00745058, 0.0055879399999999996, 0.0051222699999999999, 0.0051222699999999999, 0.0046566100000000003, 0.0041909499999999997, 0.0041909499999999997, 0.0037253, 0.0032596299999999999, 0.0027939699999999998]
 xs_sys_JES = [0.45269799999999999, 0.40563700000000003, 0.26922299999999999, 0.178814, 0.111759, 0.096857499999999999, 0.104308, 0.104522, 0.089095300000000002, 0.070780899999999994, 0.065327399999999994, 0.063649499999999998, 0.059058699999999999, 0.047415699999999998, 0.040287499999999997, 0.037224800000000002, 0.034278599999999999, 0.029636900000000001, 0.023938299999999999, 0.016763799999999999, 0.011175900000000001, 0.0083818999999999994, 0.0065192599999999998, 0.0055879399999999996, 0.0051222699999999999, 0.0046566100000000003, 0.0044237800000000004, 0.0041909499999999997, 0.0037253, 0.0034924600000000002, 0.0032596299999999999]
 xs_sys_JER = [0.43779699999999999, 0.38328499999999999, 0.24687100000000001, 0.14901200000000001, 0.104308, 0.096857499999999999, 0.104308, 0.108247, 0.081644700000000001, 0.067055600000000007, 0.06719, 0.065512200000000007, 0.055333399999999998, 0.043690399999999997, 0.039356200000000001, 0.038156099999999998, 0.033347300000000003, 0.0249803, 0.017419, 0.0130385, 0.0093132400000000004, 0.00745058, 0.0055879399999999996, 0.0051222699999999999, 0.0051222699999999999, 0.0046566100000000003, 0.0041909499999999997, 0.0039581199999999999, 0.0037253, 0.0032596299999999999, 0.0027939699999999998]
 xs_sys_btag = [0.43779699999999999, 0.38328499999999999, 0.24687100000000001, 0.14901200000000001, 0.104308, 0.096857499999999999, 0.104308, 0.104522, 0.081644700000000001, 0.063330300000000006, 0.065327399999999994, 0.065512200000000007, 0.053470700000000003, 0.041827700000000002, 0.038424899999999998, 0.038156099999999998, 0.033347300000000003, 0.0249803, 0.017419, 0.0121072, 0.0093132400000000004, 0.00745058, 0.0055879399999999996, 0.0051222699999999999, 0.0051222699999999999, 0.0046566100000000003, 0.0041909499999999997, 0.0039581199999999999, 0.0037253, 0.0032596299999999999, 0.0027939699999999998]
 xs_sys_bkg = [0.45269799999999999, 0.39073600000000003, 0.25432100000000002, 0.14901200000000001, 0.104308, 0.096857499999999999, 0.104308, 0.108247, 0.085370000000000001, 0.067055600000000007, 0.06719, 0.067374799999999999, 0.055333399999999998, 0.043690399999999997, 0.039356200000000001, 0.039087400000000001, 0.033347300000000003, 0.0249803, 0.017419, 0.0125729, 0.0093132400000000004, 0.00745058, 0.0055879399999999996, 0.0051222699999999999, 0.0051222699999999999, 0.0046566100000000003, 0.0041909499999999997, 0.0041909499999999997, 0.0037253, 0.0032596299999999999, 0.0027939699999999998]

masses = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0])

r_all = array('d')
r_lumi = array('d')
r_JES = array('d')
r_JER = array('d')
r_btag = array('d')
r_bkg = array('d')

for i in range(0,len(xs_stat)):
  r_all.append( xs_sys_all[i] / xs_stat[i] )
  r_lumi.append( xs_sys_lumi[i] / xs_stat[i] )
  r_JES.append( xs_sys_JES[i] / xs_stat[i] )
  r_JER.append( xs_sys_JER[i] / xs_stat[i] )
  r_btag.append( xs_sys_btag[i] / xs_stat[i] )
  r_bkg.append( xs_sys_bkg[i] / xs_stat[i] )

g_lumi = TGraph(len(masses),masses,r_lumi)
g_lumi.SetMarkerStyle(24)
g_lumi.SetMarkerColor(kGreen+2)
g_lumi.SetLineWidth(2)
g_lumi.SetLineStyle(2)
g_lumi.SetLineColor(kGreen+2)
g_lumi.GetXaxis().SetTitle("Resonance Mass [GeV]")
g_lumi.GetYaxis().SetTitle("Limit Ratio (Syst./Stat. only)")
g_lumi.GetYaxis().SetRangeUser(0.35,1.82)
g_lumi.GetXaxis().SetNdivisions(1005)

g_JES = TGraph(len(masses),masses,r_JES)
g_JES.SetMarkerStyle(25)
g_JES.SetMarkerColor(kBlue)
g_JES.SetLineWidth(2)
g_JES.SetLineStyle(3)
g_JES.SetLineColor(kBlue)

g_JER = TGraph(len(masses),masses,r_JER)
g_JER.SetMarkerStyle(26)
g_JER.SetMarkerColor(45)
g_JER.SetLineWidth(2)
g_JER.SetLineStyle(4)
g_JER.SetLineColor(45)

g_btag = TGraph(len(masses),masses,r_btag)
g_btag.SetMarkerStyle(27)
g_btag.SetMarkerColor(kMagenta)
g_btag.SetLineWidth(2)
g_btag.SetLineStyle(5)
g_btag.SetLineColor(kMagenta)

g_bkg = TGraph(len(masses),masses,r_bkg)
g_bkg.SetMarkerStyle(30)
g_bkg.SetMarkerColor(kRed)
g_bkg.SetLineWidth(2)
g_bkg.SetLineStyle(6)
g_bkg.SetLineColor(kRed)

g_all = TGraph(len(masses),masses,r_all)
g_all.SetMarkerStyle(20)
g_all.SetMarkerColor(kBlack)
g_all.SetLineWidth(2)
g_all.SetLineStyle(1)
g_all.SetLineColor(kBlack)

c = TCanvas("c", "",800,800)
c.cd()

g_lumi.Draw("ALP")
g_JES.Draw("LP")
g_JER.Draw("LP")
g_btag.Draw("LP")
g_bkg.Draw("LP")
g_all.Draw("LP")

legend = TLegend(.20,.65,.60,.85)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
legend.AddEntry(g_all, "All / Stat. only","lp")
legend.AddEntry(g_JES, "JES / Stat. only","lp")
legend.AddEntry(g_JER, "JER / Stat. only","lp")
legend.AddEntry(g_bkg, "Background / Stat. only","lp")
legend.AddEntry(g_btag, "b-tag SF / Stat. only","lp")
legend.AddEntry(g_lumi, "Lumi / Stat. only","lp")

legend.Draw()

l1 = TLatex()
l1.SetTextAlign(12)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.035)
l1.DrawLatex(0.70,0.35, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}")
l1.SetTextSize(0.04)
l1.DrawLatex(0.18,0.89, "RS-graviton-like, f_{b#bar{b}} = 0.5")
l1.SetTextSize(0.04)
l1.DrawLatex(0.18,0.43, "CMS Preliminary")
l1.DrawLatex(0.18,0.35, "#intLdt = 5 fb^{-1}")
l1.DrawLatex(0.19,0.30, "#sqrt{s} = 7 TeV")
l1.DrawLatex(0.18,0.25, "|#eta| < 2.5, |#Delta#eta| < 1.3")
l1.DrawLatex(0.18,0.20, "Wide Jets")
l1.SetTextSize(0.055)
if comb: l1.DrawLatex(0.54,0.205, "0, 1 and 2 b-tags")
else:    l1.DrawLatex(0.73,0.205, "2 b-tags")

if comb: c.SaveAs('CSVL_Combined_xs_limit_ratio_WideJets_RSG_fbbbar0p5.eps')
else:    c.SaveAs('CSVL_2Tag_xs_limit_ratio_WideJets_RSG_fbbbar0p5.eps')

