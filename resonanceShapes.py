#!/usr/bin/env python

from ROOT import *
from array import array

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
#gStyle.SetPadLeftMargin(0.13);
#gStyle.SetPadRightMargin(0.07);


mean_RSGToQQbar = [491.7, 686.2, 1177., 1958., 3429.]
sigma_RSGToQQbar = [39.08, 49.69, 66.8, 95.16, 143.9]

mean_RSGToBBbar = [474.7, 665., 1146., 1909., 3383.]
sigma_RSGToBBbar = [46.15, 58.62, 76.17, 121.4, 171.5]

mass = [500., 700., 1200., 2000., 3500.]

delta_mean_rel = array('d')
delta_sigma_rel = array('d')
mass_array = array('d')

for i in range(0,len(mass)):
 delta_mean_rel.append((mean_RSGToQQbar[i]-mean_RSGToBBbar[i])/mean_RSGToQQbar[i])
 delta_sigma_rel.append((sigma_RSGToBBbar[i]-sigma_RSGToQQbar[i])/sigma_RSGToQQbar[i])
 mass_array.append(mass[i])
 #print mass[i]

graph_mean = TGraph(len(mass_array),mass_array,delta_mean_rel)
graph_sigma = TGraph(len(mass_array),mass_array,delta_sigma_rel)

graph_mean.SetTitle("#Delta(Mean_{RSG#rightarrowq#bar{q}}-Mean_{RSG#rightarrowb#bar{b}})/Mean_{RSG#rightarrowq#bar{q}}")
graph_mean.SetMarkerStyle(22)
graph_mean.GetXaxis().SetTitle("Resonance Mass [GeV]")

graph_sigma.SetTitle("#Delta(Sigma_{RSG#rightarrowb#bar{b}}-Sigma_{RSG#rightarrowq#bar{q}})/Sigma_{RSG#rightarrowq#bar{q}}")
graph_sigma.SetMarkerStyle(22)
graph_sigma.GetXaxis().SetTitle("Resonance Mass [GeV]")

c = TCanvas("c", "",800,800)
c.cd()

graph_mean.Draw("ap")

c.SaveAs("delta_mean_rel.png")

graph_sigma.Draw("ap")

c.SaveAs("delta_sigma_rel.png")

## Terminate the program
#print "Press ENTER to terminate"
#wait=raw_input()
