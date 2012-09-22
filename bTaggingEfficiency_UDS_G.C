#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "tdrstyle.C"

using namespace std;


void bTagEffVsPt_comp(const string& fInputFile1, const string& fInputFile2,
                      const string& fPlot1, const string& fPlot2,
                      const Int_t nbins, const Double_t *xbins,
                      const string& fLegend1, const string& fLegend2, const string& fAlgo,
                      const string& fXAxisTitle, const string& fYAxisTitle, const string& fOutputFile,
                      const double fXmin = 0, const double fXmax = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH2D *h2_denom1 = (TH2D*)file1->Get((fPlot1 + "_All").c_str());
  TH2D *h2_num1 = (TH2D*)file1->Get((fPlot1 + "_" + fAlgo).c_str());

  TH1D *h1_denom1 = h2_denom1->ProjectionX("_px",46,55);
  TH1D *h1_num1 = h2_num1->ProjectionX("_px",46,55);

  TH1D *h1_denom1_rebinned = (TH1D*)h1_denom1->Rebin(nbins,"h1_denom1_rebinned",xbins);
  TH1D *h1_num1_rebinned = (TH1D*)h1_num1->Rebin(nbins,"h1_num1_rebinned",xbins);

  TH2D *h2_denom2 = (TH2D*)file2->Get((fPlot2 + "_B_All").c_str());
  h2_denom2->Add((TH2D*)file2->Get((fPlot2 + "_C_All").c_str()));
  h2_denom2->Add((TH2D*)file2->Get((fPlot2 + "_UDS_All").c_str()));
  h2_denom2->Add((TH2D*)file2->Get((fPlot2 + "_G_All").c_str()));
  TH2D *h2_num2 = (TH2D*)file2->Get((fPlot2 + "_B_" + fAlgo).c_str());
  h2_num2->Add((TH2D*)file2->Get((fPlot2 + "_C_" + fAlgo).c_str()));
  h2_num2->Add((TH2D*)file2->Get((fPlot2 + "_UDS_" + fAlgo).c_str()));
  h2_num2->Add((TH2D*)file2->Get((fPlot2 + "_G_" + fAlgo).c_str()));

  TH1D *h1_denom2 = h2_denom2->ProjectionX("_px",46,55);
  TH1D *h1_num2 = h2_num2->ProjectionX("_px",46,55);

  TH1D *h1_denom2_rebinned = (TH1D*)h1_denom2->Rebin(nbins,"h1_denom2_rebinned",xbins);
  TH1D *h1_num2_rebinned = (TH1D*)h1_num2->Rebin(nbins,"h1_num2_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_num1_rebinned, h1_denom1_rebinned,"cp");
//   g_efficiency1->SetMarkerSize(0.8);
  g_efficiency1->SetMarkerStyle(20);
  g_efficiency1->SetMarkerColor(kBlack);
  g_efficiency1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(fXmin!=fXmax) g_efficiency1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  g_efficiency1->GetYaxis()->SetRangeUser(0,1.);
  g_efficiency1->GetXaxis()->SetNdivisions(505);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_num2_rebinned, h1_denom2_rebinned,"cp");
//   g_efficiency2->SetMarkerSize(0.8);
  g_efficiency2->SetMarkerStyle(24);
  g_efficiency2->SetMarkerColor(kRed);
  g_efficiency2->SetLineColor(kRed);
  g_efficiency2->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency2->GetYaxis()->SetTitle(fYAxisTitle.c_str());

  g_efficiency1->Draw("AP");
  g_efficiency2->Draw("Psame");

  TLegend *legend = new TLegend(.58,.75,.83,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_efficiency1, fLegend1.c_str(),"lp");
  legend->AddEntry(g_efficiency2, fLegend2.c_str(),"lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.53,0.34, "CMS Simulation");
  l1.DrawLatex(0.54,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.53,0.24, "|#eta| < 0.5");
  if(fLegend1.find("R=")!=string::npos) l1.DrawLatex(0.53,0.19, "Anti-k_{T} PF Jets");
  else l1.DrawLatex(0.53,0.19, "Anti-k_{T} R = 0.7 PF Jets");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.19,0.19, fAlgo.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file1;
  delete file2;
}


void FlavFracVsPt(const string& fInputFile, const string& fPlot,
                  const Int_t nbins, const Double_t *xbins,
                  const string& fXAxisTitle, const string& fYAxisTitle, const string& fOutputFile,
                  const double fXmin = 0, const double fXmax = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_num1 = (TH2D*)((TH2D*)file->Get((fPlot + "_B_All").c_str()))->Clone("h2_num1");
  TH2D *h2_num2 = (TH2D*)file->Get((fPlot + "_C_All").c_str());
  TH2D *h2_num3 = (TH2D*)file->Get((fPlot + "_G_All").c_str());
  TH2D *h2_num4 = (TH2D*)file->Get((fPlot + "_UDS_All").c_str());

  TH2D *h2_denom = (TH2D*)file->Get((fPlot + "_B_All").c_str());
  h2_denom->Add((TH2D*)file->Get((fPlot + "_C_All").c_str()));
  h2_denom->Add((TH2D*)file->Get((fPlot + "_G_All").c_str()));
  h2_denom->Add((TH2D*)file->Get((fPlot + "_UDS_All").c_str()));

  TH1D *h1_denom = h2_denom->ProjectionX("_px",46,55);
  TH1D *h1_num1 = h2_num1->ProjectionX("_px",46,55);
  TH1D *h1_num2 = h2_num2->ProjectionX("_px",46,55);
  TH1D *h1_num3 = h2_num3->ProjectionX("_px",46,55);
  TH1D *h1_num4 = h2_num4->ProjectionX("_px",46,55);

  TH1D *h1_denom_rebinned = (TH1D*)h1_denom->Rebin(nbins,"h1_denom_rebinned",xbins);
  TH1D *h1_num1_rebinned = (TH1D*)h1_num1->Rebin(nbins,"h1_num1_rebinned",xbins);
  TH1D *h1_num2_rebinned = (TH1D*)h1_num2->Rebin(nbins,"h1_num2_rebinned",xbins);
  TH1D *h1_num3_rebinned = (TH1D*)h1_num3->Rebin(nbins,"h1_num3_rebinned",xbins);
  TH1D *h1_num4_rebinned = (TH1D*)h1_num4->Rebin(nbins,"h1_num4_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_flavfrac1 = new TGraphAsymmErrors(h1_num1_rebinned, h1_denom_rebinned,"cp");
//   g_flavfrac1->SetMarkerSize(0.8);
  g_flavfrac1->SetMarkerStyle(24);
  g_flavfrac1->SetMarkerColor(kRed);
  g_flavfrac1->SetLineColor(kRed);
  g_flavfrac1->SetLineStyle(4);
  g_flavfrac1->SetLineWidth(2);
  g_flavfrac1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_flavfrac1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(fXmin!=fXmax) g_flavfrac1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  g_flavfrac1->GetYaxis()->SetRangeUser(0,1.);
  g_flavfrac1->GetXaxis()->SetNdivisions(505);

  TGraphAsymmErrors *g_flavfrac2 = new TGraphAsymmErrors(h1_num2_rebinned, h1_denom_rebinned,"cp");
//   g_flavfrac2->SetMarkerSize(0.8);
  g_flavfrac2->SetMarkerStyle(25);
  g_flavfrac2->SetMarkerColor(kOrange);
  g_flavfrac2->SetLineColor(kOrange);
  g_flavfrac2->SetLineStyle(3);
  g_flavfrac2->SetLineWidth(2);

  TGraphAsymmErrors *g_flavfrac3 = new TGraphAsymmErrors(h1_num3_rebinned, h1_denom_rebinned,"cp");
//   g_flavfrac3->SetMarkerSize(0.8);
  g_flavfrac3->SetMarkerStyle(26);
  g_flavfrac3->SetMarkerColor(kGreen+1);
  g_flavfrac3->SetLineColor(kGreen+1);
  g_flavfrac3->SetLineStyle(2);
  g_flavfrac3->SetLineWidth(2);

  TGraphAsymmErrors *g_flavfrac4 = new TGraphAsymmErrors(h1_num4_rebinned, h1_denom_rebinned,"cp");
//   g_flavfrac4->SetMarkerSize(0.8);
  g_flavfrac4->SetMarkerStyle(27);
  g_flavfrac4->SetMarkerColor(kBlue);
  g_flavfrac4->SetLineColor(kBlue);
  g_flavfrac4->SetLineStyle(1);
  g_flavfrac4->SetLineWidth(2);

  g_flavfrac1->Draw("AP");
  g_flavfrac2->Draw("P");
  g_flavfrac3->Draw("P");
  g_flavfrac4->Draw("P");

  TLegend *legend = new TLegend(.58,.70,.83,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_flavfrac4, "uds","lp");
  legend->AddEntry(g_flavfrac3, "g","lp");
  legend->AddEntry(g_flavfrac2, "c","lp");
  legend->AddEntry(g_flavfrac1, "b","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.20,0.88, "G#rightarrowq#bar{q}  (q=u,d,s)");
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.53,0.34, "CMS Simulation");
  l1.DrawLatex(0.54,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.53,0.24, "|#eta| < 0.5");
  l1.DrawLatex(0.53,0.19, "Anti-k_{T} R = 0.7 PF Jets");

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void bTagEffVsPt_flav(const string& fInputFile, const string& fPlot, const string& fAlgo,
                      const Int_t nbins, const Double_t *xbins,
                      const string& fXAxisTitle, const string& fYAxisTitle, const string& fOutputFile,
                      const double fXmin = 0, const double fXmax = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_num1 = (TH2D*)file->Get((fPlot + "_B_" + fAlgo).c_str());
  TH2D *h2_num2 = (TH2D*)file->Get((fPlot + "_C_" + fAlgo).c_str());
  TH2D *h2_num3 = (TH2D*)file->Get((fPlot + "_G_" + fAlgo).c_str());
  TH2D *h2_num4 = (TH2D*)file->Get((fPlot + "_UDS_" + fAlgo).c_str());

  TH2D *h2_denom1 = (TH2D*)file->Get((fPlot + "_B_All").c_str());
  TH2D *h2_denom2 = (TH2D*)file->Get((fPlot + "_C_All").c_str());
  TH2D *h2_denom3 = (TH2D*)file->Get((fPlot + "_G_All").c_str());
  TH2D *h2_denom4 = (TH2D*)file->Get((fPlot + "_UDS_All").c_str());

  TH1D *h1_denom1 = h2_denom1->ProjectionX("_px",46,55);
  TH1D *h1_denom2 = h2_denom2->ProjectionX("_px",46,55);
  TH1D *h1_denom3 = h2_denom3->ProjectionX("_px",46,55);
  TH1D *h1_denom4 = h2_denom4->ProjectionX("_px",46,55);
  TH1D *h1_num1 = h2_num1->ProjectionX("_px",46,55);
  TH1D *h1_num2 = h2_num2->ProjectionX("_px",46,55);
  TH1D *h1_num3 = h2_num3->ProjectionX("_px",46,55);
  TH1D *h1_num4 = h2_num4->ProjectionX("_px",46,55);

  TH1D *h1_denom1_rebinned = (TH1D*)h1_denom1->Rebin(nbins,"h1_denom1_rebinned",xbins);
  TH1D *h1_denom2_rebinned = (TH1D*)h1_denom2->Rebin(nbins,"h1_denom2_rebinned",xbins);
  TH1D *h1_denom3_rebinned = (TH1D*)h1_denom3->Rebin(nbins,"h1_denom3_rebinned",xbins);
  TH1D *h1_denom4_rebinned = (TH1D*)h1_denom4->Rebin(nbins,"h1_denom4_rebinned",xbins);
  TH1D *h1_num1_rebinned = (TH1D*)h1_num1->Rebin(nbins,"h1_num1_rebinned",xbins);
  TH1D *h1_num2_rebinned = (TH1D*)h1_num2->Rebin(nbins,"h1_num2_rebinned",xbins);
  TH1D *h1_num3_rebinned = (TH1D*)h1_num3->Rebin(nbins,"h1_num3_rebinned",xbins);
  TH1D *h1_num4_rebinned = (TH1D*)h1_num4->Rebin(nbins,"h1_num4_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_num1_rebinned, h1_denom1_rebinned,"cp");
//   g_efficiency1->SetMarkerSize(0.8);
  g_efficiency1->SetMarkerStyle(24);
  g_efficiency1->SetMarkerColor(kRed);
  g_efficiency1->SetLineColor(kRed);
  g_efficiency1->SetLineStyle(4);
  g_efficiency1->SetLineWidth(2);
  g_efficiency1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(fXmin!=fXmax) g_efficiency1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  g_efficiency1->GetYaxis()->SetRangeUser(0,1.);
  g_efficiency1->GetXaxis()->SetNdivisions(505);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_num2_rebinned, h1_denom2_rebinned,"cp");
//   g_efficiency2->SetMarkerSize(0.8);
  g_efficiency2->SetMarkerStyle(25);
  g_efficiency2->SetMarkerColor(kOrange);
  g_efficiency2->SetLineColor(kOrange);
  g_efficiency2->SetLineStyle(3);
  g_efficiency2->SetLineWidth(2);

  TGraphAsymmErrors *g_efficiency3 = new TGraphAsymmErrors(h1_num3_rebinned, h1_denom3_rebinned,"cp");
//   g_efficiency3->SetMarkerSize(0.8);
  g_efficiency3->SetMarkerStyle(26);
  g_efficiency3->SetMarkerColor(kGreen+1);
  g_efficiency3->SetLineColor(kGreen+1);
  g_efficiency3->SetLineStyle(2);
  g_efficiency3->SetLineWidth(2);

  TGraphAsymmErrors *g_efficiency4 = new TGraphAsymmErrors(h1_num4_rebinned, h1_denom4_rebinned,"cp");
//   g_efficiency4->SetMarkerSize(0.8);
  g_efficiency4->SetMarkerStyle(27);
  g_efficiency4->SetMarkerColor(kBlue);
  g_efficiency4->SetLineColor(kBlue);
  g_efficiency4->SetLineStyle(1);
  g_efficiency4->SetLineWidth(2);

  g_efficiency1->Draw("AP");
  g_efficiency2->Draw("P");
  g_efficiency3->Draw("P");
  g_efficiency4->Draw("P");

  TLegend *legend = new TLegend(.58,.70,.83,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_efficiency1, "b","lp");
  legend->AddEntry(g_efficiency2, "c","lp");
  legend->AddEntry(g_efficiency3, "g","lp");
  legend->AddEntry(g_efficiency4, "uds","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.20,0.88, "G#rightarrowq#bar{q}  (q=u,d,s)");
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.53,0.34, "CMS Simulation");
  l1.DrawLatex(0.54,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.53,0.24, "|#eta| < 0.5");
  l1.DrawLatex(0.53,0.19, "Anti-k_{T} R = 0.7 PF Jets");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.19,0.19, fAlgo.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{

  Double_t PtBins[] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 2000, 3000, 4000, 5000, 6000};


  bTagEffVsPt_comp("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                   "CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching_UDS_G/Final__histograms.root",
                   "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1",
                   18, PtBins, "G#rightarrowb#bar{b}", "G#rightarrowq#bar{q}  (q=u,d,s)", "CSVL",
                   "p_{T} [GeV]", "b-tag Efficiency",
                   "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_RSGToQQbarLight_CSVL.png", 0, 2000);

  FlavFracVsPt("CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching_UDS_G/Final__histograms.root",
               "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", 18, PtBins,
               "p_{T} [GeV]", "Flavor Fraction",
               "FlavFrac_PtJ1_etalt0p5_RSGToQQbarLight.png", 0, 2000);

  bTagEffVsPt_flav("CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching_UDS_G/Final__histograms.root",
                   "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", "CSVL", 18, PtBins,
                   "p_{T} [GeV]", "b-tag Efficiency",
                   "b-tag_eff_PtJ1_etalt0p5_flav_RSGToQQbarLight.png", 0, 2000);

}
