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

  TH2D *h2_denom2 = (TH2D*)file2->Get((fPlot2 + "_All").c_str());
  TH2D *h2_num2 = (TH2D*)file2->Get((fPlot2 + "_" + fAlgo).c_str());

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


void bTagEffVsPt_comp2(const string& fInputFile1, const string& fInputFile2,
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
  h2_denom2->Add((TH2D*)file2->Get((fPlot2 + "_UDSG_All").c_str()));
  TH2D *h2_num2 = (TH2D*)file2->Get((fPlot2 + "_B_" + fAlgo).c_str());
  h2_num2->Add((TH2D*)file2->Get((fPlot2 + "_C_" + fAlgo).c_str()));
  h2_num2->Add((TH2D*)file2->Get((fPlot2 + "_UDSG_" + fAlgo).c_str()));

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


void bTagEffVsEta(const string& fInputFile, const string& fPlot,
                  const Int_t nbins, const Double_t *xbins,
                  const string& fAlgo, const string& fXAxisTitle, const string& fYAxisTitle,
                  const string& fLabel1, const string& fLabel2, const double fPtMin, const double fPtMax,
                  const string& fOutputFile)
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

  TH2D *h2_denom = (TH2D*)file->Get((fPlot + "_All").c_str());
  TH2D *h2_num = (TH2D*)file->Get((fPlot + "_" + fAlgo).c_str());

  Int_t bin1 = h2_denom->GetXaxis()->FindBin(fPtMin+5);
  Int_t bin2 = h2_denom->GetXaxis()->FindBin(fPtMax-5);
  
  TH1D *h1_denom = h2_denom->ProjectionY("_py",bin1,bin2);
  TH1D *h1_num = h2_num->ProjectionY("_py",bin1,bin2);

  TH1D *h1_denom_rebinned = (TH1D*)h1_denom->Rebin(nbins,"h1_denom_rebinned",xbins);
  TH1D *h1_num_rebinned = (TH1D*)h1_num->Rebin(nbins,"h1_num_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *h2_bkg = new TH2D("h2_bkg","",100,-2.5,2.5,100,0,1);
  h2_bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h2_bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h2_bkg->GetXaxis()->SetNdivisions(505);
  
  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_num_rebinned, h1_denom_rebinned,"cp");
//   g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(20);
  g_efficiency->SetMarkerColor(kBlack);

  h2_bkg->Draw();
  g_efficiency->Draw("P");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.19,0.86, fLabel1.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.53,0.34, "CMS Simulation");
  l1.DrawLatex(0.54,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.53,0.24, fLabel2.c_str());
  l1.DrawLatex(0.53,0.19, "Anti-k_{T} R = 0.7 PF Jets");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.19,0.19, fAlgo.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void bTagEffVsEta_2(const string& fInputFile, const string& fPlot,
                    const Int_t nbins, const Double_t *xbins,
                    const string& fAlgo, const string& fXAxisTitle, const string& fYAxisTitle,
                    const string& fLabel1, const string& fLabel2, const double fPtMin, const double fPtMax,
                    const string& fOutputFile)
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

  TH2D *h2_denom = (TH2D*)file->Get((fPlot + "_B_All").c_str());
  h2_denom->Add((TH2D*)file->Get((fPlot + "_C_All").c_str()));
  h2_denom->Add((TH2D*)file->Get((fPlot + "_UDSG_All").c_str()));
  TH2D *h2_num = (TH2D*)file->Get((fPlot + "_B_" + fAlgo).c_str());
  h2_num->Add((TH2D*)file->Get((fPlot + "_C_" + fAlgo).c_str()));
  h2_num->Add((TH2D*)file->Get((fPlot + "_UDSG_" + fAlgo).c_str()));
  

  Int_t bin1 = h2_denom->GetXaxis()->FindBin(fPtMin+5);
  Int_t bin2 = h2_denom->GetXaxis()->FindBin(fPtMax-5);

  TH1D *h1_denom = h2_denom->ProjectionY("_py",bin1,bin2);
  TH1D *h1_num = h2_num->ProjectionY("_py",bin1,bin2);

  TH1D *h1_denom_rebinned = (TH1D*)h1_denom->Rebin(nbins,"h1_denom_rebinned",xbins);
  TH1D *h1_num_rebinned = (TH1D*)h1_num->Rebin(nbins,"h1_num_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *h2_bkg = new TH2D("h2_bkg","",100,-2.5,2.5,100,0,1);
  h2_bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h2_bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h2_bkg->GetXaxis()->SetNdivisions(505);

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_num_rebinned, h1_denom_rebinned,"cp");
//   g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(20);
  g_efficiency->SetMarkerColor(kBlack);

  h2_bkg->Draw();
  g_efficiency->Draw("P");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.19,0.86, fLabel1.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.80, "CMS Simulation");
  l1.DrawLatex(0.20,0.75, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.70, fLabel2.c_str());
  l1.DrawLatex(0.19,0.65, "Anti-k_{T} R = 0.7 PF Jets");
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

  // CSVL
  bTagEffVsPt_comp("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                   "CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_AK5Jets_PartonMatching/Final__histograms.root",
                   "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B",
                   18, PtBins, "G#rightarrowb#bar{b}  (R=0.7)", "G#rightarrowb#bar{b}  (R=0.5)", "CSVL",
                   "p_{T} [GeV]", "b-tag Efficiency",
                   "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_AK5Jets_CSVL.png", 0, 2000);

  bTagEffVsPt_comp("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                   "CRAB_Jobs_bTaggingEfficiency_QCD_PartonMatching/Final__histograms.root",
                   "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "QCD_Pythia6__h2_EtaJ1_vs_PtJ1_B",
                   18, PtBins, "G#rightarrowb#bar{b}", "QCD", "CSVL",
                   "p_{T} [GeV]", "b-tag Efficiency",
                   "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_QCD_CSVL.png", 0, 2000);

  bTagEffVsPt_comp("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                   "CRAB_Jobs_bTaggingEfficiency_QCD_FCR_PartonMatching/Final__histograms.root",
                   "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "QCD_Pythia6__h2_EtaJ1_vs_PtJ1_B",
                   18, PtBins, "G#rightarrowb#bar{b}", "QCD FCR", "CSVL",
                   "p_{T} [GeV]", "b-tag Efficiency",
                   "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_QCD_FCR_CSVL.png", 0, 2000);

  bTagEffVsPt_comp2("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                    "CRAB_Jobs_bTaggingEfficiency_RSGToGG_PartonMatching/Final__histograms.root",
                    "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1",
                    18, PtBins, "G#rightarrowb#bar{b}", "G#rightarrowgg", "CSVL",
                    "p_{T} [GeV]", "b-tag Efficiency",
                    "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_RSGToGG_CSVL.png", 0, 2000);

  bTagEffVsPt_comp2("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                    "CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching/Final__histograms.root",
                    "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1",
                    18, PtBins, "G#rightarrowb#bar{b}", "G#rightarrowq#bar{q}  (q=u,d,s)", "CSVL",
                    "p_{T} [GeV]", "b-tag Efficiency",
                    "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_RSGToQQbarLight_CSVL.png", 0, 2000);

  // TCHEL
  bTagEffVsPt_comp2("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                    "CRAB_Jobs_bTaggingEfficiency_RSGToGG_PartonMatching/Final__histograms.root",
                    "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1",
                    18, PtBins, "G#rightarrowb#bar{b}", "G#rightarrowgg", "TCHEL",
                    "p_{T} [GeV]", "b-tag Efficiency",
                    "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_RSGToGG_TCHEL.png", 0, 2000);

  bTagEffVsPt_comp2("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
                    "CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching/Final__histograms.root",
                    "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", "RSGravitonJJ__h2_EtaJ1_vs_PtJ1",
                    18, PtBins, "G#rightarrowb#bar{b}", "G#rightarrowq#bar{q}  (q=u,d,s)", "TCHEL",
                    "p_{T} [GeV]", "b-tag Efficiency",
                    "b-tag_eff_PtJ1_etalt0p5_RSGToBBbar_RSGToQQbarLight_TCHEL.png", 0, 2000);

  
  Double_t EtaBins[] = {-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5};

  // CSVL
  bTagEffVsEta("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
               "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", 10, EtaBins,
               "CSVL", "#eta", "b-tag Efficiency",
               "G#rightarrowb#bar{b}", "100 < p_{T} < 300 GeV", 100, 300,
               "b-tag_eff_EtaJ1_Pt100To300_RSGToBBbar_CSVL.png");

  bTagEffVsEta("CRAB_Jobs_bTaggingEfficiency_RSGToBBbar_PartonMatching/Final__histograms.root",
               "RSGravitonJJ__h2_EtaJ1_vs_PtJ1_B", 10, EtaBins,
               "CSVL", "#eta", "b-tag Efficiency",
               "G#rightarrowb#bar{b}", "800 < p_{T} < 1000 GeV", 800, 1000,
               "b-tag_eff_EtaJ1_Pt800To1000_RSGToBBbar_CSVL.png");

  bTagEffVsEta_2("CRAB_Jobs_bTaggingEfficiency_RSGToGG_PartonMatching/Final__histograms.root",
                 "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", 10, EtaBins,
                 "CSVL", "#eta", "b-tag Efficiency",
                 "G#rightarrowgg", "100 < p_{T} < 300 GeV", 100, 300,
                 "b-tag_eff_EtaJ1_Pt100To300_RSGToGG_CSVL.png");

  bTagEffVsEta_2("CRAB_Jobs_bTaggingEfficiency_RSGToGG_PartonMatching/Final__histograms.root",
                 "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", 10, EtaBins,
                 "CSVL", "#eta", "b-tag Efficiency",
                 "G#rightarrowgg", "800 < p_{T} < 1000 GeV", 800, 1000,
                 "b-tag_eff_EtaJ1_Pt800To1000_RSGToGG_CSVL.png");

  bTagEffVsEta_2("CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching/Final__histograms.root",
                 "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", 10, EtaBins,
                 "CSVL", "#eta", "b-tag Efficiency",
                 "G#rightarrowq#bar{q}  (q=u,d,s)", "100 < p_{T} < 300 GeV", 100, 300,
                 "b-tag_eff_EtaJ1_Pt100To300_RSGToQQbarLight_CSVL.png");

  bTagEffVsEta_2("CRAB_Jobs_bTaggingEfficiency_RSGToQQbarLight_PartonMatching/Final__histograms.root",
                 "RSGravitonJJ__h2_EtaJ1_vs_PtJ1", 10, EtaBins,
                 "CSVL", "#eta", "b-tag Efficiency",
                 "G#rightarrowq#bar{q}  (q=u,d,s)", "800 < p_{T} < 1000 GeV", 800, 1000,
                 "b-tag_eff_EtaJ1_Pt800To1000_RSGToQQbarLight_CSVL.png");
    
}
