#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "tdrstyle.C"


using namespace std;

void triggerEfficiency(const string& fInputFile, const string& hNum, const string& hDenom,
                       const string& fTitle, const string& fOutputFile, const string& fLabel = "",
                       const Int_t rebin = 1, const Double_t xMin = 800, const Double_t xMax = 1300)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_num = (TH1D*)file->Get(hNum.c_str());
  TH1D *h_DijetMass_denom = (TH1D*)file->Get(hDenom.c_str());

  h_DijetMass_num->Rebin(rebin);
  h_DijetMass_denom->Rebin(rebin);

  Double_t binWidth = h_DijetMass_num->GetBinWidth(h_DijetMass_num->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int(xMin/binWidth);

  TH1D *h_Npass = new TH1D("h_Npass","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_Npass->SetBinContent(i,h_DijetMass_num->GetBinContent(i+nBinsToSkip));
    h_Npass->SetBinError(i,h_DijetMass_num->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_Ntotal = new TH1D("h_Ntotal","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_Ntotal->SetBinContent(i,h_DijetMass_denom->GetBinContent(i+nBinsToSkip));
    h_Ntotal->SetBinError(i,h_DijetMass_denom->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h_Npass, h_Ntotal,"cp");
//   g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(20);
  g_efficiency->SetMarkerColor(kRed);
  g_efficiency->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  g_efficiency->GetYaxis()->SetTitle("Trigger Efficiency");
  g_efficiency->GetXaxis()->SetRangeUser(xMin,xMax);
  g_efficiency->GetYaxis()->SetRangeUser(0.,1.05);

  g_efficiency->Draw("AP");

  TLine *line = new TLine(xMin,0.99,xMax,0.99);
  line->SetLineWidth(3.);
  line->SetLineColor(kGreen+2);
  line->SetLineStyle(2);
  line->Draw("same");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.02,0.96, fTitle.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.55,0.34, "CMS Preliminary");
  l1.DrawLatex(0.56,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.55,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.55,0.19, ("Wide Jets" + fLabel).c_str());
    
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete line;
  delete h_Npass;
  delete h_Ntotal;
  delete g_efficiency;
  delete c;
  delete file;
}


void makePlots()
{
  triggerEfficiency("CRAB_Jobs_TriggerEfficiency_CSVL_WideJets/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT650",
                    "DATA__h1_DijetMass_denom_HT650",
                    "(HLT_HT650 AND HLT_HT400)/HLT_HT400", "HLT_HT650_efficiency.png", "", 25, 650, 1150);

  triggerEfficiency("CRAB_Jobs_TriggerEfficiency_CSVL_WideJets/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT650_0tag",
                    "DATA__h1_DijetMass_denom_HT650_0tag",
                    "(HLT_HT650 AND HLT_HT400)/HLT_HT400", "HLT_HT650_efficiency_CSVL_0Tag.png", ", CSVL 0-tag", 25, 650, 1150);

  triggerEfficiency("CRAB_Jobs_TriggerEfficiency_CSVL_WideJets/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT650_1tag",
                    "DATA__h1_DijetMass_denom_HT650_1tag",
                    "(HLT_HT650 AND HLT_HT400)/HLT_HT400", "HLT_HT650_efficiency_CSVL_1Tag.png", ", CSVL 1-tag", 25, 650, 1150);

  triggerEfficiency("CRAB_Jobs_TriggerEfficiency_CSVL_WideJets/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT650_2tag",
                    "DATA__h1_DijetMass_denom_HT650_2tag",
                    "(HLT_HT650 AND HLT_HT400)/HLT_HT400", "HLT_HT650_efficiency_CSVL_2Tag.png", ", CSVL 2-tag", 25, 650, 1150);
}
