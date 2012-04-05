#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;


void overlay_DATA_MC(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax, const string& fLogy = "", const Int_t fRebin = 1,
                     const Double_t yMin = 0, const Double_t yMax = 0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH1D *h_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file->Get(hMC.c_str());

  if(fRebin>1)
  {
    h_DATA->Rebin(fRebin);
    h_MC->Rebin(fRebin);
  }

  Double_t binWidth = h_DATA->GetBinWidth(h_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DATA_range = new TH1D("h_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DATA_range->SetBinContent(i,h_DATA->GetBinContent(i+nBinsToSkip));
    h_DATA_range->SetBinError(i,h_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_MC_range = new TH1D("h_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC_range->SetBinContent(i,h_MC->GetBinContent(i+nBinsToSkip));
    h_MC_range->SetBinError(i,h_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h_DATA_range->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DATA_range->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(yMin!=yMax) h_DATA_range->GetYaxis()->SetRangeUser(yMin,yMax);
  h_DATA_range->SetMarkerStyle(20);

  h_MC_range->SetLineColor(42);
  h_MC_range->SetFillColor(42);
  h_MC_range->SetMarkerStyle(25);
//   h_MC_range->SetMarkerColor(kRed);

  h_DATA_range->Draw();
  h_MC_range->Draw("histsame");
  h_DATA_range->Draw("same");

  gPad->RedrawAxis();
  
  TLegend *legend = new TLegend(.7,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h_DATA_range,"Data","lp");
  legend->AddEntry(h_MC_range,"MC","f");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.25,0.42,"CMS Preliminary");
  l1.DrawLatex(0.25,0.34,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.25,0.29,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.25,0.24,"Wide Jets");
  else l1.DrawLatex(0.25,0.24,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.25,0.19,fLabel.c_str());

  if(fLogy=="Logy") c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void overlay_DATA_MC_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                       const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle,
                       const Double_t xMin = 1050, const Double_t xMax = 3500, const Double_t yMin = 0, const Double_t yMax = 0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file->Get(hMC.c_str());

  Double_t binWidth = h_DijetMass_DATA->GetBinWidth(h_DijetMass_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DijetMass_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DijetMass_DATA_range = new TH1D("h_DijetMass_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_DATA_range->SetBinContent(i,h_DijetMass_DATA->GetBinContent(i+nBinsToSkip));
    h_DijetMass_DATA_range->SetBinError(i,h_DijetMass_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_MC_range = new TH1D("h_DijetMass_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_MC_range->SetBinContent(i,h_DijetMass_MC->GetBinContent(i+nBinsToSkip));
    h_DijetMass_MC_range->SetBinError(i,h_DijetMass_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH1D *h_DijetMass_DATA_rebin = (TH1D*)h_DijetMass_DATA_range->Rebin(nbins,"h_DijetMass_DATA",xbins);
  TH1D *h_DijetMass_MC_rebin = (TH1D*)h_DijetMass_MC_range->Rebin(nbins,"h_DijetMass_MC",xbins);

  h_DijetMass_DATA_rebin->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DijetMass_DATA_rebin->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(yMin!=yMax) h_DijetMass_DATA_rebin->GetYaxis()->SetRangeUser(yMin,yMax);
  
  h_DijetMass_DATA_rebin->SetMarkerStyle(20);

  h_DijetMass_MC_rebin->SetLineColor(42);
  h_DijetMass_MC_rebin->SetFillColor(42);
//   h_DijetMass_MC_rebin->SetMarkerStyle(25);
//   h_DijetMass_MC_rebin->SetMarkerColor(kRed);
  
  h_DijetMass_DATA_rebin->Draw();
  h_DijetMass_MC_rebin->Draw("histsame");
  h_DijetMass_DATA_rebin->Draw("same");

  gPad->RedrawAxis();
  
  TLegend *legend = new TLegend(.7,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h_DijetMass_DATA_rebin,"Data","lp");
  legend->AddEntry(h_DijetMass_MC_rebin,"MC","f");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.3,0.82,"CMS Preliminary");
  l1.DrawLatex(0.3,0.74,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.3,0.69,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.3,0.64,"Wide Jets");
  else l1.DrawLatex(0.3,0.64,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.3,0.59,fLabel.c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_rebin;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file;
}


void data_MC_ratio(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                   const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH1D *h_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file->Get(hMC.c_str());

  Double_t binWidth = h_DATA->GetBinWidth(h_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DATA_range = new TH1D("h_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DATA_range->SetBinContent(i,h_DATA->GetBinContent(i+nBinsToSkip));
    h_DATA_range->SetBinError(i,h_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_MC_range = new TH1D("h_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC_range->SetBinContent(i,h_MC->GetBinContent(i+nBinsToSkip));
    h_MC_range->SetBinError(i,h_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h_DATA_range->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DATA_range->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h_DATA_range->SetMarkerStyle(20);

  h_MC_range->SetLineColor(kRed);
  h_MC_range->SetMarkerStyle(25);
  h_MC_range->SetMarkerColor(kRed);

  h_DATA_range->Divide(h_MC_range);
  h_DATA_range->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.42,"CMS Preliminary");
  l1.DrawLatex(0.15,0.34,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.15,0.29,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.15,0.24,"Wide Jets");
  else l1.DrawLatex(0.15,0.24,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.15,0.19,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void data_MC_ratio_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin = 1050, const Double_t xMax = 3500)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file->Get(hMC.c_str());

  Double_t binWidth = h_DijetMass_DATA->GetBinWidth(h_DijetMass_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int(xMin/binWidth);
  
  TH1D *h_DijetMass_DATA_range = new TH1D("h_DijetMass_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_DATA_range->SetBinContent(i,h_DijetMass_DATA->GetBinContent(i+nBinsToSkip));
    h_DijetMass_DATA_range->SetBinError(i,h_DijetMass_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_MC_range = new TH1D("h_DijetMass_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_MC_range->SetBinContent(i,h_DijetMass_MC->GetBinContent(i+nBinsToSkip));
    h_DijetMass_MC_range->SetBinError(i,h_DijetMass_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH1D *h_DijetMass_DATA_MC = (TH1D*)h_DijetMass_DATA_range->Rebin(nbins,"h_DijetMass_DATA_MC",xbins);
  TH1D *h_DijetMass_MC_rebin = (TH1D*)h_DijetMass_MC_range->Rebin(nbins,"h_DijetMass_MC",xbins);

  h_DijetMass_DATA_MC->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DijetMass_DATA_MC->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h_DijetMass_DATA_MC->GetYaxis()->SetRangeUser(0.,2.4);
  
  h_DijetMass_DATA_MC->Divide(h_DijetMass_MC_rebin);
  h_DijetMass_DATA_MC->Draw();
  h_DijetMass_DATA_MC->Fit("pol0");
   
  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.82,"CMS Preliminary");
  l1.DrawLatex(0.15,0.74,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.15,0.69,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.15,0.64,"Wide Jets");
  else l1.DrawLatex(0.15,0.64,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.15,0.59,fLabel.c_str());
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_MC;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file;
}


void makePlots()
{

//###################################################################################################################
//##
//## Full 2011 dataset (4.976 /fb)
//##
//###################################################################################################################

  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231,
                      2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,
                      5058, 5253, 5455, 5663, 5877, 6000};

//###################################################################################################################
//## PU reweighting applied
//###################################################################################################################

//   // DATA/MC scale factor
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_EventBins/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "data_MC_scale_factor_PUReweighted_Full2011.png", "", 43, xbins, "D