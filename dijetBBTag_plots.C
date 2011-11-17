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

  h_DATA_range->Draw();
  h_MC_range->Draw("same");

  TLegend *legend = new TLegend(.7,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h_DATA_range,"Data","lp");
  legend->AddEntry(h_MC_range,"MC","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.42,"CMS Preliminary");
  l1.DrawLatex(0.15,0.35,"#intLdt = 4.7 fb^{-1}");
  l1.DrawLatex(0.15,0.30,"#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.15,0.25,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.15,0.20,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void overlay_DATA_MC_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                       const Int_t nbins, const Double_t *xbins, const Double_t xMin = 1050, const Double_t xMax = 3500)
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

  h_DijetMass_DATA_rebin->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  h_DijetMass_DATA_rebin->GetYaxis()->SetTitle("Entries");
  h_DijetMass_DATA_rebin->SetMarkerStyle(20);

  h_DijetMass_MC_rebin->SetLineColor(kRed);
  h_DijetMass_MC_rebin->SetMarkerStyle(25);
  h_DijetMass_MC_rebin->SetMarkerColor(kRed);

  h_DijetMass_DATA_rebin->Draw();
  h_DijetMass_MC_rebin->Draw("same");

  TLegend *legend = new TLegend(.7,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h_DijetMass_DATA_rebin,"Data","lp");
  legend->AddEntry(h_DijetMass_MC_rebin,"MC","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.42,"CMS Preliminary");
  l1.DrawLatex(0.15,0.35,"#intLdt = 4.7 fb^{-1}");
  l1.DrawLatex(0.15,0.30,"#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.15,0.25,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.15,0.20,fLabel.c_str());

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

//   TLatex l1;
//   l1.SetTextAlign(12);
// //   l1.SetTextFont(42);
//   l1.SetNDC();
//   l1.SetTextSize(0.04);
//   l1.DrawLatex(0.15,0.42,"CMS Preliminary");
//   l1.DrawLatex(0.15,0.35,"#intLdt = 2.2 fb^{-1}");
//   l1.DrawLatex(0.15,0.30,"#sqrt{s} = 7 TeV");
//   l1.DrawLatex(0.15,0.25,"Anti-k_{T} R = 0.7 PFJets");
//   l1.DrawLatex(0.15,0.20,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void data_MC_ratio_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile,
                     const Int_t nbins, const Double_t *xbins, const Double_t xMin = 1050, const Double_t xMax = 3500)
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

  h_DijetMass_DATA_MC->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  h_DijetMass_DATA_MC->GetYaxis()->SetTitle("Data/MC");
  h_DijetMass_DATA_MC->GetYaxis()->SetRangeUser(0.,2.4);
  
  h_DijetMass_DATA_MC->Divide(h_DijetMass_MC_rebin);
  h_DijetMass_DATA_MC->Draw();
  h_DijetMass_DATA_MC->Fit("pol0");
   
//   TLatex l1;
//   l1.SetTextAlign(12);
// //   l1.SetTextFont(42);
//   l1.SetNDC();
//   l1.SetTextSize(0.04);
//   l1.DrawLatex(0.15,0.82,"CMS Preliminary");
//   l1.DrawLatex(0.15,0.75,"#intLdt = 2.2 fb^{-1}");
//   l1.DrawLatex(0.15,0.70,"#sqrt{s} = 7 TeV");
//   l1.DrawLatex(0.15,0.65,"Anti-k_{T} R = 0.7 PFJets");
  
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

//   Double_t xbins[] = {1050, 1110, 1170, 1230, 1300, 1370, 1440, 1510, 1590, 1670, 1750, 1830, 1920, 2010, 2100, 2200, 2300, 2400, 2500, 2610, 2720, 2830, 2950, 3070, 3190, 3320, 3500};


//   // DATA/MC scale factor  
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "data_MC_scale.png", 26, xbins, 1050, 3500);

//   // Dijet mass
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "DijetMass_data_MC.png", "", 26, xbins, 1050, 3500);

//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_TCHEM_SingleTag.png", "TCHEM SingleTag", 26, xbins, 1050, 3500);

//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_TCHEM_DoubleTag.png", "TCHEM DoubleTag", 26, xbins, 1050, 3500);

//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "data_MC_scale.png", 26, xbins, 1050, 3500);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "data_MC_scale_TCHEM_SingleTag.png", 26, xbins, 1050, 3500);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "data_MC_scale_TCHEM_DoubleTag.png", 26, xbins, 1050, 3500);
//   

//   // Muon multiplicity
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass",
//                     "MuonMultiplicity_vs_DijetMass.png", "", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                     "MuonMultiplicity_vs_DijetMass_TCHEM_SingleTag.png", "TCHEM SingleTag", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                     "MuonMultiplicity_vs_DijetMass_TCHEM_DoubleTag.png", "TCHEM DoubleTag", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass",
//                     "MuonMultiplicity_vs_DijetMass_JetEta1.3.png", "", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                     "MuonMultiplicity_vs_DijetMass_JetEta1.3_TCHEM_SingleTag.png", "TCHEM SingleTag", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                     "DATA__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                     "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                     "MuonMultiplicity_vs_DijetMass_JetEta1.3_TCHEM_DoubleTag.png", "TCHEM DoubleTag", 26, xbins, 1050, 3500);
// 
//   overlay_DATA_MC("CRAB_Jobs_DijetBBTag_MuonMultiplicity/Final__histograms.root",
//                   "DATA__cutHisto_allOtherCuts___________nMuons",
//                   "QCD_Pythia6__cutHisto_allOtherCuts___________nMuons",
//                   "MuonMultiplicity.png", "M_{jj}>1050 GeV", "Muon Multiplicity", "Events", -0.5, 3.5);

  // Muon multiplicity ratios
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass",
//                   "data_MC_ratio_MuonMultiplicity.png", 26, xbins, 1050, 3500);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                   "data_MC_ratio_MuonMultiplicity_TCHEM_SingleTag.png", 26, xbins, 1050, 3500);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                   "data_MC_ratio_MuonMultiplicity_TCHEM_DoubleTag.png", 26, xbins, 1050, 3500);
// 
//   data_MC_ratio("CRAB_Jobs_DijetBBTag_MuonMultiplicity/Final__histograms.root",
//                 "DATA__cutHisto_allOtherCuts___________nMuons",
//                 "QCD_Pythia6__cutHisto_allOtherCuts___________nMuons",
//                 "MuonMultiplicity_ratio.png", "M_{jj}>1050 GeV", "Muon Multiplicity", "DATA/MC", -0.5, 3.5);
    
//
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass",
//                   "data_MC_ratio_MuonMultiplicity_JetEta1.3.png", 26, xbins, 1050, 3500);
//    
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_SingleTag",
//                   "data_MC_ratio_MuonMultiplicity_JetEta1.3_TCHEM_SingleTag.png", 26, xbins, 1050, 3500);
//
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_MuonMultiplicity_JetEta1.3_TCHEM/Final__histograms.root",
//                   "DATA__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                   "QCD_Pythia6__h1_MuonMultiplicity_vs_DijetMass_DoubleTag",
//                   "data_MC_ratio_MuonMultiplicity_JetEta1.3_TCHEM_DoubleTag.png", 26, xbins, 1050, 3500);

// ###################################################################################################################
// ## PU and b-tag SF reweighting
// ###################################################################################################################

//   Double_t xbins[] = {1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037};//, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

  // DATA/MC scale factor
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "plots/data_MC_scale_factor_PUReweighting.png", 41, xbins, 1000, 6000);
  
  // ## TCHEM
  // Dijet mass plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "plots/DijetMass_PUSFReweighting.png", "", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "plots/DijetMass_ratio_PUSFReweighting.png", 41, xbins, 1000, 6000);
// 
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "plots/DijetMass_TCHEM_SingleTag_PUSFReweighting.png", "TCHEM SingleTag", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "plots/DijetMass_ratio_TCHEM_SingleTag_PUSFReweighting.png", 41, xbins, 1000, 6000);
// 
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "plots/DijetMass_TCHEM_DoubleTag_PUSFReweighting.png", "TCHEM DoubleTag", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "plots/DijetMass_ratio_TCHEM_DoubleTag_PUSFReweighting.png", 41, xbins, 1000, 6000);

  // Muon multiplicity plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass_pretag",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                     "plots/nMuons_vs_DijetMass_PUSFReweighting.png", "", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass_pretag",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                   "plots/nMuons_vs_DijetMass_ratio_PUSFReweighting.png", 41, xbins, 1000, 6000);
//     
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "plots/nMuons_vs_DijetMass_TCHEM_SingleTag_PUSFReweighting.png", "", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "plots/nMuons_vs_DijetMass_ratio_TCHEM_SingleTag_PUSFReweighting.png", 41, xbins, 1000, 6000);
// 
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "plots/nMuons_vs_DijetMass_TCHEM_DoubleTag_PUSFReweighting.png", "", 41, xbins, 1000, 6000);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_TCHEM_DoubleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "plots/nMuons_vs_DijetMass_ratio_TCHEM_DoubleTag_PUSFReweighting.png", 41, xbins, 1000, 6000);

  
//   overlay_DATA_MC("CRAB_Jobs_DijetBBTag_MuonMultiplicity/Final__histograms.root",
//                   "DATA__cutHisto_allOtherCuts___________nMuons",
//                   "QCD_Pythia6__cutHisto_allOtherCuts___________nMuons",
//                   "MuonMultiplicity.png", "M_{jj}>1050 GeV", "Muon Multiplicity", "Events", -0.5, 3.5);

  // ## SSVHPT DoubleTag
  // Dijet mass plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "DijetMass_JetPt520_PUReweighting.png", "", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "DijetMass_ratio_JetPt520_PUReweighting.png", 14, xbins, 1000, 2037);
// 
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_SSVHPT_DoubleTag_PUReweighting.png", "SSVHPT DoubleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_SSVHPT_DoubleTag_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_SSVHPT_DoubleTag_PUSFReweighting.png", "SSVHPT DoubleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_DoubleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_SSVHPT_DoubleTag_PUSFReweighting.png", 14, xbins, 1000, 2037);

  // ## TCHEM SingleTag
  // Dijet mass plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "DijetMass_JetPt520_PUReweighting.png", "", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "DijetMass_ratio_JetPt520_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_TCHEM_SingleTag_PUReweighting.png", "TCHEM SingleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_TCHEM_SingleTag_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_TCHEM_SingleTag_PUSFReweighting.png", "TCHEM SingleTag, b-tag SF", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_TCHEM_SingleTag_PUSFReweighting.png", 14, xbins, 1000, 2037);
// 
//   // Muon multiplicity plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass_pretag",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                     "nMuons_vs_DijetMass_JetPt520_PUReweighting.png", "", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass_pretag",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                   "nMuons_vs_DijetMass_ratio_JetPt520_PUReweighting.png", 14, xbins, 1000, 2037);
//   
//   
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_JetPt520_TCHEM_SingleTag_PUReweighting.png", "TCHEM SingleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_JetPt520_TCHEM_SingleTag_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_JetPt520_TCHEM_SingleTag_PUSFReweighting.png", "TCHEM SingleTag, b-tag SF", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_JetPt520_TCHEM_SingleTag_PUSFReweighting.png", 14, xbins, 1000, 2037);


//   // Primary vertex plots
//   overlay_DATA_MC("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices.png", "M_{jj}>1050 GeV, jet p_{T}<520 GeV", "Good Vertex Multiplicity", "Events", -0.5, 20.5);
// 
//   data_MC_ratio("CRAB_Jobs_DijetBBTag_JetPt520_TCHEM_SingleTag_PUSFReweighting/Final__histograms.root",
//                 "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                 "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                 "nGoodVertices_ratio.png", "M_{jj}>1050 GeV, jet p_{T}<520 GeV", "Good Vertex Multiplicity", "DATA/MC", -0.5, 20.5);


  // ## SSVHPT SingleTag
  // Dijet mass plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "DijetMass_JetPt520_PUReweighting.png", "", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "DijetMass_ratio_JetPt520_PUReweighting.png", 14, xbins, 1000, 2037);


//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_SSVHPT_SingleTag_PUReweighting.png", "SSVHPT SingleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_SSVHPT_SingleTag_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_JetPt520_SSVHPT_SingleTag_PUSFReweighting.png", "SSVHPT SingleTag, b-tag SF", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_JetPt520_SSVHPT_SingleTag_PUSFReweighting.png", 14, xbins, 1000, 2037);
// 
//   // Muon multiplicity plots
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_JetPt520_SSVHPT_SingleTag_PUReweighting.png", "SSVHPT SingleTag", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_JetPt520_SSVHPT_SingleTag_PUReweighting.png", 14, xbins, 1000, 2037);
// 
// 
//   overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUSFReweighting/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_JetPt520_SSVHPT_SingleTag_PUSFReweighting.png", "SSVHPT SingleTag, b-tag SF", 14, xbins, 1000, 2037);
// 
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_JetPt520_SSVHPT_SingleTag_PUSFReweighting/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_JetPt520_SSVHPT_SingleTag_PUSFReweighting.png", 14, xbins, 1000, 2037);


// ## Full 2011 dataset (4.679 /fb)

  Double_t xbins[] = {944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

// Low pile-up
  
  // DATA/MC scale factor
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "data_MC_scale_factor_PUReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);

  // Primary vertex plots
//   overlay_DATA_MC("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices_PUReweighted_Full2011_LowPU.png", "M_{jj}>944 GeV", "Good Vertex Multiplicity", "Events", -0.5, 20.5);
// 
//   data_MC_ratio("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                 "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                 "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                 "nGoodVertices_ratio_PUReweighted_Full2011_LowPU.png", "M_{jj}>944 GeV", "Good Vertex Multiplicity", "DATA/MC", -0.5, 20.5);
  
  // Dijet mass plots
  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "DijetMass_PUReweighted_Full2011_LowPU.png", "", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "DijetMass_ratio_PUReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);

  
  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_SSVHPT_SingleTag_PUReweighted_Full2011_LowPU.png", "SSVHPT SingleTag", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_SSVHPT_SingleTag_PUReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);


  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_SSVHPT_SingleTag_PUSFReweighted_Full2011_LowPU.png", "SSVHPT SingleTag, b-tag SF", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_SSVHPT_SingleTag_PUSFReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);

  // Muon multiplicity plots
  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass_pretag",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
                    "nMuons_vs_DijetMass_PUReweighted_Full2011_LowPU.png", "", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass_pretag",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
                  "nMuons_vs_DijetMass_ratio_PUReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);

  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                    "nMuons_vs_DijetMass_SSVHPT_SingleTag_PUReweighted_Full2011_LowPU.png", "SSVHPT SingleTag", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                  "nMuons_vs_DijetMass_ratio_SSVHPT_SingleTag_PUReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);


  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                    "nMuons_vs_DijetMass_SSVHPT_SingleTag_PUSFReweighted_Full2011_LowPU.png", "SSVHPT SingleTag, b-tag SF", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_LowPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                  "nMuons_vs_DijetMass_ratio_SSVHPT_SingleTag_PUSFReweighted_Full2011_LowPU.png", 42, xbins, 944, 6000);

// High pile-up

  // DATA/MC scale factor
//   data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "data_MC_scale_factor_PUReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);
  
  // Primary vertex plots
//   overlay_DATA_MC("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices_PUReweighted_Full2011_HighPU.png", "M_{jj}>944 GeV", "Good Vertex Multiplicity", "Events", -0.5, 20.5);
// 
//   data_MC_ratio("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
//                 "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                 "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                 "nGoodVertices_ratio_PUReweighted_Full2011_HighPU.png", "M_{jj}>944 GeV", "Good Vertex Multiplicity", "DATA/MC", -0.5, 20.5);

  // Dijet mass plots
  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "DijetMass_PUReweighted_Full2011_HighPU.png", "", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "DijetMass_ratio_PUReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);


  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_SSVHPT_SingleTag_PUReweighted_Full2011_HighPU.png", "SSVHPT SingleTag", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_SSVHPT_SingleTag_PUReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);


  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_SSVHPT_SingleTag_PUSFReweighted_Full2011_HighPU.png", "SSVHPT SingleTag, b-tag SF", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_SSVHPT_SingleTag_PUSFReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);

  // Muon multiplicity plots
  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass_pretag",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
                    "nMuons_vs_DijetMass_PUReweighted_Full2011_HighPU.png", "", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass_pretag",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
                  "nMuons_vs_DijetMass_ratio_PUReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);

  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                    "nMuons_vs_DijetMass_SSVHPT_SingleTag_PUReweighted_Full2011_HighPU.png", "SSVHPT SingleTag", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                  "nMuons_vs_DijetMass_ratio_SSVHPT_SingleTag_PUReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);


  overlay_DATA_MC_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                    "DATA__h1_nMuons_vs_DijetMass",
                    "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                    "nMuons_vs_DijetMass_SSVHPT_SingleTag_PUSFReweighted_Full2011_HighPU.png", "SSVHPT SingleTag, b-tag SF", 42, xbins, 944, 6000);

  data_MC_ratio_r("CRAB_Jobs_DijetBBTag_HighPU_SSVHPT_SingleTag_PUSFReweighted/Final__histograms.root",
                  "DATA__h1_nMuons_vs_DijetMass",
                  "QCD_Pythia6__h1_nMuons_vs_DijetMass",
                  "nMuons_vs_DijetMass_ratio_SSVHPT_SingleTag_PUSFReweighted_Full2011_HighPU.png", 42, xbins, 944, 6000);
  
}