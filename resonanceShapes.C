#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

void ResonanceShape(const string& fInputFile, const string& fPlot, const Int_t fRebin,
                    const string& fTitle, const string& fLabel, const string& fXAxisTitle, const string& fYAxisTitle,
                    const Double_t fXmin, const Double_t fXmax, const Double_t fFitXmin, const Double_t fFitXmax,
                    const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot = (TH1D*)file->Get(fPlot.c_str());

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot->Rebin(fRebin);
  h1_plot->SetTitle(fTitle.c_str());
  h1_plot->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot->GetXaxis()->SetRangeUser(fXmin,fXmax);
//   h1_plot->GetYaxis()->SetRangeUser(0,1.);

  TF1 *f1 = new TF1("f1","gaus",fFitXmin,fFitXmax);
  h1_plot->Fit("f1","R");
  
  h1_plot->Draw();
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}

void ResonanceShape_var_bins(const string& fInputFile, const string& fPlot1, const string& fPlot2, const Int_t fNbins, const Double_t *fXbins,
                             const string& fTitle, const string& fLegend1, const string& fLegend2, const string& fXAxisTitle, const string& fYAxisTitle,
                             const Double_t fXmin, const Double_t fXmax,
                             const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot1_ = (TH1D*)file->Get(fPlot1.c_str());
  h1_plot1_->Scale(1./h1_plot1_->Integral());
  TH1D *h1_plot1 = (TH1D*)h1_plot1_->Rebin(fNbins,"h1_DijetMass1",fXbins);

  TH1D *h1_plot2_ = (TH1D*)file->Get(fPlot2.c_str());
  h1_plot2_->Scale(1./h1_plot2_->Integral());
  TH1D *h1_plot2 = (TH1D*)h1_plot2_->Rebin(fNbins,"h1_DijetMass2",fXbins);
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot1->SetTitle(fTitle.c_str());
  h1_plot1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot1->GetXaxis()->SetRangeUser(fXmin,fXmax);

  h1_plot1->SetLineColor(kBlue);
  h1_plot1->SetLineWidth(2);

  h1_plot2->SetLineColor(kRed);
  h1_plot2->SetLineStyle(2);
  h1_plot2->SetLineWidth(2);
  
  h1_plot1->Draw("hist");
  h1_plot1->Draw("same");

  h1_plot2->Draw("histsame");
  h1_plot2->Draw("same");

  TLegend *legend = new TLegend(.6,.6,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h1_plot1,fLegend1.c_str(),"l");
  legend->AddEntry(h1_plot2,fLegend2.c_str(),"l");
  legend->Draw();

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void ResonanceShape_reshaped_var_bins(const string& fInputFile, const string& fPlot1, const string& fPlot2, const Int_t fNbins, const Double_t *fXbins,
                    const Double_t fOldMean, const Double_t fOldSigma, const Double_t fNewMean, const Double_t fNewSigma,
                    const string& fTitle, const string& fLegend1, const string& fLegend2, const string& fXAxisTitle, const string& fYAxisTitle,
                    const Double_t fXmin, const Double_t fXmax,
                    const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot1_old = (TH1D*)file->Get(fPlot1.c_str());
  h1_plot1_old->Scale(1./h1_plot1_old->Integral());

  TH1D *h1_plot1_new = new TH1D("h1_plot1_new","h1_plot1_new",h1_plot1_old->GetNbinsX(),h1_plot1_old->GetXaxis()->GetXmin(),h1_plot1_old->GetXaxis()->GetXmax());

  for(Int_t i=1; i<h1_plot1_old->GetNbinsX(); i++)
  {
    Double_t newX = fOldMean + (fNewSigma/fOldSigma)*(h1_plot1_old->GetBinCenter(i)-fOldMean) + (fNewMean-fOldMean);
    h1_plot1_new->Fill(newX,h1_plot1_old->GetBinContent(i));
  }

//   TF1 *f1 = new TF1("f1","gaus",3250,3600);
//   h1_plot1_new->Fit("f1","R");
  
  TH1D *h1_plot1 = (TH1D*)h1_plot1_new->Rebin(fNbins,"h1_DijetMass1",fXbins);

  TH1D *h1_plot2_ = (TH1D*)file->Get(fPlot2.c_str());
  h1_plot2_->Scale(1./h1_plot2_->Integral());
  TH1D *h1_plot2 = (TH1D*)h1_plot2_->Rebin(fNbins,"h1_DijetMass2",fXbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot1->SetTitle(fTitle.c_str());
  h1_plot1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot1->GetXaxis()->SetRangeUser(fXmin,fXmax);

  h1_plot1->SetLineColor(kBlue);
  h1_plot1->SetLineWidth(2);

  h1_plot2->SetLineColor(kRed);
  h1_plot2->SetLineStyle(2);
  h1_plot2->SetLineWidth(2);

  h1_plot1->Draw("hist");
  h1_plot1->Draw("same");

  h1_plot2->Draw("histsame");
  h1_plot2->Draw("same");

  TLegend *legend = new TLegend(.6,.6,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h1_plot1,fLegend1.c_str(),"l");
  legend->AddEntry(h1_plot2,fLegend2.c_str(),"l");
  legend->Draw();

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{

//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=500 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 450, 550,
//                  "ResonanceShape_RSGToQQbar_M-500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=500 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 430, 550,
//                  "ResonanceShape_RSGToBBbar_M-500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=700 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 625, 750,
//                  "ResonanceShape_RSGToQQbar_M-700.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=700 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 600, 740,
//                  "ResonanceShape_RSGToBBbar_M-700.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=1200 GeV", "" , "Dijet Mass [GeV]", "", 0, 2500, 1100, 1300,
//                  "ResonanceShape_RSGToQQbar_M-1200.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=1200 GeV", "" , "Dijet Mass [GeV]", "", 0, 2500, 1070, 1260,
//                  "ResonanceShape_RSGToBBbar_M-1200.png");
// 
//  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=2000 GeV", "" , "Dijet Mass [GeV]", "", 0, 5000, 1850, 2100,
//                  "ResonanceShape_RSGToQQbar_M-2000.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=2000 GeV", "" , "Dijet Mass [GeV]", "", 0, 5000, 1800, 2070,
//                  "ResonanceShape_RSGToBBbar_M-2000.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=3500 GeV", "" , "Dijet Mass [GeV]", "", 0, 6000, 3300, 3600,
//                  "ResonanceShape_RSGToQQbar_M-3500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=3500 GeV", "" , "Dijet Mass [GeV]", "", 0, 6000, 3250, 3600,
//                  "ResonanceShape_RSGToBBbar_M-3500.png");

  // ##########################################
  // ## With variable dijet mass binning
  // ##########################################
  
  Double_t xbins[] = {0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693,740, 788,
                      838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147,
                      3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1530,
//                           "ResonanceShape_RSG_M-500_VarBins.png");
//   
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=700 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1856,
//                           "ResonanceShape_RSG_M-700_VarBins.png");
//   
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=1200 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 2775,
//                           "ResonanceShape_RSG_M-1200_VarBins.png");
// 
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=2000 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 4509,
//                           "ResonanceShape_RSG_M-2000_VarBins.png");
// 
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=3500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 5663,
//                           "ResonanceShape_RSG_M-3500_VarBins.png");

  // ##########################################
  // ## With variable dijet mass binning
  // ##########################################

  ResonanceShape_reshaped_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
                          3429., 143.9, 3383., 171.5,
                          "M=3500 GeV", "#splitline{RSG#rightarrowq#bar{q}}{reshaped}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 5663,
                          "ResonanceShape_RSG_M-3500_reshaped_VarBins.png");
}
