#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"

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


void makePlots()
{

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
                 "RSG#rightarrowq#bar{q}, M=500 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 450, 550,
                 "ResonanceShape_RSGToQQbar_M-500.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
                 "RSG#rightarrowb#bar{b}, M=500 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 430, 550,
                 "ResonanceShape_RSGToBBbar_M-500.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
                 "RSG#rightarrowq#bar{q}, M=700 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 625, 750,
                 "ResonanceShape_RSGToQQbar_M-700.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
                 "RSG#rightarrowb#bar{b}, M=700 GeV", "" , "Dijet Mass [GeV]", "", 0, 2000, 600, 740,
                 "ResonanceShape_RSGToBBbar_M-700.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
                 "RSG#rightarrowq#bar{q}, M=1200 GeV", "" , "Dijet Mass [GeV]", "", 0, 2500, 1100, 1300,
                 "ResonanceShape_RSGToQQbar_M-1200.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
                 "RSG#rightarrowb#bar{b}, M=1200 GeV", "" , "Dijet Mass [GeV]", "", 0, 2500, 1070, 1260,
                 "ResonanceShape_RSGToBBbar_M-1200.png");

 ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
                 "RSG#rightarrowq#bar{q}, M=2000 GeV", "" , "Dijet Mass [GeV]", "", 0, 5000, 1850, 2100,
                 "ResonanceShape_RSGToQQbar_M-2000.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
                 "RSG#rightarrowb#bar{b}, M=2000 GeV", "" , "Dijet Mass [GeV]", "", 0, 5000, 1800, 2070,
                 "ResonanceShape_RSGToBBbar_M-2000.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
                 "RSG#rightarrowq#bar{q}, M=3500 GeV", "" , "Dijet Mass [GeV]", "", 0, 6000, 3300, 3600,
                 "ResonanceShape_RSGToQQbar_M-3500.png");

  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                 "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
                 "RSG#rightarrowb#bar{b}, M=3500 GeV", "" , "Dijet Mass [GeV]", "", 0, 6000, 3250, 3600,
                 "ResonanceShape_RSGToBBbar_M-3500.png");
}
