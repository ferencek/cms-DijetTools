#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

void setTDRStyle();

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
                             const string& fOutputFile, const Double_t fYmin = 0., const Double_t fYmax = 0., const string& fPlot3 = "", const string& fPlot4 = "")
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot1_ = (TH1D*)file->Get(fPlot1.c_str());
  h1_plot1_->Scale(1./h1_plot1_->Integral());
  TH1D *h1_plot1 = (TH1D*)h1_plot1_->Rebin(fNbins,"h1_DijetMass1",fXbins);

  TH1D *h1_plot2_ = (TH1D*)file->Get(fPlot2.c_str());
  if(fPlot3!="" && fPlot4!="")
  {
    h1_plot2_->Add((TH1D*)file->Get(fPlot3.c_str()));
    h1_plot2_->Add((TH1D*)file->Get(fPlot4.c_str()));
  }
  h1_plot2_->Scale(1./h1_plot2_->Integral());
  TH1D *h1_plot2 = (TH1D*)h1_plot2_->Rebin(fNbins,"h1_DijetMass2",fXbins);
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot1->SetTitle(fTitle.c_str());
  h1_plot1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  if(fYmin!=fYmax) h1_plot1->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_plot1->SetFillColor(0);
  h1_plot1->SetMarkerStyle(0);
  h1_plot1->SetLineColor(kBlue);
  h1_plot1->SetLineWidth(2);

  h1_plot2->SetFillColor(0);
  h1_plot2->SetMarkerStyle(0);
  h1_plot2->SetLineColor(kRed);
  h1_plot2->SetLineStyle(2);
  h1_plot2->SetLineWidth(2);
  
  h1_plot1->Draw("hist");
  h1_plot1->Draw("same");

  h1_plot2->Draw("histsame");
  h1_plot2->Draw("same");

  TLegend *legend = new TLegend(.5,.65,.85,.8);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
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

void ResonanceShapes_finalStates(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3, const string& fPlot1, const string& fPlot2, const string& fPlot3,
                                 const string& fTitle, const string& fLegend1, const string& fLegend2, const string& fLegend3, const string& fXAxisTitle, const string& fYAxisTitle,
                                 const string& fOutputFile, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin = 0., const Double_t fYmax = 0.)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());
  TFile *file3 = new TFile(fInputFile3.c_str());

  TH1D *h1_plot1 = (TH1D*)file1->Get(fPlot1.c_str());
  TH1D *h1_plot2 = (TH1D*)file2->Get(fPlot2.c_str());
  TH1D *h1_plot3 = (TH1D*)file3->Get(fPlot3.c_str());

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot1->SetTitle(fTitle.c_str());
  h1_plot1->SetTitleFont(42);
  h1_plot1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  if(fYmin!=fYmax) h1_plot1->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_plot1->SetLineColor(kGreen+2);
  h1_plot1->SetFillColor(0);
  h1_plot1->SetLineWidth(2);

  h1_plot2->SetLineColor(kBlue);
  h1_plot2->SetFillColor(0);
  h1_plot2->SetLineWidth(2);
  h1_plot2->SetLineStyle(7);

  h1_plot3->SetLineColor(kRed);
  h1_plot3->SetFillColor(0);
  h1_plot3->SetLineWidth(2);
  h1_plot3->SetLineStyle(4);

  h1_plot1->Draw("hist");
  h1_plot2->Draw("histsame");
  h1_plot3->Draw("histsame");

  TLegend *legend = new TLegend(.65,.65,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_plot1,fLegend1.c_str(),"l");
  legend->AddEntry(h1_plot2,fLegend2.c_str(),"l");
  legend->AddEntry(h1_plot3,fLegend3.c_str(),"l");
  legend->Draw();

  gPad->RedrawAxis();
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
}

void ResonanceShapes_masses(const string& fInputFile, const string& fPlotPrefix, const string& fTitle,
                            const string& fXAxisTitle, const string& fYAxisTitle, const string& fOutputFile,
                            const Double_t fXmin, const Double_t fXmax, const Double_t fYmin = 0., const Double_t fYmax = 0.)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  string masses[7] = {"1000", "1500", "2000", "2500", "3000" , "3500", "4000"};
  string legends[7] = {"1 TeV", "1.5 TeV", "2 TeV", "2.5 TeV", "3 TeV" , "3.5 TeV", "4 TeV"};
  
  TH1D *h1_plot[7];

  for(int i=0; i<7; ++i)
  {
    h1_plot[i] = (TH1D*)file->Get((fPlotPrefix + masses[i]).c_str());
    h1_plot[i]->SetLineColor(i+1);
    h1_plot[i]->SetFillColor(0);
    h1_plot[i]->SetLineWidth(2);

  }
  
  TCanvas *c = new TCanvas("c", "",1200,800);
  c->cd();

  h1_plot[0]->SetTitle(fTitle.c_str());
  h1_plot[0]->SetTitleFont(42);
  h1_plot[0]->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot[0]->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot[0]->GetXaxis()->SetRangeUser(fXmin,fXmax);
  if(fYmin!=fYmax) h1_plot[0]->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_plot[0]->Draw("hist");

  TLegend *legend = new TLegend(.75,.35,.95,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_plot[0],legends[0].c_str(),"l");
  
  for(int i=1; i<7; ++i)
  {
    h1_plot[i]->Draw("histsame");
    legend->AddEntry(h1_plot[i],legends[i].c_str(),"l");
  }

  legend->Draw();

  gPad->RedrawAxis();
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
}


void makePlots()
{

//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=500 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 450, 550,
//                  "ResonanceShape_RSGToQQbar_M-500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=500 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 430, 550,
//                  "ResonanceShape_RSGToBBbar_M-500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=700 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 625, 750,
//                  "ResonanceShape_RSGToQQbar_M-700.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=700 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 600, 740,
//                  "ResonanceShape_RSGToBBbar_M-700.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=1200 GeV", "", "Dijet Mass [GeV]", "", 0, 2500, 1100, 1300,
//                  "ResonanceShape_RSGToQQbar_M-1200.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=1200 GeV", "", "Dijet Mass [GeV]", "", 0, 2500, 1070, 1260,
//                  "ResonanceShape_RSGToBBbar_M-1200.png");
// 
//  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=2000 GeV", "", "Dijet Mass [GeV]", "", 0, 5000, 1850, 2100,
//                  "ResonanceShape_RSGToQQbar_M-2000.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=2000 GeV", "", "Dijet Mass [GeV]", "", 0, 5000, 1800, 2070,
//                  "ResonanceShape_RSGToBBbar_M-2000.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=3500 GeV", "", "Dijet Mass [GeV]", "", 0, 6000, 3300, 3600,
//                  "ResonanceShape_RSGToQQbar_M-3500.png");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=3500 GeV", "", "Dijet Mass [GeV]", "", 0, 6000, 3250, 3600,
//                  "ResonanceShape_RSGToBBbar_M-3500.png");

  // ##########################################
  // ## With variable dijet mass binning
  // ##########################################
  
  Double_t xbins[] = {0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693,740, 788,
                      838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147,
                      3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1530,
//                           "ResonanceShape_RSG_M-500.png");

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=700 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1856,
//                           "ResonanceShape_RSG_M-700.png");

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=1200 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 2775,
//                           "ResonanceShape_RSG_M-1200.png");

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=2000 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 4509,
//                           "ResonanceShape_RSG_M-2000.png");

//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           "M=3500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 5663,
//                           "ResonanceShape_RSG_M-3500.png");

                                            
//   // M-1200
//   // CSVL 0Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
//                           "M=1200 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 0Tag", "Dijet Mass [GeV]", "", 0, 2775,
//                           "ResonanceShape_RSG_M-1200_CSVL_0Tag.png", 0, 0.23);
//   // CSVL 1Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
//                           "M=1200 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 1Tag", "Dijet Mass [GeV]", "", 0, 2775,
//                           "ResonanceShape_RSG_M-1200_CSVL_1Tag.png", 0, 0.23);
//   // CSVL 2Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
//                           "M=1200 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 2Tag", "Dijet Mass [GeV]", "", 0, 2775,
//                           "ResonanceShape_RSG_M-1200_CSVL_2Tag.png", 0, 0.23);
//   
  // M-2000
  // RSG->bbbar
  // CSVL 0Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToBBbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 0-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToBBbar_M-2000_CSVL_0Tag.png", 0, 0.22);
  // CSVL 1Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToBBbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 1-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToBBbar_M-2000_CSVL_1Tag.png", 0, 0.22);
  // CSVL 2Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToBBbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 2-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToBBbar_M-2000_CSVL_2Tag.png", 0, 0.22);
  // RSG->qqbar
  // CSVL 0Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToQQbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_qqbarlight_0tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowq#bar{q}, CSVL 0-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToQQbar_M-2000_CSVL_0Tag.png", 0, 0.27, "myAnalyzer/h1_DijetMass_bbbar_0tag", "myAnalyzer/h1_DijetMass_ccbar_0tag");
  // CSVL 1Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToQQbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_qqbarlight_1tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowq#bar{q}, CSVL 1-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToQQbar_M-2000_CSVL_1Tag.png", 0, 0.27, "myAnalyzer/h1_DijetMass_bbbar_1tag", "myAnalyzer/h1_DijetMass_ccbar_1tag");
  // CSVL 2Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToQQbar_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_qqbarlight_2tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowq#bar{q}, CSVL 2-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToQQbar_M-2000_CSVL_2Tag.png", 0, 0.27, "myAnalyzer/h1_DijetMass_bbbar_2tag", "myAnalyzer/h1_DijetMass_ccbar_2tag");
  // RSG->gg
  // CSVL 0Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToGG_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_gg_0tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowgg", "RSG#rightarrowgg, CSVL 0-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToGG_M-2000_CSVL_0Tag.png", 0, 0.25);
  // CSVL 1Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToGG_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_gg_1tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowgg", "RSG#rightarrowgg, CSVL 1-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToGG_M-2000_CSVL_1Tag.png", 0, 0.25);
  // CSVL 2Tag
  ResonanceShape_var_bins("Resonance_shape_files/CRAB_Jobs_RSGravitonToGG_ResonanceShapes_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
                          "myAnalyzer/h1_DijetMass_gg_2tag", 80, xbins,
                          "M=2 TeV, Wide Jets", "RSG#rightarrowgg", "RSG#rightarrowgg, CSVL 2-tag", "Dijet Mass [GeV]", "", 0, 6000,
                          "ResonanceShape_RSGToGG_M-2000_CSVL_2Tag.png", 0, 0.25);
//   
//   // M-3500
//   // CSVL 0Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
//                           "M=3500 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 0Tag", "Dijet Mass [GeV]", "", 0, 6000,
//                           "ResonanceShape_RSG_M-3500_CSVL_0Tag.png", 0, 0.16);
//   // CSVL 1Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
//                           "M=3500 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 1Tag", "Dijet Mass [GeV]", "", 0, 6000,
//                           "ResonanceShape_RSG_M-3500_CSVL_1Tag.png", 0, 0.16);
//   // CSVL 2Tag
//   ResonanceShape_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar",
//                           "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
//                           "M=3500 GeV", "RSG#rightarrowb#bar{b}", "RSG#rightarrowb#bar{b}, CSVL 2Tag", "Dijet Mass [GeV]", "", 0, 6000,
//                           "ResonanceShape_RSG_M-3500_CSVL_2Tag.png", 0, 0.16);

//   // ## Resonance shapes for different final states
//   ResonanceShapes_finalStates("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_qq.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_bb.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_gg.root",
//                               "h_qq_1000", "h_bb_1000", "h_gg_1000", "M=1 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "RSG#rightarrowgg",
//                               "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-1000.png", 1, 3416);
// 
//   ResonanceShapes_finalStates("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_qq.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_bb.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_gg.root",
//                               "h_qq_2000", "h_bb_2000", "h_gg_2000", "M=2 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "RSG#rightarrowgg",
//                               "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-2000.png", 1, 3416);
// 
//   ResonanceShapes_finalStates("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_qq.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_bb.root",
//                               "Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_gg.root",
//                               "h_qq_3000", "h_bb_3000", "h_gg_3000", "M=3 TeV, Wide Jets", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "RSG#rightarrowgg",
//                               "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-3000.png", 1, 3416);

//   // ## Resonance shapes for different masses
//   ResonanceShapes_masses("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_qq.root",
//                          "h_qq_", "RSG#rightarrowq#bar{q}, Wide Jets", "Dijet Mass [GeV]", "",
//                          "ResonanceShapes_RSGToQQbar.png", 1, 5876, 0, 0.3);
//   ResonanceShapes_masses("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_bb.root",
//                          "h_bb_", "RSG#rightarrowb#bar{b}, Wide Jets", "Dijet Mass [GeV]", "",
//                          "ResonanceShapes_RSGToBBbar.png", 1, 5876, 0, 0.3);
//   ResonanceShapes_masses("Resonance_shape_files/Interpolation_code/Resonance_Shapes_WideJets_gg.root",
//                          "h_gg_", "RSG#rightarrowgg, Wide Jets", "Dijet Mass [GeV]", "",
//                          "ResonanceShapes_RSGToGG.png", 1, 5876, 0, 0.3);
  
  // ## Reshaped resonance shape

//   ResonanceShape_reshaped_var_bins("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                           "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                           3429., 143.9, 3383., 171.5,
//                           "M=3500 GeV", "#splitline{RSG#rightarrowq#bar{q}}{reshaped}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 5663,
//                           "ResonanceShape_RSG_M-3500_reshaped.png");
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(63);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

//  tdrStyle->SetEndErrorSize(0);
//   tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.07);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.07);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
