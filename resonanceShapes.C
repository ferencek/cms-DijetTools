#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tdrstyle.C"


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


void ResonanceShape_v(const string& fInputFile1, const string& fInputFile2, const string& fPlot1, const string& fPlot2, const Int_t fNbins, const Double_t *fXbins,
                      const string& fTitle, const string& fLegend1, const string& fLegend2, const string& fXAxisTitle, const string& fYAxisTitle,
                      const Double_t fXmin, const Double_t fXmax,
                      const string& fOutputFile, const Double_t fYmin = 0., const Double_t fYmax = 0., const string& fPlot3 = "", const string& fPlot4 = "")
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h1_plot1_ = (TH1D*)file1->Get(fPlot1.c_str());
  h1_plot1_->Scale(1./h1_plot1_->Integral());
  TH1D *h1_plot1 = (TH1D*)h1_plot1_->Rebin(fNbins,"h1_DijetMass1",fXbins);

  TH1D *h1_plot2_ = (TH1D*)file2->Get(fPlot2.c_str());
  if(fPlot3!="" && fPlot4!="")
  {
    h1_plot2_->Add((TH1D*)file2->Get(fPlot3.c_str()));
    h1_plot2_->Add((TH1D*)file2->Get(fPlot4.c_str()));
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
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_plot1,fLegend1.c_str(),"l");
  legend->AddEntry(h1_plot2,fLegend2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.2,0.85, fTitle.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.61,0.34, "CMS Simulation");
  l1.DrawLatex(0.62,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.61,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.61,0.19, "Wide Jets");
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file1;
  delete file2;
}


void ResonanceShapes_finalStates(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3, const string& fPlot1, const string& fPlot2, const string& fPlot3,
                                 const string& fTitle, const string& fLegend1, const string& fLegend2, const string& fLegend3, const string& fXAxisTitle, const string& fYAxisTitle,
                                 const string& fOutputFile, const Int_t fSetDiv, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin = 0., const Double_t fYmax = 0.)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());
  TFile *file3 = new TFile(fInputFile3.c_str());

  TH1D *h1_plot1 = (TH1D*)file1->Get(fPlot1.c_str());
  TH1D *h1_plot2 = (TH1D*)file2->Get(fPlot2.c_str());
  TH1D *h1_plot3 = (TH1D*)file3->Get(fPlot3.c_str());

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

//   h1_plot1->SetTitle(fTitle.c_str());
  h1_plot1->SetTitleFont(42);
  h1_plot1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot1->GetXaxis()->SetRangeUser(fXmin,fXmax);
  if(fSetDiv) h1_plot1->GetXaxis()->SetNdivisions(506);
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

  TLegend *legend = new TLegend(.58,.65,.78,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_plot1,fLegend1.c_str(),"l");
  legend->AddEntry(h1_plot2,fLegend2.c_str(),"l");
  legend->AddEntry(h1_plot3,fLegend3.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.2,0.85, fTitle.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.61,0.34, "CMS Simulation");
  l1.DrawLatex(0.62,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.61,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.61,0.19, "Wide Jets");

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
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  string masses[7] = {"1000", "1500", "2000", "2500", "3000" , "3500", "4000"};
  string legends[7] = {"M = 1 TeV", "M = 1.5 TeV", "M = 2 TeV", "M = 2.5 TeV", "M = 3 TeV" , "M = 3.5 TeV", "M = 4 TeV"};
  
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

//   h1_plot[0]->SetTitle(fTitle.c_str());
  h1_plot[0]->SetTitleFont(42);
  h1_plot[0]->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot[0]->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot[0]->SetTitleOffset(1.,"Y");
  h1_plot[0]->GetXaxis()->SetRangeUser(fXmin,fXmax);
  if(fYmin!=fYmax) h1_plot[0]->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_plot[0]->Draw("hist");

  TLegend *legend = new TLegend(.68,.45,.88,.85);
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

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.17,0.85, fTitle.c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.70,0.34, "CMS Simulation");
  l1.DrawLatex(0.71,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.70,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.70,0.19, "Wide Jets");
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
}


void makePlots()
{

//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=500 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 450, 550,
//                  "ResonanceShape_RSGToQQbar_M-500.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=500 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 430, 550,
//                  "ResonanceShape_RSGToBBbar_M-500.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=700 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 625, 750,
//                  "ResonanceShape_RSGToQQbar_M-700.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=700 GeV", "", "Dijet Mass [GeV]", "", 0, 2000, 600, 740,
//                  "ResonanceShape_RSGToBBbar_M-700.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 20,
//                  "RSG#rightarrowq#bar{q}, M=1200 GeV", "", "Dijet Mass [GeV]", "", 0, 2500, 1100, 1300,
//                  "ResonanceShape_RSGToQQbar_M-1200.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 20,
//                  "RSG#rightarrowb#bar{b}, M=1200 GeV", "", "Dijet Mass [GeV]", "", 0, 2500, 1070, 1260,
//                  "ResonanceShape_RSGToBBbar_M-1200.eps");
// 
//  ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=2000 GeV", "", "Dijet Mass [GeV]", "", 0, 5000, 1850, 2100,
//                  "ResonanceShape_RSGToQQbar_M-2000.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=2000 GeV", "", "Dijet Mass [GeV]", "", 0, 5000, 1800, 2070,
//                  "ResonanceShape_RSGToBBbar_M-2000.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar", 40,
//                  "RSG#rightarrowq#bar{q}, M=3500 GeV", "", "Dijet Mass [GeV]", "", 0, 6000, 3300, 3600,
//                  "ResonanceShape_RSGToQQbar_M-3500.eps");
// 
//   ResonanceShape("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                  "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 40,
//                  "RSG#rightarrowb#bar{b}, M=3500 GeV", "", "Dijet Mass [GeV]", "", 0, 6000, 3250, 3600,
//                  "ResonanceShape_RSGToBBbar_M-3500.eps");

  // ##########################################
  // ## With variable dijet mass binning
  // ##########################################
  
  Double_t xbins[] = {0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693,740, 788,
                      838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147,
                      3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                    "M=500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1530,
//                    "ResonanceShape_RSG_M-500.eps");
// 
//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-700_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                    "M=700 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 1856,
//                    "ResonanceShape_RSG_M-700.eps");
// 
//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                    "M=1200 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 2775,
//                    "ResonanceShape_RSG_M-1200.eps");
// 
//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                    "M=2000 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 4509,
//                    "ResonanceShape_RSG_M-2000.eps");
// 
//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiency_CSVL_PUSFReweighted/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_qqbar",
//                    "myAnalyzer/cutHisto_allPreviousCuts________DijetMass_bbbar", 80, xbins,
//                    "M=3500 GeV", "RSG#rightarrowq#bar{q}", "RSG#rightarrowb#bar{b}", "Dijet Mass [GeV]", "", 0, 5663,
//                    "ResonanceShape_RSG_M-3500.eps");
                                            
  // M-1200
  // RSG->bbbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToBBbar_M-1200_CSVL_0Tag.eps", 0, 0.24);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToBBbar_M-1200_CSVL_1Tag.eps", 0, 0.24);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToBBbar_M-1200_CSVL_2Tag.eps", 0, 0.24);
    // RSG->qqbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_0tag", 80, xbins,
                   "M=1.2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToQQbar_M-1200_CSVL_0Tag.eps", 0, 0.35);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_1tag", 80, xbins,
                   "M=1.2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToQQbar_M-1200_CSVL_1Tag.eps", 0, 0.35);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_2tag", 80, xbins,
                   "M=1.2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToQQbar_M-1200_CSVL_2Tag.eps", 0, 0.35);
  // RSG->gg
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_0tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowgg", "G#rightarrowgg, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToGG_M-1200_CSVL_0Tag.eps", 0, 0.25);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_1tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowgg", "G#rightarrowgg, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToGG_M-1200_CSVL_1Tag.eps", 0, 0.25);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-1200_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_2tag", 80, xbins,
                   "M=1.2 TeV", "G#rightarrowgg", "G#rightarrowgg, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 2775,
                   "ResonanceShape_RSGToGG_M-1200_CSVL_2Tag.eps", 0, 0.25);
                      
  // M-2000
  // RSG->bbbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-2000_CSVL_0Tag.eps", 0, 0.22);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-2000_CSVL_1Tag.eps", 0, 0.22);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-2000_CSVL_2Tag.eps", 0, 0.22);
  // RSG->qqbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_0tag", 80, xbins,
                   "M=2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-2000_CSVL_0Tag.eps", 0, 0.35);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_1tag", 80, xbins,
                   "M=2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-2000_CSVL_1Tag.eps", 0, 0.35);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_2tag", 80, xbins,
                   "M=2 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-2000_CSVL_2Tag.eps", 0, 0.35);
  // RSG->gg
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_0tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowgg", "G#rightarrowgg, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-2000_CSVL_0Tag.eps", 0, 0.25);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_1tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowgg", "G#rightarrowgg, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-2000_CSVL_1Tag.eps", 0, 0.25);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_2tag", 80, xbins,
                   "M=2 TeV", "G#rightarrowgg", "G#rightarrowgg, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-2000_CSVL_2Tag.eps", 0, 0.25);

  // M-3500
  // RSG->bbbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_0tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-3500_CSVL_0Tag.eps", 0, 0.25);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_1tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-3500_CSVL_1Tag.eps", 0, 0.25);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_bbbar",
                   "myAnalyzer/h1_DijetMass_bbbar_2tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowb#bar{b}", "G#rightarrowb#bar{b}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToBBbar_M-3500_CSVL_2Tag.eps", 0, 0.25);
  // RSG->qqbar
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_0tag", 80, xbins,
                   "M=3.5 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-3500_CSVL_0Tag.eps", 0, 0.40);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_1tag", 80, xbins,
                   "M=3.5 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-3500_CSVL_1Tag.eps", 0, 0.40);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_qqbarlight",
                   "myAnalyzer/h1_DijetMass_qqbarlight_2tag", 80, xbins,
                   "M=3.5 TeV (q=u,d,s)", "G#rightarrowq#bar{q}", "G#rightarrowq#bar{q}, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToQQbar_M-3500_CSVL_2Tag.eps", 0, 0.40);
  // RSG->gg
  // CSVL 0Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_0tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowgg", "G#rightarrowgg, 0 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-3500_CSVL_0Tag.eps", 0, 0.25);
  // CSVL 1Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_1tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowgg", "G#rightarrowgg, 1 b-tag", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-3500_CSVL_1Tag.eps", 0, 0.25);
  // CSVL 2Tag
  ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-3500_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                   "myAnalyzer/h1_DijetMass_gg",
                   "myAnalyzer/h1_DijetMass_gg_2tag", 80, xbins,
                   "M=3.5 TeV", "G#rightarrowgg", "G#rightarrowgg, 2 b-tags", "Dijet Mass [GeV]", "Probability", 0, 6000,
                   "ResonanceShape_RSGToGG_M-3500_CSVL_2Tag.eps", 0, 0.25);
  
                      
//   // ## Impact of different initial states on resonance shapes
//   // M-2000
//   // RSG->gg
//   ResonanceShape_v("CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "CRAB_Jobs_RSGraviton_ResonanceShapes_bTagEfficiencies_CSVL_PUSFReweighted_ggInitialState_WideJets/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                    "myAnalyzer/h1_DijetMass_gg",
//                    "myAnalyzer/h1_DijetMass_gg", 80, xbins,
//                    "M=2 TeV", "G#rightarrowgg", "gg#rightarrowG#rightarrowgg", "Dijet Mass [GeV]", "Probability", 0, 6000,
//                    "ResonanceShape_RSGToGG_M-2000_ggInitialState.eps", 0, 0.35);
  
  
  // ## Resonance shapes for different final states
  ResonanceShapes_finalStates("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
                              "h_qq_1000", "h_bb_1000", "h_gg_1000", "M=1 TeV", "G#rightarrowq#bar{q}  (q=u,d,s,c)", "G#rightarrowb#bar{b}", "G#rightarrowgg",
                              "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-1000.eps", 1, 1, 2600);

  ResonanceShapes_finalStates("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
                              "h_qq_2000", "h_bb_2000", "h_gg_2000", "M=2 TeV", "G#rightarrowq#bar{q}  (q=u,d,s,c)", "G#rightarrowb#bar{b}", "G#rightarrowgg",
                              "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-2000.eps", 1, 1, 4500);

  ResonanceShapes_finalStates("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
                              "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
                              "h_qq_3000", "h_bb_3000", "h_gg_3000", "M=3 TeV", "G#rightarrowq#bar{q}  (q=u,d,s,c)", "G#rightarrowb#bar{b}", "G#rightarrowgg",
                              "Dijet Mass [GeV]", "Probability", "ResonanceShapes_RSG_M-3000.eps", 0, 1, 6000);

  
  // ## Resonance shapes for different masses
  ResonanceShapes_masses("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
                         "h_qq_", "G#rightarrowq#bar{q}  (q=u,d,s,c)", "Dijet Mass [GeV]", "",
                         "ResonanceShapes_RSGToQQbar.eps", 1, 6000, 0, 0.33);
  ResonanceShapes_masses("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
                         "h_bb_", "G#rightarrowb#bar{b}", "Dijet Mass [GeV]", "Probability",
                         "ResonanceShapes_RSGToBBbar.eps", 1, 6000, 0, 0.25);
  ResonanceShapes_masses("LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
                         "h_gg_", "G#rightarrowgg", "Dijet Mass [GeV]", "",
                         "ResonanceShapes_RSGToGG.eps", 1, 6000, 0, 0.2);

}
