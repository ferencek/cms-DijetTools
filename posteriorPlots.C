#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "exoStyle.C"


using namespace std;


void posteriors(const string& fFile1, const string& fFile2, const string& fTitle,
                const string& fXTitle, const string& fYTitle, const string& fLeg1, const string& fLeg2,
                const double fXmin, const double fXmax, const double fYmin, const double fYmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  TFile *file_1 = new TFile(fFile1.c_str());
  TFile *file_2 = new TFile(fFile2.c_str());

  TGraph *g_1 = (TGraph*)file_1->Get("post_0");
  g_1->SetLineColor(kRed);
  g_1->SetLineWidth(2);
  g_1->SetLineStyle(1);
  g_1->SetMarkerStyle(20);

  TGraph *g_2 = (TGraph*)file_2->Get("post_0");
  g_2->SetLineColor(kGreen+2);
  g_2->SetLineWidth(2);
  g_2->SetLineStyle(7);
  g_2->SetMarkerStyle(20);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->GetXaxis()->SetTitle(fXTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYTitle.c_str());
  bkg->Draw();

  g_1->Draw("L");
  g_2->Draw("L");

  TLegend *legend = new TLegend(.60,.60,.88,.80);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

//   c->SetGridx();
//   c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_1;
  delete g_2;
  delete file_1;
  delete file_2;
}


void makePlots()
{
  // qq, M=1100 GeV
  posteriors("stats_1100_1_qq_bkg_1st_param_posterior_xs0.root", "stats_1100_1_qq_bkg_1st_param_posterior_xs0.4.root",
             "qq, M=1100 GeV", "n_{1}/#sigma_{1}", "Posterior", "xs=0 pb", "xs=0.4 pb", -20, 20, 0, 1.1, "bkg_1st_param_posterior_qq_M-1100.eps");

  posteriors("stats_1100_1_qq_bkg_2nd_param_posterior_xs0.root", "stats_1100_1_qq_bkg_2nd_param_posterior_xs0.4.root",
             "qq, M=1100 GeV", "n_{2}/#sigma_{2}", "Posterior", "xs=0 pb", "xs=0.4 pb", -20, 20, 0, 1.1, "bkg_2nd_param_posterior_qq_M-1100.eps");

  posteriors("stats_1100_1_qq_bkg_3rd_param_posterior_xs0.root", "stats_1100_1_qq_bkg_3rd_param_posterior_xs0.4.root",
             "qq, M=1100 GeV", "n_{3}/#sigma_{3}", "Posterior", "xs=0 pb", "xs=0.4 pb", -20, 20, 0, 1.1, "bkg_3rd_param_posterior_qq_M-1100.eps");

  posteriors("stats_1100_1_qq_bkg_4th_param_posterior_xs0.root", "stats_1100_1_qq_bkg_4th_param_posterior_xs0.4.root",
             "qq, M=1100 GeV", "n_{4}/#sigma_{4}", "Posterior", "xs=0 pb", "xs=0.4 pb", -20, 20, 0, 1.1, "bkg_4th_param_posterior_qq_M-1100.eps");

  // qq, M=3000 GeV
  posteriors("stats_3000_1_qq_bkg_1st_param_posterior_xs0.root", "stats_3000_1_qq_bkg_1st_param_posterior_xs0.004.root",
             "qq, M=3000 GeV", "n_{1}/#sigma_{1}", "Posterior", "xs=0 pb", "xs=0.004 pb", -20, 20, 0, 1.1, "bkg_1st_param_posterior_qq_M-3000.eps");

  posteriors("stats_3000_1_qq_bkg_2nd_param_posterior_xs0.root", "stats_3000_1_qq_bkg_2nd_param_posterior_xs0.004.root",
             "qq, M=3000 GeV", "n_{2}/#sigma_{2}", "Posterior", "xs=0 pb", "xs=0.004 pb", -20, 20, 0, 1.1, "bkg_2nd_param_posterior_qq_M-3000.eps");

  posteriors("stats_3000_1_qq_bkg_3rd_param_posterior_xs0.root", "stats_3000_1_qq_bkg_3rd_param_posterior_xs0.004.root",
             "qq, M=3000 GeV", "n_{3}/#sigma_{3}", "Posterior", "xs=0 pb", "xs=0.004 pb", -20, 20, 0, 1.1, "bkg_3rd_param_posterior_qq_M-3000.eps");

  posteriors("stats_3000_1_qq_bkg_4th_param_posterior_xs0.root", "stats_3000_1_qq_bkg_4th_param_posterior_xs0.004.root",
             "qq, M=3000 GeV", "n_{4}/#sigma_{4}", "Posterior", "xs=0 pb", "xs=0.004 pb", -20, 20, 0, 1.1, "bkg_4th_param_posterior_qq_M-3000.eps");

}