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

void bTaggingOptimization_2Tag_9x9(const string& fDijetMassBin, const string& fBinLabel)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLabelSize(0.035, "XYZ");

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET"; algoMap[9] = "SSVHPT";

  TFile *file1 = new TFile("CRAB_Jobs_bTaggingOptimization_bPartonMatching/Final__histograms.root");
  TFile *file2 = new TFile("CRAB_Jobs_bTaggingOptimization_bPartonMatching/Final__histograms.root");
  TFile *file3 = new TFile("CRAB_Jobs_bTaggingOptimization_RSG_bPartonMatching/Final__histograms.root");

  TH2D *h2_QCD_mistag_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_miseff_denom").c_str());
  TH2D *h2_QCD_mistag_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_miseff_num").c_str());

  TH2D *h2_QCD_eff_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_QCD_eff_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_ZprimeToBBbar_eff_denom = (TH2D*)file2->Get(("ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_ZprimeToBBbar_eff_num = (TH2D*)file2->Get(("ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_RSGravitonToBBbar_eff_denom = (TH2D*)file3->Get(("RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_RSGravitonToBBbar_eff_num = (TH2D*)file3->Get(("RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_QCD_mistag = new TH2D("h2_QCD_mistag", "h2_QCD_mistag", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_QCD_eff = new TH2D("h2_QCD_eff", "h2_QCD_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_ZprimeToBBbar_eff = new TH2D("h2_ZprimeToBBbar_eff", "h2_ZprimeToBBbar_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_RSGravitonToBBbar_eff = new TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);

  h2_QCD_mistag->SetTitle(("Mistag rate -- QCD, " + fBinLabel).c_str());
  h2_QCD_mistag->Divide(h2_QCD_mistag_num,h2_QCD_mistag_denom);

  h2_QCD_eff->SetTitle(("Efficiency -- QCD, " + fBinLabel).c_str());
  h2_QCD_eff->Divide(h2_QCD_eff_num,h2_QCD_eff_denom);

  h2_ZprimeToBBbar_eff->SetTitle(("Efficiency -- Z'#rightarrowb#bar{b}, " + fBinLabel).c_str());
  h2_ZprimeToBBbar_eff->Divide(h2_ZprimeToBBbar_eff_num,h2_ZprimeToBBbar_eff_denom);

  h2_RSGravitonToBBbar_eff->SetTitle(("Efficiency -- RSG#rightarrowb#bar{b}, " + fBinLabel).c_str());
  h2_RSGravitonToBBbar_eff->Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  for(int i=1; i<=9; i++)
  {
    h2_QCD_mistag->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_mistag->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_ZprimeToBBbar_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_ZprimeToBBbar_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
  }

  h2_QCD_mistag->Draw("colztext");
  c->SaveAs(("DoubleTag_mistag_QCD_" + fDijetMassBin + ".png").c_str());

  h2_QCD_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_QCD_" + fDijetMassBin + ".png").c_str());

  h2_ZprimeToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_ZprimeToBBbar_" + fDijetMassBin + ".png").c_str());

  h2_RSGravitonToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_RSGravitonToBBbar_" + fDijetMassBin + ".png").c_str());


  TH2D *h2_opt_QCD = new TH2D("h2_opt_QCD", "h2_opt_QCD", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_opt_ZprimeToBBbar = new TH2D("h2_opt_ZprimeToBBbar", "h2_opt_ZprimeToBBbar", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_opt_RSGravitonToBBbar = new TH2D("h2_opt_RSGravitonToBBbar", "h2_opt_RSGravitonToBBbar", 9, 0.5, 9.5, 9, 0.5, 9.5);

  for(int i=1; i<=9; i++)
  {
    h2_opt_QCD->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_QCD->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_ZprimeToBBbar->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_ZprimeToBBbar->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());

    for(int j=1; j<=9; j++)
    {
      if( h2_QCD_mistag->GetBinContent(i,j)>0 )
      {
        h2_opt_QCD->SetBinContent(i,j,(h2_QCD_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
        h2_opt_ZprimeToBBbar->SetBinContent(i,j,(h2_ZprimeToBBbar_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
        h2_opt_RSGravitonToBBbar->SetBinContent(i,j,(h2_RSGravitonToBBbar_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
      }
    }
  }

  h2_opt_QCD->SetTitle(("Optimization -- #epsilon_{b,QCD}/#sqrt{#epsilon_{non-b,QCD}}, " + fBinLabel).c_str());
  h2_opt_QCD->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_QCD_" + fDijetMassBin + ".png").c_str());

  h2_opt_ZprimeToBBbar->SetTitle(("Optimization -- #epsilon_{Z'#rightarrowb#bar{b}}/#sqrt{#epsilon_{non-b,QCD}}, " + fBinLabel).c_str());
  h2_opt_ZprimeToBBbar->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_ZprimeToBBbar_" + fDijetMassBin + ".png").c_str());

  h2_opt_RSGravitonToBBbar->SetTitle(("Optimization -- #epsilon_{RSG#rightarrowb#bar{b}}/#sqrt{#epsilon_{non-b,QCD}}, " + fBinLabel).c_str());
  h2_opt_RSGravitonToBBbar->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_RSGravitonToBBbar_" + fDijetMassBin + ".png").c_str());

  delete c;
}


void bTaggingOptimization_2Tag_9x9_QCDInclusive(const string& fDijetMassBin, const string& fBinLabel)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLabelSize(0.035, "XYZ");

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET"; algoMap[9] = "SSVHPT";

  TFile *file1 = new TFile("CRAB_Jobs_bTaggingOptimization_QCD_InclusiveTagging_bPartonMatching/Final__histograms.root");
  TFile *file2 = new TFile("CRAB_Jobs_bTaggingOptimization_bPartonMatching/Final__histograms.root");
  TFile *file3 = new TFile("CRAB_Jobs_bTaggingOptimization_RSG_bPartonMatching/Final__histograms.root");

  TH2D *h2_QCD_mistag_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_miseff_denom").c_str());
  TH2D *h2_QCD_mistag_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_miseff_num").c_str());

  TH2D *h2_QCD_eff_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_QCD_eff_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_ZprimeToBBbar_eff_denom = (TH2D*)file2->Get(("ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_ZprimeToBBbar_eff_num = (TH2D*)file2->Get(("ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());
  
  TH2D *h2_RSGravitonToBBbar_eff_denom = (TH2D*)file3->Get(("RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_RSGravitonToBBbar_eff_num = (TH2D*)file3->Get(("RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_QCD_mistag = new TH2D("h2_QCD_mistag", "h2_QCD_mistag", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_QCD_eff = new TH2D("h2_QCD_eff", "h2_QCD_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_ZprimeToBBbar_eff = new TH2D("h2_ZprimeToBBbar_eff", "h2_ZprimeToBBbar_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_RSGravitonToBBbar_eff = new TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 9, 0.5, 9.5, 9, 0.5, 9.5);

  h2_QCD_mistag->SetTitle(("Mistag rate -- QCD, " + fBinLabel).c_str());
  h2_QCD_mistag->Divide(h2_QCD_mistag_num,h2_QCD_mistag_denom);

  h2_QCD_eff->SetTitle(("Efficiency -- QCD, " + fBinLabel).c_str());
  h2_QCD_eff->Divide(h2_QCD_eff_num,h2_QCD_eff_denom);

  h2_ZprimeToBBbar_eff->SetTitle(("Efficiency -- Z'#rightarrowb#bar{b}, " + fBinLabel).c_str());
  h2_ZprimeToBBbar_eff->Divide(h2_ZprimeToBBbar_eff_num,h2_ZprimeToBBbar_eff_denom);

  h2_RSGravitonToBBbar_eff->SetTitle(("Efficiency -- RSG#rightarrowb#bar{b}, " + fBinLabel).c_str());
  h2_RSGravitonToBBbar_eff->Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom);
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  for(int i=1; i<=9; i++)
  {
    h2_QCD_mistag->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_mistag->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_ZprimeToBBbar_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_ZprimeToBBbar_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
  }
  
//   h2_QCD_mistag->Draw("colztext"); // For the inclusive case, the b-tag efficiency and the mistag rate are equal. Hence, this and the following line are commented out.
//   c->SaveAs(("DoubleTag_mistag_QCD_" + fDijetMassBin + ".png").c_str());

  h2_QCD_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_QCD_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  h2_ZprimeToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_ZprimeToBBbar_" + fDijetMassBin + ".png").c_str());

  h2_RSGravitonToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_RSGravitonToBBbar_" + fDijetMassBin + ".png").c_str());


  TH2D *h2_opt_QCD = new TH2D("h2_opt_QCD", "h2_opt_QCD", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_opt_ZprimeToBBbar = new TH2D("h2_opt_ZprimeToBBbar", "h2_opt_ZprimeToBBbar", 9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *h2_opt_RSGravitonToBBbar = new TH2D("h2_opt_RSGravitonToBBbar", "h2_opt_RSGravitonToBBbar", 9, 0.5, 9.5, 9, 0.5, 9.5);

  for(int i=1; i<=9; i++)
  {
    h2_opt_QCD->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_QCD->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_ZprimeToBBbar->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_ZprimeToBBbar->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
   
    for(int j=1; j<=9; j++)
    {
      if( h2_QCD_mistag->GetBinContent(i,j)>0 )
      {
        h2_opt_QCD->SetBinContent(i,j,(h2_QCD_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
        h2_opt_ZprimeToBBbar->SetBinContent(i,j,(h2_ZprimeToBBbar_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
        h2_opt_RSGravitonToBBbar->SetBinContent(i,j,(h2_RSGravitonToBBbar_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
      }
    }
  }

  h2_opt_QCD->SetTitle(("Optimization -- #epsilon_{b,QCD}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
  h2_opt_QCD->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_QCD_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  h2_opt_ZprimeToBBbar->SetTitle(("Optimization -- #epsilon_{Z'#rightarrowb#bar{b}}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
  h2_opt_ZprimeToBBbar->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_ZprimeToBBbar_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  h2_opt_RSGravitonToBBbar->SetTitle(("Optimization -- #epsilon_{RSG#rightarrowb#bar{b}}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
  h2_opt_RSGravitonToBBbar->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_RSGravitonToBBbar_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  delete c;
}


void bTaggingOptimization_2Tag_15x15_QCDInclusive(const string& fDijetMassBin, const string& fBinLabel)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetPaintTextFormat("1.2g");

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET";
  algoMap[9] = "SSVHPT"; algoMap[10] = "JPL"; algoMap[11] = "JPM"; algoMap[12] = "JPT"; algoMap[13] = "CSVL"; algoMap[14] = "CSVM"; algoMap[15] = "CSVT";

  TFile *file1 = new TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_QCD_InclusiveTagging_240212/Final__histograms.root");
  TFile *file2 = new TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoBBbar_240212/Final__histograms.root");

  TH2D *h2_QCD_mistag_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_miseff_denom").c_str());
  TH2D *h2_QCD_mistag_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_miseff_num").c_str());

  TH2D *h2_QCD_eff_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_QCD_eff_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_RSGravitonToBBbar_eff_denom = (TH2D*)file2->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_RSGravitonToBBbar_eff_num = (TH2D*)file2->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_QCD_mistag = new TH2D("h2_QCD_mistag", "h2_QCD_mistag", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_QCD_eff = new TH2D("h2_QCD_eff", "h2_QCD_eff", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_RSGravitonToBBbar_eff = new TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 15, 0.5, 15.5, 15, 0.5, 15.5);

  h2_QCD_mistag->SetTitle(("Mistag rate -- QCD, " + fBinLabel).c_str());
  h2_QCD_mistag->Divide(h2_QCD_mistag_num,h2_QCD_mistag_denom);

  h2_QCD_eff->SetTitle(("Efficiency -- QCD, " + fBinLabel).c_str());
  h2_QCD_eff->Divide(h2_QCD_eff_num,h2_QCD_eff_denom);

  h2_RSGravitonToBBbar_eff->SetTitle(("Efficiency -- RSG#rightarrowb#bar{b}, " + fBinLabel).c_str());
  h2_RSGravitonToBBbar_eff->Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom);

  TCanvas *c = new TCanvas("c", "",1000,1000);
  c->cd();

  for(int i=1; i<=15; i++)
  {
    h2_QCD_mistag->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_mistag->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_QCD_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToBBbar_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
  }

//   h2_QCD_mistag->Draw("colztext"); // For the inclusive case, the b-tag efficiency and the mistag rate are equal. Hence, this and the following line are commented out.
//   c->SaveAs(("DoubleTag_mistag_QCD_" + fDijetMassBin + ".png").c_str());

  h2_QCD_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_QCD_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  h2_RSGravitonToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_RSGravitonToBBbar_" + fDijetMassBin + ".png").c_str());


  TH2D *h2_opt_QCD = new TH2D("h2_opt_QCD", "h2_opt_QCD", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_opt_RSGravitonToBBbar = new TH2D("h2_opt_RSGravitonToBBbar", "h2_opt_RSGravitonToBBbar", 15, 0.5, 15.5, 15, 0.5, 15.5);

  for(int i=1; i<=15; i++)
  {
    h2_opt_QCD->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_QCD->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_opt_RSGravitonToBBbar->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());

    for(int j=1; j<=15; j++)
    {
      if( h2_QCD_mistag->GetBinContent(i,j)>0 )
      {
        h2_opt_QCD->SetBinContent(i,j,(h2_QCD_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
        h2_opt_RSGravitonToBBbar->SetBinContent(i,j,(h2_RSGravitonToBBbar_eff->GetBinContent(i,j)/sqrt(h2_QCD_mistag->GetBinContent(i,j))));
      }
    }
  }

  h2_opt_QCD->SetTitle(("Optimization -- #epsilon_{b,QCD}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
  h2_opt_QCD->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_QCD_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  h2_opt_RSGravitonToBBbar->SetTitle(("Optimization -- #epsilon_{RSG#rightarrowb#bar{b}}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
  h2_opt_RSGravitonToBBbar->Draw("colztext");

  c->SaveAs(("DoubleTag_optimization_RSGravitonToBBbar_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  delete c;
}


void makePlots()
{
  //bTaggingOptimization_2Tag_9x9("DijetMass500to1000GeV", "0.5<M_{jj}<1 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass1000to1500GeV", "1<M_{jj}<1.5 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass1500to2000GeV", "1.5<M_{jj}<2 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass2000to2500GeV", "2<M_{jj}<2.5 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass2500to3000GeV", "2.5<M_{jj}<3 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass3000to3500GeV", "3<M_{jj}<3.5 TeV");
  //bTaggingOptimization_2Tag_9x9("DijetMass3500to4000GeV", "3.5<M_{jj}<4 TeV");

  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass500to1000GeV", "0.5<M_{jj}<1 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass1000to1500GeV", "1<M_{jj}<1.5 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass1500to2000GeV", "1.5<M_{jj}<2 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass2000to2500GeV", "2<M_{jj}<2.5 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass2500to3000GeV", "2.5<M_{jj}<3 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass3000to3500GeV", "3<M_{jj}<3.5 TeV");
  //bTaggingOptimization_2Tag_9x9_QCDInclusive("DijetMass3500to4000GeV", "3.5<M_{jj}<4 TeV");

  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass500to1000GeV", "0.5<M_{jj}<1 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass1000to1500GeV", "1<M_{jj}<1.5 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass1500to2000GeV", "1.5<M_{jj}<2 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass2000to2500GeV", "2<M_{jj}<2.5 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass2500to3000GeV", "2.5<M_{jj}<3 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass3000to3500GeV", "3<M_{jj}<3.5 TeV");
  bTaggingOptimization_2Tag_15x15_QCDInclusive("DijetMass3500to4000GeV", "3.5<M_{jj}<4 TeV");
}


