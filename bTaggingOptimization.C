#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "tdrstyle.C"


void bTaggingOptimization_2Tag_15x15(const string& fDijetMassBin, const string& fBinLabel)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  //gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42, "t");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.10);
  gROOT->ForceStyle();

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET"; algoMap[9] = "SSVHPT";
  algoMap[10] = "JPL"; algoMap[11] = "JPM"; algoMap[12] = "JPT"; algoMap[13] = "CSVL"; algoMap[14] = "CSVM"; algoMap[15] = "CSVT";

  TFile *file1 = new TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_QCD_InclusiveTagging_240212/Final__histograms.root");
  TFile *file2 = new TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoBBbar_240212/Final__histograms.root");
  TFile *file3 = new TFile("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_RSGtoJJ_010312/Final__histograms.root");

  TH2D *h2_QCD_mistag_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_miseff_denom").c_str());
  TH2D *h2_QCD_mistag_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_miseff_num").c_str());

  TH2D *h2_QCD_eff_denom = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_QCD_eff_num = (TH2D*)file1->Get(("QCD_Pythia6__h2_BTaggers_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_RSGravitonToBBbar_eff_denom = (TH2D*)file2->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_RSGravitonToBBbar_eff_num = (TH2D*)file2->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_RSGravitonToGG_eff_denom = (TH2D*)file3->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_denom").c_str());
  TH2D *h2_RSGravitonToGG_eff_num = (TH2D*)file3->Get(("RSGravitonJJ__h2_BTaggers_" + fDijetMassBin + "_eff_num").c_str());

  TH2D *h2_QCD_mistag = new TH2D("h2_QCD_mistag", "h2_QCD_mistag", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_QCD_eff = new TH2D("h2_QCD_eff", "h2_QCD_eff", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_RSGravitonToBBbar_eff = new TH2D("h2_RSGravitonToBBbar_eff", "h2_RSGravitonToBBbar_eff", 15, 0.5, 15.5, 15, 0.5, 15.5);
  TH2D *h2_RSGravitonToGG_eff = new TH2D("h2_RSGravitonToGG_eff", "h2_RSGravitonToGG_eff", 15, 0.5, 15.5, 15, 0.5, 15.5);

  h2_QCD_mistag->SetTitle(("Mistag rate -- QCD, " + fBinLabel).c_str());
  h2_QCD_mistag->Divide(h2_QCD_mistag_num,h2_QCD_mistag_denom);

  h2_QCD_eff->SetTitle(("#epsilon_{QCD} -- " + fBinLabel).c_str());
  h2_QCD_eff->Divide(h2_QCD_eff_num,h2_QCD_eff_denom);

  h2_RSGravitonToBBbar_eff->SetTitle(("#epsilon_{G#rightarrowb#bar{b}} -- " + fBinLabel).c_str());
  h2_RSGravitonToBBbar_eff->Divide(h2_RSGravitonToBBbar_eff_num,h2_RSGravitonToBBbar_eff_denom);

  h2_RSGravitonToGG_eff->SetTitle(("#epsilon_{G#rightarrowgg} -- " + fBinLabel).c_str());
  h2_RSGravitonToGG_eff->Divide(h2_RSGravitonToGG_eff_num,h2_RSGravitonToGG_eff_denom);

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
    h2_RSGravitonToGG_eff->GetXaxis()->SetBinLabel(i,algoMap[i].c_str());
    h2_RSGravitonToGG_eff->GetYaxis()->SetBinLabel(i,algoMap[i].c_str());
  }

//   h2_QCD_mistag->Draw("colztext"); // For the inclusive case, the b-tag efficiency and the mistag rate are equal. Hence, this and the following line are commented out.
//   c->SaveAs(("DoubleTag_mistag_QCD_" + fDijetMassBin + ".png").c_str());

  h2_QCD_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_QCD_" + fDijetMassBin + "_Inclusive.png").c_str());

  h2_RSGravitonToBBbar_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_RSGravitonToBBbar_" + fDijetMassBin + ".png").c_str());

  h2_RSGravitonToGG_eff->Draw("colztext");
  c->SaveAs(("DoubleTag_eff_RSGravitonToGG_" + fDijetMassBin + "_Inclusive.png").c_str());

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

//   h2_opt_QCD->SetTitle(("Optimization -- #epsilon_{b,QCD}/#sqrt{#epsilon_{QCD}}, " + fBinLabel).c_str());
//   h2_opt_QCD->Draw("colztext");

//   c->SaveAs(("DoubleTag_optimization_QCD_" + fDijetMassBin + "_QCDInclusive.png").c_str());

//   h2_opt_RSGravitonToBBbar->SetTitle(("#epsilon_{G#rightarrowb#bar{b}}/#sqrt{#epsilon_{QCD}} -- " + fBinLabel).c_str());
//   gStyle->SetTitleFont(42);
  h2_opt_RSGravitonToBBbar->Draw("colztext");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.02,0.96, ("#epsilon_{G#rightarrowb#bar{b}}/#sqrt{#epsilon_{QCD}} -- " + fBinLabel).c_str());
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.87, "CMS Simulation");
  l1.DrawLatex(0.16,0.82, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.15,0.77, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.15,0.72, "Anti-k_{T} R = 0.7 PF Jets");

  c->SaveAs(("DoubleTag_optimization_RSGravitonToBBbar_" + fDijetMassBin + "_QCDInclusive.png").c_str());

  delete c;
}


void makePlots()
{

  bTaggingOptimization_2Tag_15x15("DijetMass500to1000GeV", "0.5<m_{jj}<1 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass1000to1500GeV", "1<m_{jj}<1.5 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass1500to2000GeV", "1.5<m_{jj}<2 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass2000to2500GeV", "2<m_{jj}<2.5 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass2500to3000GeV", "2.5<m_{jj}<3 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass3000to3500GeV", "3<m_{jj}<3.5 TeV");
  bTaggingOptimization_2Tag_15x15("DijetMass3500to4000GeV", "3.5<m_{jj}<4 TeV");
  
}


