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
#include "TEfficiency.h"
#include "TList.h"

void bTaggingEfficiency_2Tag(const string& fInputFile, const string& fPlotDenom, const string& fPlotNum,
                             const string& fTitle, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_denom = (TH2D*)file->Get(fPlotDenom.c_str());
  TH2D *h2_num = (TH2D*)file->Get(fPlotNum.c_str());

  TH1D *h1_bkg = new TH1D("h1_bkg","h1_bkg", 45, 0.5, 45.5);
  TH1D *h1_denom = new TH1D("h1_denom","h1_denom", 45, 0.5, 45.5);
  TH1D *h1_num = new TH1D("h1_num","h1_num", 45, 0.5, 45.5);

  h1_bkg->SetTitle(fTitle.c_str());
  h1_bkg->GetXaxis()->SetLabelSize(0.03);
  h1_bkg->GetYaxis()->SetRangeUser(0., 1.);

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET"; algoMap[9] = "SSVHPT";

  Int_t binCounter = 1;
  
  for(Int_t i=1; i<=9; ++i)
  {
    for(Int_t j=1; j<=i; ++j)
    {
      h1_bkg->GetXaxis()->SetBinLabel(binCounter,(algoMap[i]+"-"+algoMap[j]).c_str());

      h1_denom->SetBinContent(binCounter,h2_denom->GetBinContent(i,j));
      h1_num->SetBinContent(binCounter,h2_num->GetBinContent(i,j));
     
      binCounter++;
    }
  }

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_num, h1_denom,"cp");
  g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(20);
  g_efficiency->SetMarkerColor(kRed);
  
  TCanvas *c = new TCanvas("c", "",1800,800);
  c->cd();

  h1_bkg->Draw();
  g_efficiency->Draw("Psame");

  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}

void bTaggingEfficiency_2Tag_weighted(const string *fInputFiles, const Double_t *fWeights,
                                      const string& fPlotDenom, const string& fPlotNum,
                                      const string& fTitle, const string& fOutputFile
                                     )
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  map<Int_t,string> algoMap;
  algoMap[1] = "TCHEL"; algoMap[2] = "TCHEM"; algoMap[3] = "TCHET"; algoMap[4] = "TCHPL"; algoMap[5] = "TCHPM"; algoMap[6] = "TCHPT"; algoMap[7] = "SSVHEM"; algoMap[8] = "SSVHET"; algoMap[9] = "SSVHPT";

  TH1D *h1_bkg = new TH1D("h1_bkg","h1_bkg", 45, 0.5, 45.5);
  h1_bkg->SetTitle(fTitle.c_str());
  h1_bkg->GetXaxis()->SetLabelSize(0.03);
  h1_bkg->GetYaxis()->SetRangeUser(0., 1.);

  TCanvas *c = new TCanvas("c", "",1800,800);
  c->cd();

  TList *sampleEffficiencies = new TList();
  
  Int_t nFiles = sizeof(fWeights)/sizeof(*fWeights);
  
  for(Int_t n=0; n<nFiles; ++n)
  {
    TFile *file = new TFile(fInputFiles[n].c_str());

    TH2D *h2_denom = (TH2D*)file->Get(fPlotDenom.c_str());
    TH2D *h2_num = (TH2D*)file->Get(fPlotNum.c_str());

    TH1D *h1_denom = new TH1D(Form("h1_denom_%i",n),Form("h1_denom_%i",n), 45, 0.5, 45.5);
    TH1D *h1_num = new TH1D(Form("h1_num_%i",n),Form("h1_num_%i",n), 45, 0.5, 45.5);

    Int_t binCounter = 1;

    for(Int_t i=1; i<=9; ++i)
    {
      for(Int_t j=1; j<=i; ++j)
      {
        if(n==0) h1_bkg->GetXaxis()->SetBinLabel(binCounter,(algoMap[i]+"-"+algoMap[j]).c_str());

        h1_denom->SetBinContent(binCounter,h2_denom->GetBinContent(i,j));
        h1_num->SetBinContent(binCounter,h2_num->GetBinContent(i,j));

        binCounter++;
      }
    }

    TEfficiency *sampleEff = new TEfficiency(*h1_num, *h1_denom);
    sampleEff->SetWeight(fWeights[n]);

    sampleEffficiencies->Add(sampleEff);

    delete h1_denom;
    delete h1_num;
    delete file;
  }

  TEfficiency *finalEff = new TEfficiency();
  
  TGraphAsymmErrors *g_efficiency = finalEff->Combine(sampleEffficiencies,"cl=0.683",nFiles,fWeights);
  g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(20);
  g_efficiency->SetMarkerColor(kRed);

  h1_bkg->Draw();
  g_efficiency->Draw("Psame");

  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
}


void makePlots()
{
// For calculating sample weights in different dijet mass bins
  Double_t Pt170to300Eff[] = {0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt300to470Eff[] = {0.0, 2597.0, 229.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt470to600Eff[] = {13.0, 3581.0, 5138.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt600to800Eff[] = {5.0, 485.0, 7006.0, 1422.0, 11.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt800to1000Eff[] = {1.0, 42.0, 1258.0, 5920.0, 700.0, 5.0, 0.0, 0.0, 0.0};
  Double_t Pt1000to1400Eff[] = {1.0, 5.0, 103.0, 939.0, 2111.0, 346.0, 13.0, 0.0, 0.0};
  Double_t Pt1400to1800Eff[] = {0.0, 1.0, 7.0, 56.0, 347.0, 1149.0, 874.0, 82.0, 1.0};
  Double_t Pt1800Eff[] = {0.0, 0.0, 0.0, 0.0, 5.0, 23.0, 41.0, 108.0, 29.0};

  Double_t Pt170to300Mistag[] = {0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt300to470Mistag[] = {0.0, 541.0, 64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt470to600Mistag[] = {76.0, 1394.0, 901.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt600to800Mistag[] = {18.0, 525.0, 1318.0, 231.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  Double_t Pt800to1000Mistag[] = {4.0, 75.0, 521.0, 733.0, 82.0, 3.0, 0.0, 0.0, 0.0};
  Double_t Pt1000to1400Mistag[] = {2.0, 8.0, 79.0, 170.0, 214.0, 29.0, 1.0, 0.0, 0.0};
  Double_t Pt1400to1800Mistag[] = {0.0, 0.0, 3.0, 34.0, 60.0, 82.0, 54.0, 3.0, 0.0};
  Double_t Pt1800Mistag[] = {0.0, 0.0, 1.0, 1.0, 1.0, 4.0, 6.0, 3.0, 3.0};

  Double_t Pt170to300TotEff = 0, Pt300to470TotEff = 0, Pt470to600TotEff = 0, Pt600to800TotEff = 0, Pt800to1000TotEff = 0, Pt1000to1400TotEff = 0, Pt1400to1800TotEff = 0, Pt1800TotEff = 0;
  Double_t Pt170to300TotMistag = 0, Pt300to470TotMistag = 0, Pt470to600TotMistag = 0, Pt600to800TotMistag = 0, Pt800to1000TotMistag = 0, Pt1000to1400TotMistag = 0, Pt1400to1800TotMistag = 0, Pt1800TotMistag = 0;

  Int_t nBins = sizeof(Pt1800Eff)/sizeof(*Pt1800Eff);
  
  for(Int_t i=0; i<nBins; ++i)
  {
    Pt170to300TotEff+=Pt170to300Eff[i];
    Pt300to470TotEff+=Pt300to470Eff[i];
    Pt470to600TotEff+=Pt470to600Eff[i];
    Pt600to800TotEff+=Pt600to800Eff[i];
    Pt800to1000TotEff+=Pt800to1000Eff[i];
    Pt1000to1400TotEff+=Pt1000to1400Eff[i];
    Pt1400to1800TotEff+=Pt1400to1800Eff[i];
    Pt1800TotEff+=Pt1800Eff[i];

    Pt170to300TotMistag+=Pt170to300Mistag[i];
    Pt300to470TotMistag+=Pt300to470Mistag[i];
    Pt470to600TotMistag+=Pt470to600Mistag[i];
    Pt600to800TotMistag+=Pt600to800Mistag[i];
    Pt800to1000TotMistag+=Pt800to1000Mistag[i];
    Pt1000to1400TotMistag+=Pt1000to1400Mistag[i];
    Pt1400to1800TotMistag+=Pt1400to1800Mistag[i];
    Pt1800TotMistag+=Pt1800Mistag[i];
  }
 
// ###############################
// # Dijet Mass: 500to1000 GeV
// ###############################

// ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_Pt1800.png");

  string files_eff_DijetMass500to1000GeV[] = {
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass500to1000GeV[] = {
      (24300.*(Pt170to300Eff[1]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[1]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[1]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[1]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[1]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[1]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[1]/Pt1400to1800TotEff)), (3.58E-4*(Pt1800Eff[1]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass500to1000GeV, weights_eff_DijetMass500to1000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_Pt1800.png");

  string files_mistag_DijetMass500to1000GeV[] = {
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass500to1000GeV[] = {
      (24300.*(Pt170to300Mistag[1]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[1]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[1]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[1]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[1]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[1]/Pt1000to1400TotMistag))
      //, (0.0109*(Pt1400to1800Mistag[1]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[1]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass500to1000GeV, weights_mistag_DijetMass500to1000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 0.5<M_{jj}<1 TeV", "DoubleTag_Mistag_DijetMass500to1000_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 0.5<M_{jj}<1 TeV", "DoubleTag_Eff_DijetMass500to1000_RSGravitonToBBbar.png");


// ###############################
// # Dijet Mass: 1000to1500 GeV
// ###############################

// ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1800.png");

  string files_eff_DijetMass1000to1500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass1000to1500GeV[] = {
      //(24300.*(Pt170to300Eff[2]/Pt170to300TotEff)),
      (1170.*(Pt300to470Eff[2]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[2]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[2]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[2]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[2]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[2]/Pt1400to1800TotEff))
      //, (3.58E-4*(Pt1800Eff[2]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass1000to1500GeV, weights_eff_DijetMass1000to1500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_QCD_all_weighted.png");

// // #### No DeltaEta cut
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/Final__histograms.root",
//                           "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_all_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt170to300, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt170to300_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt300to470, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt300to470_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt470to600, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt470to600_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt600to800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt600to800_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt800to1000, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt800to1000_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt1000to1400, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1000to1400_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt1400to1800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1400to1800_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                           "DoubleTag Efficiency -- QCD Pt1800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_Pt1800_noDeltaEtaCut.png");
// 
//   string files_eff_DijetMass1000to1500GeV_noDelteEtaCut[] = {
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
//   };
// 
//   Double_t weights_eff_DijetMass1000to1500GeV_noDelteEtaCut[] = {
//       (24300.*(Pt170to300Eff[2]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[2]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[2]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[2]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[2]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[2]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[2]/Pt1400to1800TotEff))
//       //, (3.58E-4*(Pt1800Eff[2]/Pt1800TotEff))
//   };
// 
//   bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass1000to1500GeV_noDelteEtaCut, weights_eff_DijetMass1000to1500GeV_noDelteEtaCut,
//                                    "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
//                                    "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
//                                    "DoubleTag Efficiency -- QCD, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_QCD_all_weighted_noDelteEtaCut.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1800.png");

  string files_mistag_DijetMass1000to1500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass1000to1500GeV[] = {
      //(24300.*(Pt170to300Mistag[2]/Pt170to300TotMistag)),
      (1170.*(Pt300to470Mistag[2]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[2]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[2]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[2]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[2]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[2]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[2]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass1000to1500GeV, weights_mistag_DijetMass1000to1500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 1<M_{jj}<1.5 TeV", "DoubleTag_Mistag_DijetMass1000to1500_QCD_all_weighted.png");

// // #### No DeltaEta cut
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/Final__histograms.root",
//                           "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_all_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt170to300, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt170to300_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt300to470, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt300to470_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt470to600, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt470to600_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt600to800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt600to800_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt800to1000, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt800to1000_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt1000to1400, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1000to1400_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt1400to1800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1400to1800_noDeltaEtaCut.png");
// 
//   bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                           "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                           "DoubleTag Mistag Rate -- QCD Pt1800, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_Pt1800_noDeltaEtaCut.png");
// 
//   string files_mistag_DijetMass1000to1500GeV_noDelteEtaCut[] = {
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//       "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
//   };
// 
//   Double_t weights_mistag_DijetMass1000to1500GeV_noDelteEtaCut[] = {
//       (24300.*(Pt170to300Mistag[2]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[2]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[2]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[2]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[2]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[2]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[2]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[2]/Pt1800TotMistag))
//   };
// 
//   bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass1000to1500GeV_noDelteEtaCut, weights_mistag_DijetMass1000to1500GeV_noDelteEtaCut,
//                                    "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom",
//                                    "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num",
//                                    "DoubleTag Mistag Rate -- QCD, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Mistag_DijetMass1000to1500_QCD_all_weighted_noDelteEtaCut.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_ZprimeToBBbar.png");

// #### No DeltaEta cut

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_ZprimeToBBbar_noDeltaEtaCut.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 1<M_{jj}<1.5 TeV", "DoubleTag_Eff_DijetMass1000to1500_RSGravitonToBBbar.png");

// #### No DeltaEta cut

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_noDeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 1<M_{jj}<1.5 TeV, no #Delta#eta cut", "DoubleTag_Eff_DijetMass1000to1500_RSGravitonToBBbar_noDeltaEtaCut.png");


// ###############################
// # Dijet Mass: 1500to2000 GeV
// ###############################

// ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_Pt1800.png");

  string files_eff_DijetMass1500to2000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass1500to2000GeV[] = {
      //(24300.*(Pt170to300Eff[3]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[3]/Pt300to470TotEff)),
      (70.2*(Pt470to600Eff[3]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[3]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[3]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[3]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[3]/Pt1400to1800TotEff))
      //, (3.58E-4*(Pt1800Eff[3]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass1500to2000GeV, weights_eff_DijetMass1500to2000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_Pt1800.png");

  string files_mistag_DijetMass1500to2000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass1500to2000GeV[] = {
      //(24300.*(Pt170to300Mistag[3]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[3]/Pt300to470TotMistag)),
      (70.2*(Pt470to600Mistag[3]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[3]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[3]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[3]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[3]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[3]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass1500to2000GeV, weights_mistag_DijetMass1500to2000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 1.5<M_{jj}<2 TeV", "DoubleTag_Mistag_DijetMass1500to2000_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 1.5<M_{jj}<2 TeV", "DoubleTag_Eff_DijetMass1500to2000_RSGravitonToBBbar.png");
 

// ###############################
// # Dijet Mass: 2000to2500 GeV
// ###############################

// ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_Pt1800.png");

  string files_eff_DijetMass2000to2500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass2000to2500GeV[] = {
      //(24300.*(Pt170to300Eff[4]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[4]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[4]/Pt470to600TotEff)),
      (15.6*(Pt600to800Eff[4]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[4]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[4]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[4]/Pt1400to1800TotEff)), (3.58E-4*(Pt1800Eff[4]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass2000to2500GeV, weights_eff_DijetMass2000to2500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_Pt1800.png");

  string files_mistag_DijetMass2000to2500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass2000to2500GeV[] = {
      //(24300.*(Pt170to300Mistag[4]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[4]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[4]/Pt470to600TotMistag)),
      (15.6*(Pt600to800Mistag[4]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[4]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[4]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[4]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[4]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass2000to2500GeV, weights_mistag_DijetMass2000to2500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 2<M_{jj}<2.5 TeV", "DoubleTag_Mistag_DijetMass2000to2500_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 2<M_{jj}<2.5 TeV", "DoubleTag_Eff_DijetMass2000to2500_RSGravitonToBBbar.png");


// ###############################
// # Dijet Mass: 2500to3000 GeV
// ###############################

// ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_Pt1800.png");

  string files_eff_DijetMass2500to3000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass2500to3000GeV[] = {
      //(24300.*(Pt170to300Eff[5]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[5]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[5]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[5]/Pt600to800TotEff)),
      (1.84*(Pt800to1000Eff[5]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[5]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[5]/Pt1400to1800TotEff)), (3.58E-4*(Pt1800Eff[5]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass2500to3000GeV, weights_eff_DijetMass2500to3000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_Pt1800.png");

  string files_mistag_DijetMass2500to3000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass2500to3000GeV[] = {
      //(24300.*(Pt170to300Mistag[5]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[5]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[5]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[5]/Pt600to800TotMistag)),
      (1.84*(Pt800to1000Mistag[5]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[5]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[5]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[5]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass2500to3000GeV, weights_mistag_DijetMass2500to3000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 2.5<M_{jj}<3 TeV", "DoubleTag_Mistag_DijetMass2500to3000_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 2.5<M_{jj}<3 TeV", "DoubleTag_Eff_DijetMass2500to3000_RSGravitonToBBbar.png");


// ###############################
// # Dijet Mass: 3000to3500 GeV
// ###############################

// // ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_Pt1800.png");

  string files_eff_DijetMass3000to3500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass3000to3500GeV[] = {
      //(24300.*(Pt170to300Eff[6]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[6]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[6]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[6]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[6]/Pt800to1000TotEff)),
      (0.332*(Pt1000to1400Eff[6]/Pt1000to1400TotEff)), (0.0109*(Pt1400to1800Eff[6]/Pt1400to1800TotEff)), (3.58E-4*(Pt1800Eff[6]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass3000to3500GeV, weights_eff_DijetMass3000to3500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_Pt1800.png");

  string files_mistag_DijetMass3000to3500GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass3000to3500GeV[] = {
      //(24300.*(Pt170to300Mistag[6]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[6]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[6]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[6]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[6]/Pt800to1000TotMistag)),
      (0.332*(Pt1000to1400Mistag[6]/Pt1000to1400TotMistag)), (0.0109*(Pt1400to1800Mistag[6]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[6]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass3000to3500GeV, weights_mistag_DijetMass3000to3500GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 3<M_{jj}<3.5 TeV", "DoubleTag_Mistag_DijetMass3000to3500_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 3<M_{jj}<3.5 TeV", "DoubleTag_Eff_DijetMass3000to3500_RSGravitonToBBbar.png");


// ###############################
// # Dijet Mass: 3500to4000 GeV
// ###############################

// // ## QCD b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt170to300, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt300to470, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt470to600, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt600to800, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt800to1000, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1000to1400, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1400to1800, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- QCD Pt1800, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_Pt1800.png");

  string files_eff_DijetMass3500to4000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_eff_DijetMass3500to4000GeV[] = {
      //(24300.*(Pt170to300Eff[7]/Pt170to300TotEff)), (1170.*(Pt300to470Eff[7]/Pt300to470TotEff)), (70.2*(Pt470to600Eff[7]/Pt470to600TotEff)), (15.6*(Pt600to800Eff[7]/Pt600to800TotEff)), (1.84*(Pt800to1000Eff[7]/Pt800to1000TotEff)), (0.332*(Pt1000to1400Eff[7]/Pt1000to1400TotEff)),
      (0.0109*(Pt1400to1800Eff[7]/Pt1400to1800TotEff)), (3.58E-4*(Pt1800Eff[7]/Pt1800TotEff))
  };

  bTaggingEfficiency_2Tag_weighted(files_eff_DijetMass3500to4000GeV, weights_eff_DijetMass3500to4000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                                   "DoubleTag Efficiency -- QCD, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_QCD_all_weighted.png");

// ## QCD mistag rate

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "QCD_Pythia6__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_all.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt170to300, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt170to300.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt300to470, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt300to470.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt470to600, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt470to600.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt600to800, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt600to800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt800to1000, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt800to1000.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1000to1400, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt1000to1400.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1400to1800, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt1400to1800.png");

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                          "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                          "DoubleTag Mistag Rate -- QCD Pt1800, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_Pt1800.png");

  string files_mistag_DijetMass3500to4000GeV[] = {
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-170to300_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-300to470_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04_skim-v2__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-470to600_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-600to800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-800to1000_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      //"/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
      "/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/QCD_Pt-1800_TuneZ2_7TeV_pythia6__elhughes-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root"
  };

  Double_t weights_mistag_DijetMass3500to4000GeV[] = {
      //(24300.*(Pt170to300Mistag[7]/Pt170to300TotMistag)), (1170.*(Pt300to470Mistag[7]/Pt300to470TotMistag)), (70.2*(Pt470to600Mistag[7]/Pt470to600TotMistag)), (15.6*(Pt600to800Mistag[7]/Pt600to800TotMistag)), (1.84*(Pt800to1000Mistag[7]/Pt800to1000TotMistag)), (0.332*(Pt1000to1400Mistag[7]/Pt1000to1400TotMistag)),
      (0.0109*(Pt1400to1800Mistag[7]/Pt1400to1800TotMistag)), (3.58E-4*(Pt1800Mistag[7]/Pt1800TotMistag))
  };

  bTaggingEfficiency_2Tag_weighted(files_mistag_DijetMass3500to4000GeV, weights_mistag_DijetMass3500to4000GeV,
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom",
                                   "myAnalyzer/h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num",
                                   "DoubleTag Mistag Rate -- QCD, 3.5<M_{jj}<4 TeV", "DoubleTag_Mistag_DijetMass3500to4000_QCD_all_weighted.png");

// ## Z' b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "ZprimeBBbar__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- Z'#rightarrowb#bar{b}, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_ZprimeToBBbar.png");

// ## RSG b-tag efficiency

  bTaggingEfficiency_2Tag("/uscms_data/d3/tote/CMSSW_4_2_8/src/MyAnalysis/MyAnalyzer/test/CRAB_Jobs_bTagging_efficiency_DeltaEtaCut_030212/Final__histograms.root",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom",
                          "RSGravitonJJ__h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num",
                          "DoubleTag Efficiency -- RSG#rightarrowb#bar{b}, 3.5<M_{jj}<4 TeV", "DoubleTag_Eff_DijetMass3500to4000_RSGravitonToBBbar.png");

}
