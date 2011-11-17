#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;

void triggerEfficiency(const string& fInputFile, const string& hNum, const string& hDenom,
                       const string& fTitle, const string& fEra, const string& fOutputFile, const string& fLabel = "",
                       const Int_t rebin = 1, const Double_t xMin = 800, const Double_t xMax = 1300)
{
  gROOT->SetBatch(kTRUE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH1D *h_DijetMass_num = (TH1D*)file->Get(hNum.c_str());
  TH1D *h_DijetMass_denom = (TH1D*)file->Get(hDenom.c_str());

  h_DijetMass_num->Rebin(rebin);
  h_DijetMass_denom->Rebin(rebin);

  Double_t binWidth = h_DijetMass_num->GetBinWidth(h_DijetMass_num->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int(xMin/binWidth);
  
  TH1D *h_Npass = new TH1D("h_Npass","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_Npass->SetBinContent(i,h_DijetMass_num->GetBinContent(i+nBinsToSkip));
    h_Npass->SetBinError(i,h_DijetMass_num->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_Ntotal = new TH1D("h_Ntotal","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_Ntotal->SetBinContent(i,h_DijetMass_denom->GetBinContent(i+nBinsToSkip));
    h_Ntotal->SetBinError(i,h_DijetMass_denom->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h_Npass, h_Ntotal,"cp");
//   g_efficiency->SetMarkerSize(0.8);
  g_efficiency->SetMarkerStyle(21);
  g_efficiency->SetMarkerColor(kRed);
  g_efficiency->SetTitle(fTitle.c_str());
  g_efficiency->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  g_efficiency->GetYaxis()->SetTitle("Trigger Efficiency");
  g_efficiency->GetXaxis()->SetRangeUser(xMin,xMax);
  
  g_efficiency->Draw("AP");

  TLine *line = new TLine(xMin,0.99,xMax,0.99);
  line->SetLineWidth(3.);
  line->SetLineColor(kGreen+2);
  line->SetLineStyle(7);
  line->Draw("same");

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.5,0.40,fEra.c_str());
  l1.DrawLatex(0.5,0.35,"#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.5,0.30,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.5,0.25,fLabel.c_str());
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete line;
  delete h_Npass;
  delete h_Ntotal;
  delete g_efficiency;
  delete c;
  delete file;
}


void makePlots()
{
  // Jet Dataset
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HLTJet370_SSVHEM/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_num",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_denom",
//                     "HLT_Jet370 Efficiency (wrt HLT_Jet240)", "Run2011A", "HLT_Jet370_efficiency.png");
//   
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HLTJet370_SSVHEM/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_SingleTag",
//                     "DATA__h1_DijetMass_denom_SingleTag",
//                     "HLT_Jet370 Efficiency (wrt HLT_Jet240)", "Run2011A", "HLT_Jet370_efficiency_SSVHEM_SingleTag.png",
//                     "SSVHEM SingleTag",2);
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HLTJet370_SSVHEM/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_DoubleTag",
//                     "DATA__h1_DijetMass_denom_DoubleTag",
//                     "HLT_Jet370_efficiency_SSVHEM_DoubleTag.png", "Run2011A",
//                     "SSVHEM DoubleTag",4);
// 
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HLTJet370_TCHEM/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_SingleTag",
//                     "DATA__h1_DijetMass_denom_SingleTag",
//                     "HLT_Jet370 Efficiency (wrt HLT_Jet240)", "Run2011A", "HLT_Jet370_efficiency_TCHEM_SingleTag.png",
//                     "TCHEM SingleTag",2);
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HLTJet370_TCHEM/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_DoubleTag",
//                     "DATA__h1_DijetMass_denom_DoubleTag",
//                     "HLT_Jet370 Efficiency (wrt HLT_Jet240)", "Run2011A", "HLT_Jet370_efficiency_TCHEM_DoubleTag.png",
//                     "TCHEM DoubleTag",4);

//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_HT600",
//                     "DATA__h1_DijetMass_denom_HT600",
//                     "HLT_HT600 Efficiency (wrt HLT_Jet240)", "Run2011A", "HLT_HT600_efficiency.png","",2, 700, 1300);
// 
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_HT650",
//                     "DATA__h1_DijetMass_denom_HT650",
//                     "HLT_HT650 Efficiency (wrt HLT_Jet240)", "Run2011B", "HLT_HT650_efficiency.png","",2, 700, 1300);

//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HT750/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_HT750",
//                     "DATA__h1_DijetMass_denom_HT750",
//                     "HLT_HT750 Efficiency (wrt HLT_Jet240)", "Run2011B", "HLT_HT750_efficiency.png","",2, 700, 1300);
    
//   triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
//                     "DATA__h1_DijetMass_num_FatJet850",
//                     "DATA__h1_DijetMass_denom_FatJet850",
//                     "HLT_FatJetMass850 Efficiency (wrt HLT_Jet240)", "Run2011A+Run2011B", "HLT_FatJetMass850_efficiency.png","",2, 700, 1300);

  // HT Dataset
  triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT600",
                    "DATA__h1_DijetMass_denom_HT600",
                    "HLT_HT600 Efficiency (wrt HLT_HT400)", "", "HLT_HT600_efficiency_HTDataset.png","",20, 700, 1300);

  triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT650",
                    "DATA__h1_DijetMass_denom_HT650",
                    "HLT_HT650 Efficiency (wrt HLT_HT400)", "", "HLT_HT650_efficiency_HTDataset.png","",20, 700, 1300);

  triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
                    "DATA__h1_DijetMass_num_HT750",
                    "DATA__h1_DijetMass_denom_HT750",
                    "HLT_HT750 Efficiency (wrt HLT_HT450)", "", "HLT_HT750_efficiency_HTDataset.png","",20, 700, 1300);

  triggerEfficiency("CRAB_Jobs_DijetBBTag_TriggerEfficiency_HTFatJet/Final__histograms.root",
                    "DATA__h1_DijetMass_num_FatJet850",
                    "DATA__h1_DijetMass_denom_FatJet850",
                    "HLT_FatJetMass850 Efficiency (wrt HLT_HT450)", "", "HLT_FatJetMass850_efficiency_HTDataset.png","",20, 700, 1300);
}
