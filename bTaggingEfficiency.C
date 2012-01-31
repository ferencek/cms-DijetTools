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

using namespace std;

void bTaggingEfficiencyVsPtJ1(const string& fInputFile, const string& fSample1, const string& fSample2, const string& fAlgo,
                              const string& fTitle, const string& fLabel, const string& fOutputFile,
                              const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle)
{
  gROOT->SetBatch(kTRUE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH2D *h2_denom_Sample1 = (TH2D*)file->Get((fSample1 + "__h2_EtaJ1_vs_PtJ1_HF").c_str());
  TH2D *h2_num_Sample1 = (TH2D*)file->Get((fSample1 + "__h2_EtaJ1_vs_PtJ1_" + fAlgo).c_str());

  TH1D *h1_denom_Sample1 = h2_denom_Sample1->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample1 = h2_num_Sample1->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample1_rebinned = (TH1D*)h1_denom_Sample1->Rebin(nbins,"h1_denom_Sample1_rebinned",xbins);
  TH1D *h1_num_Sample1_rebinned = (TH1D*)h1_num_Sample1->Rebin(nbins,"h1_num_Sample1_rebinned",xbins);

  TH2D *h2_denom_Sample2 = (TH2D*)file->Get((fSample2 + "__h2_EtaJ1_vs_PtJ1_HF").c_str());
  TH2D *h2_num_Sample2 = (TH2D*)file->Get((fSample2 + "__h2_EtaJ1_vs_PtJ1_" + fAlgo).c_str());

  TH1D *h1_denom_Sample2 = h2_denom_Sample2->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample2 = h2_num_Sample2->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample2_rebinned = (TH1D*)h1_denom_Sample2->Rebin(nbins,"h1_denom_Sample2_rebinned",xbins);
  TH1D *h1_num_Sample2_rebinned = (TH1D*)h1_num_Sample2->Rebin(nbins,"h1_num_Sample2_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency_Sample1 = new TGraphAsymmErrors(h1_num_Sample1_rebinned, h1_denom_Sample1_rebinned,"cp");
//   g_efficiency_Sample1->SetMarkerSize(0.8);
  g_efficiency_Sample1->SetMarkerStyle(20);
  g_efficiency_Sample1->SetMarkerColor(kBlack);
  g_efficiency_Sample1->SetTitle(fTitle.c_str());
  g_efficiency_Sample1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample1->GetXaxis()->SetRangeUser(xMin,xMax);
  g_efficiency_Sample1->GetYaxis()->SetRangeUser(0,1.);
  
  TGraphAsymmErrors *g_efficiency_Sample2 = new TGraphAsymmErrors(h1_num_Sample2_rebinned, h1_denom_Sample2_rebinned,"cp");
//   g_efficiency_Sample2->SetMarkerSize(0.8);
  g_efficiency_Sample2->SetMarkerStyle(24);
  g_efficiency_Sample2->SetMarkerColor(kRed);
  g_efficiency_Sample2->SetLineColor(kRed);
  g_efficiency_Sample2->SetTitle(fTitle.c_str());
  g_efficiency_Sample2->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample2->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample2->GetXaxis()->SetRangeUser(xMin,xMax);

  g_efficiency_Sample1->Draw("AP");
  g_efficiency_Sample2->Draw("Psame");

  TLegend *legend = new TLegend(.6,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(g_efficiency_Sample1,fSample1.c_str(),"lp");
  legend->AddEntry(g_efficiency_Sample2,fSample2.c_str(),"lp");
  legend->Draw();
  
//   TLine *line = new TLine(xMin,0.99,xMax,0.99);
//   line->SetLineWidth(3.);
//   line->SetLineColor(kGreen+2);
//   line->SetLineStyle(7);
//   line->Draw("same");

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
//   l1.DrawLatex(0.5,0.35,"#sqrt{s} = 7 TeV");
//   l1.DrawLatex(0.5,0.30,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.2,0.8,fLabel.c_str());
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

//   delete line;
//   delete h_Npass;
//   delete h_Ntotal;
//   delete g_efficiency_Sample1;
  delete c;
  delete file;
}


void bTaggingEfficiencyVsPtJ1_ext(const string& fInputFile1, const string& fInputFile2, const string& fSample1, const string& fSample2, const string& fAlgo,
                                  const string& fTitle, const string& fLabel, const string& fLegend1, const string& fLegend2, const string& fOutputFile,
                                  const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle)
{
  gROOT->SetBatch(kTRUE);

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH2D *h2_denom_Sample1 = (TH2D*)file1->Get((fSample1 + "__h2_EtaJ1_vs_PtJ1_HF").c_str());
  TH2D *h2_num_Sample1 = (TH2D*)file1->Get((fSample1 + "__h2_EtaJ1_vs_PtJ1_" + fAlgo).c_str());

  TH1D *h1_denom_Sample1 = h2_denom_Sample1->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample1 = h2_num_Sample1->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample1_rebinned = (TH1D*)h1_denom_Sample1->Rebin(nbins,"h1_denom_Sample1_rebinned",xbins);
  TH1D *h1_num_Sample1_rebinned = (TH1D*)h1_num_Sample1->Rebin(nbins,"h1_num_Sample1_rebinned",xbins);

  TH2D *h2_denom_Sample2 = (TH2D*)file2->Get((fSample2 + "__h2_EtaJ1_vs_PtJ1_HF").c_str());
  TH2D *h2_num_Sample2 = (TH2D*)file2->Get((fSample2 + "__h2_EtaJ1_vs_PtJ1_" + fAlgo).c_str());

  TH1D *h1_denom_Sample2 = h2_denom_Sample2->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample2 = h2_num_Sample2->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample2_rebinned = (TH1D*)h1_denom_Sample2->Rebin(nbins,"h1_denom_Sample2_rebinned",xbins);
  TH1D *h1_num_Sample2_rebinned = (TH1D*)h1_num_Sample2->Rebin(nbins,"h1_num_Sample2_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency_Sample1 = new TGraphAsymmErrors(h1_num_Sample1_rebinned, h1_denom_Sample1_rebinned,"cp");
//   g_efficiency_Sample1->SetMarkerSize(0.8);
  g_efficiency_Sample1->SetMarkerStyle(20);
  g_efficiency_Sample1->SetMarkerColor(kBlack);
  g_efficiency_Sample1->SetTitle(fTitle.c_str());
  g_efficiency_Sample1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample1->GetXaxis()->SetRangeUser(xMin,xMax);
  g_efficiency_Sample1->GetYaxis()->SetRangeUser(0,1.);

  TGraphAsymmErrors *g_efficiency_Sample2 = new TGraphAsymmErrors(h1_num_Sample2_rebinned, h1_denom_Sample2_rebinned,"cp");
//   g_efficiency_Sample2->SetMarkerSize(0.8);
  g_efficiency_Sample2->SetMarkerStyle(24);
  g_efficiency_Sample2->SetMarkerColor(kRed);
  g_efficiency_Sample2->SetLineColor(kRed);
  g_efficiency_Sample2->SetTitle(fTitle.c_str());
  g_efficiency_Sample2->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample2->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample2->GetXaxis()->SetRangeUser(xMin,xMax);

  g_efficiency_Sample1->Draw("AP");
  g_efficiency_Sample2->Draw("Psame");

  TLegend *legend = new TLegend(.55,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(g_efficiency_Sample1,fLegend1.c_str(),"lp");
  legend->AddEntry(g_efficiency_Sample2,fLegend2.c_str(),"lp");
  legend->Draw();

//   TLine *line = new TLine(xMin,0.99,xMax,0.99);
//   line->SetLineWidth(3.);
//   line->SetLineColor(kGreen+2);
//   line->SetLineStyle(7);
//   line->Draw("same");

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
//   l1.DrawLatex(0.5,0.35,"#sqrt{s} = 7 TeV");
//   l1.DrawLatex(0.5,0.30,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.2,0.8,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

//   delete line;
//   delete h_Npass;
//   delete h_Ntotal;
//   delete g_efficiency_Sample1;
  delete c;
//   delete file;
}


void bTaggingEfficiencyVsPtJ2_ext(const string& fInputFile1, const string& fInputFile2, const string& fSample1, const string& fSample2, const string& fAlgo,
                                  const string& fTitle, const string& fLabel, const string& fLegend1, const string& fLegend2, const string& fOutputFile,
                                  const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle)
{
  gROOT->SetBatch(kTRUE);

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH2D *h2_denom_Sample1 = (TH2D*)file1->Get((fSample1 + "__h2_EtaJ2_vs_PtJ2_HF").c_str());
  TH2D *h2_num_Sample1 = (TH2D*)file1->Get((fSample1 + "__h2_EtaJ2_vs_PtJ2_" + fAlgo).c_str());

  TH1D *h1_denom_Sample1 = h2_denom_Sample1->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample1 = h2_num_Sample1->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample1_rebinned = (TH1D*)h1_denom_Sample1->Rebin(nbins,"h1_denom_Sample1_rebinned",xbins);
  TH1D *h1_num_Sample1_rebinned = (TH1D*)h1_num_Sample1->Rebin(nbins,"h1_num_Sample1_rebinned",xbins);

  TH2D *h2_denom_Sample2 = (TH2D*)file2->Get((fSample2 + "__h2_EtaJ2_vs_PtJ2_HF").c_str());
  TH2D *h2_num_Sample2 = (TH2D*)file2->Get((fSample2 + "__h2_EtaJ2_vs_PtJ2_" + fAlgo).c_str());

  TH1D *h1_denom_Sample2 = h2_denom_Sample2->ProjectionX("_px",48,52);
  TH1D *h1_num_Sample2 = h2_num_Sample2->ProjectionX("_px",48,52);

  TH1D *h1_denom_Sample2_rebinned = (TH1D*)h1_denom_Sample2->Rebin(nbins,"h1_denom_Sample2_rebinned",xbins);
  TH1D *h1_num_Sample2_rebinned = (TH1D*)h1_num_Sample2->Rebin(nbins,"h1_num_Sample2_rebinned",xbins);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TGraphAsymmErrors *g_efficiency_Sample1 = new TGraphAsymmErrors(h1_num_Sample1_rebinned, h1_denom_Sample1_rebinned,"cp");
//   g_efficiency_Sample1->SetMarkerSize(0.8);
  g_efficiency_Sample1->SetMarkerStyle(20);
  g_efficiency_Sample1->SetMarkerColor(kBlack);
  g_efficiency_Sample1->SetTitle(fTitle.c_str());
  g_efficiency_Sample1->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample1->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample1->GetXaxis()->SetRangeUser(xMin,xMax);
  g_efficiency_Sample1->GetYaxis()->SetRangeUser(0,1.);

  TGraphAsymmErrors *g_efficiency_Sample2 = new TGraphAsymmErrors(h1_num_Sample2_rebinned, h1_denom_Sample2_rebinned,"cp");
//   g_efficiency_Sample2->SetMarkerSize(0.8);
  g_efficiency_Sample2->SetMarkerStyle(24);
  g_efficiency_Sample2->SetMarkerColor(kRed);
  g_efficiency_Sample2->SetLineColor(kRed);
  g_efficiency_Sample2->SetTitle(fTitle.c_str());
  g_efficiency_Sample2->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  g_efficiency_Sample2->GetYaxis()->SetTitle(fYAxisTitle.c_str());
//   g_efficiency_Sample2->GetXaxis()->SetRangeUser(xMin,xMax);

  g_efficiency_Sample1->Draw("AP");
  g_efficiency_Sample2->Draw("Psame");

  TLegend *legend = new TLegend(.55,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(g_efficiency_Sample1,fLegend1.c_str(),"lp");
  legend->AddEntry(g_efficiency_Sample2,fLegend2.c_str(),"lp");
  legend->Draw();

//   TLine *line = new TLine(xMin,0.99,xMax,0.99);
//   line->SetLineWidth(3.);
//   line->SetLineColor(kGreen+2);
//   line->SetLineStyle(7);
//   line->Draw("same");

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
//   l1.DrawLatex(0.5,0.35,"#sqrt{s} = 7 TeV");
//   l1.DrawLatex(0.5,0.30,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.2,0.8,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

//   delete line;
//   delete h_Npass;
//   delete h_Ntotal;
//   delete g_efficiency_Sample1;
  delete c;
//   delete file;
}


void makePlots()
{

  Double_t PtBins[] = {0, 200, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 2000, 3000, 4000, 5000, 6000};

//   bTaggingEfficiencyVsPtJ1("CRAB_Jobs_bTaggingEfficiency_bcHadronMatching/Final__histograms.root",
//                            "QCD_Pythia6", "ZprimeBBbar" , "TCHEM",
//                            "TCHEM, b or c hadron matching", "|#eta|<0.2",
//                            "b-tag_eff_PtJ1_etalt0p2_QCD_Pythia6_ZprimeBBbar_bcHadronMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");
// 
//   bTaggingEfficiencyVsPtJ1("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
//                            "QCD_Pythia6", "ZprimeBBbar" , "TCHEM",
//                            "TCHEM, b parton matching", "|#eta|<0.2",
//                            "b-tag_eff_PtJ1_etalt0p2_QCD_Pythia6_ZprimeBBbar_bPartonMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");
//   
//   bTaggingEfficiencyVsPtJ1("CRAB_Jobs_bTaggingEfficiency_bcHadronMatching/Final__histograms.root",
//                            "QCD_Pythia6", "QCD_MadGraph" , "TCHEM",
//                            "TCHEM, b or c hadron matching", "|#eta|<0.2",
//                            "b-tag_eff_PtJ1_etalt0p2_QCD_Pythia6_MadGraph_bcHadronMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");
// 
//   bTaggingEfficiencyVsPtJ1("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
//                            "QCD_Pythia6", "QCD_MadGraph" , "TCHEM",
//                            "TCHEM, b parton matching", "|#eta|<0.2",
//                            "b-tag_eff_PtJ1_etalt0p2_QCD_Pythia6_MadGraph_bPartonMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");

//   bTaggingEfficiencyVsPtJ1_ext("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
//                                "CRAB_Jobs_bTaggingEfficiency_QCD_FCR_bPartonMatching/Final__histograms.root",
//                                "QCD_Pythia6", "QCD_Pythia6" , "TCHEM",
//                                "TCHEM, b parton matching", "|#eta|<0.2", "QCD Pythia6", "QCD Pythia6 FCR",
//                                "b-tag_eff_PtJ1_etalt0p2_QCD_FCR_Pythia6_bPartonMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");
// 
//   bTaggingEfficiencyVsPtJ1_ext("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
//                                "CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
//                                "ZprimeBBbar", "RSGravitonJJ" , "TCHEM",
//                                "TCHEM, b parton matching", "|#eta|<0.2", "Z'#rightarrowb#bar{b}", "RS Graviton",
//                                "b-tag_eff_PtJ1_etalt0p2_ZprimeBBbar_RSGraviton_bPartonMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");

  bTaggingEfficiencyVsPtJ1_ext("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
                               "CRAB_Jobs_bTaggingEfficiency_QCD_GSP_bPartonMatching/Final__histograms.root",
                               "QCD_Pythia6", "QCD_Pythia6" , "TCHEM",
                               "TCHEM, b parton matching", "|#eta|<0.2", "QCD Pythia6", "QCD Pythia6 GSP",
                               "b-tag_eff_PtJ1_etalt0p2_QCD_GSP_Pythia6_bPartonMatching.png", 15, PtBins, "p_{T,1} [GeV]", "b-tag efficiency");

  bTaggingEfficiencyVsPtJ2_ext("CRAB_Jobs_bTaggingEfficiency_bPartonMatching/Final__histograms.root",
                               "CRAB_Jobs_bTaggingEfficiency_QCD_GSP_bPartonMatching/Final__histograms.root",
                               "QCD_Pythia6", "QCD_Pythia6" , "TCHEM",
                               "TCHEM, b parton matching", "|#eta|<0.2", "QCD Pythia6", "QCD Pythia6 GSP",
                               "b-tag_eff_PtJ2_etalt0p2_QCD_GSP_Pythia6_bPartonMatching.png", 15, PtBins, "p_{T,2} [GeV]", "b-tag efficiency");
  
}
