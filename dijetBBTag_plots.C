#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tdrstyle.C"


using namespace std;


void DATA_MC_scale_factor(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC,
                          const Int_t nbins, const Double_t *xbins, const Double_t xMin = 890, const Double_t xMax = 6000)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file1->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file2->Get(hMC.c_str());

  Double_t binWidth = h_DijetMass_DATA->GetBinWidth(h_DijetMass_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int(xMin/binWidth);

  TH1D *h_DijetMass_DATA_range = new TH1D("h_DijetMass_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_DATA_range->SetBinContent(i,h_DijetMass_DATA->GetBinContent(i+nBinsToSkip));
    h_DijetMass_DATA_range->SetBinError(i,h_DijetMass_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_MC_range = new TH1D("h_DijetMass_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_MC_range->SetBinContent(i,h_DijetMass_MC->GetBinContent(i+nBinsToSkip));
    h_DijetMass_MC_range->SetBinError(i,h_DijetMass_MC->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_DATA_rebin = (TH1D*)h_DijetMass_DATA_range->Rebin(nbins,"h_DijetMass_DATA",xbins);
  TH1D *h_DijetMass_MC_rebin = (TH1D*)h_DijetMass_MC_range->Rebin(nbins,"h_DijetMass_MC",xbins);

  double data = h_DijetMass_DATA_rebin->Integral(1, h_DijetMass_DATA_rebin->GetNbinsX()+1);
  double mc = h_DijetMass_MC_rebin->Integral(1, h_DijetMass_MC_rebin->GetNbinsX()+1);
  cout <<"Data = " << data << endl;
  cout <<"MC = " << mc << endl;
  
  double scaleFactor = data / mc;
  cout <<"DATA/MC scale factor = " << scaleFactor << endl;

  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete file1;
  delete file2;
}


void overlay_DATA_MC(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax, const Double_t yMin = 0, const Double_t yMax = 0,
                     const string& fLogy = "", const Int_t fRebin = 1)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_DATA = (TH1D*)file1->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file2->Get(hMC.c_str());

  if(fRebin>1)
  {
    h_DATA->Rebin(fRebin);
    h_MC->Rebin(fRebin);
  }

  Double_t binWidth = h_DATA->GetBinWidth(h_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DATA_range = new TH1D("h_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DATA_range->SetBinContent(i,h_DATA->GetBinContent(i+nBinsToSkip));
    h_DATA_range->SetBinError(i,h_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_MC_range = new TH1D("h_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC_range->SetBinContent(i,h_MC->GetBinContent(i+nBinsToSkip));
    h_MC_range->SetBinError(i,h_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h_DATA_range->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DATA_range->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(yMin!=yMax) h_DATA_range->GetYaxis()->SetRangeUser(yMin,yMax);
  h_DATA_range->SetMarkerStyle(20);

  h_MC_range->SetLineColor(42);
  h_MC_range->SetFillColor(42);
  h_MC_range->SetMarkerStyle(25);
//   h_MC_range->SetMarkerColor(kRed);

  h_DATA_range->Draw();
  h_MC_range->Draw("histsame");
  h_DATA_range->Draw("same");

  gPad->RedrawAxis();
  
  TLegend *legend = new TLegend(.35,.74,.7,.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_DATA_range,"Data","lp");
  legend->AddEntry(h_MC_range,"QCD PYTHIA6 MC (#times 1.16)","f");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.42, "CMS Preliminary");
  l1.DrawLatex(0.19,0.34, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.19, "Wide Jets, M_{jj}>890 GeV");
  else l1.DrawLatex(0.19,0.19, "Anti-k_{T} R = 0.7 PFJets, M_{jj}>890 GeV");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.74,0.19, fLabel.c_str());

  if(fLogy=="Logy") c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file1;
  delete file2;
}


void overlay_MC(const string& fInputFile1, const string& fInputFile2, const string& hMC1, const string& hMC2, const string& fOutputFile, const string& fLabel, const string& fLegend,
                     const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax, const Double_t yMin = 0, const Double_t yMax = 0,
                     const string& fLogy = "", const Int_t fRebin = 1)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_MC1 = (TH1D*)file1->Get(hMC1.c_str());
  TH1D *h_MC2 = (TH1D*)file2->Get(hMC2.c_str());

  if(fRebin>1)
  {
    h_MC1->Rebin(fRebin);
    h_MC2->Rebin(fRebin);
  }

  Double_t binWidth = h_MC1->GetBinWidth(h_MC1->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_MC1->GetBinLowEdge(1))/binWidth);

  TH1D *h_MC1_range = new TH1D("h_MC1_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC1_range->SetBinContent(i,h_MC1->GetBinContent(i+nBinsToSkip));
    h_MC1_range->SetBinError(i,h_MC1->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_MC2_range = new TH1D("h_MC2_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC2_range->SetBinContent(i,h_MC2->GetBinContent(i+nBinsToSkip));
    h_MC2_range->SetBinError(i,h_MC2->GetBinError(i+nBinsToSkip));
  }

  h_MC1_range->Scale(1./h_MC1_range->Integral());
  h_MC2_range->Scale(1./h_MC2_range->Integral());
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h_MC1_range->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_MC1_range->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(yMin!=yMax) h_MC1_range->GetYaxis()->SetRangeUser(yMin,yMax);
  
  h_MC1_range->SetLineColor(kRed);
  h_MC1_range->SetMarkerStyle(24);
  h_MC1_range->SetMarkerColor(kRed);
  
  h_MC2_range->SetLineColor(kBlue);
  h_MC2_range->SetMarkerStyle(26);
  h_MC2_range->SetMarkerColor(kBlue);

  h_MC1_range->Draw();
  h_MC2_range->Draw("same");

  gPad->RedrawAxis();

  TLegend *legend = new TLegend(.28,.74,.63,.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_MC1_range,"QCD PYTHIA6 MC, M_{jj}>890 GeV","lp");
  legend->AddEntry(h_MC2_range,fLegend.c_str(),"lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.34, "CMS Simulation");
  l1.DrawLatex(0.20,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.19, "Wide Jets");
  else l1.DrawLatex(0.19,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.74,0.19, fLabel.c_str());

  if(fLogy=="Logy") c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_MC1_range;
  delete h_MC2_range;
  delete c;
  delete file1;
  delete file2;
}


void overlay_DATA_MC_v(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                       const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle,
                       const Double_t xMin = 1050, const Double_t xMax = 3500, const Double_t yMin = 0, const Double_t yMax = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file1->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file2->Get(hMC.c_str());

  Double_t binWidth = h_DijetMass_DATA->GetBinWidth(h_DijetMass_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DijetMass_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DijetMass_DATA_range = new TH1D("h_DijetMass_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_DATA_range->SetBinContent(i,h_DijetMass_DATA->GetBinContent(i+nBinsToSkip));
    h_DijetMass_DATA_range->SetBinError(i,h_DijetMass_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_MC_range = new TH1D("h_DijetMass_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_MC_range->SetBinContent(i,h_DijetMass_MC->GetBinContent(i+nBinsToSkip));
    h_DijetMass_MC_range->SetBinError(i,h_DijetMass_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH1D *h_DijetMass_DATA_rebin = (TH1D*)h_DijetMass_DATA_range->Rebin(nbins,"h_DijetMass_DATA",xbins);
  TH1D *h_DijetMass_MC_rebin = (TH1D*)h_DijetMass_MC_range->Rebin(nbins,"h_DijetMass_MC",xbins);

  h_DijetMass_DATA_rebin->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DijetMass_DATA_rebin->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  if(yMin!=yMax) h_DijetMass_DATA_rebin->GetYaxis()->SetRangeUser(yMin,yMax);
  
  h_DijetMass_DATA_rebin->SetMarkerStyle(20);

  h_DijetMass_MC_rebin->SetLineColor(42);
  h_DijetMass_MC_rebin->SetFillColor(42);
//   h_DijetMass_MC_rebin->SetMarkerStyle(25);
//   h_DijetMass_MC_rebin->SetMarkerColor(kRed);
  
  h_DijetMass_DATA_rebin->Draw();
  h_DijetMass_MC_rebin->Draw("histsame");
  h_DijetMass_DATA_rebin->Draw("same");

  gPad->RedrawAxis();
  
  TLegend *legend = new TLegend(.35,.74,.7,.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_DijetMass_DATA_rebin,"Data","lp");
  legend->AddEntry(h_DijetMass_MC_rebin,"QCD PYTHIA6 MC (#times 1.16)","f");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.42, "CMS Preliminary");
  l1.DrawLatex(0.19,0.34, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.29, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.24, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.19, "Wide Jets");
  else l1.DrawLatex(0.19,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.74,0.19, fLabel.c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_rebin;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file1;
  delete file2;
}


void data_MC_ratio(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                   const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_DATA = (TH1D*)file1->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file2->Get(hMC.c_str());

  Double_t binWidth = h_DATA->GetBinWidth(h_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int((xMin-h_DATA->GetBinLowEdge(1))/binWidth);

  TH1D *h_DATA_range = new TH1D("h_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DATA_range->SetBinContent(i,h_DATA->GetBinContent(i+nBinsToSkip));
    h_DATA_range->SetBinError(i,h_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_MC_range = new TH1D("h_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_MC_range->SetBinContent(i,h_MC->GetBinContent(i+nBinsToSkip));
    h_MC_range->SetBinError(i,h_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h_DATA_range->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DATA_range->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h_DATA_range->SetMarkerStyle(20);

  h_MC_range->SetLineColor(kRed);
  h_MC_range->SetMarkerStyle(25);
  h_MC_range->SetMarkerColor(kRed);

  h_DATA_range->Divide(h_MC_range);
  h_DATA_range->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.87, "CMS Preliminary");
  l1.DrawLatex(0.19,0.79, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.74, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.69, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.64, "Wide Jets");
  else l1.DrawLatex(0.19,0.64, "Anti-k_{T} R = 0.7 PFJets");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.74,0.19, fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file1;
  delete file2;
}


void data_MC_ratio_v(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin = 1050, const Double_t xMax = 3500,
                     const Double_t fUp = 0., const Double_t fDown = 0.)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.93);
  gROOT->ForceStyle();

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file1->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file2->Get(hMC.c_str());

  Double_t binWidth = h_DijetMass_DATA->GetBinWidth(h_DijetMass_DATA->GetXaxis()->FindBin(xMin));
  Int_t nBins = int((xMax-xMin)/binWidth);
  Int_t nBinsToSkip = int(xMin/binWidth);
  
  TH1D *h_DijetMass_DATA_range = new TH1D("h_DijetMass_DATA_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_DATA_range->SetBinContent(i,h_DijetMass_DATA->GetBinContent(i+nBinsToSkip));
    h_DijetMass_DATA_range->SetBinError(i,h_DijetMass_DATA->GetBinError(i+nBinsToSkip));
  }

  TH1D *h_DijetMass_MC_range = new TH1D("h_DijetMass_MC_range","",nBins,xMin,xMax);
  for(Int_t i=1; i<=nBins; ++i)
  {
    h_DijetMass_MC_range->SetBinContent(i,h_DijetMass_MC->GetBinContent(i+nBinsToSkip));
    h_DijetMass_MC_range->SetBinError(i,h_DijetMass_MC->GetBinError(i+nBinsToSkip));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH1D *h_DijetMass_DATA_MC = (TH1D*)h_DijetMass_DATA_range->Rebin(nbins,"h_DijetMass_DATA_MC",xbins);
  TH1D *h_DijetMass_MC_rebin = (TH1D*)h_DijetMass_MC_range->Rebin(nbins,"h_DijetMass_MC",xbins);

  h_DijetMass_DATA_MC->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h_DijetMass_DATA_MC->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h_DijetMass_DATA_MC->GetYaxis()->SetRangeUser(0.,2.4);
  
  h_DijetMass_DATA_MC->Divide(h_DijetMass_MC_rebin);

  Double_t mass_array[5] = {xMin, xMax, xMax, xMin};
  Double_t f_array[5] = {fDown, fDown, fUp, fUp};

  TGraph *g_f = new TGraph(4, mass_array, f_array);
  g_f->SetFillColor(kYellow);

  TF1 *fit = new TF1("fit","pol0",xMin,xMax);
  fit->SetParameter(0, 1.);
  fit->SetLineWidth(2);
  fit->SetLineColor(kRed);
  
  h_DijetMass_DATA_MC->Fit("fit","R");

  h_DijetMass_DATA_MC->Draw();
  if(fUp!=fDown)
  {
    g_f->Draw("Fsame");
    h_DijetMass_DATA_MC->Draw("same");
  }
   
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.87, "CMS Preliminary");
  l1.DrawLatex(0.19,0.79, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.74, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.69, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.64, "Wide Jets");
  else l1.DrawLatex(0.19,0.64, "Anti-k_{T} R = 0.7 PFJets");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.74,0.19, fLabel.c_str());
  
  gPad->RedrawAxis();
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_MC;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file1;
  delete file2;
}


void makePlots()
{
  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231,
                      2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,
                      5058, 5253, 5455, 5663, 5877, 6000};

//###################################################################################################################
//##
//## Full 2011 dataset (4.976 /fb)
//##
//###################################################################################################################
                      
//###################################################################################################################
//##
//## Wide jets
//##
//###################################################################################################################
  
//###################################################################################################################
//## PU reweighting applied
//###################################################################################################################

//   // DATA/MC scale factor
//   DATA_MC_scale_factor("CRAB_Jobs_MainAnalysis_CSVL_PUReweighted_PartonMatching_WideJets_EventBins/Final__histograms.root",
//                        "CRAB_Jobs_MainAnalysis_CSVL_PUReweighted_PartonMatching_WideJets_EventBins/Final__histograms.root",
//                        "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                        "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                        43, xbins, 890, 6000);

  // Dijet mass plots
  overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                    "DijetMass_PUReweighted_WideJets_Full2011.eps", "", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
                  "DijetMass_ratio_PUReweighted_WideJets_Full2011.eps", "", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//   // Primary vertex plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                   "nGoodVertices_PUReweighted_Full2011.eps", "", "Good Vertex Multiplicity", "Events", -0.5, 22.5, 0, 1.5e+05);

//   data_MC_ratio("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                 "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                 "DATA__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                 "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                 "nGoodVertices_ratio_PUReweighted_Full2011.eps", "", "Good Vertex Multiplicity", "DATA/MC", -0.5, 22.5);

//   // MET/SumET
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET_pretag",
//                   "METoSumET_PUReweighted_Full2011.eps", "", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 0.1, 3e6, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//                   "DeltaPhiJ1J2_PUReweighted_Full2011.eps", "", "|#Delta#phi|", "Events", 0, 3.15, 0.1, 6e5, "Logy");
// 
//   // |DeltaEta(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//                   "DeltaEtaJ1J2_PUReweighted_Full2011.eps", "", "|#Delta#eta|", "Events", 0, 1.5, 0.1, 8e4);
// 
//   // Pt plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ1_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1_pretag",
//                   "PtJ1_PUReweighted_Full2011.eps", "", "p_{T,1} [GeV]", "Events", 0, 2500, 0.1, 2.0e6, "Logy", 40);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2_pretag",
//                   "PtJ2_PUReweighted_Full2011.eps", "", "p_{T,2} [GeV]", "Events", 0, 2500, 0.1, 2.0e6, "Logy", 40);
// 
//   // eta plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ1_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ1_pretag",
//                   "EtaJ1_PUReweighted_Full2011.eps", "", "#eta_{1}", "Events", -3, 3, 0, 6e4);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ2_pretag",
//                   "EtaJ2_PUReweighted_Full2011.eps", "", "#eta_{2}", "Events", -3, 3, 0, 6e4);
// 
//   // phi plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ1_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ1_pretag",
//                   "PhiJ1_PUReweighted_Full2011.eps", "", "#phi_{1}", "Events", -3.15, 3.15, 1e4, 6e4, "", 10);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ2_pretag",
//                   "PhiJ2_PUReweighted_Full2011.eps", "", "#phi_{2}", "Events", -3.15, 3.15, 1e4, 6e4, "", 10);
  
//   // Muon multiplicity plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nMuons_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons_pretag",
//                   "nMuons_PUReweighted_Full2011.eps", "", "Muon Multiplicity", "Events", -0.5, 5.5, 0.1, 2e6, "Logy", 1);
// 
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass_pretag",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                     "nMuons_vs_DijetMass_PUReweighted_Full2011.eps", "", 43, xbins, "Dijet Mass [GeV]", "Entries", 890, 6000, 0.01, 3e4);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass_pretag",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                   "nMuons_vs_DijetMass_ratio_PUReweighted_Full2011.eps", "", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//###################################################################################################################
//## PU and k-factor reweighting applied
//###################################################################################################################

//   // Dijet mass plots
//   // CSVL 0Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_CSVL_0Tag_PUkFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUkFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 1Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_CSVL_1Tag_PUkFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUkFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 2Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_CSVL_2Tag_PUkFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUkFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//###################################################################################################################
//## PU and b-tag SF reweighting applied
//###################################################################################################################

//   // Dijet mass plots
//   // CSVL 0Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_CSVL_0Tag_PUSFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 1Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_CSVL_1Tag_PUSFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 2Tag
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_CSVL_2Tag_PUSFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
  
//###################################################################################################################
//## PU, b-tag SF, and k-factor reweighting applied
//###################################################################################################################

  // Dijet mass plots
  // CSVL 0Tag
  overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_CSVL_0Tag_PUSFkFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.03451, 0.984602);
//   // SFb Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_0tag",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

  // CSVL 1Tag
  overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_CSVL_1Tag_PUSFkFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.05427, 0.947269);
//   // SFb Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_1tag",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

  // CSVL 2Tag
  overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                    "DijetMass_CSVL_2Tag_PUSFkFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass",
                  "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
                  "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.13408, 0.875629);
//   // SFb Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__h1_DijetMass_2tag",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

  
//   // CSVL 0Tag
//   // Primary vertex plot
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "Good Vertex Multiplicity", "Events", -0.5, 22.5, 0, 1.1e+05);
//   // MET/SumET
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 0.1, 2.2e6, "Logy", 2);
//   // |DeltaPhi(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "DeltaPhiJ1J2_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "|#Delta#phi|", "Events", 0, 3.15, 0.1, 4.4e5, "Logy");
//   // |DeltaEta(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "DeltaEtaJ1J2_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "|#Delta#eta|", "Events", 0, 1.5, 0.1, 5.8e4);
//   // Pt plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                   "PtJ1_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "p_{T,1} [GeV]", "Events", 0, 2500, 0.1, 1.5e6, "Logy", 40);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                   "PtJ2_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "p_{T,2} [GeV]", "Events", 0, 2500, 0.1, 1.5e6, "Logy", 40);
//   // eta plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ1",
//                   "EtaJ1_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "#eta_{1}", "Events", -3, 3, 0, 4.4e4);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ2",
//                   "EtaJ2_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "#eta_{2}", "Events", -3, 3, 0, 4.4e4);
//   // phi plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ1",
//                   "PhiJ1_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "#phi_{1}", "Events", -3.15, 3.15, 7.3e3, 4.4e4, "", 10);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ2",
//                   "PhiJ2_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "#phi_{2}", "Events", -3.15, 3.15, 7.3e3, 4.4e4, "", 10);
// 
//   // CSVL 1Tag
//   // Primary vertex plot
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "Good Vertex Multiplicity", "Events", -0.5, 22.5, 0, 4.0e+04);
//   // MET/SumET
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 0.1, 8e5, "Logy", 2);
//   // |DeltaPhi(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "DeltaPhiJ1J2_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "|#Delta#phi|", "Events", 0, 3.15, 0.1, 1.6e5, "Logy");
//   // |DeltaEta(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "DeltaEtaJ1J2_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "|#Delta#eta|", "Events", 0, 1.5, 0.1, 2.1e4);
//   // Pt plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                   "PtJ1_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "p_{T,1} [GeV]", "Events", 0, 2500, 0.1, 5.3e5, "Logy", 40);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                   "PtJ2_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "p_{T,2} [GeV]", "Events", 0, 2500, 0.1, 5.3e5, "Logy", 40);
//   // eta plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ1",
//                   "EtaJ1_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "#eta_{1}", "Events", -3, 3, 0, 1.6e4);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ2",
//                   "EtaJ2_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "#eta_{2}", "Events", -3, 3, 0, 1.6e4);
//   // phi plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ1",
//                   "PhiJ1_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "#phi_{1}", "Events", -3.15, 3.15, 2.6e3, 1.6e4, "", 10);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ2",
//                   "PhiJ2_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "#phi_{2}", "Events", -3.15, 3.15, 2.6e3, 1.6e4, "", 10);
// 
//   // CSVL 2Tag
//   // Primary vertex plot
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices",
//                   "nGoodVertices_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "Good Vertex Multiplicity", "Events", -0.5, 22.5, 0, 4.0e+03);
//   // MET/SumET
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 0.1, 8e4, "Logy", 2);
//   // |DeltaPhi(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//                   "DeltaPhiJ1J2_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "|#Delta#phi|", "Events", 0, 3.15, 0.1, 1.6e4, "Logy");
//   // |DeltaEta(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//                   "DeltaEtaJ1J2_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "|#Delta#eta|", "Events", 0, 1.5, 0.1, 2.1e3);
//   // Pt plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                   "PtJ1_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "p_{T,1} [GeV]", "Events", 0, 2500, 0.1, 5.3e4, "Logy", 40);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PtJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                   "PtJ2_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "p_{T,2} [GeV]", "Events", 0, 2500, 0.1, 5.3e4, "Logy", 40);
//   // eta plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ1",
//                   "EtaJ1_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "#eta_{1}", "Events", -3, 3, 0, 1.6e3);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________EtaJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ2",
//                   "EtaJ2_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "#eta_{2}", "Events", -3, 3, 0, 1.6e3);
//   // phi plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ1",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ1",
//                   "PhiJ1_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "#phi_{1}", "Events", -3.15, 3.15, 2.6e2, 1.7e3, "", 10);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________PhiJ2",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ2",
//                   "PhiJ2_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "#phi_{2}", "Events", -3.15, 3.15, 2.6e2, 1.7e3, "", 10);


//   // Muon multiplicity plots
//   // CSVL 0Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nMuons",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                   "nMuons_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", "Muon Multiplicity", "Events", -0.5, 5.5, 0.1, 2e6, "Logy", 1);
// 
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Entries", 890, 6000, 0.01, 3e4);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_Full2011.eps", "0 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // CSVL 1Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nMuons",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                   "nMuons_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", "Muon Multiplicity", "Events", -0.5, 5.5, 0.1, 2e6, "Logy", 1);
// 
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Entries", 890, 6000, 0.01, 1e4);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_Full2011.eps", "1 b-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 1.1e4);
//   // CSVL 2Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nMuons",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                   "nMuons_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", "Muon Multiplicity", "Events", -0.5, 5.5, 0.1, 2e6, "Logy", 1);
// 
//   overlay_DATA_MC_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Entries", 890, 6000, 0.01, 1.1e3);
// 
//   data_MC_ratio_v("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_Full2011.eps", "2 b-tags", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//###################################################################################################################
//##
//## Signal and Background MC
//##
//###################################################################################################################

//###################################################################################################################
//## M=2 TeV, G->bbbar
//###################################################################################################################

//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET_pretag",
//              "METoSumET_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "DeltaPhiJ1J2_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "DeltaEtaJ1J2_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1_pretag",
//              "PtJ1_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2_pretag",
//              "PtJ2_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "EtaJ1_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "EtaJ2_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "PhiJ1_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "PhiJ2_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "", "G->b#bar{b}, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   // CSVL 0Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "0 b-tags", "G->b#bar{b}, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 1Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "1 b-tag", "G->b#bar{b}, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 2Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToBBbar/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToBBbar.eps", "2 b-tags", "G->b#bar{b}, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
// //###################################################################################################################
// //## M=2 TeV, G->qqbar (q=u,d,s)
// //###################################################################################################################
// 
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET_pretag",
//              "METoSumET_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "DeltaPhiJ1J2_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "DeltaEtaJ1J2_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1_pretag",
//              "PtJ1_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2_pretag",
//              "PtJ2_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "EtaJ1_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "EtaJ2_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "PhiJ1_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "PhiJ2_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   // CSVL 0Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "0 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 1Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "1 b-tag", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 2Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToQQbarLight/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToQQbarLight.eps", "2 b-tags", "G->q#bar{q}  (q=u,d,s), M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   
// //###################################################################################################################
// //## M=2 TeV, G->gg
// //###################################################################################################################
// 
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET_pretag",
//              "METoSumET_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//              "DeltaPhiJ1J2_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//              "DeltaEtaJ1J2_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1_pretag",
//              "PtJ1_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2_pretag",
//              "PtJ2_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1_pretag",
//              "EtaJ1_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2_pretag",
//              "EtaJ2_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1_pretag",
//              "PhiJ1_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2_pretag",
//              "PhiJ2_PUSFReweighted_MC_QCD_RSGToGG.eps", "", "G->gg", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   // CSVL 0Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_0Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "0 b-tags", "G->gg, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 1Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_1Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "1 b-tag", "G->gg, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
//   // CSVL 2Tag
//   // MET/SumET
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________METoSumET",
//              "myAnalyzer/cutHisto_allPreviousCuts________METoSumET",
//              "METoSumET_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaPhi(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaPhiJ1J2",
//              "DeltaPhiJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "|#Delta#phi|", "Events", 0, 3.15, 1.e-05, 1, "Logy", 2);
// 
//   // |DeltaEta(J1,J2)|
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "myAnalyzer/cutHisto_allPreviousCuts________DeltaEtaJ1J2",
//              "DeltaEtaJ1J2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "|#Delta#eta|", "Events", 0, 1.5, 0, 0.08);
// 
//   // Pt plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ1",
//              "PtJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "p_{T,1} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PtJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PtJ2",
//              "PtJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "p_{T,2} [GeV]", "Events", 0, 2500, 1.e-05, 1., "Logy", 40);
// 
//   // eta plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ1",
//              "EtaJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "#eta_{1}", "Events", -3, 3, 0, 0.12);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________EtaJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________EtaJ2",
//              "EtaJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "#eta_{2}", "Events", -3, 3, 0, 0.12);
// 
//   // phi plots
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ1",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ1",
//              "PhiJ1_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "#phi_{1}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
// 
//   overlay_MC("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFReweighted_PartonMatching_WideJets_RSGToGG/RSGravitonToJJ_M-2000_TuneZ2_7TeV_pythia6__ferencek-Summer11-PU_S4_START42_V11-v1_EDMTuple_V00-00-04__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________PhiJ2",
//              "myAnalyzer/cutHisto_allPreviousCuts________PhiJ2",
//              "PhiJ2_CSVL_2Tag_PUSFReweighted_MC_QCD_RSGToGG.eps", "2 b-tags", "G->gg, M=2 TeV", "#phi_{2}", "Events", -3.15, 3.15, 0, 0.1, "", 10);
  
}