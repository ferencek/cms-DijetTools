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


double DATA_MC_scale_factor(const string& fInputFile, const string& hDATA, const string& hMC,
                          const Int_t nbins, const Double_t *xbins, const Double_t xMin = 890, const Double_t xMax = 6000)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file->Get(hMC.c_str());

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

  double scaleFactor = h_DijetMass_DATA_rebin->Integral(1, h_DijetMass_DATA_rebin->GetNbinsX()+1) / h_DijetMass_MC_rebin->Integral(1, h_DijetMass_MC_rebin->GetNbinsX()+1);
  cout <<"DATA/MC scale factor = " << scaleFactor << endl;

  return scaleFactor;
}


void overlay_DATA_MC(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax, const Double_t yMin = 0, const Double_t yMax = 0,
                     const string& fLogy = "", const Int_t fRebin = 1)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH1D *h_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file->Get(hMC.c_str());

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
  
  TLegend *legend = new TLegend(.7,.7,.85,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h_DATA_range,"Data","lp");
  legend->AddEntry(h_MC_range,"MC","f");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.25,0.42,"CMS Preliminary");
  l1.DrawLatex(0.25,0.34,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.25,0.29,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.25,0.24,"Wide Jets");
  else l1.DrawLatex(0.25,0.24,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.25,0.19,fLabel.c_str());

  if(fLogy=="Logy") c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void overlay_DATA_MC_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                       const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle,
                       const Double_t xMin = 1050, const Double_t xMax = 3500, const Double_t yMin = 0, const Double_t yMax = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file->Get(hMC.c_str());

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
  
  TLegend *legend = new TLegend(.35,.72,.7,.87);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_DijetMass_DATA_rebin,"Data","lp");
  legend->AddEntry(h_DijetMass_MC_rebin,"QCD PYTHIA6 MC (#times 1.19)","f");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.40, "CMS Preliminary");
  l1.DrawLatex(0.19,0.32, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.27, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.22, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.17, ("Wide Jets" + fLabel).c_str());
  else l1.DrawLatex(0.19,0.17, ("Anti-k_{T} R = 0.7 PFJets" + fLabel).c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_rebin;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file;
}


void overlay_DATA_MC_r_2(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                         const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle,
                         const Double_t xMin = 1050, const Double_t xMax = 3500, const Double_t yMin = 0, const Double_t yMax = 0)
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

  TLegend *legend = new TLegend(.35,.72,.7,.87);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_DijetMass_DATA_rebin,"Data","lp");
  legend->AddEntry(h_DijetMass_MC_rebin,"QCD PYTHIA6 MC (#times 1.19)","f");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.40, "CMS Preliminary");
  l1.DrawLatex(0.19,0.32, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.27, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.22, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.17, ("Wide Jets" + fLabel).c_str());
  else l1.DrawLatex(0.19,0.17, ("Anti-k_{T} R = 0.7 PFJets" + fLabel).c_str());

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


void data_MC_ratio(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                   const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin, const Double_t xMax)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);

  TFile *file = new TFile(fInputFile.c_str());

//   Int_t nKeys = file->GetListOfKeys()->GetEntries();
//   for(Int_t i=0; i<nKeys; ++i) cout << file->GetListOfKeys()->At(i).GetName() << endl;

  TH1D *h_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_MC = (TH1D*)file->Get(hMC.c_str());

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
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.15,0.42,"CMS Preliminary");
  l1.DrawLatex(0.15,0.34,"#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.15,0.29,"#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.15,0.24,"Wide Jets");
  else l1.DrawLatex(0.15,0.24,"Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.15,0.19,fLabel.c_str());

//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DATA_range;
  delete h_MC_range;
  delete c;
  delete file;
}


void data_MC_ratio_r(const string& fInputFile, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                     const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin = 1050, const Double_t xMax = 3500,
                     const Double_t fUp = 0., const Double_t fDown = 0.)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.92);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h_DijetMass_DATA = (TH1D*)file->Get(hDATA.c_str());
  TH1D *h_DijetMass_MC = (TH1D*)file->Get(hMC.c_str());

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
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.64, ("Wide Jets" + fLabel).c_str());
  else l1.DrawLatex(0.19,0.64, ("Anti-k_{T} R = 0.7 PFJets" + fLabel).c_str());

  gPad->RedrawAxis();
  
//   c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete h_DijetMass_DATA_MC;
  delete h_DijetMass_MC_rebin;
  delete h_DijetMass_DATA_range;
  delete h_DijetMass_MC_range;
  delete c;
  delete file;
}


void data_MC_ratio_r_2(const string& fInputFile1, const string& fInputFile2, const string& hDATA, const string& hMC, const string& fOutputFile, const string& fLabel,
                       const Int_t nbins, const Double_t *xbins, const string& fXAxisTitle, const string& fYAxisTitle, const Double_t xMin = 1050, const Double_t xMax = 3500)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1111);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.92);
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
  h_DijetMass_DATA_MC->Draw();

  TF1 *fit = new TF1("fit","pol0",xMin,xMax);
  fit->SetParameter(0, 1.);
  fit->SetLineWidth(2);
  fit->SetLineColor(kRed);

  h_DijetMass_DATA_MC->Fit("fit","R");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.87, "CMS Preliminary");
  l1.DrawLatex(0.19,0.79, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.74, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.69, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  if(fInputFile1.find("WideJets")!=string::npos) l1.DrawLatex(0.19,0.64, ("Wide Jets" + fLabel).c_str());
  else l1.DrawLatex(0.19,0.64, ("Anti-k_{T} R = 0.7 PFJets" + fLabel).c_str());

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

//###################################################################################################################
//##
//## Full 2011 dataset (4.976 /fb)
//##
//###################################################################################################################

  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231,
                      2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,
                      5058, 5253, 5455, 5663, 5877, 6000};

//###################################################################################################################
//##
//## Wide jets
//##
//###################################################################################################################
  
//###################################################################################################################
//## PU reweighting applied
//###################################################################################################################

  // DATA/MC scale factor
//   DATA_MC_scale_factor("CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_EventBins/Final__histograms.root",
//                        "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                        "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                        43, xbins, 890, 6000);

  // Primary vertex plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                   "nGoodVertices_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "Good Vertex Multiplicity", "Events", -0.5, 22.5, 0, 1.2e+05);
// 
//   data_MC_ratio("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                 "DATA__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                 "QCD_Pythia6__cutHisto_allPreviousCuts________nGoodVertices_pretag",
//                 "nGoodVertices_ratio_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "Good Vertex Multiplicity", "DATA/MC", -0.5, 22.5);


//   // Dijet mass plots
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                     "DijetMass_PUReweighted_WideJets_Full2011.png", "", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass_pretag",
//                   "DijetMass_ratio_PUReweighted_WideJets_Full2011.png", "", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//                       
//   // CSVL 0Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_0Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_0Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);  
//   // CSVL 1Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_1Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_1Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // CSVL 2Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_2Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_2Tag_PUReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//   // MET/SumET
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET_pretag",
//                   "METoSumET_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, 0.1, 6e5, "Logy");
//   // TCHEL 0Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_TCHEL_0Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 0Tag", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, "Logy");
//   // TCHEL 1Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_TCHEL_1Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 1Tag", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, "Logy");
//   // TCHEL 2Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________METoSumET",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________METoSumET",
//                   "METoSumET_TCHEL_2Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 2Tag", "#slash{E}_{T}/#SigmaE_{T}", "Events", 0, 1, "Logy");

//   // |DeltaPhi(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaPhiJ1J2_pretag",
//                   "DeltaPhiJ1J2_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "|#Delta#phi(j_{1},j_{2})|", "Events", 0, 3.15, 0.1, 6e5, "Logy");
// 
//   // |DeltaEta(J1,J2)|
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DeltaEtaJ1J2_pretag",
//                   "DeltaEtaJ1J2_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "|#Delta#eta(j_{1},j_{2})|", "Events", 0, 1.5, 0.1, 8e4);
// 
// 
//   // Pt plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ1_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1_pretag",
//                    "PtJ1_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "p_{T,1}", "Events", 0, 2200, 0.1, 6e5, "Logy", 20);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ2_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2_pretag",
//                    "PtJ2_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "p_{T,2}", "Events", 0, 2200, 0.1, 6e5, "Logy", 20);
//   // TCHEL 0Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ1",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                    "PtJ1_TCHEL_0Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 0Tag", "p_{T,1}", "Events", 0, 2200, "Logy", 20);
//
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ2",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                    "PtJ2_TCHEL_0Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 0Tag", "p_{T,2}", "Events", 0, 2200, "Logy", 20);
//   // TCHEL 1Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ1",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                    "PtJ1_TCHEL_1Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 1Tag", "p_{T,1}", "Events", 0, 2200, "Logy", 20);
//
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ2",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                    "PtJ2_TCHEL_1Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 1Tag", "p_{T,2}", "Events", 0, 2200, "Logy", 20);
//   // TCHEL 2Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ1",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ1",
//                    "PtJ1_TCHEL_2Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 2Tag", "p_{T,1}", "Events", 0, 2200, "Logy", 20);
//
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PtJ2",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PtJ2",
//                    "PtJ2_TCHEL_2Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 2Tag", "p_{T,2}", "Events", 0, 2200, "Logy", 20);

//   // eta plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________EtaJ1_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ1_pretag",
//                    "EtaJ1_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "#eta_{1}", "Events", -3, 3, 0, 5e4);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________EtaJ2_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________EtaJ2_pretag",
//                    "EtaJ2_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "#eta_{2}", "Events", -3, 3, 0, 5e4);
// 
//   // phi plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PhiJ1_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ1_pretag",
//                    "PhiJ1_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "#phi_{1}", "Events", -3.15, 3.15, 1e4, 6e4, "", 10);
// 
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUReweighted_PartonMatching_WideJets/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________PhiJ2_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________PhiJ2_pretag",
//                    "PhiJ2_PUReweighted_Full2011.png", "M_{jj}>890 GeV", "#phi_{2}", "Events", -3.15, 3.15, 1e4, 6e4, "", 10);

  // Muon multiplicity plots
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________nMuons_pretag",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons_pretag",
//                    "nMuons_PUReweighted_Full2011.png", "M_{jj}>944 GeV", "Muon Multiplicity", "Events", -0.5, 5.5, "Logy", 1, 0.1, 1e6);
//
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass_pretag",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                     "nMuons_vs_DijetMass_PUReweighted_Full2011.png", "", 42, xbins, "Dijet Mass [GeV]", "Entries", 944, 6000, 0.01, 20000);
//
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass_pretag",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass_pretag",
//                   "nMuons_vs_DijetMass_ratio_PUReweighted_Full2011.png", "", 42, xbins, "Dijet Mass [GeV]", "Data/MC", 944, 6000);
//   // TCHEL 0Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________nMuons",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                    "nMuons_TCHEL_0Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 0Tag", "Muon Multiplicity", "Events", -0.5, 5.5, "Logy", 1, 0.1, 1e6);
//
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_TCHEL_0Tag_PUReweighted_Full2011.png", "TCHEL 0Tag", 42, xbins, "Dijet Mass [GeV]", "Entries", 944, 6000, 0.01, 2000);
//
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_TCHEL_0Tag_PUReweighted_Full2011.png", "TCHEL 0Tag", 42, xbins, "Dijet Mass [GeV]", "Data/MC", 944, 6000);
//   // TCHEL 1Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________nMuons",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                    "nMuons_TCHEL_1Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 1Tag", "Muon Multiplicity", "Events", -0.5, 5.5, "Logy", 1, 0.1, 1e6);
//
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_TCHEL_1Tag_PUReweighted_Full2011.png", "TCHEL 1Tag", 42, xbins, "Dijet Mass [GeV]", "Entries", 944, 6000, 0.01, 1e4);
//
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_TCHEL_1Tag_PUReweighted_Full2011.png", "TCHEL 1Tag", 42, xbins, "Dijet Mass [GeV]", "Data/MC", 944, 6000);
//   // TCHEL 2Tag
//   overlay_DATA_MC("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                    "DATA__cutHisto_allPreviousCuts________nMuons",
//                    "QCD_Pythia6__cutHisto_allPreviousCuts________nMuons",
//                    "nMuons_TCHEL_2Tag_PUReweighted_Full2011.png", "M_{jj}>944 GeV, TCHEL 2Tag", "Muon Multiplicity", "Events", -0.5, 5.5, "Logy", 1, 0.1, 1e6);
//
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                     "DATA__h1_nMuons_vs_DijetMass",
//                     "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                     "nMuons_vs_DijetMass_TCHEL_2Tag_PUReweighted_Full2011.png", "TCHEL 2Tag", 42, xbins, "Dijet Mass [GeV]", "Entries", 944, 6000, 0.01, 5000);
//
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
//                   "DATA__h1_nMuons_vs_DijetMass",
//                   "QCD_Pythia6__h1_nMuons_vs_DijetMass",
//                   "nMuons_vs_DijetMass_ratio_TCHEL_2Tag_PUReweighted_Full2011.png", "TCHEL 2Tag", 42, xbins, "Dijet Mass [GeV]", "Data/MC", 944, 6000);

//###################################################################################################################
//## PU and k-factor reweighting applied
//###################################################################################################################

//   // CSVL 0Tag
//   overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                       "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                       "DATA__cutHisto_allPreviousCuts________DijetMass",
//                       "QCD_Pythia6__h1_DijetMass_0tag",
//                       "DijetMass_CSVL_0Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_ratio_CSVL_0Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 1Tag
//   overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                       "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                       "DATA__cutHisto_allPreviousCuts________DijetMass",
//                       "QCD_Pythia6__h1_DijetMass_1tag",
//                       "DijetMass_CSVL_1Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_ratio_CSVL_1Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
// 
//   // CSVL 2Tag
//   overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                       "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                       "DATA__cutHisto_allPreviousCuts________DijetMass",
//                       "QCD_Pythia6__h1_DijetMass_2tag",
//                       "DijetMass_CSVL_2Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUkFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_ratio_CSVL_2Tag_PUkFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//###################################################################################################################
//## PU and b-tag SF reweighting applied
//###################################################################################################################

  // CSVL 0Tag
  overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                      "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                      "DATA__cutHisto_allPreviousCuts________DijetMass",
                      "QCD_Pythia6__h1_DijetMass_0tag",
                      "DijetMass_CSVL_0Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__h1_DijetMass_0tag",
                    "DijetMass_ratio_CSVL_0Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

  // CSVL 1Tag
  overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                      "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                      "DATA__cutHisto_allPreviousCuts________DijetMass",
                      "QCD_Pythia6__h1_DijetMass_1tag",
                      "DijetMass_CSVL_1Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__h1_DijetMass_1tag",
                    "DijetMass_ratio_CSVL_1Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

  // CSVL 2Tag
  overlay_DATA_MC_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                      "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                      "DATA__cutHisto_allPreviousCuts________DijetMass",
                      "QCD_Pythia6__h1_DijetMass_2tag",
                      "DijetMass_CSVL_2Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);

  data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                    "CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_QCD/Final__histograms.root",
                    "DATA__cutHisto_allPreviousCuts________DijetMass",
                    "QCD_Pythia6__h1_DijetMass_2tag",
                    "DijetMass_ratio_CSVL_2Tag_PUSFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
  
//###################################################################################################################
//## PU, b-tag SF, and k-factor reweighting applied
//###################################################################################################################

//   // CSVL 0Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_0Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.02623, 0.977853);
//   // SFb Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_0tag",
//                     "DijetMass_ratio_CSVL_0Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//   // CSVL 1Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_1Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 1-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.04482, 0.939841);
//   // SFb Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_1tag",
//                     "DijetMass_ratio_CSVL_1Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);

//   // CSVL 2Tag
//   overlay_DATA_MC_r("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                     "DijetMass_CSVL_2Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Events", 890, 6000, 0.01, 5e5);
// 
//   data_MC_ratio_r("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass",
//                   "QCD_Pythia6__cutHisto_allPreviousCuts________DijetMass",
//                   "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_WideJets_Full2011.png", ", CSVL 2-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000, 1.12124, 0.86858);
//   // SFb Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFbUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFb Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFbDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFbDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Up
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlUp/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFlUp_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
//   // SFl Down
//   data_MC_ratio_r_2("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                     "CRAB_Jobs_MainAnalysis_CSVL_PUSFkFReweighted_PartonMatching_WideJets_QCD_SFlDown/Final__histograms.root",
//                     "DATA__cutHisto_allPreviousCuts________DijetMass",
//                     "QCD_Pythia6__h1_DijetMass_2tag",
//                     "DijetMass_ratio_CSVL_2Tag_PUSFkFReweighted_SFlDown_WideJets_Full2011.png", ", CSVL 0-tag", 43, xbins, "Dijet Mass [GeV]", "Data/MC", 890, 6000);
}