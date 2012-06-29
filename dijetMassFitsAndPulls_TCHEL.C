#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLine.h"
#include "tdrstyle.C"
#include "ResonancesCrossSection.h"


using namespace std;

// TCHEL 0- and 2-tag efficienies and efficieny erros for heavy and light flavor final states
double masses_eff[5] = {500.0, 700.0, 1200.0, 2000.0, 3500.0};
double eff0_h[5] = {0.04695263596288976, 0.048874698629358393, 0.072071168683632919, 0.11364769382708144, 0.12268614779446146};
double eff2_h[5] = {0.62224739946394214, 0.59487876997338773, 0.52113194683979203, 0.44451634223700387, 0.4203374639490588};
double eff0_l[5] = {0.52971242981698297, 0.438814755302947, 0.32625361910162476, 0.25817068642932711, 0.19069406258687238};
double eff2_l[5] = {0.077383024858261651, 0.11406735842487882, 0.18287425028050511, 0.2428131772568019, 0.3185558099993932};

TGraph *g_eff0_h = new TGraph(5, masses_eff, eff0_h);
TGraph *g_eff2_h = new TGraph(5, masses_eff, eff2_h);
TGraph *g_eff0_l = new TGraph(5, masses_eff, eff0_l);
TGraph *g_eff2_l = new TGraph(5, masses_eff, eff2_l);


// Default fit function
Double_t fitQCD( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

// Alternate fit function A (4 parameter fit function)
Double_t fitQCD_A( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(x,p[2]);
}

// Alternate fit function B (3 parameter fit function)
Double_t fitQCD_B( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x,p[1])/pow(x,p[2]);
}


void performFit(const string& fInputFile, const string& fPlot,
                const string& fbbFile, const string& fqqFile, const string& fggFile,
                const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile, const Double_t fP0=1e-04, const Double_t fP1=1e+01, const Double_t fP2=4e+00, const Double_t fP3=-0.1e-01, const string& fixP3 = "")
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.97);
//   gStyle->SetStatH(0);  // uncomment for PAS
//   gStyle->SetStatW(0);  // uncomment for PAS
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02);
  gROOT->ForceStyle();

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(10000);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot = (TH1D*)file->Get(fPlot.c_str());
  TH1D *h1_plot_r = (TH1D*)h1_plot->Rebin(fNbins,"h1_plot_r",fBins);

  TH1D *h1_pulls = (TH1D*)h1_plot_r->Clone("h1_pulls");
  h1_pulls->Reset();

  const double alpha = 1 - 0.6827;
  double vx[1000],vy[1000],vexl[1000],vexh[1000],veyl[1000],veyh[1000];
  
//   int bins_to_skip = 0;
//   
//   for(int i=h1_plot_r->GetNbinsX(); i>0; --i)
//   {
//       if( h1_plot_r->GetBinContent(i)!=0. ) break;
//       ++bins_to_skip;
//   }
  
  int bins_to_skip = 10; // to have the same ndof in all three b-tag multiplicity bins
  int bins_to_process = (h1_plot_r->GetNbinsX()-bins_to_skip);
  
  for(int i=0; i<bins_to_process; ++i)
  {
      double n    = h1_plot_r->GetBinContent(i+1);
      double dm   = h1_plot_r->GetBinWidth(i+1);
//       double mass = h1_plot_r->GetBinCenter(i+1);
      double xl   = h1_plot_r->GetBinLowEdge(i+1);
      double xh   = xl+dm;
      vx[i]   = (xl+xh)/2.;
      vexl[i] = dm/2.;
      vexh[i] = dm/2.;
      vy[i]   = n / (dm*fLumi);

      double l = 0.5*TMath::ChisquareQuantile(alpha/2,2*n);
      double h = (n==0) ? ( 0.5*TMath::ChisquareQuantile(1-alpha,2*(n+1)) ) : ( 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(n+1)) );

      veyl[i] = (n-l)/(fLumi*dm);
      veyh[i] = (h-n)/(fLumi*dm);
  }

  // data in the graph format
  TGraphAsymmErrors *g = new TGraphAsymmErrors(bins_to_process,vx,vy,vexl,vexh,veyl,veyh);

  // Fit to data
  TF1 *fit = new TF1("fit",fitQCD,fFitXmin,fFitXmax,4); // 4 Par. Fit
  gStyle->SetOptFit(1111);
  fit->SetParameter(0,fP0);
  fit->SetParameter(1,fP1);
  fit->SetParameter(2,fP2);
  fit->SetParameter(3,fP3);
//   fit->SetParLimits(0, 0., 1.);
//   fit->SetParLimits(1, 0., 100.);
//   fit->SetParLimits(2, 0., 100.);
//   fit->SetParLimits(3,-5., 5.);
  if(fixP3=="FixP3") fit->FixParameter(3,fP3);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);

  cout << "*********************************************************"<<endl;
  TFitResultPtr s = g->Fit("fit","SR");
  TString status_default = gMinuit->fCstatu.Data();
  // Results of the fit
  cout << "*********************************************************"<<endl;
  Double_t chi_fit = fit->GetChisquare();
  Double_t ndf_fit = fit->GetNDF();
  cout << "Chi2/ndf: " << chi_fit << "/" << ndf_fit << " = " << chi_fit/ndf_fit << endl;
  cout << "Status: "<<status_default<<endl;
  cout << "*********************************************************"<<endl;

  // Print fit results
  s->Print("V");
  
  for(int i=0; i<bins_to_process; ++i)
    {
      double data = vy[i];
      double error = veyh[i];
      double m = vx[i];
      double fit_default = fit->Eval(m);
      if(data>=fit_default) error = veyl[i];
      if(data<fit_default)  error = veyh[i];

      if(error!=0.)
      {
        //cout << "m = " << m << std::endl;
        h1_pulls->SetBinContent(i+1,(data-fit_default)/error);
        h1_pulls->SetBinError(i+1,1.);
      }
    }
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->Divide(1,2,0,0,0);

  // Begin top part
  c->cd(1);

  TPad *p_1 = (TPad*)c->GetPad(1);
  p_1->SetPad(0.,0.25,1.,1.);
  p_1->SetLogy();
  p_1->SetLeftMargin(0.13);
  p_1->SetRightMargin(0.03);
  p_1->SetTopMargin(0.03);

  TH1F *vFrame = p_1->DrawFrame(700.0,4.0e-08,4337.0,10.0);
  vFrame->SetTitle("");
  vFrame->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  //vFrame->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  vFrame->GetYaxis()->SetTitle("(1/#scale[0.5]{#int}Ldt)(dN/dm) [pb/GeV]");
  vFrame->GetYaxis()->SetTitleOffset(1.05);

  g->Draw("P");
 
  TLegend *legend = new TLegend(.4,.7,.6,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->AddEntry(g, "Data","lp");
  legend->AddEntry(fit, "Fit","l");
  legend->Draw();


  // resonance masses
  double rsg1 = 1400;
  double rsg2 = 2000;
  double zprime1 = 1700;
  double zprime2 = 2300;

  unsigned rsgBin1 = (rsg1 / 100) - 5;
  unsigned rsgBin2 = (rsg2 / 100) - 5;
  unsigned zprimeBin1 = (zprime1 / 100) - 5;
  unsigned zprimeBin2 = (zprime2 / 100) - 5;
  
  // fbb
  double fbb_rsg = 0.1;
  double fbb_zprime = 0.216;

  // get resonance shapes
  TFile *file_bb = new TFile(fbbFile.c_str());
  TFile *file_qq = new TFile(fqqFile.c_str());
  TFile *file_gg = new TFile(fggFile.c_str());

  TH1D *h1_bb_rsg1 = (TH1D*)file_bb->Get(Form("h_bb_%.0f", rsg1));
  TH1D *h1_gg_rsg1 = (TH1D*)file_gg->Get(Form("h_gg_%.0f", rsg1));
  TH1D *h1_bb_rsg2 = (TH1D*)file_bb->Get(Form("h_bb_%.0f", rsg2));
  TH1D *h1_gg_rsg2 = (TH1D*)file_gg->Get(Form("h_gg_%.0f", rsg2));
  TH1D *h1_bb_zprime1 = (TH1D*)file_bb->Get(Form("h_bb_%.0f", zprime1));
  TH1D *h1_qq_zprime1 = (TH1D*)file_qq->Get(Form("h_qq_%.0f", zprime1));
  TH1D *h1_bb_zprime2 = (TH1D*)file_bb->Get(Form("h_bb_%.0f", zprime2));
  TH1D *h1_qq_zprime2 = (TH1D*)file_qq->Get(Form("h_qq_%.0f", zprime2));

  // calculate b-tag efficiencies
  double eff_h_rsg1, eff_l_rsg1, eff_h_rsg2, eff_l_rsg2, eff_h_zprime1, eff_l_zprime1, eff_h_zprime2, eff_l_zprime2 = 0.;
  
  if( fLabel.find("0 b-tags")!=string::npos )
  {
    eff_h_rsg1 = g_eff0_h->Eval(rsg1);
    eff_l_rsg1 = g_eff0_l->Eval(rsg1);
    eff_h_rsg2 = g_eff0_h->Eval(rsg2);
    eff_l_rsg2 = g_eff0_l->Eval(rsg2);
    eff_h_zprime1 = g_eff0_h->Eval(zprime1);
    eff_l_zprime1 = g_eff0_l->Eval(zprime1);
    eff_h_zprime2 = g_eff0_h->Eval(zprime2);
    eff_l_zprime2 = g_eff0_l->Eval(zprime2);
  }
  else if( fLabel.find("1 b-tag")!=string::npos )
  {
    eff_h_rsg1 = 1 - g_eff0_h->Eval(rsg1) - g_eff2_h->Eval(rsg1);
    eff_l_rsg1 = 1 - g_eff0_l->Eval(rsg1) - g_eff2_l->Eval(rsg1);
    eff_h_rsg2 = 1 - g_eff0_h->Eval(rsg2) - g_eff2_h->Eval(rsg2);
    eff_l_rsg2 = 1 - g_eff0_l->Eval(rsg2) - g_eff2_l->Eval(rsg2);
    eff_h_zprime1 = 1 - g_eff0_h->Eval(zprime1) - g_eff2_h->Eval(zprime1);
    eff_l_zprime1 = 1 - g_eff0_l->Eval(zprime1) - g_eff2_l->Eval(zprime1);
    eff_h_zprime2 = 1 - g_eff0_h->Eval(zprime2) - g_eff2_h->Eval(zprime2);
    eff_l_zprime2 = 1 - g_eff0_l->Eval(zprime2) - g_eff2_l->Eval(zprime2);
  }
  else
  {
    eff_h_rsg1 = g_eff2_h->Eval(rsg1);
    eff_l_rsg1 = g_eff2_l->Eval(rsg1);
    eff_h_rsg2 = g_eff2_h->Eval(rsg2);
    eff_l_rsg2 = g_eff2_l->Eval(rsg2);
    eff_h_zprime1 = g_eff2_h->Eval(zprime1);
    eff_l_zprime1 = g_eff2_l->Eval(zprime1);
    eff_h_zprime2 = g_eff2_h->Eval(zprime2);
    eff_l_zprime2 = g_eff2_l->Eval(zprime2);
  }

  std::vector<double> v_rsg1, v_rsg1_m, v_rsg2, v_rsg2_m, v_zprime1, v_zprime1_m, v_zprime2, v_zprime2_m;

  for(int i=0;i<h1_plot_r->GetNbinsX();i++)
  {
    double dm = h1_plot_r->GetBinWidth(i+1);
    double mass = h1_plot_r->GetBinCenter(i+1);
    
    if(mass>0.4*rsg1 && mass<1.3*rsg1)
    {
      double prob = eff_h_rsg1*fbb_rsg*h1_bb_rsg1->GetBinContent(h1_bb_rsg1->GetXaxis()->FindBin(mass)) + eff_l_rsg1*(1-fbb_rsg)*h1_gg_rsg1->GetBinContent(h1_gg_rsg1->GetXaxis()->FindBin(mass));
//       cout << "mass = " << mass << " prob = " << prob << endl;
      if( prob > 0 )
      {
        v_rsg1.push_back(prob * rsgraviton[rsgBin1] / dm);
        v_rsg1_m.push_back(mass);
      }
    }

    if(mass>0.4*rsg2 && mass<1.3*rsg2)
    {
      double prob = eff_h_rsg2*fbb_rsg*h1_bb_rsg2->GetBinContent(h1_bb_rsg2->GetXaxis()->FindBin(mass)) + eff_l_rsg2*(1-fbb_rsg)*h1_gg_rsg2->GetBinContent(h1_gg_rsg2->GetXaxis()->FindBin(mass));
//       cout << "mass = " << mass << " prob = " << prob << endl;
      if( prob > 0 )
      {
        v_rsg2.push_back(prob * rsgraviton[rsgBin2] / dm);
        v_rsg2_m.push_back(mass);
      }
    }

    if(mass>0.4*zprime1 && mass<1.3*zprime1)
    {
      double prob = eff_h_zprime1*fbb_zprime*h1_bb_zprime1->GetBinContent(h1_bb_zprime1->GetXaxis()->FindBin(mass)) + eff_l_zprime1*(1-fbb_zprime)*h1_qq_zprime1->GetBinContent(h1_qq_zprime1->GetXaxis()->FindBin(mass));
//       cout << "mass = " << mass << " prob = " << prob << endl;
      if( prob > 0 )
      {
        v_zprime1.push_back(prob * zprime[zprimeBin1] / dm);
        v_zprime1_m.push_back(mass);
      }
    }

    if(mass>0.4*zprime2 && mass<1.3*zprime2)
    {
      double prob = eff_h_zprime2*fbb_zprime*h1_bb_zprime2->GetBinContent(h1_bb_zprime2->GetXaxis()->FindBin(mass)) + eff_l_zprime2*(1-fbb_zprime)*h1_qq_zprime2->GetBinContent(h1_qq_zprime2->GetXaxis()->FindBin(mass));
//       cout << "mass = " << mass << " prob = " << prob << endl;
      if( prob > 0 )
      {
        v_zprime2.push_back(prob * zprime[zprimeBin2] / dm);
        v_zprime2_m.push_back(mass);
      }
    }
  }

  const unsigned size1 = v_rsg1.size();
  double xs1[size1], mass1[size1];
  for(unsigned i=0; i<size1; i++)
  {
    mass1[i] = v_rsg1_m[i];
    xs1[i] = v_rsg1[i];
//     cout << "mass1[" << i << "] = " << mass1[i] << std::endl;
//     cout << "xs1[" << i << "] = " << xs1[i] << std::endl;
  }

  const unsigned size2 = v_rsg2.size();
  double xs2[size2], mass2[size2];
  for(unsigned i=0; i<size2; i++)
  {
    mass2[i] = v_rsg2_m[i];
    xs2[i] = v_rsg2[i];
//     cout << "mass2[" << i << "] = " << mass2[i] << std::endl;
//     cout << "xs2[" << i << "] = " << xs2[i] << std::endl;
  }

  const unsigned size3 = v_zprime1.size();
  double xs3[size3], mass3[size3];
  for(unsigned i=0; i<size3; i++)
  {
    mass3[i] = v_zprime1_m[i];
    xs3[i] = v_zprime1[i];
//     cout << "mass3[" << i << "] = " << mass3[i] << std::endl;
//     cout << "xs3[" << i << "] = " << xs3[i] << std::endl;
  }

  const unsigned size4 = v_zprime2.size();
  double xs4[size4], mass4[size4];
  for(unsigned i=0; i<size4; i++)
  {
    mass4[i] = v_zprime2_m[i];
    xs4[i] = v_zprime2[i];
//     cout << "mass4[" << i << "] = " << mass4[i] << std::endl;
//     cout << "xs4[" << i << "] = " << xs4[i] << std::endl;
  }
  
  TGraph* gr_rsg1 = new TGraph(size1, mass1, xs1);
  TGraph* gr_rsg2 = new TGraph(size2, mass2, xs2);
  TGraph* gr_zprime1 = new TGraph(size3, mass3, xs3);
  TGraph* gr_zprime2 = new TGraph(size4, mass4, xs4);

  gr_rsg1->SetLineColor(kGreen+2);
  gr_rsg1->SetLineStyle(7);
  gr_rsg1->SetLineWidth(2);
  gr_rsg2->SetLineColor(kGreen+2);
  gr_rsg2->SetLineStyle(7);
  gr_rsg2->SetLineWidth(2);
  gr_zprime1->SetLineColor(kRed);
  gr_zprime1->SetLineStyle(5);
  gr_zprime1->SetLineWidth(2);
  gr_zprime2->SetLineColor(kRed);
  gr_zprime2->SetLineStyle(5);
  gr_zprime2->SetLineWidth(2);

  gr_rsg1->Draw("sameC");
  gr_rsg2->Draw("sameC");
  gr_zprime1->Draw("sameC");
  gr_zprime2->Draw("sameC");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.70,0.57, "CMS Preliminary");
  l1.DrawLatex(0.70,0.49, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.71,0.44, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.70,0.39, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.70,0.34, "Wide Jets");
  l1.SetTextColor(kGreen+2);
  if( fLabel.find("0 b-tags")!=string::npos )
  {
    l1.DrawLatex(0.20,0.45, "G (1.4 TeV)");
    l1.DrawLatex(0.35,0.30, "G (2 TeV)");
  }
  else if( fLabel.find("1 b-tag")!=string::npos )
  {
    l1.DrawLatex(0.20,0.40, "G (1.4 TeV)");
    l1.DrawLatex(0.35,0.26, "G (2 TeV)");
  }
  else
  {
    l1.DrawLatex(0.20,0.32, "G (1.4 TeV)");
    l1.DrawLatex(0.35,0.16, "G (2 TeV)");
  }
  l1.SetTextColor(kRed);
  if( fLabel.find("0 b-tags")!=string::npos )
  {
    l1.DrawLatex(0.35,0.42, "Z' (1.7 TeV)");
    l1.DrawLatex(0.50,0.29, "Z' (2.3 TeV)");
  }
  else if( fLabel.find("1 b-tag")!=string::npos )
  {
    l1.DrawLatex(0.35,0.38, "Z' (1.7 TeV)");
    l1.DrawLatex(0.50,0.25, "Z' (2.3 TeV)");
  }
  else
  {
    l1.DrawLatex(0.35,0.30, "Z' (1.7 TeV)");
    l1.DrawLatex(0.50,0.16, "Z' (2.3 TeV)");
  }
  l1.SetTextColor(kBlack);
  l1.SetTextSize(0.06);
  l1.DrawLatex(0.17,0.89, fLabel.c_str());
  
  // End top part
  
  // Begin bottom part
  c->cd(2);
  TPad *p_2 = (TPad*)c->GetPad(2);
  p_2->SetPad(0.,0.,1.,0.2535);
  p_2->SetBottomMargin(0.35);
  p_2->SetLeftMargin(0.13);
  p_2->SetRightMargin(0.03);
  p_2->SetGridx();

  TH1F *vFrame2 = p_2->DrawFrame(700.0, -3.3, 4337.0, 3.3);
  vFrame2->SetTitle("");
  vFrame2->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  vFrame2->GetYaxis()->SetTitle("(Data-Fit)/Error");
  vFrame2->GetYaxis()->SetTitleSize(0.12);
  vFrame2->GetYaxis()->SetLabelSize(0.10);
  vFrame2->GetYaxis()->SetTitleOffset(0.5);
  vFrame2->GetXaxis()->SetTitleOffset(0.90);
  vFrame2->GetXaxis()->SetTitleSize(0.175);
  vFrame2->GetXaxis()->SetLabelSize(0.155);

  h1_pulls->SetLineWidth(0);
  h1_pulls->SetLineColor(kRed);
  h1_pulls->SetFillColor(kRed);

  h1_pulls->Draw("samehist");

  double line_y1 = h1_pulls->GetBinContent(1);
  if( line_y1>0. ) line_y1 = 0.;
  TLine *line2 = new TLine(890.,line_y1,890.,-3.3);
  line2->SetLineColor(0);
  line2->SetLineWidth(2);
  line2->Draw("");

  TLine *line = new TLine(700.,0,4337,0);
  line->Draw("");

  gPad->RedrawAxis();
  
  c->SaveAs(fOutputFile.c_str());

  // End bottom part

//   delete c;
//   delete file;
}


void makePlots()
{
  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856,
                      1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558,
                      3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};
                      
  performFit("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUSFkFReweighted_PartonMatching_WideJets_DATA/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
             43, xbins, 4976, 890, 4200, "0 b-tags",
             "DijetMass_fit_pulls_TCHEL_0Tag_WideJets.eps", 9.16513e-07, 5.42799e+00, 7.12254e+00, 2.55957e-01);

  performFit("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUSFkFReweighted_PartonMatching_WideJets_DATA/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
             43, xbins, 4976, 890, 4200, "1 b-tag",
             "DijetMass_fit_pulls_TCHEL_1Tag_WideJets.eps", 1.71330e-05, 7.86344e+00, 5.89232e+00, 1.58363e-01);

  performFit("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUSFkFReweighted_PartonMatching_WideJets_DATA/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_bb.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_qq.root",
             "LimitCode/Data_and_ResonanceShapes/Resonance_Shapes_WideJets_gg.root",
             43, xbins, 4976, 890, 4200, "2 b-tags",
             "DijetMass_fit_pulls_TCHEL_2Tag_WideJets.eps", 3.11395e-05, 9.02689e+00, 5.34280e+00, 1.98697e-01);
}
