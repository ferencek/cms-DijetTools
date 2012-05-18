#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
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
#include "TPolyLine.h"
#include "TColor.h"
#include "TMatrixD.h"
#include "tdrstyle.C"

using namespace std;


Int_t nbtags = 0;

// Default fit function
Double_t fitQCD( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

// Diagonalized default fit function
Double_t fitQCD_diag(Double_t *m, Double_t *p)
{
    double data_0[] = {
       1.328e-05,  -3.847e-05,  -0.0001867,          -1,
          0.9989,    -0.03026,     0.03594,    7.72e-06,
          0.0191,      0.9604,      0.2779,  -8.857e-05,
         0.04293,      0.2769,     -0.9599,   0.0001691
    };
    double data_1[] = {
       4.444e-06,  -1.266e-05,  -6.138e-05,          -1,
          0.9989,     -0.0301,     0.03689,   2.555e-06,
         0.01866,      0.9603,      0.2783,  -2.916e-05,
         0.04381,      0.2773,     -0.9598,    5.56e-05
    };
    double data_2[] = {
       8.507e-08,  -2.392e-07,  -1.157e-06,          -1,
          0.9989,    -0.02876,     0.03758,   4.836e-08,
         0.01719,      0.9604,      0.2781,  -5.502e-07,
         0.04409,      0.2771,     -0.9598,   1.048e-06
    };

    TMatrixD d = TMatrixD(4,4);
    if(nbtags==0)      d.SetMatrixArray(data_0);
    else if(nbtags==1) d.SetMatrixArray(data_1);
    else               d.SetMatrixArray(data_2);

    double p0 = d(0,0)*p[0] + d(0,1)*p[1] + d(0,2)*p[2] + d(0,3)*p[3];
    double p1 = d(1,0)*p[0] + d(1,1)*p[1] + d(1,2)*p[2] + d(1,3)*p[3];
    double p2 = d(2,0)*p[0] + d(2,1)*p[1] + d(2,2)*p[2] + d(2,3)*p[3];
    double p3 = d(3,0)*p[0] + d(3,1)*p[1] + d(3,2)*p[2] + d(3,3)*p[3];

    double x=m[0]/7000.;
    return p0*pow(1.-x,p1)/pow(x,p2+p3*log(x));
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
                const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile,
                const Double_t fPD0, const Double_t fPD1, const Double_t fPD2, const Double_t fPD3, const Double_t fPD0e,  const Double_t fPD1e, const Double_t fPD2e, const Double_t fPD3e,
                const Double_t fPA0, const Double_t fPA1, const Double_t fPA2, const Double_t fPA3,
                const Double_t fPB0, const Double_t fPB1, const Double_t fPB2,
                const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3, const Int_t fRestrict = 0)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetStatH(0.2);
  //gStyle->SetStatW(0.2);
  //gStyle->SetStatX(0.97);
  //gStyle->SetStatY(0.97);
  gStyle->SetStatH(0);
  gStyle->SetStatW(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02);
  gROOT->ForceStyle();

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(30000);

  if(fLabel.find("0")!=string::npos)      nbtags=0;
  else if(fLabel.find("1")!=string::npos) nbtags=1;
  else                                    nbtags=2;

  Int_t MyPalette[31];
  Double_t red[]   = {0., 0., 0., 1., 1.};
  Double_t green[] = {0., 1., 1., 1., 0.};
  Double_t blue[]  = {1., 1., 0., 0., 0.};
  Double_t stop[]  = {0., 0.25, .50, 0.75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(5, stop, red, green, blue, 31);
  for (int i=0; i<31; ++i) MyPalette[i] = FI+i;
  
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
  g->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  g->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  g->GetXaxis()->SetNdivisions(1005);
  g->GetYaxis()->SetRangeUser(4.0e-08,10.0);
  g->GetYaxis()->SetTitleOffset(1.05);
    
  // Diagonalized default fit
  TF1 *fit_default_diag = new TF1("fit_default_diag",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag->FixParameter(0,fPD0);
  fit_default_diag->FixParameter(1,fPD1);
  fit_default_diag->FixParameter(2,fPD2);
  fit_default_diag->FixParameter(3,fPD3);
  fit_default_diag->SetLineWidth(2);
  fit_default_diag->SetLineColor(kBlack);

  // Diagonalized default fit + P0
  TF1 *fit_default_diag_0p = new TF1("fit_default_diag_0p",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_0p->FixParameter(0,fPD0+fPD0e);
  fit_default_diag_0p->FixParameter(1,fPD1);
  fit_default_diag_0p->FixParameter(2,fPD2);
  fit_default_diag_0p->FixParameter(3,fPD3);
  fit_default_diag_0p->SetLineWidth(1);
  fit_default_diag_0p->SetLineStyle(2);
  fit_default_diag_0p->SetLineColor(kViolet+7);

    // Diagonalized default fit - P0
  TF1 *fit_default_diag_0m = new TF1("fit_default_diag_0m",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_0m->FixParameter(0,fPD0-fPD0e);
  fit_default_diag_0m->FixParameter(1,fPD1);
  fit_default_diag_0m->FixParameter(2,fPD2);
  fit_default_diag_0m->FixParameter(3,fPD3);
  fit_default_diag_0m->SetLineWidth(1);
  fit_default_diag_0m->SetLineStyle(2);
  fit_default_diag_0m->SetLineColor(kRed);

  // Diagonalized default fit + P1
  TF1 *fit_default_diag_1p = new TF1("fit_default_diag_1p",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_1p->FixParameter(0,fPD0);
  fit_default_diag_1p->FixParameter(1,fPD1+fPD1e);
  fit_default_diag_1p->FixParameter(2,fPD2);
  fit_default_diag_1p->FixParameter(3,fPD3);

  // Diagonalized default fit - P1
  TF1 *fit_default_diag_1m = new TF1("fit_default_diag_1m",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_1m->FixParameter(0,fPD0);
  fit_default_diag_1m->FixParameter(1,fPD1-fPD1e);
  fit_default_diag_1m->FixParameter(2,fPD2);
  fit_default_diag_1m->FixParameter(3,fPD3);

  // Diagonalized default fit + P2
  TF1 *fit_default_diag_2p = new TF1("fit_default_diag_2p",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_2p->FixParameter(0,fPD0);
  fit_default_diag_2p->FixParameter(1,fPD1);
  fit_default_diag_2p->FixParameter(2,fPD2+fPD2e);
  fit_default_diag_2p->FixParameter(3,fPD3);

  // Diagonalized default fit - P2
  TF1 *fit_default_diag_2m = new TF1("fit_default_diag_2m",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_2m->FixParameter(0,fPD0);
  fit_default_diag_2m->FixParameter(1,fPD1);
  fit_default_diag_2m->FixParameter(2,fPD2-fPD2e);
  fit_default_diag_2m->FixParameter(3,fPD3);

  // Diagonalized default fit + P3
  TF1 *fit_default_diag_3p = new TF1("fit_default_diag_3p",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_3p->FixParameter(0,fPD0);
  fit_default_diag_3p->FixParameter(1,fPD1);
  fit_default_diag_3p->FixParameter(2,fPD2);
  fit_default_diag_3p->FixParameter(3,fPD3+fPD3e);

  // Diagonalized default fit - P3
  TF1 *fit_default_diag_3m = new TF1("fit_default_diag_3m",fitQCD_diag,fFitXmin,fFitXmax,4);
  fit_default_diag_3m->FixParameter(0,fPD0);
  fit_default_diag_3m->FixParameter(1,fPD1);
  fit_default_diag_3m->FixParameter(2,fPD2);
  fit_default_diag_3m->FixParameter(3,fPD3-fPD3e);
  
  // Alternate fit function A (4 parameter fit function)
  TF1 *fit_A = new TF1("fit_A",fitQCD_A,fFitXmin,fFitXmax,4);
  fit_A->SetParameter(0,fPA0);
  fit_A->SetParameter(1,fPA1);
  fit_A->SetParameter(2,fPA2);
  fit_A->SetParameter(3,fPA3);
  fit_A->SetLineWidth(2);
  fit_A->SetLineStyle(1);
  fit_A->SetLineColor(kMagenta);

  // Results of the fit
  cout << "*********************************************************"<<endl;
  g->Fit("fit_A","R");
  cout << "*********************************************************"<<endl;
  cout << "Chi2/ndf: " << fit_A->GetChisquare() << "/" << fit_A->GetNDF() << " = " << fit_A->GetChisquare()/fit_A->GetNDF() << endl;
  cout << "Status: " << gMinuit->fCstatu.Data() << endl;
  cout << "*********************************************************"<<endl;

  // Alternate fit function B (3 parameter fit function)
  TF1 *fit_B = new TF1("fit_B",fitQCD_B,fFitXmin,fFitXmax,3);
  fit_B->SetParameter(0,fPB0);
  fit_B->SetParameter(1,fPB1);
  fit_B->SetParameter(2,fPB2);
  fit_B->SetLineWidth(2);
  fit_B->SetLineStyle(1);
  fit_B->SetLineColor(kCyan+3);

  // Results of the fit
  cout << "*********************************************************"<<endl;
  g->Fit("fit_B","R");
  cout << "*********************************************************"<<endl;
  cout << "Chi2/ndf: " << fit_B->GetChisquare() << "/" << fit_B->GetNDF() << " = " << fit_B->GetChisquare()/fit_B->GetNDF() << endl;
  cout << "Status: " << gMinuit->fCstatu.Data() << endl;
  cout << "*********************************************************"<<endl;

  int npoints = 1000;
  double step = (fFitXmax-fFitXmin)/(npoints-1);
  double vdx[npoints], vdvar[npoints], vAy[npoints], vBy[npoints];

  for(int i=0; i<npoints; ++i)
  {
    double mass = fFitXmin + i*step;
    vdx[i] = mass;
    
    vdvar[i] = sqrt( pow((fabs(fit_default_diag_0p->Eval(mass)-fit_default_diag->Eval(mass)) + fabs(fit_default_diag_0m->Eval(mass)-fit_default_diag->Eval(mass)))/2.,2) +
                     pow((fabs(fit_default_diag_1p->Eval(mass)-fit_default_diag->Eval(mass)) + fabs(fit_default_diag_1m->Eval(mass)-fit_default_diag->Eval(mass)))/2.,2) +
                     pow((fabs(fit_default_diag_2p->Eval(mass)-fit_default_diag->Eval(mass)) + fabs(fit_default_diag_2m->Eval(mass)-fit_default_diag->Eval(mass)))/2.,2) +
                     pow((fabs(fit_default_diag_3p->Eval(mass)-fit_default_diag->Eval(mass)) + fabs(fit_default_diag_3m->Eval(mass)-fit_default_diag->Eval(mass)))/2.,2) );

//    // alternate definition of the variation in the default fit
//    vdvar[i] = sqrt( pow(max(fabs(fit_default_diag_0p->Eval(mass)-fit_default_diag->Eval(mass)),fabs(fit_default_diag_0m->Eval(mass)-fit_default_diag->Eval(mass))),2) +
//                      pow(max(fabs(fit_default_diag_1p->Eval(mass)-fit_default_diag->Eval(mass)),fabs(fit_default_diag_1m->Eval(mass)-fit_default_diag->Eval(mass))),2) +
//                      pow(max(fabs(fit_default_diag_2p->Eval(mass)-fit_default_diag->Eval(mass)),fabs(fit_default_diag_2m->Eval(mass)-fit_default_diag->Eval(mass))),2) +
//                      pow(max(fabs(fit_default_diag_3p->Eval(mass)-fit_default_diag->Eval(mass)),fabs(fit_default_diag_3m->Eval(mass)-fit_default_diag->Eval(mass))),2) );

    vAy[i] = ( fit_A->Eval(mass)-fit_default_diag->Eval(mass) )/vdvar[i];
    vBy[i] = ( fit_B->Eval(mass)-fit_default_diag->Eval(mass) )/vdvar[i];
  }

  double vdvarx[] = {fFitXmin, fFitXmax, fFitXmax, fFitXmin};
  double vdvary[] = {-1., -1., 1., 1.};
  
  TPolyLine *l_default_var = new TPolyLine(4,vdvarx,vdvary);
  l_default_var->SetFillColor(kGray);

  TLine *l_default = new TLine(fFitXmin, 0, fFitXmax, 0);
  l_default->SetLineStyle(1);
  l_default->SetLineWidth(2);
  l_default->SetLineColor(kBlack);

  // not plotted, used for legend only
  TGraph *g_default = new TGraph(npoints,vdx,vdvar);
  g_default->SetLineWidth(2);
  g_default->SetLineStyle(1);
  g_default->SetFillColor(kGray);
  g_default->SetLineColor(kBlack);
  
  TGraph *g_A = new TGraph(npoints,vdx,vAy);
  g_A->SetLineWidth(2);
  g_A->SetLineStyle(1);
  g_A->SetLineColor(kMagenta);

  TGraph *g_B = new TGraph(npoints,vdx,vBy);
  g_B->SetLineWidth(2);
  g_B->SetLineStyle(1);
  g_B->SetLineColor(kCyan+3);

  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  c1->cd();

  g->Draw("AP");
  fit_default_diag->Draw("same");
  fit_default_diag_1p->Draw("same");
  fit_default_diag_1m->Draw("same");
  fit_A->Draw("same");
  fit_B->Draw("same");

  TCanvas *c2 = new TCanvas("c2", "", 800, 800);
  c2->cd();
  
  TH2D *bkg = new TH2D("bkg", "", 100, fFitXmin, fFitXmax, 100, -6, 6);
  bkg->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  bkg->GetYaxis()->SetTitle("(Fit-Default fit)/(Default fit variation)");
  bkg->GetXaxis()->SetNdivisions(1005);
  bkg->GetYaxis()->SetTitleOffset(1.05);

  bkg->Draw();
  l_default_var->Draw("Fsame");
  gPad->RedrawAxis();
  l_default->Draw("same");
  g_A->Draw("Lsame");
  g_B->Draw("Lsame");

  
  // Excluded window fits
  TF1 *fit_ew[bins_to_process-2];
  TGraph *g_ew[bins_to_process-2];

  for(int i=0; i<(bins_to_process-2); ++i)
  //for(int i=0; i<1; ++i)
  {
    double vx_ew[1000],vy_ew[1000],vexl_ew[1000],vexh_ew[1000],veyl_ew[1000],veyh_ew[1000];

    int k=0;
    for(int j=0; j<bins_to_process; ++j)
    {
      if( j>=i && j<(i+3) ) continue; // skip a 3-bin-wide window

      vx_ew[k]=vx[j];
      vy_ew[k]=vy[j];
      vexl_ew[k]=vexl[j];
      vexh_ew[k]=vexh[j];
      veyl_ew[k]=veyl[j];
      veyh_ew[k]=veyh[j];

      ++k;
    }

    TGraphAsymmErrors *ga_ew = new TGraphAsymmErrors((bins_to_process-3),vx_ew,vy_ew,vexl_ew,vexh_ew,veyl_ew,veyh_ew);

    // Default fit function
    fit_ew[i] = new TF1(Form("fit_ew%i",i),fitQCD,fFitXmin,fFitXmax,4);
    fit_ew[i]->SetParameter(0,fP0);
    fit_ew[i]->SetParameter(1,fP1);
    fit_ew[i]->SetParameter(2,fP2);
    fit_ew[i]->SetParameter(3,fP3);
    fit_ew[i]->SetLineWidth(1);
    fit_ew[i]->SetLineStyle(2);
    fit_ew[i]->SetLineColor(MyPalette[i]);

    // Results of the fit
    cout << "************" << Form(" Excluded window fit: %i ",(i+1)) << "********************"<<endl;
    ga_ew->Fit(Form("fit_ew%i",i),"R");
    cout << "*********************************************************"<<endl;
    cout << "Chi2/ndf: " << fit_ew[i]->GetChisquare() << "/" << fit_ew[i]->GetNDF() << " = " << fit_ew[i]->GetChisquare()/fit_ew[i]->GetNDF() << endl;
    cout << "Status: " << gMinuit->fCstatu.Data() << endl;
    cout << "*********************************************************"<<endl;

    c1->cd();
    //fit_ew[i]->Draw("same");

    double vgx_ew[npoints], vgy_ew[npoints];
    
    int m=0;
    for(int j=0; j<npoints; ++j)
    {
      if( fRestrict && !(vdx[j]>=h1_plot_r->GetBinLowEdge(i+1) && vdx[j]<h1_plot_r->GetBinLowEdge(i+4)) ) continue;
     
      vgx_ew[m] = vdx[j];
      vgy_ew[m] = ( fit_ew[i]->Eval(vdx[j])-fit_default_diag->Eval(vdx[j]) )/vdvar[j];

      ++m;
    }

    g_ew[i] = new TGraph(m,vgx_ew,vgy_ew);
    g_ew[i]->SetLineWidth(1);
    g_ew[i]->SetLineStyle(2);
    g_ew[i]->SetLineColor(MyPalette[i]);
    c2->cd();
    g_ew[i]->Draw("Lsame");

    delete ga_ew;
  }

  c1->cd();
  c1->SetLogy(1);
  //c1->SaveAs(("Fits_" + fOutputFile).c_str());

  c2->cd();
  
  TLegend *legend = new TLegend(.56,.75,.98,.95);
  //legend->SetBorderSize(0);
  legend->SetShadowColor(0);
  legend->SetFillColor(0);
  //legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(g_default, "Default fit with variation","fl");
  legend->AddEntry(g_A, "Alternate fit A (4 par.)","l");
  legend->AddEntry(g_B, "Alternate fit B (3 par.)","l");
  legend->AddEntry(g_ew[0],  "Excl. window fits (low)","l");
  legend->AddEntry(g_ew[14], "Excl. window fits (middle)","l");
  legend->AddEntry(g_ew[29], "Excl. window fits (high)","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.19,0.41, "CMS Preliminary");
  l1.DrawLatex(0.19,0.33, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.20,0.28, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.19,0.23, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.19,0.18, "Wide Jets");
  l1.SetTextSize(0.055);
  l1.DrawLatex(0.77,0.19, fLabel.c_str());

  c2->SetLogy(0);
  c2->SaveAs(fOutputFile.c_str());
}


void makePlots()
{
  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856,
                      1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558,
                      3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

  performFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "0 b-tags",
             "DijetMass_fit_variation_CSVL_0Tag_WideJets.eps",
             8.621, 4.951, 1.762, -0.0004518, 0.0462292, 0.00888221, 0.00191325, 5.87739e-08,
             6.30667e-05, 8.74030e+00, 5.17865e+00, -1.23820e-02,
             6.53929e-05, 8.84105e+00, 5.16671e+00,
             4.68740e-05, 8.52546e+00, 5.40954e+00, 4.96618e-02);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "1 b-tag",
             "DijetMass_fit_variation_CSVL_1Tag_WideJets.eps",
             7.956, 4.969, 1.747, -0.0001502, 0.0755236, 0.0147203, 0.00313092, 3.21208e-08,
             2.40482e-05, 8.48056e+00, 5.12759e+00, 3.85425e-02,
             2.15507e-05, 8.17614e+00, 5.16377e+00,
             1.54555e-05, 7.86212e+00, 5.40672e+00, 4.97943e-02);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "2 b-tags",
             "DijetMass_fit_variation_CSVL_2Tag_WideJets.eps",
             5.576, 6.026, 1.771, -3.309e-06, 0.241067, 0.0479228, 0.0101133, 1.97153e-09,
             7.86958e-07, 5.64247e+00, 5.47277e+00, -1.99716e-01,
             1.19464e-06, 6.76675e+00, 5.33372e+00,
             2.91837e-07, 5.46300e+00, 6.37524e+00, 2.15882e-01);

  // fit curves plotted in the excluded window only
  performFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "0 b-tags",
             "DijetMass_fit_variation_CSVL_0Tag_WideJets_restricted.eps",
             8.621, 4.951, 1.762, -0.0004518, 0.0462292, 0.00888221, 0.00191325, 5.87739e-08,
             6.30667e-05, 8.74030e+00, 5.17865e+00, -1.23820e-02,
             6.53929e-05, 8.84105e+00, 5.16671e+00,
             4.68740e-05, 8.52546e+00, 5.40954e+00, 4.96618e-02, 1);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "1 b-tag",
             "DijetMass_fit_variation_CSVL_1Tag_WideJets_restricted.eps",
             7.956, 4.969, 1.747, -0.0001502, 0.0755236, 0.0147203, 0.00313092, 3.21208e-08,
             2.40482e-05, 8.48056e+00, 5.12759e+00, 3.85425e-02,
             2.15507e-05, 8.17614e+00, 5.16377e+00,
             1.54555e-05, 7.86212e+00, 5.40672e+00, 4.97943e-02, 1);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass",
             43, xbins, 4976, 890, 4200, "2 b-tags",
             "DijetMass_fit_variation_CSVL_2Tag_WideJets_restricted.eps",
             5.576, 6.026, 1.771, -3.309e-06, 0.241067, 0.0479228, 0.0101133, 1.97153e-09,
             7.86958e-07, 5.64247e+00, 5.47277e+00, -1.99716e-01,
             1.19464e-06, 6.76675e+00, 5.33372e+00,
             2.91837e-07, 5.46300e+00, 6.37524e+00, 2.15882e-01, 1);
}

