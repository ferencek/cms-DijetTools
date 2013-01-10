#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "exoStyle.C"

using namespace std;

Double_t fitQCD(Double_t *m, Double_t *p)
{
    double x=m[0]/7000.;
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}


void drawFit(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
             const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
             const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot = (TH1D*)file->Get(fPlot.c_str());

  TH1D *h1_plot_r = (TH1D*)h1_plot->Rebin(fNbins,"h1_plot_r",fBins);

  TH1D *h1_plot_diff = (TH1D*)h1_plot_r->Clone();
  h1_plot_diff->Reset();

  for(Int_t i=1; i<=h1_plot_r->GetNbinsX(); ++i)
  {
    Double_t n   = h1_plot_r->GetBinContent(i);
    Double_t err = h1_plot_r->GetBinError(i);
    Double_t dm  = h1_plot_r->GetBinWidth(i);

    h1_plot_diff->SetBinContent(i, n/(dm*fLumi));
    h1_plot_diff->SetBinError(i, err/(dm*fLumi));
  }

  h1_plot_diff->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  h1_plot_diff->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  h1_plot_diff->GetXaxis()->SetRangeUser(fFitXmin, 5000);
  h1_plot_diff->GetXaxis()->SetNdivisions(1005);
  h1_plot_diff->SetMarkerStyle(20);
  h1_plot_diff->SetMarkerSize(0.8);
  h1_plot_diff->SetTitleOffset(1.1,"Y");

  // Fit
  TF1 *fit = new TF1("fit",fitQCD,fFitXmin,fFitXmax,4); // 4 Par. Fit
  fit->FixParameter(0,fP0/fLumi);
  fit->FixParameter(1,fP1);
  fit->FixParameter(2,fP2);
  fit->FixParameter(3,fP3);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot_diff->Draw();
  fit->Draw("same");

  TLegend *legend = new TLegend(.55,.55,.85,.7);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_plot_diff, "Data","lp");
  legend->AddEntry(fit, "Background Fit","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.03);
  l1.DrawLatex(0.17,0.36, "CMS Preliminary");
  l1.DrawLatex(0.17,0.30, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.17,0.26, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.17,0.22,"Wide Jets");
  l1.DrawLatex(0.17,0.18, fLabel.c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void drawFitVariations(const Double_t fLumi, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                       const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3, const Int_t fSigma, const Double_t fPlus, const Double_t fMinus,
                       const string& fLeg1, const string& fLeg2, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gROOT->ForceStyle();

  TMatrixDSym covMatrix = TMatrixDSym(4);
  TVectorD eigenValues = TVectorD(4);
  TMatrixD eigenVectors = TMatrixD(4,4);

  double cov_matrix[] = {
      2.157e-05,  1.392e-04, -1.291e-05,  5.015e-06,
      1.392e-04,  1.479e-03,  2.486e-05,  6.315e-05,
     -1.291e-05,  2.486e-05,  5.229e-05,  1.532e-05,
      5.015e-06,  6.315e-05,  1.532e-05,  9.371e-06
  };

  covMatrix.SetMatrixArray(cov_matrix);

  const TMatrixDSymEigen eigen_data(covMatrix);
  eigenValues = eigen_data.GetEigenValues();
  eigenValues.Sqrt();
  //eigenValues.Print();
  eigenVectors = eigen_data.GetEigenVectors();

  double P0_p = fP0, P1_p = fP1, P2_p = fP2, P3_p = fP3, P0_m = fP0, P1_m = fP1, P2_m = fP2, P3_m = fP3;

  double d_p[4] = {0.}, d_m[4] = {0.};
  for(int v=0; v<4; ++v)
  {
    if( v!=(fSigma-1) ) continue;
    for(int k=0; k<4; ++k)
    {
      d_p[k]=fPlus*eigenValues(v)*eigenVectors[k][v];
      d_m[k]=fMinus*eigenValues(v)*eigenVectors[k][v];
    }
    P0_p += d_p[0];
    P1_p += d_p[1];
    P2_p += d_p[2];
    P3_p += d_p[3];
    P0_m += d_m[0];
    P1_m += d_m[1];
    P2_m += d_m[2];
    P3_m += d_m[3];
  }

  TF1 *fit = new TF1("fit",fitQCD,fXmin,fXmax,4);
  fit->FixParameter(0,fP0/fLumi);
  fit->FixParameter(1,fP1);
  fit->FixParameter(2,fP2);
  fit->FixParameter(3,fP3);

  TF1 *fit_p = new TF1("fit",fitQCD,fXmin,fXmax,4);
  fit_p->FixParameter(0,P0_p/fLumi);
  fit_p->FixParameter(1,P1_p);
  fit_p->FixParameter(2,P2_p);
  fit_p->FixParameter(3,P3_p);

  TF1 *fit_m = new TF1("fit",fitQCD,fXmin,fXmax,4);
  fit_m->FixParameter(0,P0_m/fLumi);
  fit_m->FixParameter(1,P1_m);
  fit_m->FixParameter(2,P2_m);
  fit_m->FixParameter(3,P3_m);
  
  int npoints = 1000;
  double step = (fXmax-fXmin)/(npoints-1);
  double vx[npoints], vy[npoints], vy_p[npoints], vy_m[npoints];

  for(int i=0; i<npoints; ++i)
  {
    double mass = fXmin + i*step;
    vx[i] = mass;

    vy[i] = 1.;
    vy_p[i] = fit_p->Eval(mass) / fit->Eval(mass);
    vy_m[i] = fit_m->Eval(mass) / fit->Eval(mass);
  }

  TGraph *g = new TGraph(npoints,vx,vy);
  g->SetLineWidth(2);
  g->SetLineStyle(1);
  g->SetLineColor(kBlack);

  TGraph *g_p = new TGraph(npoints,vx,vy_p);
  g_p->SetLineWidth(2);
  g_p->SetLineStyle(7);
  g_p->SetLineColor(kGreen+2);

  TGraph *g_m = new TGraph(npoints,vx,vy_m);
  g_m->SetLineWidth(2);
  g_m->SetLineStyle(3);
  g_m->SetLineColor(kRed);

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->SetTitleOffset(1.2,"Y");
  bkg->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  bkg->GetYaxis()->SetTitle("(Modified Fit)/(Default Fit)");
  bkg->Draw();

  g->Draw("L");
  g_p->Draw("L");
  g_m->Draw("L");

  TLegend *legend = new TLegend(.2,.65,.35,.8);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g, "Default","l");
  legend->AddEntry(g_p, fLeg1.c_str(),"l");
  legend->AddEntry(g_m, fLeg2.c_str(),"l");
  legend->Draw();

  //c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
}


void makePlots()
{

  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856,
                      1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558,
                      3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

  // Draw fit
  drawFit("Data_and_ResonanceShapes/Final__histograms_CSVL_0Tag_WideJets.root",
          "DATA__cutHisto_allPreviousCuts________DijetMass_pretag", 43, xbins,
          4976, 890, 4200, "M_{jj}>890 GeV", 3.27759e-01, 8.33793e+00, 5.37127e+00, 4.06333e-02, "DijetMass_Fit.eps");

  // Draw fit variations
  drawFitVariations(4976, 890, 4200, 0.97, 1.03, 3.27759e-01, 8.33793e+00, 5.37127e+00, 4.06333e-02, 1, 1, -1, "+1#sigma_{1}", "-1#sigma_{1}", "DijetMass_Fit_Variation_sigma1.eps");

  drawFitVariations(4976, 890, 4200, 0.995, 1.005, 3.27759e-01, 8.33793e+00, 5.37127e+00, 4.06333e-02, 2, 1, -1, "+1#sigma_{2}", "-1#sigma_{2}", "DijetMass_Fit_Variation_sigma2.eps");

  drawFitVariations(4976, 890, 4200, 0.99, 1.01, 3.27759e-01, 8.33793e+00, 5.37127e+00, 4.06333e-02, 3, 1, -1, "+1#sigma_{3}", "-1#sigma_{3}", "DijetMass_Fit_Variation_sigma3.eps");

  drawFitVariations(4976, 890, 4200, 0.995, 1.005, 3.27759e-01, 8.33793e+00, 5.37127e+00, 4.06333e-02, 4, 1, -1, "+1#sigma_{4}", "-1#sigma_{4}", "DijetMass_Fit_Variation_sigma4.eps");
}

