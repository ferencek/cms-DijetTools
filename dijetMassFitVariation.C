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
#include "TGraph.h"
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


void drawFits(const Double_t fLumi, const Int_t fNpoints, const Double_t fFitXmin, const Double_t fFitXmax)
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


  Int_t MyPalette[31];
  Double_t red[]   = {0., 0., 0., 1., 1.};
  Double_t green[] = {0., 1., 1., 1., 0.};
  Double_t blue[]  = {1., 1., 0., 0., 0.};
  Double_t stop[]  = {0., 0.25, .50, 0.75, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(5, stop, red, green, blue, 31);
  for (int i=0; i<31; ++i) MyPalette[i] = FI+i;

  double p0_0[] = {  3.51794e-01,  5.34627e-01,  2.27731e-01,  2.28260e-01,  2.28252e-01,  2.28216e-01,  2.28122e-01,  2.28341e-01,  2.28146e-01,  2.28178e-01, 1.85143e-01, 2.07922e-01,  2.32111e-01,  2.31017e-01,  2.92651e-01,  4.33266e-01,  4.57946e-01,  3.21107e-01,  2.28579e-01,  2.28465e-01,  2.28483e-01,  2.28414e-01,  2.28399e-01,  2.28395e-01,  2.28375e-01,  2.28368e-01,  2.28362e-01,  2.28356e-01,  2.28347e-01,  2.28325e-01,  2.28318e-01  };
  double p1_0[] = {  8.90758e+00,  9.26177e+00,  8.50514e+00,  8.51137e+00,  8.51139e+00,  8.51123e+00,  8.51046e+00,  8.51553e+00,  8.51073e+00,  8.51110e+00, 8.36396e+00, 8.49504e+00,  8.57168e+00,  8.55201e+00,  8.79841e+00,  9.20916e+00,  9.25090e+00,  8.86881e+00,  8.51261e+00,  8.51226e+00,  8.51299e+00,  8.51246e+00,  8.51231e+00,  8.51229e+00,  8.51203e+00,  8.51192e+00,  8.51186e+00,  8.51182e+00,  8.51177e+00,  8.51159e+00,  8.51162e+00  };
  double p2_0[] = {  5.10321e+00,  4.77275e+00,  5.42108e+00,  5.42169e+00,  5.42174e+00,  5.42185e+00,  5.42198e+00,  5.42259e+00,  5.42199e+00,  5.42199e+00, 5.59246e+00, 5.51631e+00,  5.42627e+00,  5.42419e+00,  5.26067e+00,  4.99024e+00,  4.94588e+00,  5.18655e+00,  5.42064e+00,  5.42104e+00,  5.42129e+00,  5.42140e+00,  5.42141e+00,  5.42142e+00,  5.42141e+00,  5.42139e+00,  5.42140e+00,  5.42140e+00,  5.42142e+00,  5.42147e+00,  5.42151e+00  };
  double p3_0[] = { -1.24481e-02, -8.62189e-02,  5.17362e-02,  5.23870e-02,  5.24005e-02,  5.24205e-02,  5.24120e-02,  5.27720e-02,  5.24323e-02,  5.24534e-02, 9.07554e-02, 7.69416e-02,  5.66691e-02,  5.51616e-02,  2.35968e-02, -2.84159e-02, -3.82604e-02,  7.18425e-03,  5.21640e-02,  5.22548e-02,  5.23705e-02,  5.23721e-02,  5.23632e-02,  5.23665e-02,  5.23464e-02,  5.23359e-02,  5.23327e-02,  5.23312e-02,  5.23335e-02,  5.23420e-02,  5.23494e-02  };
  double p0_1[] = {  2.29245e-01,  4.48654e-01,  1.40928e-01,  1.41463e-01,  1.41450e-01,  1.41406e-01,  1.41371e-01,  1.41336e-01,  1.41381e-01,  1.41334e-01, 1.00755e-01, 1.10531e-01,  1.45285e-01,  1.44162e-01,  1.89573e-01,  2.93423e-01,  3.18963e-01,  2.12208e-01,  1.41851e-01,  1.41754e-01,  1.41652e-01,  1.41601e-01,  1.41613e-01,  1.41620e-01,  1.41602e-01,  1.41587e-01,  1.41568e-01,  1.41547e-01,  1.41530e-01,  1.41514e-01,  1.41617e-01  };
  double p1_1[] = {  8.92325e+00,  9.47911e+00,  8.46680e+00,  8.47716e+00,  8.47713e+00,  8.47675e+00,  8.47668e+00,  8.47972e+00,  8.47665e+00,  8.47586e+00, 8.23630e+00, 8.35643e+00,  8.57311e+00,  8.54165e+00,  8.83529e+00,  9.30929e+00,  9.37119e+00,  8.91470e+00,  8.47969e+00,  8.47907e+00,  8.47802e+00,  8.47774e+00,  8.47790e+00,  8.47827e+00,  8.47801e+00,  8.47788e+00,  8.47779e+00,  8.47769e+00,  8.47762e+00,  8.47754e+00,  8.47983e+00  };
  double p2_1[] = {  4.62123e+00,  4.08494e+00,  4.97078e+00,  4.97186e+00,  4.97194e+00,  4.97213e+00,  4.97240e+00,  4.97306e+00,  4.97231e+00,  4.97233e+00, 5.24751e+00, 5.19292e+00,  4.97890e+00,  4.97587e+00,  4.78976e+00,  4.49525e+00,  4.42812e+00,  4.69750e+00,  4.96984e+00,  4.97034e+00,  4.97070e+00,  4.97100e+00,  4.97097e+00,  4.97109e+00,  4.97111e+00,  4.97117e+00,  4.97128e+00,  4.97140e+00,  4.97151e+00,  4.97161e+00,  4.97181e+00  };
  double p3_1[] = { -1.05731e-01, -2.26709e-01, -3.71594e-02, -3.60717e-02, -3.60533e-02, -3.60211e-02, -3.59493e-02, -3.57946e-02, -3.59723e-02, -3.60210e-02, 2.56332e-02, 1.71321e-02, -2.93715e-02, -3.16805e-02, -6.69307e-02, -1.22173e-01, -1.37121e-01, -8.77410e-02, -3.64913e-02, -3.63885e-02, -3.63521e-02, -3.62813e-02, -3.62793e-02, -3.62226e-02, -3.62332e-02, -3.62254e-02, -3.61990e-02, -3.61745e-02, -3.61468e-02, -3.61187e-02, -3.59192e-02  };
  double p0_2[] = {  9.95357e-03,  1.56099e-01,  4.04380e-03,  4.83679e-03,  4.86055e-03,  4.96731e-03,  4.78769e-03,  4.89263e-03,  4.62509e-03,  4.22097e-03, 2.20968e-03, 2.32598e-03,  4.34732e-03,  4.26892e-03,  5.68114e-03,  9.14732e-03,  9.93709e-03,  6.40626e-03,  4.16562e-03,  4.12695e-03,  4.11463e-03,  4.11287e-03,  4.11939e-03,  4.12164e-03,  4.12334e-03,  4.12405e-03,  4.12343e-03,  4.12243e-03,  4.12129e-03,  4.11957e-03,  4.11812e-03  };
  double p1_2[] = {  7.20298e+00,  9.38678e+00,  6.33364e+00,  6.52945e+00,  6.53340e+00,  6.55384e+00,  6.51978e+00,  6.55070e+00,  6.48786e+00,  6.40357e+00, 5.97593e+00, 6.09321e+00,  6.58068e+00,  6.51583e+00,  6.86323e+00,  7.43503e+00,  7.45635e+00,  6.89574e+00,  6.39128e+00,  6.38289e+00,  6.37716e+00,  6.37762e+00,  6.37965e+00,  6.38083e+00,  6.38129e+00,  6.38155e+00,  6.38148e+00,  6.38136e+00,  6.38120e+00,  6.38063e+00,  6.38067e+00  };
  double p2_2[] = {  4.97690e+00,  2.72146e+00,  5.58156e+00,  5.46635e+00,  5.46246e+00,  5.44650e+00,  5.47376e+00,  5.45937e+00,  5.49935e+00,  5.56699e+00, 6.10173e+00, 6.09192e+00,  5.59992e+00,  5.59708e+00,  5.41885e+00,  5.11946e+00,  5.03767e+00,  5.30127e+00,  5.57671e+00,  5.58361e+00,  5.58450e+00,  5.58526e+00,  5.58435e+00,  5.58427e+00,  5.58399e+00,  5.58388e+00,  5.58403e+00,  5.58424e+00,  5.58447e+00,  5.58479e+00,  5.58518e+00  };
  double p3_2[] = { -5.68285e-02, -5.78912e-01,  4.50036e-02,  2.49305e-02,  2.40715e-02,  2.07789e-02,  2.64348e-02,  2.35231e-02,  3.17502e-02,  4.57711e-02, 1.66715e-01, 1.70387e-01,  6.31545e-02,  5.95588e-02,  2.93106e-02, -2.20244e-02, -4.29895e-02, -6.12648e-04,  4.77759e-02,  4.91979e-02,  4.91020e-02,  4.93570e-02,  4.92251e-02,  4.92838e-02,  4.92247e-02,  4.92054e-02,  4.92462e-02,  4.92911e-02,  4.93466e-02,  4.94216e-02,  4.95263e-02  };

  // B-only fit -- 0-tag
  TF1 *fit_B_0Tag = new TF1("fit_B_0Tag",fitQCD,fFitXmin,fFitXmax,4);
  fit_B_0Tag->FixParameter(0,2.28253e-01/fLumi);
  fit_B_0Tag->FixParameter(1,8.51128e+00);
  fit_B_0Tag->FixParameter(2,5.42169e+00);
  fit_B_0Tag->FixParameter(3,5.23788e-02);

  // B-only fit -- 1-tag
  TF1 *fit_B_1Tag = new TF1("fit_B_1Tag",fitQCD,fFitXmin,fFitXmax,4);
  fit_B_1Tag->FixParameter(0,1.41458e-01/fLumi);
  fit_B_1Tag->FixParameter(1,8.47705e+00);
  fit_B_1Tag->FixParameter(2,4.97186e+00);
  fit_B_1Tag->FixParameter(3,-3.60799e-02);

  // B-only fit -- 2-tag
  TF1 *fit_B_2Tag = new TF1("fit_B_2Tag",fitQCD,fFitXmin,fFitXmax,4);
  fit_B_2Tag->FixParameter(0,4.11527e-03/fLumi);
  fit_B_2Tag->FixParameter(1,6.37892e+00);
  fit_B_2Tag->FixParameter(2,5.58519e+00);
  fit_B_2Tag->FixParameter(3,4.94202e-02);

  TCanvas *c_0Tag = new TCanvas("c_0Tag", "", 800, 800);
  TCanvas *c_1Tag = new TCanvas("c_1Tag", "", 800, 800);
  TCanvas *c_2Tag = new TCanvas("c_2Tag", "", 800, 800);

  TLegend *legend = new TLegend(.2,.17,.5,.47);
  legend->SetBorderSize(0);
  legend->SetShadowColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);

  for(int i=0; i<31; i+=5)
  {
    // S+B fit -- 0-tag
    TF1 *fit_SB_0Tag = new TF1("fit_SB_0Tag",fitQCD,fFitXmin,fFitXmax,4);
    fit_SB_0Tag->FixParameter(0,p0_0[i]/fLumi);
    fit_SB_0Tag->FixParameter(1,p1_0[i]);
    fit_SB_0Tag->FixParameter(2,p2_0[i]);
    fit_SB_0Tag->FixParameter(3,p3_0[i]);

    // S+B fit -- 1-tag
    TF1 *fit_SB_1Tag = new TF1("fit_SB_1Tag",fitQCD,fFitXmin,fFitXmax,4);
    fit_SB_1Tag->FixParameter(0,p0_1[i]/fLumi);
    fit_SB_1Tag->FixParameter(1,p1_1[i]);
    fit_SB_1Tag->FixParameter(2,p2_1[i]);
    fit_SB_1Tag->FixParameter(3,p3_1[i]);

    // S+B fit -- 2-tag
    TF1 *fit_SB_2Tag = new TF1("fit_SB_2Tag",fitQCD,fFitXmin,fFitXmax,4);
    fit_SB_2Tag->FixParameter(0,p0_2[i]/fLumi);
    fit_SB_2Tag->FixParameter(1,p1_2[i]);
    fit_SB_2Tag->FixParameter(2,p2_2[i]);
    fit_SB_2Tag->FixParameter(3,p3_2[i]);

    double m[fNpoints];
    double SB_B_ratio_0Tag[fNpoints];
    double SB_B_ratio_1Tag[fNpoints];
    double SB_B_ratio_2Tag[fNpoints];

    double step = (fFitXmax-fFitXmin)/(fNpoints-1);

    for(int j=0; j<fNpoints; ++j)
    {
      double mass = fFitXmin + step*j;

      m[j] = mass;
      SB_B_ratio_0Tag[j] = fit_SB_0Tag->Eval(mass) / fit_B_0Tag->Eval(mass);
      SB_B_ratio_1Tag[j] = fit_SB_1Tag->Eval(mass) / fit_B_1Tag->Eval(mass);
      SB_B_ratio_2Tag[j] = fit_SB_2Tag->Eval(mass) / fit_B_2Tag->Eval(mass);
    }

    c_0Tag->cd();
    TGraph *g_SB_B_ratio_0Tag = new TGraph(fNpoints,m,SB_B_ratio_0Tag);
    g_SB_B_ratio_0Tag->SetLineWidth(3);
    g_SB_B_ratio_0Tag->SetLineStyle(int(i/5)+1);
    g_SB_B_ratio_0Tag->SetLineColor(MyPalette[i]);
    if(i==0)
    {
      g_SB_B_ratio_0Tag->GetXaxis()->SetTitle("Dijet Mass [GeV]");
      g_SB_B_ratio_0Tag->GetYaxis()->SetTitle("B Fit (S+B)/B Fit (B-only)");
      g_SB_B_ratio_0Tag->GetXaxis()->SetNdivisions(1005);
      g_SB_B_ratio_0Tag->GetYaxis()->SetRangeUser(0.9,1.05);
      g_SB_B_ratio_0Tag->GetYaxis()->SetTitleOffset(1.2);
      g_SB_B_ratio_0Tag->Draw("AL");
    }
    else g_SB_B_ratio_0Tag->Draw("L");

    if(i==0)  legend->AddEntry(g_SB_B_ratio_0Tag, "M=1 TeV","l");
    if(i==5)  legend->AddEntry(g_SB_B_ratio_0Tag, "M=1.5 TeV","l");
    if(i==10) legend->AddEntry(g_SB_B_ratio_0Tag, "M=2 TeV","l");
    if(i==15) legend->AddEntry(g_SB_B_ratio_0Tag, "M=2.5 TeV","l");
    if(i==20) legend->AddEntry(g_SB_B_ratio_0Tag, "M=3 TeV","l");
    if(i==25) legend->AddEntry(g_SB_B_ratio_0Tag, "M=3.5 TeV","l");
    if(i==30) legend->AddEntry(g_SB_B_ratio_0Tag, "M=4 TeV","l");

    c_1Tag->cd();
    TGraph *g_SB_B_ratio_1Tag = new TGraph(fNpoints,m,SB_B_ratio_1Tag);
    g_SB_B_ratio_1Tag->SetLineWidth(3);
    g_SB_B_ratio_1Tag->SetLineStyle(int(i/5)+1);
    g_SB_B_ratio_1Tag->SetLineColor(MyPalette[i]);
    if(i==0)
    {
      g_SB_B_ratio_1Tag->GetXaxis()->SetTitle("Dijet Mass [GeV]");
      g_SB_B_ratio_1Tag->GetYaxis()->SetTitle("B Fit (S+B)/B Fit (B-only)");
      g_SB_B_ratio_1Tag->GetXaxis()->SetNdivisions(1005);
      g_SB_B_ratio_1Tag->GetYaxis()->SetRangeUser(0.9,1.05);
      g_SB_B_ratio_1Tag->GetYaxis()->SetTitleOffset(1.2);
      g_SB_B_ratio_1Tag->Draw("AL");
    }
    else g_SB_B_ratio_1Tag->Draw("L");

    c_2Tag->cd();
    TGraph *g_SB_B_ratio_2Tag = new TGraph(fNpoints,m,SB_B_ratio_2Tag);
    g_SB_B_ratio_2Tag->SetLineWidth(3);
    g_SB_B_ratio_2Tag->SetLineStyle(int(i/5)+1);
    g_SB_B_ratio_2Tag->SetLineColor(MyPalette[i]);
    if(i==0)
    {
      g_SB_B_ratio_2Tag->GetXaxis()->SetTitle("Dijet Mass [GeV]");
      g_SB_B_ratio_2Tag->GetYaxis()->SetTitle("B Fit (S+B)/B Fit (B-only)");
      g_SB_B_ratio_2Tag->GetXaxis()->SetNdivisions(1005);
      g_SB_B_ratio_2Tag->GetYaxis()->SetRangeUser(0.9,1.05);
      g_SB_B_ratio_2Tag->GetYaxis()->SetTitleOffset(1.2);
      g_SB_B_ratio_2Tag->Draw("AL");
    }
    else g_SB_B_ratio_2Tag->Draw("L");

    delete fit_SB_0Tag;
    delete fit_SB_1Tag;
    delete fit_SB_2Tag;
  }

  TLine *line = new TLine(fFitXmin-200,1.,fFitXmax+200,1.);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(1);

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();

  c_0Tag->cd();
  legend->Draw();
  line->Draw("same");
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.2,0.88, "qq/bb, f_{b#bar{b}}=0.5");
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.73,0.22, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.75,0.88, "0 b-tags");

  c_0Tag->SaveAs("Fit_variations_0Tag.eps");


  c_1Tag->cd();
  legend->Draw();
  line->Draw("same");
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.2,0.88, "qq/bb, f_{b#bar{b}}=0.5");
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.73,0.22, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.75,0.88, "1 b-tag");
  
  c_1Tag->SaveAs("Fit_variations_1Tag.eps");


  c_2Tag->cd();
  legend->Draw();
  line->Draw("same");
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.2,0.88, "qq/bb, f_{b#bar{b}}=0.5");
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.73,0.22, "f_{b#bar{b}} = #frac{BR(X#rightarrowb#bar{b})}{BR(X#rightarrowjj)}");
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.75,0.88, "2 b-tags");

  c_2Tag->SaveAs("Fit_variations_2Tag.eps");

}


void makePlots()
{

  drawFits(4976, 100, 890, 4200);

}
