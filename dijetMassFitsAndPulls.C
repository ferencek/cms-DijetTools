//Dijet Mass spectrum is compared to pythia QCD prediction, smooth fit and dijet resonance signals.
//Author: Sertac Ozturk (Cukurova Univ. & FNAL) //sertac@fnal.gov , sertac.ozturk@cern.ch32
//Modified by chiyoung Jeong -- chiyoung.jeong@gmail.com

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
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;

void setTDRStyle();

// QCD fit function -- alternate 4 parameter fit function -- also used for QCD fit.
Double_t fitQCD( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(m[0],p[2]);
}

// Default fit
Double_t fitQCD1( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

// QCD fit function -- alternate 3 parameter fit function -- also used for QCD fit.
Double_t fitQCD2( Double_t *m, Double_t *p)
{
  double x=m[0]/7000.;
  return p[0]*pow(1.-x,p[1])/pow(m[0],p[2]);
}

Double_t fitQCD3( Double_t *m, Double_t *p)
{
  return p[0]/pow(m[0]+p[1],p[2]);
}


void performFit(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile, const Double_t fP0=1e-04, const Double_t fP1=1e+01, const Double_t fP2=4e+00, const Double_t fP3=-0.1e-01)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
//   gROOT->ForceStyle();

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(30000);

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot = (TH1D*)file->Get(fPlot.c_str());
  TH1D *h1_plot_r = (TH1D*)h1_plot->Rebin(fNbins,"h1_plot_r",fBins);

  TH1D *h1_pulls = (TH1D*)h1_plot_r->Clone("h1_pulls");
  h1_pulls->Reset();

  double a = 0.3173/2;
  double n, nl, nh, dm, mass, xl, xh;
  double vx[1000],vy[1000],vexl[1000],vexh[1000],veyl[1000],veyh[1000];
  
  int bins_to_skip = 0;
  
  for(int i=h1_plot_r->GetNbinsX(); i>0; --i)
  {
      if( h1_plot_r->GetBinContent(i)!=0. ) break;
      ++bins_to_skip;
  }

  int bins_to_process = (h1_plot_r->GetNbinsX()-bins_to_skip);
  
  for(int i=0; i<bins_to_process; ++i)
  {
      n    = h1_plot_r->GetBinContent(i+1);
      dm   = h1_plot_r->GetBinWidth(i+1);
      mass = h1_plot_r->GetBinCenter(i+1);
      xl   = h1_plot_r->GetBinLowEdge(i+1);
      xh   = xl+dm;
      vx[i]   = (xl+xh)/2.;
      vexl[i] = dm/2.;
      vexh[i] = dm/2.;
      vy[i]   = n / (dm*fLumi);


      if (n<25)
      {
          nl = n-0.5*TMath::ChisquareQuantile(a,2*n);
          nh = 0.5*TMath::ChisquareQuantile(1-a,2*(n+1))-n;
          veyl[i] = nl/(fLumi*dm);
          veyh[i] = nh/(fLumi*dm);
      }
      else if (n>=25)
      {
          veyl[i] = sqrt(n)/(fLumi*dm);
          veyh[i] = sqrt(n)/(fLumi*dm);
      }
  }

  // data in the graph format
  TGraphAsymmErrors *g = new TGraphAsymmErrors(bins_to_process,vx,vy,vexl,vexh,veyl,veyh);

  // Fit to data
  TF1 *fit = new TF1("fit",fitQCD1,fFitXmin,fFitXmax,4); // 4 Par. Fit
  gStyle->SetOptFit(1111);
  fit->SetParameter(0,fP0);
  fit->SetParameter(1,fP1);
  fit->SetParameter(2,fP2);
  fit->SetParameter(3,fP3);
//   fit->SetParLimits(0, 0., 1.);
//   fit->SetParLimits(1, 0., 100.);
//   fit->SetParLimits(2, 0., 100.);
//   fit->SetParLimits(3,-5., 5.);
//   fit->FixParameter(3,0);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);
  g->Fit("fit","R");
  
  TString status_default = gMinuit->fCstatu.Data();
  // Results of the fit
  cout << "*********************************************************"<<endl;
  Double_t chi_fit = fit->GetChisquare();
  Double_t ndf_fit = fit->GetNDF();
  cout << "Chi2/ndf: " << chi_fit << "/" << ndf_fit << " = " << chi_fit/ndf_fit << endl;
  cout << "Status: "<<status_default<<endl;
  cout << "*********************************************************"<<endl;

  for(int i=0; i<bins_to_process; ++i)
    {
      double data = vy[i];
      double error = veyh[i];
      double m = vx[i];
      double fit_default = fit->Eval(m);

      if(error!=0.)
      {
        //cout << "m = " << m << std::endl;
        h1_pulls->SetBinContent(i+1,(data-fit_default)/error);
        h1_pulls->SetBinError(i+1,1.);
      }
    }
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->Divide(1,2,0,0,0);
  c->cd(1);

  TPad *p_1 = (TPad*)c->GetPad(1);
  p_1->SetPad(0.01,0.22,0.99,0.99);
  p_1->SetLogy();
  p_1->SetRightMargin(0.05);
  p_1->SetTopMargin(0.05);

  TH1F *vFrame = p_1->DrawFrame(700.0,3e-07,4337.0,10.0);
  vFrame->SetTitle("");
  vFrame->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  vFrame->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");

  g->Draw("P");
 
  TLegend *legend = new TLegend(.4,.7,.6,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g, "Data","lp");
  legend->AddEntry(fit, "Fit","l");
  legend->Draw();
// 
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.18,0.33, "CMS Preliminary");
  l1.DrawLatex(0.18,0.25, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.18,0.20, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.18,0.15, "|#eta| < 2.5, |#Delta#eta| < 1.3");
  l1.DrawLatex(0.18,0.10, ("Wide Jets"+ fLabel).c_str());

  c->cd(2);
  TPad *p_2 = (TPad*)c->GetPad(2);
  p_2->SetPad(0.01,0.02,0.99,0.22);
  p_2->SetBottomMargin(0.35);
  p_2->SetRightMargin(0.05);
  p_2->SetGridx();
  
  TH1F *vFrame2 = p_2->DrawFrame(700.0, -3.3, 4337.0, 3.3);
  vFrame2->SetTitle("");
  vFrame2->GetXaxis()->SetTitle("Dijet Mass (GeV)");
  vFrame2->GetYaxis()->SetTitle("Significance");
  vFrame2->GetYaxis()->SetTitleSize(0.12);
  vFrame2->GetYaxis()->SetLabelSize(0.10);
  vFrame2->GetYaxis()->SetTitleOffset(0.50);
  vFrame2->GetXaxis()->SetTitleOffset(0.90);
  vFrame2->GetXaxis()->SetTitleSize(0.18);
  vFrame2->GetXaxis()->SetLabelSize(0.18);

  h1_pulls->SetFillColor(kRed);

  h1_pulls->Draw("samehist");

  TLine *line = new TLine(700.,0,4337,0);
  line->Draw("");
  
  c->SaveAs(fOutputFile.c_str());

//   delete c;
//   delete file;
}


void makePlots()
{
  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856,
                      1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558,
                      3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};
                      
  performFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
             4976, 890, 4200, ", CSVL 0-tag", "DijetMass_fit_pulls_CSVL_0Tag_WideJets.png", 1e-03, 1e+01, 5e+00, 1e-01);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
             4976, 890, 3779, ", CSVL 1-tag", "DijetMass_fit_pulls_CSVL_1Tag_WideJets.png", 1e-03, 1e+01, 5e+00, 1e-01);

  performFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
             4976, 890, 3631, ", CSVL 2-tag", "DijetMass_fit_pulls_CSVL_2Tag_WideJets.png", 1e-05, 1e+00, 1e+01, 1e+00);


}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(63);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

//  tdrStyle->SetEndErrorSize(0);
//   tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.07);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.07);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}