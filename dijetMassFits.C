#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
// #include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
// #include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;

Double_t fitQCD1(Double_t *m, Double_t *p)
{
    double x=m[0]/7000.;
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

void performFit(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(5000);

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
  h1_plot_diff->SetMarkerStyle(20);
  h1_plot_diff->SetMarkerSize(0.8);
  h1_plot_diff->SetTitleOffset(1.4,"Y");

  // Fit to data
  TF1 *fit = new TF1("fit",fitQCD1,fFitXmin,fFitXmax,4); // 4 Par. Fit
//   gStyle->SetOptFit(1111);
  fit->SetParameter(0,1e-04);
  fit->SetParameter(1,1e+01);
  fit->SetParameter(2,4e+00);
  fit->SetParameter(3,-0.1e-01);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);
  h1_plot_diff->Fit("fit","R");
  TString status_default = gMinuit->fCstatu.Data();
  // Results of the fit
  cout << "*********************************************************"<<endl;
  Double_t chi_fit = fit->GetChisquare();
  Double_t ndf_fit = fit->GetNDF();
  cout << "Chi2/ndf: " << chi_fit << "/" << ndf_fit << " = " << chi_fit/ndf_fit << endl;
  cout << "Status: "<<status_default<<endl;
  cout << "*********************************************************"<<endl;
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot_diff->Draw();

  TLegend *legend = new TLegend(.7,.55,.85,.7);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(h1_plot_diff, "Data","lp");
  legend->AddEntry(fit, "Fit","l");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
//   l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.03);
  l1.DrawLatex(0.17,0.33, "CMS Preliminary");
  l1.DrawLatex(0.17,0.27, "#intLdt = 4.7 fb^{-1}");
  l1.DrawLatex(0.17,0.23, "#sqrt{s} = 7 TeV");
  l1.DrawLatex(0.17,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.17,0.15, fLabel.c_str());
  
  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());



  delete c;
  delete file;
}


void makePlots()
{

  Double_t xbins[] = {944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

  
  performFit("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass_pretag", 42, xbins,
             4679, 944, 4000, "M_{jj}>944 GeV", "DijetMass_fit.png");

  performFit("CRAB_Jobs_MainAnalysis_TCHEL_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 4000, "M_{jj}>944 GeV, TCHEL 0Tag", "DijetMass_fit_TCHEL_0Tag.png");

  performFit("CRAB_Jobs_MainAnalysis_TCHEL_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 4000, "M_{jj}>944 GeV, TCHEL 1Tag", "DijetMass_fit_TCHEL_1Tag.png");

  performFit("CRAB_Jobs_MainAnalysis_TCHEL_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 3600, "M_{jj}>944 GeV, TCHEL 2Tag", "DijetMass_fit_TCHEL_2Tag.png");

  performFit("CRAB_Jobs_MainAnalysis_SSVHPT_0Tag_PUReweighted_bPartonMatching/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 4000, "M_{jj}>944 GeV, SSVHPT 0Tag", "DijetMass_fit_SSVHPT_0Tag.png");

  performFit("CRAB_Jobs_MainAnalysis_SSVHPT_1Tag_PUReweighted_bPartonMatching_EventBins/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 3400, "M_{jj}>944 GeV, SSVHPT 1Tag", "DijetMass_fit_SSVHPT_1Tag.png");

  performFit("CRAB_Jobs_MainAnalysis_SSVHPT_2Tag_PUReweighted_bPartonMatching/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 42, xbins,
             4679, 944, 2600, "M_{jj}>944 GeV, SSVHPT 2Tag", "DijetMass_fit_SSVHPT_2Tag.png");
}

