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
#include "TFitResult.h"
#include "TMatrixD.h"

using namespace std;

Double_t fitQCD1(Double_t *m, Double_t *p)
{
    double x=m[0]/7000.;
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitQCD1_diag(Double_t *m, Double_t *p)
{
    double data[] = {
//         -0.507488, -0.508077,  0.350645,  0.601134,
//         -0.613174, -0.199281, -0.715634, -0.268651,
//         -0.218923,  0.744478, -0.238875,  0.583749,
//         -0.564399,  0.384574,  0.554847, -0.475081

        4.32e-05,  -0.0001246,  -0.0006045,          -1,
        0.9989,    -0.03008,     0.03618,   2.503e-05,
        0.01886,      0.9604,       0.278,  -0.0002869,
        0.04311,       0.277,     -0.9599,   0.0005476

//          4.32e-05,      0.9989,     0.01886,     0.04311,
//        -0.0001246,   -0.03008,      0.9604,       0.277,
//        -0.0006045,     0.03618,       0.278,     -0.9599,
//                -1,   2.503e-05,  -0.0002869,   0.0005476
        
    };
    
    TMatrixD d = TMatrixD(4,4);
    d.SetMatrixArray(data);
 
    double p0 = d(0,0)*p[0] + d(0,1)*p[1] + d(0,2)*p[2] + d(0,3)*p[3];
    double p1 = d(1,0)*p[0] + d(1,1)*p[1] + d(1,2)*p[2] + d(1,3)*p[3];
    double p2 = d(2,0)*p[0] + d(2,1)*p[1] + d(2,2)*p[2] + d(2,3)*p[3];
    double p3 = d(3,0)*p[0] + d(3,1)*p[1] + d(3,2)*p[2] + d(3,3)*p[3];
    
    double x=m[0]/7000.;
    return p0*pow(1.-x,p1)/pow(x,p2+p3*log(x));
}


void performFit(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile, const Double_t fP0=1e-04, const Double_t fP1=1e+01, const Double_t fP2=4e+00, const Double_t fP3=-0.1e-01)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(30000);

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
  fit->SetParameter(0,fP0);
  fit->SetParameter(1,fP1);
  fit->SetParameter(2,fP2);
  fit->SetParameter(3,fP3);
//   fit->SetParLimits(0, 0., 1.);
//   fit->SetParLimits(1, 0., 100.);
//   fit->SetParLimits(2, 0., 100.);
//   fit->SetParLimits(3,-5., 5.);
  fit->FixParameter(3,fP3);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);
  cout << "*********************************************************"<<endl;
  TFitResultPtr s = h1_plot_diff->Fit("fit","SR");
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
  l1.DrawLatex(0.17,0.27, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.17,0.23, "#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.17,0.19,"Wide Jets");
  else l1.DrawLatex(0.17,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.17,0.15, fLabel.c_str());
  
  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void performFit_diag(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
                     const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                     const string& fOutputFile, const Double_t fP0=1e-04, const Double_t fP1=1e+01, const Double_t fP2=4e+00, const Double_t fP3=-0.1e-01)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

//   cout << "Default number of iterations: " << TVirtualFitter::GetMaxIterations() << endl;
  TVirtualFitter::SetMaxIterations(30000);

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
  TF1 *fit = new TF1("fit",fitQCD1_diag,fFitXmin,fFitXmax,4); // 4 Par. Fit
//   gStyle->SetOptFit(1111);
  fit->SetParameter(0,fP0);
  fit->SetParameter(1,fP1);
  fit->SetParameter(2,fP2);
  fit->SetParameter(3,fP3);
//   fit->FixParameter(0,fP0);
//   fit->FixParameter(1,fP1);
//   fit->FixParameter(2,fP2);
//   fit->FixParameter(3,fP3);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);
  cout << "*********************************************************"<<endl;
  TFitResultPtr s = h1_plot_diff->Fit("fit","SR");
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
  l1.DrawLatex(0.17,0.27, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.17,0.23, "#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.17,0.19,"Wide Jets");
  else l1.DrawLatex(0.17,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.17,0.15, fLabel.c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void drawFit(const string& fInputFile, const string& fPlot, const Int_t fNbins, const Double_t *fBins,
                const Double_t fLumi, const Double_t fFitXmin, const Double_t fFitXmax, const string& fLabel,
                const string& fOutputFile, const Double_t fP0=1e-04, const Double_t fP1=1e+01, const Double_t fP2=4e+00, const Double_t fP3=-0.1e-01)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.07);

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

  // Fit
  TF1 *fit = new TF1("fit",fitQCD1,fFitXmin,fFitXmax,4); // 4 Par. Fit
//   gStyle->SetOptFit(1111);
//   fit->SetParameter(0,fP0/fLumi);
//   fit->SetParameter(1,fP1);
//   fit->SetParameter(2,fP2);
//   fit->SetParameter(3,fP3);
  fit->FixParameter(0,fP0/fLumi);
  fit->FixParameter(1,fP1);
  fit->FixParameter(2,fP2);
  fit->FixParameter(3,fP3);
  fit->SetLineWidth(2);
  fit->SetLineColor(kBlue);
//   h1_plot_diff->Fit("fit","R");
//   TString status_default = gMinuit->fCstatu.Data();
//   // Results of the fit
//   cout << "*********************************************************"<<endl;
//   Double_t chi_fit = fit->GetChisquare();
//   Double_t ndf_fit = fit->GetNDF();
//   cout << "Chi2/ndf: " << chi_fit << "/" << ndf_fit << " = " << chi_fit/ndf_fit << endl;
//   cout << "Status: "<<status_default<<endl;
//   cout << "*********************************************************"<<endl;

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_plot_diff->Draw();
  fit->Draw("same");

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
  l1.DrawLatex(0.17,0.27, "#intLdt = 5 fb^{-1}");
  l1.DrawLatex(0.17,0.23, "#sqrt{s} = 7 TeV");
  if(fInputFile.find("WideJets")!=string::npos) l1.DrawLatex(0.17,0.19,"Wide Jets");
  else l1.DrawLatex(0.17,0.19, "Anti-k_{T} R = 0.7 PFJets");
  l1.DrawLatex(0.17,0.15, fLabel.c_str());

  c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{

  Double_t xbins[] = {890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856,
                      1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558,
                      3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6000};

// #############
// ## Wide jets
// #############

// ## CSVL
                      
//   performFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass_pretag", 43, xbins,
//              4976, 890, 4200, "M_{jj}>890 GeV", "DijetMass_fit_WideJets.png", 1e-03, 1e+01, 5e+00, 1e-01);
  // CSVL 0Tag
  performFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
             "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
             4976, 890, 4200, "M_{jj}>890 GeV, CSVL 0Tag", "DijetMass_fit_CSVL_0Tag_WideJets.png", 1.51777e-04, 9.70316e+00, 4.57886e+00, -1.14765e-01);

  performFit_diag("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
                  "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
                  4976, 890, 4200, "M_{jj}>890 GeV, CSVL 0Tag", "DijetMass_fit_diag_CSVL_0Tag_WideJets.png", 9.774, 4.074, 1.734, -0.001285);
//   performFit_diag("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//                   "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//                   4976, 890, 4200, "M_{jj}>890 GeV, CSVL 0Tag", "DijetMass_fit_diag_CSVL_0Tag_WideJets.png",  0.1108, -0.1261, 10.59, -1.708);
  
//   // CSVL 1Tag
//   performFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//              4976, 890, 4000, "M_{jj}>890 GeV, CSVL 1Tag", "DijetMass_fit_CSVL_1Tag_WideJets.png", 1.83940e-05, 8.04984e+00, 5.28857e+00, 2.74169e-02);
//   // CSVL 2Tag
//   performFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//              4976, 890, 4000, "M_{jj}>890 GeV, CSVL 2Tag", "DijetMass_fit_CSVL_2Tag_WideJets.png", 6.23477e-10, -4.97863e-01, 1.07937e+01, 1.10468e+00);
                      
//   // Draw fit
//   //CSVL 0Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, CSVL 0Tag", "DijetMass_drawFit_CSVL_0Tag_WideJets.png", 2.28123e-01, 8.50964e+00, 5.42146e+00, 5.21746e-02);
//   // CSVL 1Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4000, "M_{jj}>890 GeV, CSVL 1Tag", "DijetMass_drawFit_CSVL_1Tag_WideJets.png", 1.41441e-01, 8.47680e+00, 4.97182e+00, -3.61227e-02);
//   // CSVL 2Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, CSVL 2Tag", "DijetMass_drawFit_CSVL_2Tag_WideJets.png", 4.11085e-03, 6.37674e+00, 5.58509e+00, 4.91493e-02);

// ## CSVM

//   // CSVM 0Tag
//   performFit("CRAB_Jobs_MainAnalysis_CSVM_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//              4976, 890, 4200, "M_{jj}>890 GeV, CSVM 0Tag", "DijetMass_fit_CSVM_0Tag_WideJets.png", 2.76841e-04, 9.93281e+00, 4.28966e+00, -1.70391e-01);
//   // CSVM 1Tag
//   performFit("CRAB_Jobs_MainAnalysis_CSVM_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//              4976, 890, 4000, "M_{jj}>890 GeV, CSVM 1Tag", "DijetMass_fit_CSVM_1Tag_WideJets.png", 4.22853e-06, 7.43758e+00, 5.62396e+00, 6.27649e-02);
//   // CSVM 2Tag
//   performFit("CRAB_Jobs_MainAnalysis_CSVM_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//              "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//              4976, 890, 4000, "M_{jj}>890 GeV, CSVM 2Tag", "DijetMass_fit_CSVM_2Tag_WideJets.png", 8.42408e-07, 1.07030e+01, 4.87934e+00, 0);

//   // Draw fit
//   //CSVM 0Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVM_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, CSVM 0Tag", "DijetMass_drawFit_CSVM_0Tag_WideJets.png", 4.18653e-01, 8.74889e+00, 5.12920e+00, -3.77561e-03);
//   // CSVM 1Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVM_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4000, "M_{jj}>890 GeV, CSVM 1Tag", "DijetMass_drawFit_CSVM_1Tag_WideJets.png", 9.63612e-03, 6.58225e+00, 6.14230e+00, 1.59159e-01);
//   // CSVM 2Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVM_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, CSVM 2Tag", "DijetMass_drawFit_CSVM_2Tag_WideJets.png", 1.53383e-03, 8.74605e+00, 5.23789e+00, 0);
                                            
}

