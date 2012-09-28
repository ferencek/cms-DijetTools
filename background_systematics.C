#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"

using namespace std;

Double_t fitQCD(Double_t *m, Double_t *p)
{
    double x=m[0]/7000.;
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

TF1 func("fit",fitQCD,900.,5000.,4);

double evalFit(const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3,
               const Double_t fN0, const Double_t fN1, const Double_t fN2, const Double_t fN3,
	       const TVectorD& fEigenValues, const TMatrixD& fEigenVectors,
	       const Double_t fLumi, const Double_t fMass)
{
  Double_t p[4] = {0.};
  p[0]=fP0;
  p[1]=fP1;
  p[2]=fP2;
  p[3]=fP3;

  Double_t n[4] = {0.};
  n[0]=fN0;
  n[1]=fN1;
  n[2]=fN2;
  n[3]=fN3;

  Double_t g[4] = {0.};

  for(Int_t v=0; v<4; ++v)
  {
    for(Int_t k=0; k<4; ++k) g[k]=n[v]*fEigenValues(v)*fEigenVectors[k][v];
    p[0] += g[0];
    p[1] += g[1];
    p[2] += g[2];
    p[3] += g[3];
  }

  func.FixParameter(0,p[0]/fLumi);
  func.FixParameter(1,p[1]);
  func.FixParameter(2,p[2]);
  func.FixParameter(3,p[3]);

  return func.Eval(fMass);
}


void evalSyst_M(const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3,
	      const TVectorD& fEigenValues, const TMatrixD& fEigenVectors,
	      const Double_t fLumi, const Double_t fMass, Double_t& fRelSysUp, Double_t& fRelSysDown)
{
  fRelSysUp = -99.;
  fRelSysDown= 99.;
  Double_t initial = evalFit(fP0, fP1, fP2, fP3, 0, 0, 0, 0, fEigenValues, fEigenVectors, fLumi, fMass);

  for(Int_t i0=-1; i0<=1; ++i0)
  {
    for(Int_t i1=-1; i1<=1; ++i1)
    {
      for(Int_t i2=-1; i2<=1; ++i2)
      {
        for(Int_t i3=-1; i3<=1; ++i3)
        {
	  Double_t changed = evalFit(fP0, fP1, fP2, fP3, i0, i1, i2, i3, fEigenValues, fEigenVectors, fLumi, fMass);
	  Double_t variation = (changed - initial)/initial;
	  if(variation>0 && variation>fRelSysUp) fRelSysUp = variation;
	  if(variation<0 && variation<fRelSysDown) fRelSysDown = variation;
        }
      }
    }
  }
}


void evalSyst(const Double_t fP0, const Double_t fP1, const Double_t fP2, const Double_t fP3,
	      const TVectorD& fEigenValues, const TMatrixD& fEigenVectors,
	      const Double_t fLumi, const string& fTitle, const Int_t fSize, const Double_t *fMasses)
{
  cout << fTitle << endl << endl;
  for(Int_t i=0; i<fSize; ++i)
  {
    cout << "For m=" << fMasses[i] << " GeV" << endl;
    Double_t rel_syst_up=0, rel_syst_down=0;
    evalSyst_M(fP0, fP1, fP2, fP3, fEigenValues, fEigenVectors, fLumi, fMasses[i], rel_syst_up, rel_syst_down);
    cout << "rel_syst_up: " << rel_syst_up << endl;
    cout << "rel_syst_down: " << rel_syst_down << endl << endl;
  }
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
  TF1 *fit = new TF1("fit",fitQCD,fFitXmin,fFitXmax,4); // 4 Par. Fit
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

//   // Inclusive
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass_pretag", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, Inclusive", "DijetMass_drawFit_Inclusive_WideJets.png", 3.27721e-01, 8.33760e+00, 5.37123e+00, 4.06015e-02);
// 
//   // CSVL 0Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_0Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4200, "M_{jj}>890 GeV, CSVL 0Tag", "DijetMass_drawFit_CSVL_0Tag_WideJets.png", 2.28252e-01, 8.51127e+00, 5.42169e+00, 5.23794e-02);
// 
//   // CSVL 1Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_1Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4000, "M_{jj}>890 GeV, CSVL 1Tag", "DijetMass_drawFit_CSVL_1Tag_WideJets.png", 1.41458e-01, 8.47705e+00, 4.97186e+00, -3.60807e-02);
// 
//   // CSVL 2Tag
//   drawFit("CRAB_Jobs_MainAnalysis_CSVL_2Tag_PUSFkFReweighted_PartonMatching_WideJets/Final__histograms.root",
//           "DATA__cutHisto_allPreviousCuts________DijetMass", 43, xbins,
//           4976, 890, 4000, "M_{jj}>890 GeV, CSVL 2Tag", "DijetMass_drawFit_CSVL_2Tag_WideJets.png", 4.13861e-03, 6.38444e+00, 5.58112e+00, 4.86048e-02);


  TMatrixDSym covMatrix_incl = TMatrixDSym(4);
  TVectorD eigenValues_incl = TVectorD(4);
  TMatrixD eigenVectors_incl = TMatrixD(4,4);

  double cov_matrix_incl[] = {
    2.155e-05,  1.391e-04, -1.287e-05,  5.023e-06,
    1.391e-04,  1.478e-03,  2.533e-05,  6.331e-05,
   -1.287e-05,  2.533e-05,  5.235e-05,  1.536e-05,
    5.023e-06,  6.331e-05,  1.536e-05,  9.391e-06
  };

  covMatrix_incl.SetMatrixArray(cov_matrix_incl);
  const TMatrixDSymEigen eigen_incl(covMatrix_incl);
  eigenValues_incl = eigen_incl.GetEigenValues();
  eigenValues_incl.Sqrt();
  eigenVectors_incl = eigen_incl.GetEigenVectors();


  TMatrixDSym covMatrix_0tag = TMatrixDSym(4);
  TVectorD eigenValues_0tag = TVectorD(4);
  TMatrixD eigenVectors_0tag = TMatrixD(4,4);

  double cov_matrix_0tag[] = {
    1.478e-05,  1.381e-04, -1.239e-05,  5.050e-06,
    1.381e-04,  2.119e-03,  3.826e-05,  9.090e-05,
   -1.239e-05,  3.826e-05,  7.346e-05,  2.167e-05,
    5.050e-06,  9.090e-05,  2.167e-05,  1.331e-05
  };

  covMatrix_0tag.SetMatrixArray(cov_matrix_0tag);
  const TMatrixDSymEigen eigen_0tag(covMatrix_0tag);
  eigenValues_0tag = eigen_0tag.GetEigenValues();
  eigenValues_0tag.Sqrt();
  eigenVectors_0tag = eigen_0tag.GetEigenVectors();


  TMatrixDSym covMatrix_1tag = TMatrixDSym(4);
  TVectorD eigenValues_1tag = TVectorD(4);
  TMatrixD eigenVectors_1tag = TMatrixD(4,4);

  double cov_matrix_1tag[] = {
    1.536e-05,  2.287e-04, -2.128e-05,  8.333e-06,
    2.287e-04,  5.567e-03,  9.569e-05,  2.434e-04,
   -2.128e-05,  9.569e-05,  2.011e-04,  5.922e-05,
    8.333e-06,  2.434e-04,  5.922e-05,  3.621e-05
  };

  covMatrix_1tag.SetMatrixArray(cov_matrix_1tag);
  const TMatrixDSymEigen eigen_1tag(covMatrix_1tag);
  eigenValues_1tag = eigen_1tag.GetEigenValues();
  eigenValues_1tag.Sqrt();
  eigenVectors_1tag = eigen_1tag.GetEigenVectors();


  TMatrixDSym covMatrix_2tag = TMatrixDSym(4);
  TVectorD eigenValues_2tag = TVectorD(4);
  TMatrixD eigenVectors_2tag = TMatrixD(4,4);

  double cov_matrix_2tag[] = {
    1.339e-07,  6.682e-05, -6.544e-06,  2.432e-06,
    6.682e-05,  5.432e-02,  8.737e-04,  2.441e-03,
   -6.544e-06,  8.737e-04,  2.086e-03,  6.132e-04,
    2.432e-06,  2.441e-03,  6.132e-04,  3.723e-04
  };

  covMatrix_2tag.SetMatrixArray(cov_matrix_2tag);
  const TMatrixDSymEigen eigen_2tag(covMatrix_2tag);
  eigenValues_2tag = eigen_2tag.GetEigenValues();
  eigenValues_2tag.Sqrt();
  eigenVectors_2tag = eigen_2tag.GetEigenVectors();


  Double_t masses[] = {1000, 2000, 3000, 3500};

  evalSyst(3.27721e-01, 8.33760e+00, 5.37123e+00, 4.06015e-02, eigenValues_incl, eigenVectors_incl, 4976, ">> Incl", 4, masses);
  evalSyst(2.28252e-01, 8.51127e+00, 5.42169e+00, 5.23794e-02, eigenValues_0tag, eigenVectors_0tag, 4976, ">> 0 tag", 4, masses);
  evalSyst(1.41458e-01, 8.47705e+00, 4.97186e+00, -3.60807e-02, eigenValues_1tag, eigenVectors_1tag, 4976, ">> 1 tag", 4, masses);
  evalSyst(4.13861e-03, 6.38444e+00, 5.58112e+00, 4.86048e-02, eigenValues_2tag, eigenVectors_2tag, 4976, ">> 2 tag", 4, masses);
}

