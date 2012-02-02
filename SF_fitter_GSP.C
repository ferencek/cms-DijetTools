/*###############################################################
  Bin number is a four digit number of the following form:  ABCD

  A: muon bin
    1: 0 muons
    2: >=1 muons

  B: max(|eta_J1|,|eta_J2|) bin
    1: |eta|<1.2
    2: 1.2=<|eta|<=2.5

  C: dijet mass bin
    1: 944<=mass<1000 GeV
    2: 1000<=mass<1200 GeV
    3: 1200<=mass<1800 GeV
    4: 1800<=mass<2500 GeV
    5: 2500<=mass<6000 GeV
    
  D: primary vertex bin
    1: 0<=nPV<=5
    2: 6<=nPV<=8
    3: 9<=nPV<=40
  ###############################################################
*/

#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TMinuit.h"

using namespace std;

bool DEBUG=false;
TFile *file_tchpt;
TFile *file_ssvhpt;
double MCn_GSP=0.;
double MCn_noGSP=0.;

void INIT(void)
{
  file_tchpt=new TFile("CRAB_Jobs_DijetBBTag_TCHPT_SingleTag_PUReweighted_bPartonMatching_GSP_EventBins/Final__histograms.root");
  file_ssvhpt=new TFile("CRAB_Jobs_DijetBBTag_SSVHPT_SingleTag_PUReweighted_bPartonMatching_GSP_EventBins/Final__histograms.root");
  return;
}

Double_t getNumber(const string& fNumber, const string& fGSP_label, const Int_t fBin, const string& fBTagger)
{
  gROOT->SetBatch(kTRUE);

  string allowedNumbers = "n0, n1, n2, N000, N001, N002, N100, N101, N110, N111, N200, N210, N220";
  if(allowedNumbers.find(fNumber)==string::npos)
  {
    cout<<"Allowed numbers are: n0, n1, n2, N000, N001, N002, N100, N101, N110, N111, N200, N210, and N220."<<endl;
    return -1;
  }

  if(fGSP_label!="GSP" && fGSP_label!="noGSP")
  {
    cout<<"Allowed gluon splitting labels are GSP and noGSP."<<endl;
    return -1;
  }
  
  Int_t bin, muBin, etaBin, massBin, pvBin;

  bin = fBin;
  pvBin = bin%10;
  bin = bin/10;
  massBin = bin%10;
  bin = bin/10;
  etaBin = bin%10;
  bin = bin/10;
  muBin = bin%10;

  if( (muBin<1 || muBin>2) || (etaBin<1 || etaBin>2) || (massBin<1 || massBin>5) || (pvBin<1 || pvBin>3) )
  {
    cout<<"Invalid bin number. Valid bin numbers are of the following form: ABCD, where A,B=1,2; C=1,2,3,4,5; D=1,2,3."<<endl;
    assert(0);
  }

  if(fBTagger!="TCHPT" && fBTagger!="SSVHPT")
  {
    cout<<"Allowed b-taggers are TCHPT and SSVHPT."<<endl;
    return -1;
  }
  
  TFile* file=0;
  if(fBTagger=="TCHPT") file=file_tchpt;
  else if(fBTagger=="SSVHPT") file=file_ssvhpt;
  file->cd();

  string muString, etaString;

  if(muBin==1) muString = "0mu";
  else if(muBin==2) muString = "ge1mu";

  if(etaBin==1) etaString = "maxEta_lt_1p2";
  else if(etaBin==2) etaString = "maxEta_ge_1p2";

  string histoName = (fNumber.find("n")!=string::npos ? "DATA__h2_" + fNumber + "_" + muString + "_" + etaString : "QCD_Pythia6__h2_" + fNumber + "_" + fGSP_label + "_" + muString + "_" + etaString);
  
  TH2D *h2 = (TH2D*)file->Get(histoName.c_str());

  Int_t binsX[] = {944, 1000, 1200, 1800, 2500, 6000};
//   Int_t binsX[] = {944, 980, 1500, 2000, 2500, 6000};
  Int_t binsY[] = {0,6,9,41};
  
  Double_t Integral = h2->Integral(binsX[massBin-1]+1,binsX[massBin],binsY[pvBin-1]+1,binsY[pvBin]);

  return Integral;
}

void INIT_HFFraction(const string& fBTagger)
{
  for(int muBin=1; muBin<=2; muBin++) {
    for(int etaBin=1; etaBin<=2; etaBin++) {
      for(int massBin=1; massBin<=5; massBin++) {
        for(int pvBin=1; pvBin<=3; pvBin++) {
          int fBin=pvBin+massBin*10+etaBin*100+muBin*1000;

          double N000_GSP=getNumber("N000", "GSP", fBin, fBTagger);
          double N001_GSP=getNumber("N001", "GSP", fBin, fBTagger);
          double N002_GSP=getNumber("N002", "GSP", fBin, fBTagger);
          double N100_GSP=getNumber("N100", "GSP", fBin, fBTagger);
          double N101_GSP=getNumber("N101", "GSP", fBin, fBTagger);
          double N110_GSP=getNumber("N110", "GSP", fBin, fBTagger);
          double N111_GSP=getNumber("N111", "GSP", fBin, fBTagger);
          double N200_GSP=getNumber("N200", "GSP", fBin, fBTagger);
          double N210_GSP=getNumber("N210", "GSP", fBin, fBTagger);
          double N220_GSP=getNumber("N220", "GSP", fBin, fBTagger);

          double N000_noGSP=getNumber("N000", "noGSP", fBin, fBTagger);
          double N001_noGSP=getNumber("N001", "noGSP", fBin, fBTagger);
          double N002_noGSP=getNumber("N002", "noGSP", fBin, fBTagger);
          double N100_noGSP=getNumber("N100", "noGSP", fBin, fBTagger);
          double N101_noGSP=getNumber("N101", "noGSP", fBin, fBTagger);
          double N110_noGSP=getNumber("N110", "noGSP", fBin, fBTagger);
          double N111_noGSP=getNumber("N111", "noGSP", fBin, fBTagger);
          double N200_noGSP=getNumber("N200", "noGSP", fBin, fBTagger);
          double N210_noGSP=getNumber("N210", "noGSP", fBin, fBTagger);
          double N220_noGSP=getNumber("N220", "noGSP", fBin, fBTagger);

          MCn_GSP+=(N000_GSP+N001_GSP+N002_GSP+N100_GSP+N101_GSP+N110_GSP+N111_GSP+N200_GSP+N210_GSP+N220_GSP);
          MCn_noGSP+=(N000_noGSP+N001_noGSP+N002_noGSP+N100_noGSP+N101_noGSP+N110_noGSP+N111_noGSP+N200_noGSP+N210_noGSP+N220_noGSP);

        }
      }
    }
  }
}

// nll calculation for a particular tagger
double nll(const string& fBTagger, double SFh, double SFl, double K, double muSF)
{
  if(DEBUG) cout << "SFh=" << SFh << "; SFl=" << SFl << "; K=" << K << "; muSF=" << muSF << endl;

  double nll=0;
  // scan over each bin of the likelihood for the given tagger
  for(int muBin=1; muBin<=2; muBin++) {
    for(int etaBin=1; etaBin<=2; etaBin++) {
      for(int massBin=1; massBin<=5; massBin++) {
	for(int pvBin=1; pvBin<=3; pvBin++) {
	  int fBin=pvBin+massBin*10+etaBin*100+muBin*1000;
	  double n0=getNumber("n0", "GSP", fBin, fBTagger);
	  double n1=getNumber("n1", "GSP", fBin, fBTagger);
	  double n2=getNumber("n2", "GSP", fBin, fBTagger);
          
          double N000_GSP=getNumber("N000", "GSP", fBin, fBTagger);
          double N001_GSP=getNumber("N001", "GSP", fBin, fBTagger);
          double N002_GSP=getNumber("N002", "GSP", fBin, fBTagger);
          double N100_GSP=getNumber("N100", "GSP", fBin, fBTagger);
          double N101_GSP=getNumber("N101", "GSP", fBin, fBTagger);
          double N110_GSP=getNumber("N110", "GSP", fBin, fBTagger);
          double N111_GSP=getNumber("N111", "GSP", fBin, fBTagger);
          double N200_GSP=getNumber("N200", "GSP", fBin, fBTagger);
          double N210_GSP=getNumber("N210", "GSP", fBin, fBTagger);
          double N220_GSP=getNumber("N220", "GSP", fBin, fBTagger);

          double N000_noGSP=getNumber("N000", "noGSP", fBin, fBTagger);
          double N001_noGSP=getNumber("N001", "noGSP", fBin, fBTagger);
          double N002_noGSP=getNumber("N002", "noGSP", fBin, fBTagger);
          double N100_noGSP=getNumber("N100", "noGSP", fBin, fBTagger);
          double N101_noGSP=getNumber("N101", "noGSP", fBin, fBTagger);
          double N110_noGSP=getNumber("N110", "noGSP", fBin, fBTagger);
          double N111_noGSP=getNumber("N111", "noGSP", fBin, fBTagger);
          double N200_noGSP=getNumber("N200", "noGSP", fBin, fBTagger);
          double N210_noGSP=getNumber("N210", "noGSP", fBin, fBTagger);
          double N220_noGSP=getNumber("N220", "noGSP", fBin, fBTagger);

          double r=1+(1-K)*(MCn_GSP/MCn_noGSP);

          double MCn0=(K*N220_GSP+r*N220_noGSP)*(1-SFh)*(1-SFh)+(K*N111_GSP+r*N111_noGSP)*(1-SFh)*(1-SFl)+
                      (K*N002_GSP+r*N002_noGSP)*(1-SFl)*(1-SFl)+
                      ((K*N210_GSP+r*N210_noGSP)+(K*N110_GSP+r*N110_noGSP))*(1-SFh)+
                      ((K*N101_GSP+r*N101_noGSP)+(K*N001_GSP+r*N001_noGSP))*(1-SFl)+
                      K*(N200_GSP+N100_GSP+N000_GSP)+r*(N200_noGSP+N100_noGSP+N000_noGSP);
          double MCn1=(K*N220_GSP+r*N220_noGSP)*2*SFh*(1-SFh)+
                      (K*N111_GSP+r*N111_noGSP)*(SFh*(1-SFl)+(1-SFh)*SFl)+
                      (K*N002_GSP+r*N002_noGSP)*2*SFl*(1-SFl)+
                      ((K*N210_GSP+r*N210_noGSP)+(K*N110_GSP+r*N110_noGSP))*SFh+
                      ((K*N101_GSP+r*N101_noGSP)+(K*N001_GSP+r*N001_noGSP))*SFl;
          double MCn2=(K*N220_GSP+r*N220_noGSP)*SFh*SFh+(K*N111_GSP+r*N111_noGSP)*SFh*SFl+
                      (K*N002_GSP+r*N002_noGSP)*SFl*SFl;

	  if(muBin==2) { MCn0*=muSF; MCn1*=muSF; MCn2*=muSF; }

          if(DEBUG) {
            cout << "muBin=" << muBin << "; massBin=" << massBin << "; pvBin=" << pvBin << "; etaBin=" << etaBin << "; n0=" << n0 << "; n1=" << n1 << "; n2=" << n2 << "; MCn0=" << MCn0 << "; MCn1=" << MCn1 << "; MCn2=" << MCn2 << endl;
          }
	  
	  if(MCn0>0.) nll -= (n0*TMath::Log(MCn0) - TMath::LnGamma(n0+1) - MCn0);
	  if(MCn1>0.) nll -= (n1*TMath::Log(MCn1) - TMath::LnGamma(n1+1) - MCn1);
	  if(MCn2>0.) nll -= (n2*TMath::Log(MCn2) - TMath::LnGamma(n2+1) - MCn2);
	}
      }
    }
  }
  
  if(DEBUG) std::cout << "NLL=" << nll << std::endl;
  
  return nll;
}

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double SFh1=par[0];
  double SFl1=par[1];
  double SFh2=par[2];
  double SFl2=par[3];
  double K=par[4];
  double muSF=par[5];
  f=nll("TCHPT", SFh1, SFl1, K, muSF)+nll("SSVHPT", SFh2, SFl2, K, muSF);
  return;
}

void fcn_TCOnly(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double SFh1=par[0];
  double SFl1=par[1];
  double K=par[4];
  double muSF=par[5];
  f=nll("TCHPT", SFh1, SFl1, K, muSF);
  return;
}

void fcn_SSVOnly(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double SFh2=par[2];
  double SFl2=par[3];
  double K=par[4];
  double muSF=par[5];
  f=nll("SSVHPT", SFh2, SFl2, K, muSF);
  return;
}

void minimize(void)
{
  INIT();
  INIT_HFFraction("SSVHPT"); // here it does not matter which b-tagger is used (it can be either "SSVHPT" or "TCHPT")
  
  cout << "GSP HF fraction=" << (MCn_GSP/(MCn_GSP+MCn_noGSP)) << endl;
  
  // make some parameter choices
  //  muSF=0.9; // muon SF
  int whichTagger=0; // 0 = combined, 1=TC only, 2=SSV only
  int printlevel=2; // -1 = suppressed, 0 = normal, 1 = verbose
  bool fixMuSF=true; // fix the SF to unity

  TMinuit t;
  t.SetPrintLevel(printlevel);
  if(whichTagger==0) t.SetFCN(fcn);
  if(whichTagger==1) t.SetFCN(fcn_TCOnly);
  if(whichTagger==2) t.SetFCN(fcn_SSVOnly);
  t.DefineParameter(0, "TCHPT SFh", 1.0, 0.1, 0.0, 2.0);
  t.DefineParameter(1, "TCHPT SFl", 1.0, 0.1, 0.0, 2.0);
  t.DefineParameter(2, "SSVHPT SFh", 1.0, 0.1, 0.0, 2.0);
  t.DefineParameter(3, "SSVHPT SFl", 1.0, 0.1, 0.0, 2.0);
  t.DefineParameter(4, "HF K factor", 1.0, 0.1, 0.0, 5.0);
  t.DefineParameter(5, "muon SF", 1.0, 0.1, 0.0, 2.0);
  
  if(whichTagger==1) { t.FixParameter(2); t.FixParameter(3); }
  if(whichTagger==2) { t.FixParameter(0); t.FixParameter(1); }
  if(fixMuSF) t.FixParameter(5);
  
  t.Migrad();
//   t.mnmnos();
}
