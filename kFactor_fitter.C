/*###############################################################
  Bin number is a four digit number of the following form:  ABCD

  A: muon bin
    1: 0 muons
    2: >=1 muons

  B: max(|eta_J1|,|eta_J2|) bin
    1: |eta|<1.2
    2: 1.2=<|eta|<=2.5

  C: dijet mass bin
    1: 890<=mass<950 GeV
    2: 950<=mass<1200 GeV
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
TFile *file_csvl;
TFile *file_csvm;
double MCnHtotal=0.;
double MCnLtotal=0.;

void INIT(void)
{
  file_csvl=new TFile("CRAB_Jobs_MainAnalysis_CSVL_PUSFReweighted_PartonMatching_WideJets_EventBins/Final__histograms.root");
  file_csvm=new TFile("CRAB_Jobs_MainAnalysis_CSVM_PUSFReweighted_PartonMatching_WideJets_EventBins/Final__histograms.root");
  return;
}

Double_t getNumber(const string& fNumber, const Int_t fBin, const string& fBTagger)
{
  gROOT->SetBatch(kTRUE);

  string allowedNumbers = "n0, n1, n2, N00, N01, N02, N10, N11, N12, N20, N21, N22";
  if(allowedNumbers.find(fNumber)==string::npos)
  {
    cout<<"Allowed numbers are: n0, n1, n2, N00, N01, N02, N10, N11, N12, N20, N21, and N22."<<endl;
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

  if(fBTagger!="CSVL" && fBTagger!="CSVM")
  {
    cout<<"Allowed b-taggers are CSVL and CSVM."<<endl;
    return -1;
  }

  TFile* file=0;
  if(fBTagger=="CSVL") file=file_csvl;
  else if(fBTagger=="CSVM") file=file_csvm;
  file->cd();

  string muString, etaString;

  if(muBin==1) muString = "0mu";
  else if(muBin==2) muString = "ge1mu";

  if(etaBin==1) etaString = "maxEta_lt_1p2";
  else if(etaBin==2) etaString = "maxEta_ge_1p2";

  string histoName = (fNumber.find("n")!=string::npos ? "DATA__h2_" : "QCD_Pythia6__h2_") + fNumber + "_" + muString + "_" + etaString;

  TH2D *h2 = (TH2D*)file->Get(histoName.c_str());

  Int_t binsX[] = {890, 950, 1200, 1800, 2500, 6000};
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

          double N00=getNumber("N00", fBin, fBTagger);
          double N01=getNumber("N01", fBin, fBTagger);
          double N02=getNumber("N02", fBin, fBTagger);
          double N10=getNumber("N10", fBin, fBTagger);
          double N11=getNumber("N11", fBin, fBTagger);
          double N12=getNumber("N12", fBin, fBTagger);
          double N20=getNumber("N20", fBin, fBTagger);
          double N21=getNumber("N21", fBin, fBTagger);
          double N22=getNumber("N22", fBin, fBTagger);

          MCnLtotal+=(N00+N01+N02);
          MCnHtotal+=(N10+N11+N12+N20+N21+N22);

        }
      }
    }
  }
}

// nll calculation for a particular tagger
double nll(const string& fBTagger, double K, double muSF)
{
  if(DEBUG) cout << "K=" << K << "; muSF=" << muSF << endl;

  double nll=0;
  // scan over each bin of the likelihood for the given tagger
  for(int muBin=1; muBin<=2; muBin++) {
    for(int etaBin=1; etaBin<=2; etaBin++) {
      for(int massBin=1; massBin<=5; massBin++) {
        for(int pvBin=1; pvBin<=3; pvBin++) {
          int fBin=pvBin+massBin*10+etaBin*100+muBin*1000;
          double n0=getNumber("n0", fBin, fBTagger);
          double n1=getNumber("n1", fBin, fBTagger);
          double n2=getNumber("n2", fBin, fBTagger);
          double N00=getNumber("N00", fBin, fBTagger);
          double N01=getNumber("N01", fBin, fBTagger);
          double N02=getNumber("N02", fBin, fBTagger);
          double N10=getNumber("N10", fBin, fBTagger);
          double N11=getNumber("N11", fBin, fBTagger);
          double N12=getNumber("N12", fBin, fBTagger);
          double N20=getNumber("N20", fBin, fBTagger);
          double N21=getNumber("N21", fBin, fBTagger);
          double N22=getNumber("N22", fBin, fBTagger);

          double r=1+(1-K)*(MCnHtotal/MCnLtotal);

          double MCn0=K*N20+K*N10+r*N00;
          double MCn1=K*N21+K*N11+r*N01;
          double MCn2=K*N22+K*N12+r*N02;

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
  double K=par[0];
  double muSF=par[1];
  f=nll("CSVL", K, muSF)+nll("CSVM", K, muSF);
  return;
}

void fcn_CSVLOnly(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double K=par[0];
  double muSF=par[1];
  f=nll("CSVL", K, muSF);
  return;
}

void fcn_CSVMOnly(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double K=par[0];
  double muSF=par[1];
  f=nll("CSVM", K, muSF);
  return;
}

void minimize(void)
{
  INIT();
  INIT_HFFraction("CSVL"); // here it does not matter which b-tagger is used (it can be either "CSVM" or "CSVL")

  cout << "MCnTotal=" << (MCnLtotal+MCnHtotal) << endl;
  cout << "HF fraction=" << (MCnHtotal/(MCnLtotal+MCnHtotal)) << endl;

  // make some parameter choices
  //  muSF=0.9; // muon SF
  int whichTagger=1; // 0 = combined, 1=CSVL only, 2=CSVM only
  int printlevel=2; // -1 = suppressed, 0 = normal, 1 = verbose
  bool fixMuSF=true; // fix the SF to unity

  TMinuit t;
  t.SetPrintLevel(printlevel);
  if(whichTagger==0) t.SetFCN(fcn);
  if(whichTagger==1) t.SetFCN(fcn_CSVLOnly);
  if(whichTagger==2) t.SetFCN(fcn_CSVMOnly);
  t.DefineParameter(0, "HF K factor", 1.0, 0.1, 0.0, 5.0);
  t.DefineParameter(1, "muon SF", 1.0, 0.1, 0.0, 2.0);

  if(fixMuSF) t.FixParameter(1);

  t.Migrad();
//   t.mnmnos();
}
