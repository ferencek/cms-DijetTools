#include <iostream>
#include <TFile.h>
#include <TH1D.h>


void convertResonanceShapes(const string& fFileIn, const string& fFileOut, Int_t fNbins, Double_t* fBinBoundaries, const string& fResonanceType)
{

   TFile *file_in = new TFile(fFileIn.c_str());
   TFile *file_out = new TFile(fFileOut.c_str(),"RECREATE");

   Int_t nKeys = file_in->GetListOfKeys()->GetEntries();

   for(Int_t i=0; i<nKeys; ++i)
   {
     string histname = file_in->GetListOfKeys()->At(i)->GetName();
     if(histname.find("cdf")!=string::npos || histname.find("Efficiency")!=string::npos) continue;
     std::cout << "Converting " << ("h_" + fResonanceType + "_" + histname.substr(8)) << std::endl;

     TH1D* h_shape = new TH1D(("h_" + fResonanceType + "_" + histname.substr(8)).c_str(), (fResonanceType + " Resonance Shape").c_str(), fNbins, fBinBoundaries);
     TH1D* h_pdf = new TH1D(Form("h_pdf_%i",i), "", 14000, 0, 14000);
     TH1D* h_cdf = new TH1D(("h_" + fResonanceType + "_" + histname.substr(8) + "_cdf").c_str(), (fResonanceType + " Resonance Shape CDF").c_str(), 14000, 0, 14000);

     TH1D* h_shape_in = (TH1D*)file_in->Get(histname.c_str());

     h_shape_in->Scale(5.); // Maxime's resonance shapes are normalized to 0.2 so they need to be scaled up by a factor of 5

     for(Int_t j=1; j<=h_shape_in->GetNbinsX(); ++j)
     {
       Int_t bin_min = h_pdf->GetXaxis()->FindBin(h_shape_in->GetXaxis()->GetBinLowEdge(j)+0.5);
       Int_t bin_max = h_pdf->GetXaxis()->FindBin(h_shape_in->GetXaxis()->GetBinUpEdge(j)-0.5);
       Double_t bin_content = h_shape_in->GetBinContent(j)/Double_t(bin_max-bin_min+1);
       for(Int_t b=bin_min; b<=bin_max; ++b)
         h_pdf->SetBinContent(b, bin_content);
     }

     for(Int_t j=1; j<=h_cdf->GetNbinsX(); ++j)
     {
       Int_t bin_min = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinLowEdge(j)+0.5);
       Int_t bin_max = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinUpEdge(j)-0.5);

       Double_t curr = 0.;
       for(Int_t b=bin_min; b<=bin_max; ++b)
         curr += h_pdf->GetBinContent(b);

       Double_t prev = h_cdf->GetBinContent(j-1);

       h_cdf->SetBinContent(j, prev+curr);
     }

     for(Int_t j=1; j<=h_shape->GetNbinsX(); ++j)
     {
       Int_t bin_min = h_pdf->GetXaxis()->FindBin(h_shape->GetXaxis()->GetBinLowEdge(j)+0.5);
       Int_t bin_max = h_pdf->GetXaxis()->FindBin(h_shape->GetXaxis()->GetBinUpEdge(j)-0.5);
       Double_t bin_content = h_pdf->Integral(bin_min,bin_max);
       h_shape->SetBinContent(j, bin_content);
     }

     file_out->cd();
     h_shape->Write();
     h_cdf->Write();

     delete h_pdf;
  }

  file_out->Close();
  file_in->Close();

  return;
}


void runConversion()
{
  const Int_t nBins = 103;

  Double_t binBoundaries[nBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
  354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
  1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
  4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430,
  10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};


  convertResonanceShapes("Resonance_Shapes_RSGraviton_2012_D6T_ak5_QQtoQQ_fat30.root", "Resonance_Shapes_qq.root", nBins, binBoundaries, "qq");
  convertResonanceShapes("Resonance_Shapes_RSGraviton_2012_D6T_ak5_GGtoGG_fat30.root", "Resonance_Shapes_gg.root", nBins, binBoundaries, "gg");
}
