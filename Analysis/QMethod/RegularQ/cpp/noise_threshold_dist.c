#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );

//char root_file_name[128] = "run2c350MeV.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
char root_file_name1[128] = "run2D_thresh100_noisecorrection_multiplier_10.root";
char root_file_name2[128] = "run2D_thresh150_noisecorrection_multiplier_10.root";
char root_file_name3[128] = "run2D_thresh200_noisecorrection_multiplier_10.root";
char root_file_name4[128] = "run2D_thresh250_noisecorrection_multiplier_10.root";
char root_file_name5[128] = "run2D_thresh300_noisecorrection_multiplier_10.root";
char root_file_name6[128] = "run2D_thresh350_noisecorrection_multiplier_10.root";
char root_file_name7[128] = "run2D_thresh400_noisecorrection_multiplier_10.root";
char root_file_name8[128] = "run2D_thresh450_noisecorrection_multiplier_10.root";
char root_file_name9[128] = "run2D_thresh500_noisecorrection_multiplier_10.root";
char root_file_name10[128] = "run2D_thresh550_noisecorrection_multiplier_10.root";
char root_file_name11[128] = "run2D_thresh600_noisecorrection_multiplier_10.root";
char root_file_name12[128] = "run2D_thresh650_noisecorrection_multiplier_10.root";
char root_file_name13[128] = "run2D_thresh700_noisecorrection_multiplier_10.root";
char root_file_name14[128] = "run2D_thresh750_noisecorrection_multiplier_10.root";
char root_file_name15[128] = "run2D_thresh800_noisecorrection_multiplier_10.root";

/*
char root_file_name1[128] = "run2D_thresh300_noisecorrection_multiplier_5.root";
char root_file_name2[128] = "run2D_thresh300_noisecorrection_multiplier_6.root";
char root_file_name3[128] = "run2D_thresh300_noisecorrection_multiplier_7.root";
char root_file_name4[128] = "run2D_thresh300_noisecorrection_multiplier_8.root";
char root_file_name5[128] = "run2D_thresh300_noisecorrection_multiplier_9.root";
char root_file_name6[128] = "run2D_thresh300_noisecorrection_multiplier_10.root";
char root_file_name7[128] = "run2D_thresh300_noisecorrection_multiplier_11.root";
char root_file_name8[128] = "run2D_thresh300_noisecorrection_multiplier_13.root";
char root_file_name9[128] = "run2D_thresh300_noisecorrection_multiplier_15.root";
*/

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_res[25], *hthresh_[25], *hnoisy_[25];
//TH1 *hm;
TH1F *hlm, h0;
TCanvas *c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10;
TF1 *fit_func, *fit_func9, *fit_func13, *fit_func17, *fit_func21, *fit_func22, *fit_func26, *fit_func29;
TH1D *h_res9, *h_res13, *h_res17, *h_res21, *h_res22, *h_res26, *h_res29;
TH1  *hm9, *hm13, *hm17, *hm21, *hm22, *hm26, *hm29, *hm[25];
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;
Int_t m=0;
Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25],chi2[25],Entries[25];



Double_t fprec5(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym * cos(omega_a*time + phi));

}



void noise_threshold_dist()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;


 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");
  

     _file[1]=TFile::Open(root_file_name1);
      hthresh_[1]= (TH1D*)(_file[1]->FindObjectAny("hthreshold_"));
       hnoisy_[1]= (TH1D*)(_file[1]->FindObjectAny("HnoisyCaloXtal"));
      
      
     
      _file[2]=TFile::Open(root_file_name2);
      hthresh_[2]= (TH1D*)(_file[2]->FindObjectAny("hthreshold_"));
       hnoisy_[2]= (TH1D*)(_file[2]->FindObjectAny("HnoisyCaloXtal"));
     

      _file[3]=TFile::Open(root_file_name3);
      hthresh_[3]= (TH1D*)(_file[3]->FindObjectAny("hthreshold_"));
       hnoisy_[3]= (TH1D*)(_file[3]->FindObjectAny("HnoisyCaloXtal"));

     _file[4]=TFile::Open(root_file_name4);
      hthresh_[4]= (TH1D*)(_file[4]->FindObjectAny("hthreshold_"));
       hnoisy_[4]= (TH1D*)(_file[4]->FindObjectAny("HnoisyCaloXtal"));
     
     _file[5]=TFile::Open(root_file_name5);
      hthresh_[5]= (TH1D*)(_file[5]->FindObjectAny("hthreshold_"));
       hnoisy_[5]= (TH1D*)(_file[5]->FindObjectAny("HnoisyCaloXtal"));
     
     _file[6]=TFile::Open(root_file_name6);
      hthresh_[6]= (TH1D*)(_file[6]->FindObjectAny("hthreshold_"));
       hnoisy_[6]= (TH1D*)(_file[6]->FindObjectAny("HnoisyCaloXtal"));

     
        _file[7]=TFile::Open(root_file_name7);
      hthresh_[7]= (TH1D*)(_file[7]->FindObjectAny("hthreshold_"));
       hnoisy_[7]= (TH1D*)(_file[7]->FindObjectAny("HnoisyCaloXtal"));
     
        _file[8]=TFile::Open(root_file_name8);
      hthresh_[8]= (TH1D*)(_file[8]->FindObjectAny("hthreshold_"));
       hnoisy_[8]= (TH1D*)(_file[8]->FindObjectAny("HnoisyCaloXtal"));

     
        _file[9]=TFile::Open(root_file_name9);
      hthresh_[9]= (TH1D*)(_file[9]->FindObjectAny("hthreshold_"));
       hnoisy_[9]= (TH1D*)(_file[9]->FindObjectAny("HnoisyCaloXtal"));

     
           _file[10]=TFile::Open(root_file_name10);
      hthresh_[10]= (TH1D*)(_file[10]->FindObjectAny("hthreshold_"));
       hnoisy_[10]= (TH1D*)(_file[10]->FindObjectAny("HnoisyCaloXtal"));

     
        _file[11]=TFile::Open(root_file_name11);
      hthresh_[11]= (TH1D*)(_file[11]->FindObjectAny("hthreshold_"));
       hnoisy_[11]= (TH1D*)(_file[11]->FindObjectAny("HnoisyCaloXtal"));

     
        _file[12]=TFile::Open(root_file_name12);
      hthresh_[12]= (TH1D*)(_file[12]->FindObjectAny("hthreshold_"));
       hnoisy_[12]= (TH1D*)(_file[12]->FindObjectAny("HnoisyCaloXtal"));
     
        _file[13]=TFile::Open(root_file_name13);
      hthresh_[13]= (TH1D*)(_file[13]->FindObjectAny("hthreshold_"));
       hnoisy_[13]= (TH1D*)(_file[13]->FindObjectAny("HnoisyCaloXtal"));

     
        _file[14]=TFile::Open(root_file_name14);
      hthresh_[14]= (TH1D*)(_file[14]->FindObjectAny("hthreshold_"));
       hnoisy_[14]= (TH1D*)(_file[14]->FindObjectAny("HnoisyCaloXtal"));

          _file[15]=TFile::Open(root_file_name15);
      hthresh_[15]= (TH1D*)(_file[15]->FindObjectAny("hthreshold_"));
       hnoisy_[15]= (TH1D*)(_file[15]->FindObjectAny("HnoisyCaloXtal"));

       

   
    
  c2 = new TCanvas("c2","per multiplier residuals");
  c2->Divide(5,3);
  c3 = new TCanvas("c3","per multiplier fft");
  c3->Divide(5,3);
   
 for(int i=1;i<=15;i++)
   {
     c2->cd(i);
     hthresh_[i]->GetXaxis()->SetTitle("Energy [MeV]");
     hthresh_[i]->GetYaxis()->SetTitle("Counts [arb]");
     hthresh_[i]->Draw();
     c3->cd(i);
     hnoisy_[i]->GetXaxis()->SetTitle("xtal");
     hnoisy_[i]->GetYaxis()->SetTitle("calo");
     hnoisy_[i]->Draw("colz");
   }


     
 }

