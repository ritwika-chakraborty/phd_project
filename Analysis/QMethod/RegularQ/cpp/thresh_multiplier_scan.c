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
char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
char root_file_name1[128] = "run2D_thresh300_noisecorrection_multiplier_5.root";
char root_file_name2[128] = "run2D_thresh300_noisecorrection_multiplier_6.root";
char root_file_name3[128] = "run2D_thresh300_noisecorrection_multiplier_7.root";
char root_file_name4[128] = "run2D_thresh300_noisecorrection_multiplier_8.root";
char root_file_name5[128] = "run2D_thresh300_noisecorrection_multiplier_9.root";
char root_file_name6[128] = "run2D_thresh300_noisecorrection_multiplier_10.root";
char root_file_name7[128] = "run2D_thresh300_noisecorrection_multiplier_11.root";
char root_file_name8[128] = "run2D_thresh300_noisecorrection_multiplier_13.root";
char root_file_name9[128] = "run2D_thresh300_noisecorrection_multiplier_15.root";

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_res[25], *h_multi_1[25],*h_multi_2[25],*h_multi_3[25],*h_multi_4[25],*h_multi_5[25],*h_multi_6[25],*h_multi_7[25], *h_multi_8[25], *h_multi_9[25];
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
Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25],chi2[25], Entries[25];



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

Double_t fprec9(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    return norm * Ncbox * exp(-time/life) * (1 + asym * cos(omega_a*time + phi));

}
Double_t fprec13(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    return norm * Ncbox * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + phi + phicbox));

}


Double_t fprec17(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double A_vbo = par[13];

    double Tau_vbo = par[14];

    double omega_vbo = par[15];

    double phi_vbo = par[16];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    return norm * Ncbox * Nvoy1 * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + phi + phicbox));

}


Double_t fprec21(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double A_vbo = par[13];

    double Tau_vbo = par[14];

    double omega_vbo = par[15];

    double phi_vbo = par[16];

    double A_vbo2 = par[17];

    double Tau_vbo2 = par[18];

    double omega_vbo2 = par[19];

    double phi_vbo2 = par[20];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    return norm * Ncbox * Nvoy1 * Nvoy2 * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + phi + phicbox));

}


Double_t fprec22(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double A_vbo = par[13];

    double Tau_vbo = par[14];

    double omega_vbo = par[15];

    double phi_vbo = par[16];

    double A_vbo2 = par[17];

    double Tau_vbo2 = par[18];

    double omega_vbo2 = par[19];

    double phi_vbo2 = par[20];

    double lost_muon_amp = par[21];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));
    
    return norm * Ncbox * Nvoy1 * Nvoy2 * muloss * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + phi + phicbox));

}


Double_t fprec26(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double A_vbo = par[13];

    double Tau_vbo = par[14];

    double omega_vbo = par[15];

    double phi_vbo = par[16];

    double A_vbo2 = par[17];

    double Tau_vbo2 = par[18];

    double omega_vbo2 = par[19];

    double phi_vbo2 = par[20];

    double lost_muon_amp = par[21];

    double A_2cbo = par[22];

    double Tau_2cbo = par[23];

    double omega_2cbo = par[24];

    double phi_2cbo = par[25];
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double N2cboX =  (1 + A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time+phi_2cbo)));

    return norm * Ncbox * Nvoy1 * Nvoy2 * muloss * N2cboX * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + (phi + phicbox)));

}

Double_t fprec29(Double_t *x, Double_t *par)
{
     double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double A_vbo = par[13];

    double Tau_vbo = par[14];

    double omega_vbo = par[15];

    double phi_vbo = par[16];

    double A_vbo2 = par[17];

    double Tau_vbo2 = par[18];

    double omega_vbo2 = par[19];

    double phi_vbo2 = par[20];

    double lost_muon_amp = par[21];

    double A_pr = par[22];

    double Tau_pr = par[23];

    double omega_pr = par[24];

    double phi_pr = par[25];

    double A_vd = par[26];

    double tau_vd = par[27];

    double C_vd = par[28];
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double pedring =  (1 - A_pr * exp(-time/Tau_pr) * (cos(omega_pr*time+phi_pr)));

    double vertdrift =  (1 + (A_vd * exp(-time/tau_vd)) + C_vd);

    //  return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)))  * (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) * (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2))) * (1-A_pr * exp(-time/Tau_pr) * (cos(omega_pr*time+phi_pr))) * (1+(A_vd * exp(-time/tau_vd)) + C_vd);

   return norm  * Ncbox * Nvoy1 * Nvoy2 * muloss * pedring * vertdrift * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + (phi + phicbox)));

   // return norm * Ncbox *  exp(-time/life) * (1 + asym * cos(omega_a*time + phi));

    }

void thresh_multiplier_scan()

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
     for(int i=1;i<=24;i++)
      {
      h_multi_1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

     
      _file[2]=TFile::Open(root_file_name2);
     for(int i=1;i<=24;i++)
      {
      h_multi_2[i]= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }
     

      _file[3]=TFile::Open(root_file_name3);
     for(int i=1;i<=24;i++)
      {
      h_multi_3[i]= (TH1D*)(_file[3]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }


     _file[4]=TFile::Open(root_file_name4);
     for(int i=1;i<=24;i++)
      {
      h_multi_4[i]= (TH1D*)(_file[4]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

     
     _file[5]=TFile::Open(root_file_name5);
     for(int i=1;i<=24;i++)
      {
      h_multi_5[i]= (TH1D*)(_file[5]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

     
     _file[6]=TFile::Open(root_file_name6);
     for(int i=1;i<=24;i++)
      {
      h_multi_6[i]= (TH1D*)(_file[6]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

       _file[7]=TFile::Open(root_file_name7);
     for(int i=1;i<=24;i++)
      {
      h_multi_7[i]= (TH1D*)(_file[7]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

         _file[8]=TFile::Open(root_file_name8);
     for(int i=1;i<=24;i++)
      {
      h_multi_8[i]= (TH1D*)(_file[8]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

         _file[9]=TFile::Open(root_file_name9);
     for(int i=1;i<=24;i++)
      {
      h_multi_9[i]= (TH1D*)(_file[9]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }

   
    h[1]= new TH1D("calo histogram sum multiplier 1", "h_sum_multi_1", h_multi_1[1]->GetNbinsX(), 100001, 352001);
    h[1]->Sumw2(kTRUE);

    h[2]= new TH1D("calo histogram sum multiplier 2", "h_sum_multi_2", h_multi_2[1]->GetNbinsX(), 100001, 352001);
    h[2]->Sumw2(kTRUE);

    h[3]= new TH1D("calo histogram sum multiplier 3", "h_sum_multi_3", h_multi_3[1]->GetNbinsX(), 100001, 352001);
    h[3]->Sumw2(kTRUE);

    h[4]= new TH1D("calo histogram sum multiplier 4", "h_sum_multi_4", h_multi_4[1]->GetNbinsX(), 100001, 352001);
    h[4]->Sumw2(kTRUE);

    h[5]= new TH1D("calo histogram sum multiplier 5", "h_sum_multi_5", h_multi_5[1]->GetNbinsX(), 100001, 352001);
    h[5]->Sumw2(kTRUE);

    h[6]= new TH1D("calo histogram sum multiplier 6", "h_sum_multi_6", h_multi_6[1]->GetNbinsX(), 100001, 352001);
    h[6]->Sumw2(kTRUE);

    h[7]= new TH1D("calo histogram sum multiplier 7", "h_sum_multi_7", h_multi_7[1]->GetNbinsX(), 100001, 352001);
    h[7]->Sumw2(kTRUE);

    h[8]= new TH1D("calo histogram sum multiplier 8", "h_sum_multi_8", h_multi_8[1]->GetNbinsX(), 100001, 352001);
    h[8]->Sumw2(kTRUE);

    h[9]= new TH1D("calo histogram sum multiplier 9", "h_sum_multi_9", h_multi_9[1]->GetNbinsX(), 100001, 352001);
    h[9]->Sumw2(kTRUE);

 
   for(Int_t i=1; i<=24; i++)
    {
     h[1]->Add(h_multi_1[i],1);
     h[2]->Add(h_multi_2[i],1);
     h[3]->Add(h_multi_3[i],1);
     h[4]->Add(h_multi_4[i],1);
     h[5]->Add(h_multi_5[i],1);
     h[6]->Add(h_multi_6[i],1);
     h[7]->Add(h_multi_7[i],1);
     h[8]->Add(h_multi_8[i],1);
     h[9]->Add(h_multi_9[i],1);
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);


   cout<<h[1]->GetBinContent(500)<<" "<<endl;

  double binwidth=h[1]->GetBinWidth(1);
   int nbins=h[1]->GetNbinsX();
   double binlowedge=h[1]->GetBinLowEdge(1);
   double binhighedge=h[1]->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   //   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   for(int i=1; i<=9; i++)
     {
      h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }

   cout<<h[7]->GetBinWidth(1)<<" "<<h[7]->GetNbinsX()<<" "<<h[7]->GetBinLowEdge(1)<<" "<<(h[7]->GetBinLowEdge(h[7]->GetNbinsX())+h[7]->GetBinWidth(1))<<" "<<endl;

 
 


 fit_func= new TF1("fprec5", fprec5,  30000,  309000, 5);
 fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi");


   
  fit_func->SetNpx(1000000);


  c1 = new TCanvas("c1","hslice wiggle fit");



      fit_start=30000;
      fit_stop=300000;

  

  fit_func= new TF1("fprec5", fprec5,  30000,  309000, 5);
  fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi");
  fit_func->SetNpx(1000000);

  fit_func9= new TF1("fprec9", fprec9,  30000,  309000, 9);
  fit_func9->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func9->SetNpx(10000000);

  fit_func13= new TF1("fprec13", fprec13,  30000,  309000, 13);
  fit_func13->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func13->SetParName(9, "A2_cbo");
  fit_func13->SetParName(10, "phi2_cbo");
  fit_func13->SetParName(11, "A3_cbo");
  fit_func13->SetParName(12, "phi3_cbo");
  fit_func13->SetNpx(10000000);


  fit_func17= new TF1("fprec17", fprec17,  30000,  309000, 17);
  fit_func17->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func17->SetParName(9, "A2_cbo");
  fit_func17->SetParName(10, "phi2_cbo");
  fit_func17->SetParName(11, "A3_cbo");
  fit_func17->SetParName(12, "phi3_cbo");
  fit_func17->SetParName(13, "A_vw");
  fit_func17->SetParName(14, "Tau_vw");
  fit_func17->SetParName(15, "omega_vw");
  fit_func17->SetParName(16, "phi_vw");
  fit_func17->SetNpx(10000000);

  fit_func21= new TF1("fprec21", fprec21,  30000,  309000, 21);
  fit_func21->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func21->SetParName(9, "A2_cbo");
  fit_func21->SetParName(10, "phi2_cbo");
  fit_func21->SetParName(11, "A3_cbo");
  fit_func21->SetParName(12, "phi3_cbo");
  fit_func21->SetParName(13, "A_vw");
  fit_func21->SetParName(14, "Tau_vw");
  fit_func21->SetParName(15, "omega_vw");
  fit_func21->SetParName(16, "phi_vw");
  fit_func21->SetParName(17, "A_y");
  fit_func21->SetParName(18, "Tau_y");
  fit_func21->SetParName(19, "omega_y");
  fit_func21->SetParName(20, "phi_y");
  fit_func21->SetNpx(10000000);

  fit_func22= new TF1("fprec22", fprec22,  30000,  309000, 22);
  fit_func22->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func22->SetParName(9, "A2_cbo");
  fit_func22->SetParName(10, "phi2_cbo");
  fit_func22->SetParName(11, "A3_cbo");
  fit_func22->SetParName(12, "phi3_cbo");
  fit_func22->SetParName(13, "A_vw");
  fit_func22->SetParName(14, "Tau_vw");
  fit_func22->SetParName(15, "omega_vw");
  fit_func22->SetParName(16, "phi_vw");
  fit_func22->SetParName(17, "A_y");
  fit_func22->SetParName(18, "Tau_y");
  fit_func22->SetParName(19, "omega_y");
  fit_func22->SetParName(20, "phi_y");
  fit_func22->SetParName(21, "muloss_amp");
  fit_func22->SetNpx(10000000);


  fit_func26= new TF1("fprec26", fprec26,  30000,  309000, 26);
  fit_func26->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func26->SetParName(9, "A2_cbo");
  fit_func26->SetParName(10, "phi2_cbo");
  fit_func26->SetParName(11, "A3_cbo");
  fit_func26->SetParName(12, "phi3_cbo");
  fit_func26->SetParName(13, "A_2cbo");
  fit_func26->SetParName(14, "Tau_2cbo");
  fit_func26->SetParName(15, "omega_2cbo");
  fit_func26->SetParName(16, "phi_2cbo");
  fit_func26->SetParName(17, "A_vw");
  fit_func26->SetParName(18, "Tau_vw");
  fit_func26->SetParName(19, "omega_vw");
  fit_func26->SetParName(20, "phi_vw");
  fit_func26->SetParName(21, "muloss_amp");
  fit_func26->SetParName(22, "A_vbo");
  fit_func26->SetParName(23, "Tau_vbo");
  fit_func26->SetParName(24, "omega_vbo");
  fit_func26->SetParName(25, "phi_vbo");
  fit_func26->SetNpx(10000000);

      fit_func29= new TF1("fprec29", fprec29,  30000,  309000, 29);
      fit_func29->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
      fit_func29->SetParName(9, "A2_cbo");
      fit_func29->SetParName(10, "phi2_cbo");
      fit_func29->SetParName(11, "A3_cbo");
      fit_func29->SetParName(12, "phi3_cbo");
      fit_func29->SetParName(13, "A_vbo");
      fit_func29->SetParName(14, "Tau_vbo");
      fit_func29->SetParName(15, "omega_vbo");
      fit_func29->SetParName(16, "phi_vbo");
      fit_func29->SetParName(17, "A_vbo2");
      fit_func29->SetParName(18, "Tau_vbo2");
      fit_func29->SetParName(19, "omega_vbo2");
      fit_func29->SetParName(20, "phi_vbo2");
      fit_func29->SetParName(21, "muloss_amp");
      fit_func29->SetParName(22, "A_pr");
      fit_func29->SetParName(23, "Tau_pr");
      fit_func29->SetParName(24, "omega_pr");
      fit_func29->SetParName(25, "phi_pr");
      fit_func29->SetParName(26, "A_vd");
      fit_func29->SetParName(27, "Tau_vd");
      fit_func29->SetParName(28, "constant_vd");
      fit_func29->SetNpx(10000000);





 for(int i=1; i<=9; i++)
 {


   cout<<"multiplier "<<i<<endl;
   
 
  h[i]->Rebin(8);
  fit_func->SetParameters(h[i]->GetBinContent(83), 64440, 0.221, 0, 2.22);


  h[i]->GetYaxis()->SetRangeUser(-100., 40000000000.);
  h[i]->GetYaxis()->SetTitle("ADC counts");
  h[i]->GetXaxis()->SetTitle("time [ns]");
  h[i]->Fit(fit_func,"R","",fit_start,fit_stop);

 

  /*  fit_func9->SetParameter(0,fit_func->GetParameter(0));
  fit_func9->SetParameter(1,fit_func->GetParameter(1));
  fit_func9->SetParameter(2,fit_func->GetParameter(2));
  fit_func9->SetParameter(3,fit_func->GetParameter(3));
  fit_func9->SetParameter(4,fit_func->GetParameter(4));
  fit_func9->SetParameter(5,0.02);
  fit_func9->SetParameter(6,230000);
  fit_func9->SetParameter(7,0.00234);
  fit_func9->SetParameter(8,4.1);

  h[i]->Fit(fit_func9,"R","",fit_start,fit_stop);


   fit_func13->SetParameter(0,fit_func9->GetParameter(0));
  fit_func13->SetParameter(1,fit_func9->GetParameter(1));
  fit_func13->SetParameter(2,fit_func9->GetParameter(2));
  fit_func13->SetParameter(3,fit_func9->GetParameter(3));
  fit_func13->SetParameter(4,fit_func9->GetParameter(4));
  fit_func13->SetParameter(5,fit_func9->GetParameter(5));
  fit_func13->SetParameter(6,fit_func9->GetParameter(6));
  fit_func13->SetParameter(7,fit_func9->GetParameter(7));
  fit_func13->SetParameter(8,fit_func9->GetParameter(8));
  fit_func13->SetParameter(9, 0.023);
  fit_func13->SetParameter(10, 2.1);
  fit_func13->SetParameter(11, -0.003);
  fit_func13->SetParameter(12, 1.7);

  h[i]->Fit(fit_func13,"R","",fit_start,fit_stop);


      fit_func17->SetParameter(0,fit_func13->GetParameter(0));
      fit_func17->SetParameter(1,fit_func13->GetParameter(1));
      fit_func17->SetParameter(2,fit_func13->GetParameter(2));
      fit_func17->SetParameter(3,fit_func13->GetParameter(3));
      fit_func17->SetParameter(4,fit_func13->GetParameter(4));
      fit_func17->SetParameter(5,fit_func13->GetParameter(5));
      fit_func17->SetParameter(6,fit_func13->GetParameter(6));
      fit_func17->SetParameter(7,fit_func13->GetParameter(7));
      fit_func17->SetParameter(8,fit_func13->GetParameter(8));
      fit_func17->SetParameter(9,fit_func13->GetParameter(9));
      fit_func17->SetParameter(10,fit_func13->GetParameter(10));
      fit_func17->SetParameter(11,fit_func13->GetParameter(11));
      fit_func17->SetParameter(12,fit_func13->GetParameter(12));
      fit_func17->SetParameter(13,0.005);
      fit_func17->SetParameter(14,24000);
      fit_func17->SetParameter(15, 0.01404);
      fit_func17->SetParameter(16, 2.6);

  h[i]->Fit(fit_func17,"R","",fit_start,fit_stop);

      fit_func21->SetParameter(0,fit_func17->GetParameter(0));
      fit_func21->SetParameter(1,fit_func17->GetParameter(1));
      fit_func21->SetParameter(2,fit_func17->GetParameter(2));
      fit_func21->SetParameter(3,fit_func17->GetParameter(3));
      fit_func21->SetParameter(4,fit_func17->GetParameter(4));
      fit_func21->SetParameter(5,fit_func17->GetParameter(5));
      fit_func21->SetParameter(6,fit_func17->GetParameter(6));
      fit_func21->SetParameter(7,fit_func17->GetParameter(7));
      fit_func21->SetParameter(8,fit_func17->GetParameter(8));
      fit_func21->SetParameter(9,fit_func17->GetParameter(9));
      fit_func21->SetParameter(10,fit_func17->GetParameter(10));
      fit_func21->SetParameter(11,fit_func17->GetParameter(11));
      fit_func21->SetParameter(12,fit_func17->GetParameter(12));
      fit_func21->SetParameter(13,fit_func17->GetParameter(13));
      fit_func21->SetParameter(14,fit_func17->GetParameter(14));
      fit_func21->SetParameter(15,fit_func17->GetParameter(15));
      fit_func21->SetParameter(16,fit_func17->GetParameter(16));
      fit_func21->SetParameter(17, 0.00008);
      fit_func21->SetParameter(18, 600000);
      fit_func21->SetParameter(19, 0.01393);
      fit_func21->SetParameter(20, -0.3);

      h[i]->Fit(fit_func21,"R","",fit_start,fit_stop);


         
      fit_func22->SetParameter(0,fit_func21->GetParameter(0));
      fit_func22->SetParameter(1,fit_func21->GetParameter(1));
      fit_func22->SetParameter(2,fit_func21->GetParameter(2));
      fit_func22->SetParameter(3,fit_func21->GetParameter(3));
      fit_func22->SetParameter(4,fit_func21->GetParameter(4));
      fit_func22->SetParameter(5,fit_func21->GetParameter(5));
      fit_func22->SetParameter(6,fit_func21->GetParameter(6));
      fit_func22->SetParameter(7,fit_func21->GetParameter(7));
      fit_func22->SetParameter(8,fit_func21->GetParameter(8));
      fit_func22->SetParameter(9,fit_func21->GetParameter(9));
      fit_func22->SetParameter(10,fit_func21->GetParameter(10));
      fit_func22->SetParameter(11,fit_func21->GetParameter(11));
      fit_func22->SetParameter(12,fit_func21->GetParameter(12));
      fit_func22->SetParameter(13,fit_func21->GetParameter(13));
      fit_func22->SetParameter(14,fit_func21->GetParameter(14));
      fit_func22->SetParameter(15,fit_func21->GetParameter(15));
      fit_func22->SetParameter(16,fit_func21->GetParameter(16));
      fit_func22->SetParameter(17,fit_func21->GetParameter(17));
      fit_func22->SetParameter(18,fit_func21->GetParameter(18));
      fit_func22->SetParameter(19,fit_func21->GetParameter(19));
      fit_func22->SetParameter(20,fit_func21->GetParameter(20));
      fit_func22->SetParameter(21,0);

      h[i]->Fit(fit_func22,"R","",fit_start,fit_stop);



       fit_func26->SetParameter(0,fit_func22->GetParameter(0));
      fit_func26->SetParameter(1,fit_func22->GetParameter(1));
      fit_func26->SetParameter(2,fit_func22->GetParameter(2));
      fit_func26->SetParameter(3,fit_func22->GetParameter(3));
      fit_func26->SetParameter(4,fit_func22->GetParameter(4));
      fit_func26->SetParameter(5,fit_func22->GetParameter(5));
      fit_func26->SetParameter(6,fit_func22->GetParameter(6));
      fit_func26->SetParameter(7,fit_func22->GetParameter(7));
      fit_func26->SetParameter(8,fit_func22->GetParameter(8));
      fit_func26->SetParameter(9,fit_func22->GetParameter(9));
      fit_func26->SetParameter(10,fit_func22->GetParameter(10));
      fit_func26->SetParameter(11,fit_func22->GetParameter(11));
      fit_func26->SetParameter(12,fit_func22->GetParameter(12));
      fit_func26->SetParameter(13,fit_func22->GetParameter(13));
      fit_func26->SetParameter(14,fit_func22->GetParameter(14));
      fit_func26->SetParameter(15,fit_func22->GetParameter(15));
      fit_func26->SetParameter(16,fit_func22->GetParameter(16));
      fit_func26->SetParameter(17,fit_func22->GetParameter(17));
      fit_func26->SetParameter(18,fit_func22->GetParameter(18));
      fit_func26->SetParameter(19,fit_func22->GetParameter(19));
      fit_func26->SetParameter(20,fit_func22->GetParameter(20));
      fit_func26->SetParameter(21,fit_func22->GetParameter(21));
      fit_func26->SetParameter(22, 0.001);
      fit_func26->SetParameter(23, 300000);
      fit_func26->SetParameter(24, 0.01392);
      fit_func26->SetParameter(25, 2.5);

      
      
      h[i]->Fit(fit_func26,"R","",fit_start,fit_stop);


       fit_func29->SetNpx(10000000);
      fit_func29->SetParameter(0,fit_func26->GetParameter(0));
      fit_func29->SetParameter(1,fit_func26->GetParameter(1));
      fit_func29->SetParameter(2,fit_func26->GetParameter(2));
      fit_func29->SetParameter(3,fit_func26->GetParameter(3));
      fit_func29->SetParameter(4,fit_func26->GetParameter(4));
      fit_func29->SetParameter(5,fit_func26->GetParameter(5));
      fit_func29->SetParameter(6,fit_func26->GetParameter(6));
      fit_func29->SetParameter(7,fit_func26->GetParameter(7));
      fit_func29->SetParameter(8,fit_func26->GetParameter(8));
      fit_func29->SetParameter(9,fit_func26->GetParameter(9));
      fit_func29->SetParameter(10,fit_func26->GetParameter(10));
      fit_func29->SetParameter(11,fit_func26->GetParameter(11));
      fit_func29->SetParameter(12,fit_func26->GetParameter(12));
      fit_func29->SetParameter(13,fit_func26->GetParameter(13));
      fit_func29->SetParameter(14,fit_func26->GetParameter(14));
      fit_func29->SetParameter(15,fit_func26->GetParameter(15));
      fit_func29->SetParameter(16,fit_func26->GetParameter(16));
      fit_func29->SetParameter(17,fit_func26->GetParameter(17));
      fit_func29->SetParameter(18,fit_func26->GetParameter(18));
      fit_func29->SetParameter(19,fit_func26->GetParameter(19));
      fit_func29->SetParameter(20,fit_func26->GetParameter(20));
      fit_func29->SetParameter(21,fit_func26->GetParameter(21));
      fit_func29->SetParameter(22,fit_func26->GetParameter(22));
      fit_func29->SetParameter(23,fit_func26->GetParameter(23));
      fit_func29->SetParameter(24,fit_func26->GetParameter(24));
      fit_func29->SetParameter(25,fit_func26->GetParameter(25));
      fit_func29->SetParameter(26,-0.00842);
      fit_func29->SetParameter(27, 248023);
      fit_func29->SetParameter(28, 0.00569);

      h[i]->Fit(fit_func29,"R","",fit_start,fit_stop);
   
  */
      gStyle->SetOptFit(1111);

      blindR[m]=fit_func->GetParameter(3);
      //   phase[m]=fit_func22->GetParameter(4);
      chi2[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      dblindR[m]=fit_func->GetParError(3);
      //dphase[m]=fit_func22->GetParError(4);
      Entries[m]=h[i]->GetEntries();
      
      n[m]=m+4;
      if(n[m]==12){n[m]=13;}
      else if(n[m]==13){n[m]=15;}
      m=m+1;


      h_res[i] = new TH1D("residual hist","h_res",h[i]->GetNbinsX(),h[i]->GetBinLowEdge(1), (h[i]->GetBinLowEdge(h[i]->GetNbinsX())+h[i]->GetBinWidth(1)));
      
      for (int ibin = ((fit_start + 1- h[i]->GetBinLowEdge(1))/h[i]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h[i]->GetBinLowEdge(1))/h[i]->GetBinWidth(1))+1; ++ibin)
      {
       double res =  (h[i]->GetBinContent(ibin)- fit_func->Eval( h[i]->GetBinCenter(ibin) ) );
       if(h[i]->GetBinError(ibin)!=0){res=(res/h[i]->GetBinError(ibin));}
       h_res[i]->SetBinContent(ibin, (res)  );
      }


      h_res[i]->GetYaxis()->SetTitle("ADC counts");
      h_res[i]->GetXaxis()->SetTitle("time [ns]");


     hm[i] = h_res[i]->FFT(hm[i], "MAG");
     hm[i]->SetLineColor(kBlack);
     hm[i]->SetBins(h_res[i]->GetNbinsX(), 0, 1000/h_res[i]->GetBinWidth(1));
     hm[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
     hm[i]->GetXaxis()->SetTitle("Freq [MHz]");
     hm[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");


   
 }
 
  c2 = new TCanvas("c2","per multiplier residuals");
  c2->Divide(3,3);
  c3 = new TCanvas("c3","per multiplier fft");
  c3->Divide(3,3);
   
 for(int i=1;i<=9;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hm[i]->Draw();
   }

 c4=new TCanvas("c4","Blind R vs caloriemeter index");
 c4->Divide(1,3);
 c4->cd(1);
    gr1=new TGraphErrors(m,n,blindR,0,dblindR);
    gr1->SetTitle("R vs threshold multiplier");
    gr1->GetXaxis()->SetTitle("thresh multiplier");
    gr1->GetYaxis()->SetTitle("R [ppm]");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->RemovePoint(0);
    gr1->GetXaxis()->SetRangeUser(4,16);
    gr1->GetXaxis()->SetLabelSize(0.05);
    gr1->GetYaxis()->SetLabelSize(0.05);
    // gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr1->Draw();

    //   c5=new TCanvas("c5","phase vs caloriemeter index");
    c4->cd(2);
    gr2=new TGraphErrors(m,n,dblindR,0,0);
    gr2->SetTitle("dR vs threshold multiplier");
    gr2->GetXaxis()->SetTitle("thresh multiplier");
    gr2->GetYaxis()->SetTitle("dR [ppm]");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->RemovePoint(0);
    gr2->GetXaxis()->SetRangeUser(4,16);
     gr2->GetXaxis()->SetLabelSize(0.05);
    gr2->GetYaxis()->SetLabelSize(0.05);
    //  gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr2->Draw();


     c4->cd(3);
    gr3=new TGraphErrors(m,n,chi2,0,0);
    gr3->SetTitle("chi-squared vs thresh multiplier");
    gr3->GetXaxis()->SetTitle("chi-squared");
    gr3->GetYaxis()->SetTitle("thresh multiplier");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->RemovePoint(0);
    gr3->GetXaxis()->SetRangeUser(4,16);
     gr3->GetXaxis()->SetLabelSize(0.05);
    gr3->GetYaxis()->SetLabelSize(0.05);
    //  gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr3->Draw();

    c5=new TCanvas("c5","number of entries vs thresh multiplier");
    gr4=new TGraphErrors(m,n,Entries,0,0);
    gr4->SetTitle("Entries vs threshold multiplier");
    gr4->GetXaxis()->SetTitle("thresh multiplier");
    gr4->GetYaxis()->SetTitle("Entries");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->RemovePoint(0);
    gr4->GetXaxis()->SetRangeUser(4,16);
    gr4->GetXaxis()->SetLabelSize(0.05);
    gr4->GetYaxis()->SetLabelSize(0.05);
    // gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr4->Draw();


     
 }

