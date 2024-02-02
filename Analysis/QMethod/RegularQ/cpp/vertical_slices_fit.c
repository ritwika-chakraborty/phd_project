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

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_res[25];
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
Double_t blindR[25],dblindR[25],n[25];



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

    double A_pr = par[22];

    double Tau_pr = par[23];

    double omega_pr = par[24];

    double phi_pr = par[25];
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double pedring =  (1 - A_pr * exp(-time/Tau_pr) * (cos(omega_pr*time+phi_pr)));

    return norm * Ncbox * Nvoy1 * Nvoy2 * muloss * pedring * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + (phi + phicbox)));

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

void vertical_slices_fit()

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

  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzerDB",dir[1]);
  dir[1]->GetObject("qHist1D_sig_vslice_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzerDB",dir[2]);
  dir[2]->GetObject("qHist1D_sig_vslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzerDB",dir[3]);
  dir[3]->GetObject("qHist1D_sig_vslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzerDB",dir[4]);
  dir[4]->GetObject("qHist1D_sig_vslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzerDB",dir[5]);
  dir[5]->GetObject("qHist1D_sig_vslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

    _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzerDB",dir[6]);
  dir[6]->GetObject("qHist1D_sig_vslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzerDB",dir[7]);
  dir[7]->GetObject("qHist1D_sig_vslice_6",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzerDB",dir[8]);
  dir[8]->GetObject("qHist1D_sig_vslice_7",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzerDB",dir[9]);
  dir[9]->GetObject("qHist1D_sig_vslice_8",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   Double_t sum=0;
  
 for(Int_t k=1; k<=9; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=9; i++)
    {
     h_sum->Add(h[i],1);
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);


   cout<<h_sum->GetBinContent(500)<<" "<<endl;

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   for(int i=1; i<=9; i++)
     {
      h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;



   hcalo = new TH1D("calosum histogram for fr", "hcalo", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));

   for(Int_t k=1; k<=h_sum->GetNbinsX(); k++)
     {
       hcalo->SetBinContent(k,h_sum->GetBinContent(k));
       hcalo->SetBinError(k,h_sum->GetBinError(k));
     }
   hcalo->Rebin(4);
   


 
 


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
  fit_func21->SetParName(17, "A_vbo");
  fit_func21->SetParName(18, "Tau_vbo");
  fit_func21->SetParName(19, "omega_vbo");
  fit_func21->SetParName(20, "phi_vbo");
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
  fit_func22->SetParName(17, "A_vbo");
  fit_func22->SetParName(18, "Tau_vbo");
  fit_func22->SetParName(19, "omega_vbo");
  fit_func22->SetParName(20, "phi_vbo");
  fit_func22->SetParName(21, "muloss_amp");
  fit_func22->SetNpx(10000000);

  fit_func26= new TF1("fprec26", fprec26,  30000,  309000, 26);
  fit_func26->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func26->SetParName(9, "A2_cbo");
  fit_func26->SetParName(10, "phi2_cbo");
  fit_func26->SetParName(11, "A3_cbo");
  fit_func26->SetParName(12, "phi3_cbo");
  fit_func26->SetParName(13, "A_vw");
  fit_func26->SetParName(14, "Tau_vw");
  fit_func26->SetParName(15, "omega_vw");
  fit_func26->SetParName(16, "phi_vw");
  fit_func26->SetParName(17, "A_vbo");
  fit_func26->SetParName(18, "Tau_vbo");
  fit_func26->SetParName(19, "omega_vbo");
  fit_func26->SetParName(20, "phi_vbo");
  fit_func26->SetParName(21, "muloss_amp");
  fit_func26->SetParName(22, "A_pr");
  fit_func26->SetParName(23, "Tau_pr");
  fit_func26->SetParName(24, "omega_pr");
  fit_func26->SetParName(25, "phi_pr");
  fit_func26->SetNpx(10000000);








  c1 = new TCanvas("c1","hslice wiggle fit");



      fit_start=30000;
      fit_stop=300000;

   h_sum->Rebin(8);


 for(int i=1; i<=9; i++)
 {


   cout<<"V-Slice "<<i<<endl;
   


  h[i]->Rebin(8);
  fit_func->SetParameters(h[i]->GetBinContent(83), 64440, 0.221, 0, 2.22);


  h[i]->GetYaxis()->SetRangeUser(-100., 2500000000.);
  h[i]->GetYaxis()->SetTitle("ADC counts");
  h[i]->GetXaxis()->SetTitle("time [ns]");
  h[i]->Fit(fit_func,"R","",fit_start,fit_stop);

  fit_func9->SetParameter(0,fit_func->GetParameter(0));
  fit_func9->SetParameter(1,fit_func->GetParameter(1));
  fit_func9->SetParameter(2,fit_func->GetParameter(2));
  fit_func9->SetParameter(3,fit_func->GetParameter(3));
  fit_func9->SetParameter(4,fit_func->GetParameter(4));
  fit_func9->SetParameter(5,0.001);
  fit_func9->SetParameter(6,350000);
  fit_func9->SetParameter(7,0.00234);
  fit_func9->SetParameter(8,0.5);

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
  fit_func13->SetParameter(10, 1.67);
  fit_func13->SetParameter(11, -0.0039);
  fit_func13->SetParameter(12, 1.29);

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
      fit_func17->SetParameter(13,0.0006);
      fit_func17->SetParameter(14,275000);
      //   fit_func17->SetParLimits(14,50000,900000);
      fit_func17->SetParameter(15, 0.01391);
      fit_func17->SetParameter(16, 0.3);

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
      fit_func21->SetParameter(17, 0.003);
      fit_func21->SetParameter(18, 39000);
      //  fit_func21->SetParLimits(18,5000,90000);
      fit_func21->SetParameter(19, 0.01404);
      fit_func21->SetParameter(20, 2.5);

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
      //   fit_func22->SetParLimits(14,50000,900000);
      // fit_func22->SetParLimits(18,5000,90000);

      h[i]->Fit(fit_func22,"R","",fit_start,fit_stop);



      /*      fit_func26->SetParameter(0,fit_func22->GetParameter(0));
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
      fit_func26->SetParameter(22, 0.0003);
      fit_func26->SetParameter(23, 77000);
      fit_func26->SetParameter(24, 0.0107);
      fit_func26->SetParameter(25, 8.3);

      
      
      h[i]->Fit(fit_func26,"R","",fit_start,fit_stop);*/
  
      gStyle->SetOptFit(1111);

      blindR[m]=fit_func22->GetParameter(4);
      // blindR[m]=fit_func22->GetChisquare()/fit_func22->GetNDF();
      dblindR[m]=fit_func22->GetParError(4);
      n[m]=m;
      m=m+1;


      h_res[i] = new TH1D("residual hist","h_res",h[i]->GetNbinsX(),h[i]->GetBinLowEdge(1), (h[i]->GetBinLowEdge(h[i]->GetNbinsX())+h[i]->GetBinWidth(1)));
      
      for (int ibin = ((fit_start + 1- h[i]->GetBinLowEdge(1))/h[i]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h[i]->GetBinLowEdge(1))/h[i]->GetBinWidth(1))+1; ++ibin)
      {
       double res =  (h[i]->GetBinContent(ibin)- fit_func22->Eval( h[i]->GetBinCenter(ibin) ) );
       if(h[i]->GetBinError(ibin)!=0){res=(res/h[i]->GetBinError(ibin));}
       h_res[i]->SetBinContent(ibin, (res)  );
      }

      h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
      h_res[i]->GetXaxis()->SetTitle("time [ns]");
      h_res[i]->GetXaxis()->SetRangeUser(fit_start, fit_stop);
      h_res[i]->GetXaxis()->SetLabelSize(0.09);
      h_res[i]->GetYaxis()->SetLabelSize(0.09);


     hm[i] = h_res[i]->FFT(hm[i], "MAG");
     hm[i]->SetLineColor(kBlack);
     hm[i]->SetBins(h_res[i]->GetNbinsX(), 0, 1000/h_res[i]->GetBinWidth(1));
     hm[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
     hm[i]->GetXaxis()->SetTitle("Freq [MHz]");
     hm[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
     hm[i]->GetXaxis()->SetLabelSize(0.09);
     hm[i]->GetYaxis()->SetLabelSize(0.09);



   
 }

  c2 = new TCanvas("c2","per calo residuals");
  c2->Divide(1,9);
  c3 = new TCanvas("c3","per calo fft");
  c3->Divide(1,9);
   
 for(int i=1;i<=9;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hm[i]->Draw();
   }

 c4=new TCanvas("c1","phase vs vertical slice");
    gr1=new TGraphErrors(m,n,blindR,0,dblindR);
    gr1->SetTitle("g-2 phase vs vertical slice number");
    gr1->GetXaxis()->SetTitle("vslice #");
    gr1->GetYaxis()->SetTitle("phi");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->RemovePoint(0);
    // gr1->GetYaxis()->SetRangeUser(-55,-30);
    gr1->Draw();


   /*      fit_func29= new TF1("fprec29", fprec29,  30000,  309000, 29);
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

      




      c9 = new TCanvas("c9","29 parameter fit");
      c9->Divide(1,3);
      c9->cd(1);
      
      		h_sum->Fit(fit_func29,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res29= new TH1D("residual histogram 29", "h_res 29", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c9->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func29->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res29->SetBinContent(ibin, (res)  );
     
      }

	h_res29->GetYaxis()->SetTitle("ADC counts");
	h_res29->GetXaxis()->SetTitle("time [ns]");	
        h_res29->Draw();

  c9->cd(3);
     hm29 = h_res29->FFT(hm29, "MAG");
   hm29->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm29->SetBins(h_res29->GetNbinsX(), 0, 1/h_res29->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm29->GetXaxis()->SetRangeUser(0,1/(2*h_res29->GetBinWidth(1)));
   hm29->GetXaxis()->SetTitle("Freq [GHz]");
   hm29->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm29->Draw();     
   */
      
 }
