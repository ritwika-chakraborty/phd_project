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
char root_file_name[128] = "run2_final.root";

char muonloss_rootfile[128]="LostMuonsSpectra.root";

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo;
//TH1 *hm;
TH1F *hlm, h0;
TCanvas *c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10;
TF1 *fit_func, *fit_func9, *fit_func13, *fit_func17, *fit_func21, *fit_func25, *fit_func26, *fit_func29 ,*fit_func30, *fit_func31,*fit_func32,*fit_func33;
TH1D  *h_res[50], *h_res9, *h_res13, *h_res17, *h_res21, *h_res25, *h_res26, *h_res29,*h_res30, *h_res31, *h_res32, *h_res33, *h_res_final;
TH1  *hm[50], *hm9, *hm13, *hm17, *hm21, *hm25, *hm26,  *hm29, *hm30, *hm31, *hm32, *hm33, *hFFTsq, *hfft_final, *hAuto;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;
Double_t blindR[25], dblindR[25], A[25], dA[25];



Double_t fprec5(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym * cos(omega_a*time - phi));

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

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    return norm * Ncbox * exp(-time/life) * (1 + asym * cos(omega_a*time - phi));

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

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    return norm * Ncbox * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

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

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    return norm * Ncbox * Nvoy1 * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

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

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    return norm * Ncbox * Nvoy1 * Nvoy2 * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

}

Double_t fprec25(Double_t *x, Double_t *par)
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

    double A_2cbo = par[21];

    double Tau_2cbo = par[22];

    double omega_2cbo = par[23];

    double phi_2cbo = par[24];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    return norm * Ncbox * Nvoy1 * Nvoy2 * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

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

    double A_2cbo = par[21];

    double Tau_2cbo = par[22];

    double omega_2cbo = par[23];

    double phi_2cbo = par[24];

    double A_new = par[25];

    double Tau_new = par[26];

    double omega_new = par[27];

    double phi_new = par[28];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    double new_freq =  (1 + A_new * exp(-time/Tau_new) * (cos(omega_new*time - phi_new)));
    
    return norm * Ncbox * Nvoy1 * Nvoy2 * new_freq * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

}


Double_t fprec31(Double_t *x, Double_t *par)
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

    double A_2cbo = par[21];

    double Tau_2cbo = par[22];

    double omega_2cbo = par[23];

    double phi_2cbo = par[24];

    double A_new = par[25];

    double Tau_new = par[26];

    double omega_new = par[27];

    double phi_new = par[28];

    double A_VD = par[29];

    double Tau_VD = par[30];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    double new_freq =  (1 + A_new * exp(-time/Tau_new) * (cos(omega_new*time - phi_new)));

    double vertical_drift =  (1 + A_VD * exp (-time/Tau_VD));

    return norm * Ncbox * Nvoy1 * Nvoy2 * vertical_drift * new_freq * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

}


Double_t fprec32(Double_t *x, Double_t *par)
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

    double A_2cbo = par[21];

    double Tau_2cbo = par[22];

    double omega_2cbo = par[23];

    double phi_2cbo = par[24];

    double A_new = par[25];

    double Tau_new = par[26];

    double omega_new = par[27];

    double phi_new = par[28];

    double A_VD = par[29];

    double Tau_VD = par[30];

    double lost_muon_amp = par[31];

    double tau_rlx = par[32];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    double new_freq =  (1 + A_new * exp(-time/Tau_new) * (cos(omega_new*time - phi_new)));

    double vertical_drift =  (1 + A_VD * exp (-time/Tau_VD));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time )));

    return norm * Ncbox * Nvoy1 * Nvoy2 * vertical_drift * muloss * new_freq * exp(-time/life) * (1 + asym * exp(-time/tau_rlx) * Acbox * cos(omega_a*time - phi - phicbox));

}



void starttime_scan()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

   _file[0]=TFile::Open(muonloss_rootfile);
  _file[0]->GetObject("Run2All",dir[0]);
  dir[0]->GetObject("triple_losses_spectra_integral",hlm);
  hlm->SetName("hlm");

  hlm->Scale(1./hlm->GetBinContent(hlm->GetNbinsX()));


    _file[1]=TFile::Open(root_file_name);
    // ############# My root file
    for(int i=1;i<=24;i++)
    {
     h[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=24; i++)
    {
     h_sum->Add(h[i],1);
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_sum->Draw();
   return;
   cout<<h_sum->GetBinContent(500)<<" "<<endl;

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   h_sum->Draw();

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


   fit_start=30000;
   fit_stop=300000;

   h_sum->Rebin(8);
   fit_func->SetParameters(h_sum->GetBinContent(83), 64439.1, 0.2337, 0, 2.22);



   h_sum->GetYaxis()->SetRangeUser(-100., 159000000000.);
	h_sum->GetYaxis()->SetTitle("ADC counts");
	h_sum->GetXaxis()->SetTitle("time [ns]");
	//fit_func->SetParLimits(8, -3.14, 3.14);
	 // fit_func->SetParLimits(10, -3.14, 3.14);
	//	h_sum->Draw();
        //fit_func->Draw("same");
    
  

      fit_func32= new TF1("fprec32", fprec32,  30000,  309000, 33);
      fit_func32->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
      fit_func32->SetParName(9, "A2_cbo");
      fit_func32->SetParName(10, "phi2_cbo");
      fit_func32->SetParName(11, "A3_cbo");
      fit_func32->SetParName(12, "phi3_cbo");
      fit_func32->SetParName(13, "A_vw");
      fit_func32->SetParName(14, "Tau_vw");
      fit_func32->SetParName(15, "omega_vw");
      fit_func32->SetParName(16, "phi_vw");
      fit_func32->SetParName(17, "A_vbo");
      fit_func32->SetParName(18, "Tau_vbo");
      fit_func32->SetParName(19, "omega_vbo");
      fit_func32->SetParName(20, "phi_vbo");
      fit_func32->SetParName(21, "A_2cbo");
      fit_func32->SetParName(22, "Tau_2cbo");
      fit_func32->SetParName(23, "omega_2cbo");
      fit_func32->SetParName(24, "phi_2cbo");
      fit_func32->SetParName(25, "A_new");
      fit_func32->SetParName(26, "Tau_new");
      fit_func32->SetParName(27, "omega_new");
      fit_func32->SetParName(28, "phi_new");
      fit_func32->SetParName(29, "A_VD");
      fit_func32->SetParName(30, "Tau_VD");
      fit_func32->SetParName(31, "lm_amp");
      fit_func32->SetParName(32, "Tau_rlx");


      fit_func32->SetNpx(10000000);
      fit_func32->SetParameter(0,1.34134e+11);
      fit_func32->SetParameter(1,6.44488e+04);
      fit_func32->SetParameter(2,2.29226e-01);
      fit_func32->SetParameter(3,0.0);
      fit_func32->SetParameter(4,4.06836e+00);
      fit_func32->SetParameter(5,2.48247e-03);
      fit_func32->FixParameter(6,2.43986e+05);
      fit_func32->FixParameter(7,2.34050e-03);
      fit_func32->SetParameter(8,-4.83424e-01);
      fit_func32->SetParameter(9,7.25802e-04);
      fit_func32->SetParameter(10,-3.91582e-01);
      fit_func32->SetParameter(11,-8.32201e-05);
      fit_func32->SetParameter(12,1.93346e+00);
      fit_func32->SetParameter(13,3.07163e-03);
      fit_func32->FixParameter(14,2.75847e+04);
      fit_func32->FixParameter(15,1.40374e-02);
      fit_func32->SetParameter(16,3.30456e+00);
      fit_func32->SetParameter(17,2.75015e-04);
      fit_func32->FixParameter(18,1.20477e+05);
      fit_func32->FixParameter(19,1.39319e-02);
      fit_func32->SetParameter(20,2.53853e-01);
      fit_func32->SetParameter(21,9.32236e-05);
      fit_func32->FixParameter(22,fit_func32->GetParameter(6)/2);
      fit_func32->FixParameter(23,2*fit_func32->GetParameter(7));
      fit_func32->SetParameter(24,3.45349e+00);
      fit_func32->SetParameter(25,1.42371e-04);
      fit_func32->FixParameter(26,5.14311e+04);
      fit_func32->FixParameter(27,1.19733e-02);
      fit_func32->SetParameter(28,4.22322e+00);
      fit_func32->SetParameter(29,-1.61204e-03);
      fit_func32->SetParameter(30,4.40814e+04);
      fit_func32->SetParameter(31,0.0);
      fit_func32->SetParameter(32,290000000);

      

      fit_func32->SetNpx(10000000);

      //fit_func32->SetParameter(0,139000000000);

      //fit_func32->Draw("same");

      // h_sum->Fit(fit_func32,"R","",fit_start,fit_stop);
      //return;

      for(int i=1; i<=20; i++)
      {


      cout<<"start time is "<<fit_start<<" ns"<<endl;
   
      fit_func32->SetParameter(0,fit_func32->GetParameter(0));
      fit_func32->SetParameter(1,fit_func32->GetParameter(1));
      fit_func32->SetParameter(2,fit_func32->GetParameter(2));
      fit_func32->SetParameter(3,fit_func32->GetParameter(3));
      fit_func32->SetParameter(4,fit_func32->GetParameter(4));
      fit_func32->SetParameter(5,fit_func32->GetParameter(5));
      fit_func32->SetParameter(6,fit_func32->GetParameter(6));
      fit_func32->SetParameter(7,fit_func32->GetParameter(7));
      fit_func32->SetParameter(8,fit_func32->GetParameter(8));
      fit_func32->SetParameter(9,fit_func32->GetParameter(9));
      fit_func32->SetParameter(10,fit_func32->GetParameter(10));
      fit_func32->SetParameter(11,fit_func32->GetParameter(11));
      fit_func32->SetParameter(12,fit_func32->GetParameter(12));
      fit_func32->SetParameter(13,fit_func32->GetParameter(13));
      fit_func32->SetParameter(14,fit_func32->GetParameter(14));
      fit_func32->SetParameter(15,fit_func32->GetParameter(15));
      fit_func32->SetParameter(16,fit_func32->GetParameter(16));
      fit_func32->SetParameter(17,fit_func32->GetParameter(17));
      fit_func32->SetParameter(18,fit_func32->GetParameter(18));
      fit_func32->SetParameter(19,fit_func32->GetParameter(19));
      fit_func32->SetParameter(20,fit_func32->GetParameter(20));
      fit_func32->SetParameter(21,fit_func32->GetParameter(21));
      fit_func32->SetParameter(22,fit_func32->GetParameter(22));
      fit_func32->SetParameter(23,fit_func32->GetParameter(23));
      fit_func32->SetParameter(24,fit_func32->GetParameter(24));
      fit_func32->SetParameter(25,fit_func32->GetParameter(25));
      fit_func32->SetParameter(26,fit_func32->GetParameter(26));
      fit_func32->SetParameter(27,fit_func32->GetParameter(27));
      fit_func32->SetParameter(28,fit_func32->GetParameter(28));
      fit_func32->SetParameter(29,fit_func32->GetParameter(29));
      fit_func32->SetParameter(30,fit_func32->GetParameter(30));
      fit_func32->SetParameter(31,fit_func32->GetParameter(31));
      fit_func32->SetParameter(32,fit_func32->GetParameter(32));

      




 

       h_sum->Fit(fit_func32,"R","",fit_start,fit_stop);

       gStyle->SetOptFit(1111);

       blindR[m]=fit_func32->GetParameter(3);
       dblindR[m]=fit_func32->GetParError(3);
       A[m]=fit_func32->GetParameter(2);
       dA[m]=fit_func32->GetParError(2);
       n[m]=fit_start;
       m=m+1;
 

      h_res[i] = new TH1D("residual hist","h_res",h_sum->GetNbinsX(),h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));
      
      for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
      {
       double res =  (h_sum->GetBinContent(ibin)- fit_func32->Eval( h_sum->GetBinCenter(ibin) ) );
       if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
       h_res[i]->SetBinContent(ibin, (res)  );
      }


      h_res[i]->GetYaxis()->SetTitle("ADC counts");
      h_res[i]->GetXaxis()->SetTitle("time [ns]");


     hm[i] = h_res[i]->FFT(hm[i], "MAG");
     hm[i]->SetLineColor(kBlack);
     hm[i]->SetBins(h_res[i]->GetNbinsX(), 0, 1/h_res[i]->GetBinWidth(1));
     hm[i]->GetXaxis()->SetRangeUser(0,1/(2*h_res[i]->GetBinWidth(1)));
     hm[i]->GetXaxis()->SetTitle("Freq [GHz]");
     hm[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");

     fit_start=fit_start+1000;

   
 }
     
  c2 = new TCanvas("c2","starttime residuals");
  c2->Divide(4,5);
  c3 = new TCanvas("c3","starttime fft");
  c3->Divide(4,5);
   
 for(int i=1;i<=20;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hm[i]->Draw();
   }

    c1=new TCanvas("c1","blind R vs fit start time");
    gr1=new TGraphErrors(m,n,blindR,0,dblindR);
    gr1->SetTitle("blindR vs fit start time");
    gr1->GetXaxis()->SetTitle("Start time [ns]");
    gr1->GetYaxis()->SetTitle("Blind R [ppm]");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->RemovePoint(0);
    gr1->Draw();

    c4=new TCanvas("c4","Asymmetry vs fit start time");
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle("A vs fit start time");
    gr2->GetXaxis()->SetTitle("Start time [ns]");
    gr2->GetYaxis()->SetTitle("Asymmetry [ppm]");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->RemovePoint(0);
    gr2->Draw();


 }
