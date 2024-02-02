#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TString.h"
#include <sys/time.h>
#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );

//char root_file_name[128] = "r2_w8_calos.root";
//char root_file_name[128] = "r3no_w8_xtalreject.root";
/*char root_file_name_0[128] = "run3B_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_1[128] = "run3C_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_2[128] = "run3D_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_3[128] = "run3E_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_4[128] = "run3F_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_5[128] = "run3G_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_6[128] = "run3I_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_7[128] = "run3J_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_8[128] = "run3K_thresh300_fbfDQC_wndw_8_new.root";
char root_file_name_9[128] = "run3L_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_10[128] = "run3M_thresh300_nofbfDQC_wndw_8_new.root";
char root_file_name_11[128] = "run3N_thresh300_nofbfDQC_wndw_8_new_2.root";
char root_file_name_12[128] = "run3O_thresh300_nofbfDQC_wndw_8.root";
*/

char root_file_name_0[128] = "r3b_w8_xtalreject.root";
char root_file_name_1[128] = "r3c_w8_xtalreject.root";
char root_file_name_2[128] = "r3d_w8_xtalreject.root";
char root_file_name_3[128] = "r3e_w8_xtalreject.root";
char root_file_name_4[128] = "r3f_w8_xtalreject.root";
char root_file_name_5[128] = "r3g_w8_xtalreject.root";
char root_file_name_6[128] = "r3ij_w8_xtalreject.root";
char root_file_name_7[128] = "r3ij_w8_xtalreject.root";
char root_file_name_8[128] = "r3k_w8_xtalreject.root";
char root_file_name_9[128] = "r3l_w8_xtalreject.root";
char root_file_name_10[128] = "r3m_w8_xtalreject.root";
char root_file_name_11[128] = "r3n_w8_xtalreject.root";
char root_file_name_12[128] = "r3o_w8_xtalreject.root";

char muonloss_rootfile[128]="LostMuonsSpectra.root";

//char root_file_name[128] = "run2_thresh300_fbfDQC_wndw_8.root";
char *histname = new char[10];
char *histname2 = new char[10];

TFile *_file[25];
TDirectoryFile *dir[25];

TVirtualFitter *gFitter;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *hrm1[25], *hrm2[25], *hrm3[25], *hrm4[25], *h_calosum, *hr_shift, *hr_noshift;
//TH1 *hm;
TH1F *hlm, h0,*hloss[20];
TCanvas *c1,*c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cfit, *cauto;
TF1 *fit_func, *fit_func9, *fit_func13, *fit_func17, *fit_func21, *fit_func25, *fit_func26, *fit_func29;
TH1D *h_res, *h_res9, *h_res13, *h_res17, *h_res21, *h_res25, *h_res26, *h_res29, *h_res_final;
TH1 *hm, *hm9, *hm13, *hm17, *hm21, *hm25, *hm26, *hm29, *hFFTsq, *hfft_final, *hAuto;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;

Double_t chi[40],R[40],dR[40],A[40],dA[40],phi[40],dphi[40];
Double_t n[40];
Int_t m=1;
Double_t N;

Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start=29926.250;
Double_t fit_stop=299926.25;
Double_t corr_coeff;
Int_t countfcn=0;

Int_t idataset;

Double_t fit_par[20][50];
Double_t dfit_par[20][50];

int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
TMatrixD cov(mdim,mdim);
//,cov2(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim)),data2((mdim)*(mdim));

TArrayD ratio_cov_matrix(TH1D *hr_sum)
{
  
 double covar;

 for(Int_t k=hr_sum->FindBin(10000); k<=hr_sum->GetNbinsX(); k++)
    {
      if(hr_sum->GetBinError(k)!=0)
	{ i0=k;
	  break;}
    }
  cout<<"i0 is"<<i0<<" "<<endl;
  
   for (Int_t i = 0; i < (mdim)*(mdim); i++)
      {
	const Int_t ir = (i)/(mdim);
	const Int_t ic = (i)%(mdim);

        
	
        data[i] = 0.0;
 
        if(ir==ic)
	  {
	    //  data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir))/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
           data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	  }

      	  if(ir==ic-1)
	  {
	   data[i]=(corr_coeff)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	  }


	  if(ir==ic+1)
	  {
	   data[i]=(corr_coeff)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
      	  }
	

	
      }
  cout<<"flgs "<<flg<<" "<<endl;

  return data;

}




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


Double_t fprec26(Double_t *x, Double_t *par)
{
  //hlm=(TH1F*)hloss[idataset]->Clone();
    
    //hlm->Scale(1./hlm->GetBinContent(hlm->GetNbinsX()));
  
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

    double lost_muon_amp = par[25];

    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time )));
    
    return norm * Ncbox * Nvoy1 * Nvoy2 * muloss * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

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

    double lost_muon_amp = par[25];

    double A_vd = par[26];

    double tau_vd = par[27];

    double C_vd = par[28];
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi3_cbo));

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time - phi_cbo)));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time - phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time - phi_vbo2)));

    double N2cbox = (A_2cbo * exp(-time/Tau_2cbo) * (cos(omega_2cbo*time - phi_2cbo)));

    Ncbox=Ncbox+N2cbox;

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time )));

    double vertdrift =  (1 + (A_vd * exp(-time/tau_vd)) + C_vd);

     
    return norm * Ncbox * Nvoy1 * Nvoy2 * muloss * vertdrift * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time - phi - phicbox));

}

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
  TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(h_sum->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
   for(Int_t i=h_sum->FindBin(fit_start); i<=h_sum->FindBin(fit_stop); i++)
    {
      for(Int_t j=h_sum->FindBin(fit_start); j<=h_sum->FindBin(fit_stop); j++)
	{
	    if(j>=i-5 && j<=i+5)
	    {
	      ch=ch+((h_sum->GetBinContent(i))-(fuser->Eval(h_sum->GetBinCenter(i))))*cov[i][j]*((h_sum->GetBinContent(j))-(fuser->Eval(h_sum->GetBinCenter(j))));
	    }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  //cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}


void run3_regular_q_runscan()

 {


 c1=new TCanvas("c1","9 parameter wiggle_fit");  



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


  fit_func25= new TF1("fprec25", fprec25,  30000,  309000, 25);
  fit_func25->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func25->SetParName(9, "A2_cbo");
  fit_func25->SetParName(10, "phi2_cbo");
  fit_func25->SetParName(11, "A3_cbo");
  fit_func25->SetParName(12, "phi3_cbo");
  fit_func25->SetParName(13, "A_vw");
  fit_func25->SetParName(14, "Tau_vw");
  fit_func25->SetParName(15, "omega_vw");
  fit_func25->SetParName(16, "phi_vw");
  fit_func25->SetParName(17, "A_vbo");
  fit_func25->SetParName(18, "Tau_vbo");
  fit_func25->SetParName(19, "omega_vbo");
  fit_func25->SetParName(20, "phi_vbo");
  fit_func25->SetParName(21, "A_2cbo");
  fit_func25->SetParName(22, "tau_2cbo");
  fit_func25->SetParName(23, "omega_2cbo");
  fit_func25->SetParName(24, "phi_2cbo");
  fit_func25->SetNpx(10000000);

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
  fit_func26->SetParName(21, "A_2cbo");
  fit_func26->SetParName(22, "tau_2cbo");
  fit_func26->SetParName(23, "omega_2cbo");
  fit_func26->SetParName(24, "phi_2cbo");
  fit_func26->SetParName(25, "mu_loss");
  fit_func26->SetNpx(10000000);


  fit_func29= new TF1("fprec29", fprec29,  30000,  309000, 29);
  fit_func29->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  fit_func29->SetParName(9, "A2_cbo");
  fit_func29->SetParName(10, "phi2_cbo");
  fit_func29->SetParName(11, "A3_cbo");
  fit_func29->SetParName(12, "phi3_cbo");
  fit_func29->SetParName(13, "A_vw");
  fit_func29->SetParName(14, "Tau_vw");
  fit_func29->SetParName(15, "omega_vw");
  fit_func29->SetParName(16, "phi_vw");
  fit_func29->SetParName(17, "A_vbo");
  fit_func29->SetParName(18, "Tau_vbo");
  fit_func29->SetParName(19, "omega_vbo");
  fit_func29->SetParName(20, "phi_vbo");
  fit_func29->SetParName(21, "A_2cbo");
  fit_func29->SetParName(22, "tau_2cbo");
  fit_func29->SetParName(23, "omega_2cbo");
  fit_func29->SetParName(24, "phi_2cbo");
  fit_func29->SetParName(25, "mu_loss");
  fit_func29->SetParName(26, "A_vd");
  fit_func29->SetParName(27, "Tau_vd");
  fit_func29->SetParName(28, "constant_vd");
  fit_func29->SetNpx(10000000);
  

 



  /* _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");
  */
    _file[1]=TFile::Open(root_file_name_0);
    _file[2]=TFile::Open(root_file_name_1);
    _file[3]=TFile::Open(root_file_name_2);
    _file[4]=TFile::Open(root_file_name_3);
    _file[5]=TFile::Open(root_file_name_4);
    _file[6]=TFile::Open(root_file_name_5);
    _file[7]=TFile::Open(root_file_name_6);
    _file[8]=TFile::Open(root_file_name_7);
    _file[9]=TFile::Open(root_file_name_8);
    _file[10]=TFile::Open(root_file_name_9);
    _file[11]=TFile::Open(root_file_name_10);
    _file[12]=TFile::Open(root_file_name_11);
    _file[13]=TFile::Open(root_file_name_12);



   _file[14]=TFile::Open(muonloss_rootfile);
   
   _file[14]->GetObject("Run3B",dir[1]);
   dir[1]->GetObject("triple_losses_spectra_integral",hloss[1]);
   
    _file[14]->GetObject("Run3C",dir[2]);
   dir[2]->GetObject("triple_losses_spectra_integral",hloss[2]);

    _file[14]->GetObject("Run3D",dir[3]);
   dir[3]->GetObject("triple_losses_spectra_integral",hloss[3]);

    _file[14]->GetObject("Run3E",dir[4]);
   dir[4]->GetObject("triple_losses_spectra_integral",hloss[4]);

    _file[14]->GetObject("Run3F",dir[5]);
   dir[5]->GetObject("triple_losses_spectra_integral",hloss[5]);

    _file[14]->GetObject("Run3G",dir[6]);
   dir[6]->GetObject("triple_losses_spectra_integral",hloss[6]);

     _file[14]->GetObject("Run3I",dir[7]);
   dir[7]->GetObject("triple_losses_spectra_integral",hloss[7]);
   
    _file[14]->GetObject("Run3J",dir[8]);
   dir[8]->GetObject("triple_losses_spectra_integral",hloss[8]);

    _file[14]->GetObject("Run3K",dir[9]);
   dir[9]->GetObject("triple_losses_spectra_integral",hloss[9]);

    _file[14]->GetObject("Run3L",dir[10]);
   dir[10]->GetObject("triple_losses_spectra_integral",hloss[10]);

    _file[14]->GetObject("Run3M",dir[11]);
   dir[11]->GetObject("triple_losses_spectra_integral",hloss[11]);

   /*    _file[14]->GetObject("Run3N",dir[12]);
   dir[12]->GetObject("triple_losses_spectra_integral",hloss[12]);

     _file[14]->GetObject("Run3O",dir[13]);
   dir[13]->GetObject("triple_losses_spectra_integral",hloss[13]);
   */


   

   
  for(idataset=1;idataset<=13;idataset++)
  {
   
    //   hlm=(TH1F*)hloss[idataset]->Clone();
    
    // hlm->Scale(1./hlm->GetBinContent(hlm->GetNbinsX()));

    // hlm->Draw();
    //return;
    
     // ############# My root file
    
    /*   for(int i=1;i<=24;i++)
    {
     h[i]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }
    
     cout<<" test !!!!!!!!!!!!!!! dataset "<<idataset<<", binwidth of h[1] "<<h[1]->GetBinWidth(100)<<", binwidth of h[2] "<<h[2]->GetBinWidth(100)<<endl;
    */
    // ############### Tim's root file



    if(idataset==12 || idataset==13)
      {
       hrm1[1]= (TH1D*)(_file[idataset]->FindObjectAny("qHist1D_sig_rm_sum_0_0"));
       hrm2[1]= (TH1D*)(_file[idataset]->FindObjectAny("qHist1D_sig_rm_sum_1_0"));
       hrm3[1]= (TH1D*)(_file[idataset]->FindObjectAny("qHist1D_sig_rm_sum_2_0"));
       hrm4[1]= (TH1D*)(_file[idataset]->FindObjectAny("qHist1D_sig_rm_sum_3_0"));
      
      
       for(int i=1;i<=1;i++)
       {
        sprintf(histname2,"calo_hist_sum_%d",i);	
        h[i]=new TH1D(histname2, histname2, hrm1[1]->GetNbinsX(), 100001, 352001);
        h[i]->Sumw2(kTRUE);
        h[i]->Add(hrm1[i],1);
        h[i]->Add(hrm2[i],1);
        h[i]->Add(hrm3[i],1);
        h[i]->Add(hrm4[i],1);
       }
      }
    else
      {
         for(int i=1;i<=24;i++)
         {
          hrm1[i]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_0_0", i)));
          hrm2[i]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_1_0", i)));
          hrm3[i]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_2_0", i)));
          hrm4[i]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_3_0", i)));
         }

	  for(int i=1;i<=24;i++)
          {
           sprintf(histname2,"calo_hist_sum_%d",i);	
           h[i]=new TH1D(histname2, histname2, hrm1[1]->GetNbinsX(), 100001, 352001);
           h[i]->Sumw2(kTRUE);
           h[i]->Add(hrm1[i],1);
           h[i]->Add(hrm2[i],1);
           h[i]->Add(hrm3[i],1);
           h[i]->Add(hrm4[i],1);
           }
    
      }
    
  /*h[24]->Fit("gaus","","",104750,104850);
  TF1 *fit = h[24]->GetFunction("gaus");
  inject_time=fit->GetParameter(1);
  */

  Double_t sum=0;
  if(idataset==12 || idataset==13)
   {
    for(Int_t k=1; k<=1; k++)
    {
     sum=sum+h[k]->GetBinContent(500);
    }
   }
  else
   {
    for(Int_t k=1; k<=24; k++)
    {
     sum=sum+h[k]->GetBinContent(500);
    }
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
  // h_calosum=(TH1D*)h[10]->Clone();
  
  h_calosum= new TH1D("calo histogram sum", "h_calosum", h[1]->GetNbinsX(), 100001, 352001);
  h_calosum->Sumw2(kTRUE);

  // cout<<" test !!!!!!!!!!!!!!! dataset "<<idataset<<", binwidth of h[1] "<<h[1]->GetBinWidth(100)<<", binwidth of h[2] "<<h[2]->GetBinWidth(100)<<endl;
   if(idataset==12 || idataset==13)
   {
    for(Int_t i=1; i<=1; i++)
    {
     h_calosum->Add(h[i],1);
    }
   }
   else
     {
      for(Int_t i=1; i<=24; i++)
      {
       h_calosum->Add(h[i],1);
      }
     }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_calosum->Draw();

   cout<<h_calosum->GetBinContent(500)<<" "<<endl;
    double binwidth=h_calosum->GetBinWidth(1);
    int nbins=h_calosum->GetNbinsX();
    double binlowedge=h_calosum->GetBinLowEdge(1);
    double binhighedge=h_calosum->GetBinLowEdge(nbins)+binwidth;

    cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

    h_calosum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   
   cout<<h_calosum->GetBinWidth(1)<<" "<<h_calosum->GetNbinsX()<<" "<<h_calosum->GetBinLowEdge(1)<<" "<<(h_calosum->GetBinLowEdge(h_calosum->GetNbinsX())+h_calosum->GetBinWidth(1))<<" "<<endl;

   h_calosum->Draw();


   
   h_calosum->Rebin(4);
   h_calosum->Scale(0.25);



    hr_noshift=(TH1D*)h_calosum->Clone();

    hr_shift=(TH1D*)hr_noshift->Clone();
    hr_shift->Reset();
   

    for(int ibin=1;ibin<=hr_shift->GetNbinsX();ibin++)
	{
	  hr_shift->SetBinContent(ibin,hr_noshift->GetBinContent(ibin+1));
	  hr_shift->SetBinError(ibin,hr_noshift->GetBinError(ibin+1));
      	}

	h_sum=(TH1D*)h_calosum->Clone();
	h_sum->Reset();

	for(int ibin=1;ibin<=h_calosum->GetNbinsX();ibin++)
	  {
	    h_sum->SetBinContent(ibin,0.5*(hr_noshift->GetBinContent(ibin)+hr_shift->GetBinContent(ibin)));
	    h_sum->SetBinError(ibin,0.5*sqrt((hr_noshift->GetBinError(ibin)*hr_noshift->GetBinError(ibin))+(hr_shift->GetBinError(ibin)*hr_shift->GetBinError(ibin))));
	  }

	
	h_sum->Rebin(2);
	h_sum->Scale(0.5);

   

 c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);
 
 





   fit_func->SetParameters(h_sum->GetBinContent(83), 64439.1, 0.2337, 0, 4.01);


           h_sum->GetYaxis()->SetRangeUser(-100., 35000000000.);
	h_sum->GetYaxis()->SetTitle("Energy [MeV]");
	h_sum->GetXaxis()->SetTitle("time [ns]");
	//fit_func->SetParLimits(8, -3.14, 3.14);
	 // fit_func->SetParLimits(10, -3.14, 3.14);
	//	h_sum->Draw();
        //fit_func->Draw("same");
			h_sum->Fit(fit_func,"R","",fit_start,fit_stop);
     
		h_res= new TH1D("residual histogram", "h_res", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("Energy [MeV]");
	h_res->GetXaxis()->SetTitle("time [ns]");
	h_res->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res->GetXaxis()->SetLabelSize(0.04);
	h_res->GetYaxis()->SetLabelSize(0.04);
        h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm->SetBins(h_res->GetNbinsX(), 0, 1000/h_res->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm->GetXaxis()->SetRangeUser(0,1000/(2*h_res->GetBinWidth(1)));
   hm->GetXaxis()->SetTitle("Freq [MHz]");
   hm->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm->GetXaxis()->SetLabelSize(0.04);
   hm->GetYaxis()->SetLabelSize(0.04);
   hm->Draw();     
   
	
      fit_func9->SetParameter(0,fit_func->GetParameter(0));
      fit_func9->SetParameter(1,fit_func->GetParameter(1));
      fit_func9->SetParameter(2,fit_func->GetParameter(2));
      fit_func9->SetParameter(3,fit_func->GetParameter(3));
      fit_func9->SetParameter(4,fit_func->GetParameter(4));
      fit_func9->SetParameter(5,0.0025);
      fit_func9->SetParameter(6,228960);
      fit_func9->SetParameter(7,0.002329);
      fit_func9->SetParameter(8,0.465);

      c3 = new TCanvas("c3","9 parameter fit");
      c3->Divide(1,3);
      c3->cd(1);
      
      		h_sum->Fit(fit_func9,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res9= new TH1D("residual histogram 9", "h_res 9", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c3->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func9->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res9->SetBinContent(ibin, (res)  );
     
      }

	h_res9->GetYaxis()->SetTitle("ADC counts");
	h_res9->GetXaxis()->SetTitle("time [ns]");	
        h_res9->Draw();

  c3->cd(3);
     hm9 = h_res9->FFT(hm9, "MAG");
   hm9->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm9->SetBins(h_res9->GetNbinsX(), 0, 1/h_res9->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm9->GetXaxis()->SetRangeUser(0,1/(2*h_res9->GetBinWidth(1)));
   hm9->GetXaxis()->SetTitle("Freq [GHz]");
   hm9->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm9->Draw();     
   //  return;
   
      fit_func13->SetParameter(0,fit_func9->GetParameter(0));
      fit_func13->SetParameter(1,fit_func9->GetParameter(1));
      fit_func13->SetParameter(2,fit_func9->GetParameter(2));
      fit_func13->SetParameter(3,fit_func9->GetParameter(3));
      fit_func13->SetParameter(4,fit_func9->GetParameter(4));
      fit_func13->SetParameter(5,fit_func9->GetParameter(5));
      fit_func13->SetParameter(6,fit_func9->GetParameter(6));
      fit_func13->SetParameter(7,fit_func9->GetParameter(7));
      fit_func13->SetParameter(8,fit_func9->GetParameter(8));
      fit_func13->SetParameter(9, 0.00067);
      fit_func13->SetParameter(10, 0.41);
      fit_func13->SetParameter(11, 0.00014);
      fit_func13->SetParameter(12, 2.77);


      c4 = new TCanvas("c4","13 parameter fit");
      c4->Divide(1,3);
      c4->cd(1);
      
      		h_sum->Fit(fit_func13,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res13= new TH1D("residual histogram 13", "h_res 13", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c4->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func13->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res13->SetBinContent(ibin, (res)  );
     
      }

	h_res13->GetYaxis()->SetTitle("ADC counts");
	h_res13->GetXaxis()->SetTitle("time [ns]");	
        h_res13->Draw();

  c4->cd(3);
     hm13 = h_res13->FFT(hm13, "MAG");
   hm13->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm13->SetBins(h_res13->GetNbinsX(), 0, 1/h_res13->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm13->GetXaxis()->SetRangeUser(0,1/(2*h_res13->GetBinWidth(1)));
   hm13->GetXaxis()->SetTitle("Freq [GHz]");
   hm13->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm13->Draw();     
   

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
      fit_func17->SetParameter(13, 0.0075);
      fit_func17->SetParameter(14,38384);
      fit_func17->SetParameter(15, 0.01405);
      fit_func17->SetParameter(16, 2.77);



      c5 = new TCanvas("c5","17 parameter fit");
      c5->Divide(1,3);
      c5->cd(1);
      
      		h_sum->Fit(fit_func17,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res17= new TH1D("residual histogram 17", "h_res 17", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c5->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func17->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res17->SetBinContent(ibin, (res)  );
     
      }

	h_res17->GetYaxis()->SetTitle("ADC counts");
	h_res17->GetXaxis()->SetTitle("time [ns]");	
        h_res17->Draw();

  c5->cd(3);
     hm17 = h_res17->FFT(hm17, "MAG");
   hm17->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm17->SetBins(h_res13->GetNbinsX(), 0, 1/h_res17->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm17->GetXaxis()->SetRangeUser(0,1/(2*h_res17->GetBinWidth(1)));
   hm17->GetXaxis()->SetTitle("Freq [GHz]");
   hm17->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm17->Draw();     
   
   
   
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
      fit_func21->FixParameter(14,fit_func17->GetParameter(14));
      fit_func21->FixParameter(15,fit_func17->GetParameter(15));
      fit_func21->SetParameter(16,fit_func17->GetParameter(16));
      fit_func21->SetParameter(17, 0.000208);
      fit_func21->SetParameter(18,102568);
      fit_func21->SetParameter(19, 0.0139);
      fit_func21->SetParameter(20, 0.0167);




      c6 = new TCanvas("c6","21 parameter fit");
      c6->Divide(1,3);
      c6->cd(1);
      
      		h_sum->Fit(fit_func21,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res21= new TH1D("residual histogram 21", "h_res 21", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c6->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func21->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res21->SetBinContent(ibin, (res)  );
     
      }

	h_res21->GetYaxis()->SetTitle("ADC counts");
	h_res21->GetXaxis()->SetTitle("time [ns]");	
        h_res21->Draw();

  c6->cd(3);
     hm21 = h_res21->FFT(hm21, "MAG");
   hm21->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm21->SetBins(h_res21->GetNbinsX(), 0, 1000/h_res21->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm21->GetXaxis()->SetRangeUser(0,1000/(2*h_res21->GetBinWidth(1)));
   hm21->GetXaxis()->SetTitle("Freq [MHz]");
   hm21->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm21->Draw();
   
   //  return;
   
      fit_func25->SetParameter(0,fit_func21->GetParameter(0));
      fit_func25->SetParameter(1,fit_func21->GetParameter(1));
      fit_func25->SetParameter(2,fit_func21->GetParameter(2));
      fit_func25->SetParameter(3,fit_func21->GetParameter(3));
      fit_func25->SetParameter(4,fit_func21->GetParameter(4));
      fit_func25->SetParameter(5,fit_func21->GetParameter(5));
      fit_func25->SetParameter(6,fit_func21->GetParameter(6));
      fit_func25->SetParameter(7,fit_func21->GetParameter(7));
      fit_func25->SetParameter(8,fit_func21->GetParameter(8));
      fit_func25->SetParameter(9,fit_func21->GetParameter(9));
      fit_func25->SetParameter(10,fit_func21->GetParameter(10));
      fit_func25->SetParameter(11,fit_func21->GetParameter(11));
      fit_func25->SetParameter(12,fit_func21->GetParameter(12));
      fit_func25->SetParameter(13,fit_func21->GetParameter(13));
      fit_func25->FixParameter(14,fit_func21->GetParameter(14));
      fit_func25->FixParameter(15,fit_func21->GetParameter(15));
      fit_func25->SetParameter(16,fit_func21->GetParameter(16));
      fit_func25->SetParameter(17,fit_func21->GetParameter(17));
      fit_func25->FixParameter(18,fit_func21->GetParameter(18));
      fit_func25->FixParameter(19,fit_func21->GetParameter(19));
      fit_func25->SetParameter(20,fit_func21->GetParameter(20));
      fit_func25->SetParameter(21,0.001);
      fit_func25->FixParameter(22,fit_func21->GetParameter(6)/2);
      fit_func25->FixParameter(23,2*fit_func21->GetParameter(7));
      fit_func25->SetParameter(24,3.5);
      




      c7 = new TCanvas("c7","25 parameter fit");
      c7->Divide(1,3);
      c7->cd(1);
      
      		h_sum->Fit(fit_func25,"R","",fit_start,fit_stop);
                h_sum->Draw();
		h_res25= new TH1D("residual histogram 25", "h_res 25", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c7->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func25->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res25->SetBinContent(ibin, (res)  );
     
      }

        h_res25->GetYaxis()->SetTitle("Energy [MeV]");
	h_res25->GetXaxis()->SetTitle("time [ns]");
	h_res25->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res25->GetXaxis()->SetLabelSize(0.04);
	h_res25->GetYaxis()->SetLabelSize(0.04);
        h_res25->Draw();

  c7->cd(3);
     hm25 = h_res25->FFT(hm25, "MAG");
   hm25->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm25->SetBins(h_res25->GetNbinsX(), 0, 1000/h_res25->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm25->GetXaxis()->SetRangeUser(0,1000/(2*h_res25->GetBinWidth(1)));
   hm25->GetXaxis()->SetTitle("Freq [MHz]");
   hm25->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm25->GetXaxis()->SetLabelSize(0.04);
   hm25->GetYaxis()->SetLabelSize(0.04);
   hm25->Draw();     

   hm25->GetXaxis()->SetTitle("Freq [GHz]");
   hm25->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm25->Draw();

   



  


   /*    fit_func26->SetParameter(0,fit_func25->GetParameter(0));
      fit_func26->SetParameter(1,fit_func25->GetParameter(1));
      fit_func26->SetParameter(2,fit_func25->GetParameter(2));
      fit_func26->SetParameter(3,fit_func25->GetParameter(3));
      fit_func26->SetParameter(4,fit_func25->GetParameter(4));
      fit_func26->SetParameter(5,fit_func25->GetParameter(5));
      fit_func26->SetParameter(6,fit_func25->GetParameter(6));
      fit_func26->SetParameter(7,fit_func25->GetParameter(7));
      fit_func26->SetParameter(8,fit_func25->GetParameter(8));
      fit_func26->SetParameter(9,fit_func25->GetParameter(9));
      fit_func26->SetParameter(10,fit_func25->GetParameter(10));
      fit_func26->SetParameter(11,fit_func25->GetParameter(11));
      fit_func26->SetParameter(12,fit_func25->GetParameter(12));
      fit_func26->SetParameter(13,fit_func25->GetParameter(13));
      fit_func26->FixParameter(14,fit_func25->GetParameter(14));
      fit_func26->FixParameter(15,fit_func25->GetParameter(15));
      fit_func26->SetParameter(16,fit_func25->GetParameter(16));
      fit_func26->SetParameter(17,fit_func25->GetParameter(17));
      fit_func26->FixParameter(18,fit_func25->GetParameter(18));
      fit_func26->FixParameter(19,fit_func25->GetParameter(19));
      fit_func26->SetParameter(20,fit_func25->GetParameter(20));
      fit_func26->SetParameter(21,fit_func25->GetParameter(21));
      fit_func26->FixParameter(22,fit_func25->GetParameter(6)/2);
      fit_func26->FixParameter(23,2*fit_func25->GetParameter(7));
      fit_func26->SetParameter(24,fit_func25->GetParameter(24));
      fit_func26->SetParameter(25,0);

       c8 = new TCanvas("c8","26 parameter fit");
      c8->Divide(1,3);
      c8->cd(1);
      
     	h_sum->Fit(fit_func26,"R","",fit_start,fit_stop);
        h_sum->Draw();
	//fit_func26->SetLineColor(kGreen);
	//fit_func26->Draw("same");
	h_res26= new TH1D("residual histogram 26", "h_res 26", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c8->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func26->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res26->SetBinContent(ibin, (res)  );
     
      }

	h_res26->GetYaxis()->SetTitle("Energy [MeV]");
	h_res26->GetXaxis()->SetTitle("time [ns]");
	h_res26->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res26->GetXaxis()->SetLabelSize(0.04);
	h_res26->GetYaxis()->SetLabelSize(0.04);
        h_res26->Draw();

  c8->cd(3);
     hm26 = h_res26->FFT(hm26, "MAG");
   hm26->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm26->SetBins(h_res26->GetNbinsX(), 0, 1000/h_res26->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm26->GetXaxis()->SetRangeUser(0,1000/(2*h_res26->GetBinWidth(1)));
   hm26->GetXaxis()->SetTitle("Freq [MHz]");
   hm26->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm26->GetXaxis()->SetLabelSize(0.04);
   hm26->GetYaxis()->SetLabelSize(0.04);
   hm26->Draw();     
 
   

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
      fit_func29->FixParameter(14,fit_func26->GetParameter(14));
      fit_func29->FixParameter(15,fit_func26->GetParameter(15));
      fit_func29->SetParameter(16,fit_func26->GetParameter(16));
      fit_func29->SetParameter(17,fit_func26->GetParameter(17));
      fit_func29->FixParameter(18,fit_func26->GetParameter(18));
      fit_func29->FixParameter(19,fit_func26->GetParameter(19));
      fit_func29->SetParameter(20,fit_func26->GetParameter(20));
      fit_func29->SetParameter(21,fit_func26->GetParameter(21));
      fit_func29->FixParameter(22,fit_func26->GetParameter(6)/2);
      fit_func29->FixParameter(23,2*fit_func26->GetParameter(7));
      fit_func29->SetParameter(24,fit_func26->GetParameter(24));
      fit_func29->SetParameter(25,fit_func26->GetParameter(25));
      fit_func29->SetParameter(26,-0.00842);
      fit_func29->SetParameter(27, 12000);
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
	h_res29->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res29->GetXaxis()->SetLabelSize(0.04);
	h_res29->GetYaxis()->SetLabelSize(0.04);
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

    hFFTsq = (TH1D*)hm25->Clone();
       hFFTsq->Reset();
       for (int ib = 1; ib <= hFFTsq->GetNbinsX(); ib++)
        {
         hFFTsq->SetBinContent( ib, hm25->GetBinContent(ib)*hm25->GetBinContent(ib));
        }

       cauto=new TCanvas("cauto","catuo");
       
       hAuto = hFFTsq->FFT(hAuto, "MAG");
       hAuto->GetXaxis()->SetLabelSize(0.04);
       hAuto->GetXaxis()->SetTitle("bins");
       hAuto->GetXaxis()->SetRangeUser(0,h_res25->GetNbinsX()/2);
       hAuto->SetLineColor(kRed);
       hAuto->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hAuto->Draw("HIST");

       corr_coeff=hAuto->GetBinContent(2)/hAuto->GetBinContent(1);


       data2=ratio_cov_matrix(h_sum);

       TMatrixD cov2(mdim,mdim);

       cov2.SetMatrixArray(data2.GetArray());

       cov2.ResizeTo(h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), -1);

       cov2.SetTol(1.e-23);
       Double_t det1;
       cov2.Invert(&det1);

       cov.ResizeTo(h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), -1);

;
       for(int ix=h_sum->FindBin(fit_start);ix<=cov.GetRowUpb();ix++)
	 {
	   for(int jx=h_sum->FindBin(fit_start);jx<=cov.GetColUpb();jx++)
	     {
	       cov[ix][jx]=cov2[ix][jx];
	     }
	 }
   

       fit_func25->SetParameter(0, fit_func25->GetParameter(0));
       fit_func25->SetParameter(1, fit_func25->GetParameter(1));
       fit_func25->SetParameter(2, fit_func25->GetParameter(2));
       fit_func25->SetParameter(3, fit_func25->GetParameter(3));
       fit_func25->SetParameter(4, fit_func25->GetParameter(4));
       fit_func25->SetParameter(5, fit_func25->GetParameter(5));
       fit_func25->SetParameter(6, fit_func25->GetParameter(6));
       fit_func25->SetParameter(7, fit_func25->GetParameter(7));
       fit_func25->SetParameter(8, fit_func25->GetParameter(8));
       fit_func25->SetParameter(9, fit_func25->GetParameter(9));
       fit_func25->SetParameter(10, fit_func25->GetParameter(10));
       fit_func25->SetParameter(11, fit_func25->GetParameter(11));
       fit_func25->SetParameter(12, fit_func25->GetParameter(12));
       fit_func25->SetParameter(13, fit_func25->GetParameter(13));
       fit_func25->FixParameter(14, fit_func25->GetParameter(14));
       fit_func25->FixParameter(15, fit_func25->GetParameter(15));
       fit_func25->SetParameter(16, fit_func25->GetParameter(16));
       fit_func25->SetParameter(17, fit_func25->GetParameter(17));
       fit_func25->FixParameter(18, fit_func25->GetParameter(18));
       fit_func25->FixParameter(19, fit_func25->GetParameter(19));
       fit_func25->SetParameter(20, fit_func25->GetParameter(20));
       fit_func25->SetParameter(21, fit_func25->GetParameter(21));
       fit_func25->FixParameter(22, fit_func25->GetParameter(22));
       fit_func25->FixParameter(23, fit_func25->GetParameter(23));
       fit_func25->SetParameter(24, fit_func25->GetParameter(24));

       

  
       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(h_sum);
       gFitter->SetFCN(chi2);

       cfit=new TCanvas("cfit","cfit");
       cfit->Divide(1,3);
       cfit->cd(1);

       h_sum->Fit("fprec25","RUE","",fit_start,fit_stop);

       countfcn=0;
       gStyle->SetOptFit(1111);
 
       gettimeofday(&t_end, NULL);
       printf("QRatio fit duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));


         h_res_final= new TH1D("residual histogram final", "h_res_final", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


    
	 for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (h_sum->GetBinContent(ibin)- fit_func25->Eval( h_sum->GetBinCenter(ibin) ) );
          if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
          h_res_final->SetBinContent(ibin, (res)  );
     
        }

	h_res_final->GetYaxis()->SetTitle("Energy [MeV]");
	h_res_final->GetXaxis()->SetTitle("time [ns]");
	h_res_final->GetXaxis()->SetLabelSize(0.05);
        h_res_final->GetYaxis()->SetLabelSize(0.05);
        cfit->cd(2);
	h_res_final->Draw();
      

 
       hfft_final = h_res_final->FFT(hfft_final, "MAG");
       hfft_final->SetLineColor(kBlack);
       hfft_final->SetBins(h_sum->GetNbinsX(),0,1000/h_res_final->GetBinWidth(1));
       hfft_final->GetXaxis()->SetRangeUser(0,1000/(2*h_res_final->GetBinWidth(1)));
       hfft_final->GetXaxis()->SetTitle("Freq [MHz]");
       hfft_final->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft_final->GetXaxis()->SetLabelSize(0.05);
       hfft_final->GetYaxis()->SetLabelSize(0.05);
       cfit->cd(3);
       hfft_final->Draw();
       
   if(idataset==12 || idataset==13)
   {
    for(int i=1;i<=1;i++)
    {
      h[i]->Delete();
    }
   }
   else
   {
    for(int i=1;i<=24;i++)
    {
      h[i]->Delete();
    }
   }
    h_calosum->Reset();

    for(int ipar=0;ipar<=28;ipar++)
     {
      fit_par[idataset][ipar]=fit_func25->GetParameter(ipar);
      dfit_par[idataset][ipar]=fit_func25->GetParError(ipar);
     }

  }
   
           

 }
