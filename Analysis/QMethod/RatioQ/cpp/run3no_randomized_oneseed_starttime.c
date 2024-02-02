#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include <sys/time.h>
#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );


//char root_file_name[128] = "highstat_wigsim_output_statfluc_binintegral_e16.root";
char root_file_name[128] = "r3no_w8_xtalreject.root";
//char root_file_name[128] = "run2C_thresh300_calosum_ratio.root";
//char root_file_name[128] = "run2_thresh300_fbfDQC.root";
TVirtualFitter *gFitter;
TH1 *hfft[50], *hFFTsq[50], *hAuto[50];
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*hr1[25],*hr2[25],*hr3[25],*hr4[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[50], *hr_sum, *hsum1,*hsum2,*hsum3,*hsum4,*hosum1, *hosum2, *hosum3, *hosum4, *hr_shift1, *hr_shift2, *hr_shift3, *hr_shift4,*hr_noshift1, *hr_noshift2, *hr_noshift3, *hr_noshift4, *hcomp_sum1, *hcomp_sum2, *hcomp_sum3, *hcomp_sum4;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cauto,*cfit;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28;
TGraph *gkband1p,*gkband1m,*gkband2p,*gkband2m,*gkband3p,*gkband3m,*gkband4p,*gkband4m,*gkband5p,*gkband5m,*gkband6p,*gkband6m,*gkband7p,*gkband7m,*gkband8p,*gkband8m,*gkband9p,*gkband9m,*gkband10p,*gkband10m,*gkband11p,*gkband11m,*gkband12p,*gkband12m,*gkband13p,*gkband13m,*gkband14p,*gkband14m,*gkband15p,*gkband15m,*gkband16p,*gkband16m,*gkband17p,*gkband17m,*gkband18p,*gkband18m,*gkband19p,*gkband19m,*gkband20p,*gkband20m,*gkband21p,*gkband21m,*gkband22p,*gkband22m,*gkband23p,*gkband23m;
//TH1 *hfft, *hFFTsq, *hAuto;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=30000;
Double_t fit_stop=300000;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
bool useFR=true;
Int_t m=0;
Int_t corr_coeff;
Int_t countfcn=0;


Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi_0[50], dphi_0[50], A_cbo_N[50], tau_cbo[50], omega_cbo[50], phi_cbo_N[50], A_cbo_A[50], phi_cbo_A[50], A_cbo_phi[50], phi_cbo_phi[50], A_vw[50], tau_vw[50], omega_vw[50], phi_vw[50], A_y[50], tau_y[50], omega_y[50], phi_y[50], A_2cbo[50], tau_2cbo[50], omega_2cbo[50], phi_2cbo[50],  dA_cbo_N[50], dtau_cbo[50], domega_cbo[50], dphi_cbo_N[50], dA_cbo_A[50], dphi_cbo_A[50], dA_cbo_phi[50], dphi_cbo_phi[50], dA_vw[50], dtau_vw[50], domega_vw[50], dphi_vw[50], dA_y[50], dtau_y[50], domega_y[50], dphi_y[50], dA_2cbo[50], dtau_2cbo[50], domega_2cbo[50], dphi_2cbo[50], n[50],kband1p[50],kband1m[50],kband2p[50],kband2m[50],kband3p[50],kband3m[50],kband4p[50],kband4m[50],kband5p[50],kband5m[50],kband6p[50],kband6m[50],kband7p[50],kband7m[50],kband8p[50],kband8m[50],kband9p[50],kband9m[50],kband10p[50],kband10m[50],kband11p[50],kband11m[50],kband12p[50],kband12m[50],kband13p[50],kband13m[50],kband14p[50],kband14m[50],kband15p[50],kband15m[50],kband16p[50],kband16m[50],kband17p[50],kband17m[50],kband18p[50],kband18m[50],kband19p[50],kband19m[50],kband20p[50],kband20m[50],kband21p[50],kband21m[50],kband22p[50],kband22m[50],kband23p[50],kband23m[50];

//Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25];
TFile *foutput;
TH2D *hcov;
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
TMatrixD cov(mdim,mdim);
//,cov2(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim)),data2((mdim)*(mdim));



TArrayD ratio_cov_matrix(TH1D *hr_sum){

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



Double_t fprec3(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double  f=(1+ asym*cos(omega_a*time - phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) - phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) - phi));

    //  return 0.22*cos(omega_a*time + phi);
    return (2*f - ff - fb)/(2*f + ff + fb);

}


Double_t fprec7(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo));

    double  f=(1+ asym*cos(omega_a*time - phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) - phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) - phi));

    return (2*f*Ncbo - ff*Ncbof - fb*Ncbob)/(2*f*Ncbo + ff*Ncbof + fb*Ncbob);

}


Double_t fprec11(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_cbo_A = par[7];

    double phi_cbo_A= par[8];

    double A_cbo_phi= par[9];

    double phi_cbo_phi=par[10];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_phi));

    double  f=(1+ asym*Acbo*cos(omega_a*time - phi - phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) - phi - phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) - phi - phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo - ff*Ncbof - fb*Ncbob)/(2*f*Ncbo + ff*Ncbof + fb*Ncbob);

}


Double_t fprec15(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_cbo_A = par[7];

    double phi_cbo_A= par[8];

    double A_cbo_phi= par[9];

    double phi_cbo_phi=par[10];

    double asym_vw= par[11];

    double tau_vw = par[12];

    double omega_vw = par[13];

    double phi_vw = par[14];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time - phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) - phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) - phi_vw));

    double  f=(1+ asym*Acbo*cos(omega_a*time - phi - phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) - phi - phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) - phi - phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw - ff*Ncbof*Nvwf - fb*Ncbob*Nvwb)/(2*f*Ncbo*Nvw + ff*Ncbof*Nvwf + fb*Ncbob*Nvwb);

}


Double_t fprec19(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_cbo_A = par[7];

    double phi_cbo_A= par[8];

    double A_cbo_phi= par[9];

    double phi_cbo_phi=par[10];

    double asym_vw= par[11];

    double tau_vw = par[12];

    double omega_vw = par[13];

    double phi_vw = par[14];

    double asym_vbo= par[15];

    double tau_vbo = par[16];

    double omega_vbo = par[17];

    double phi_vbo = par[18];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time - phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) - phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) - phi_vw));

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time - phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) - phi_vbo));

    double Nvbob=(1+ asym_vbo*exp(-(time - T_a/2)/tau_vbo)*cos(omega_vbo*(time - T_a/2) - phi_vbo));

    double  f=(1+ asym*Acbo*cos(omega_a*time - phi - phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) - phi - phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) - phi - phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof - fb*Ncbob*Nvwb*Nvbob)/(2*f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof + fb*Ncbob*Nvwb*Nvbob);

}

Double_t fprec23(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_cbo_A = par[7];

    double phi_cbo_A= par[8];

    double A_cbo_phi= par[9];

    double phi_cbo_phi=par[10];

    double asym_vw= par[11];

    double tau_vw = par[12];

    double omega_vw = par[13];

    double phi_vw = par[14];

    double asym_vbo= par[15];

    double tau_vbo = par[16];

    double omega_vbo = par[17];

    double phi_vbo = par[18];

    double asym_2cbo= par[19];

    double tau_2cbo = par[20];

    double omega_2cbo = par[21];

    double phi_2cbo = par[22];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1 + asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time - phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) - phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) - phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time - phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) - phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) - phi_vw));

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time - phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) - phi_vbo));

    double Nvbob=(1+ asym_vbo*exp(-(time - T_a/2)/tau_vbo)*cos(omega_vbo*(time - T_a/2) - phi_vbo));

    double N2cbo=(asym_2cbo*exp(-time/tau_2cbo)*cos(omega_2cbo*time - phi_2cbo));

    double N2cbof=(asym_2cbo*exp(-(time + T_a/2)/tau_2cbo)*cos(omega_2cbo*(time + T_a/2) - phi_2cbo));

    double N2cbob=(asym_2cbo*exp(-(time - T_a/2)/tau_2cbo)*cos(omega_2cbo*(time - T_a/2) - phi_2cbo));

    Ncbo=Ncbo+N2cbo;
    Ncbof=Ncbof+N2cbof;
    Ncbob=Ncbob+N2cbob;

    double  f=(1+ asym*Acbo*cos(omega_a*time - phi - phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) - phi - phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) - phi - phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof - fb*Ncbob*Nvwb*Nvbob)/(2*f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof + fb*Ncbob*Nvwb*Nvbob);

}

TH1D* construct_rhist_rand(TH1D *hsum1, TH1D *hsum2, TH1D *hsum3, TH1D *hsum4,TH1D *hosum1, TH1D *hosum2, TH1D *hosum3, TH1D *hosum4)
{
 


   //create the component histograms of the ratio histograms   
   h1=(TH1D*)hsum1->Clone();
   h2=(TH1D*)hsum2->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

      nbinshift=(0.5*T_a_true)/h1->GetBinWidth(1);
      cout<<"!!!!!!! NBINSHIFT is "<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=hsum3->GetNbinsX(); ibin++)             
     {
       hp->SetBinContent( ibin, hsum3->GetBinContent( ibin + nbinshift));
       hp->SetBinError(ibin, hsum3->GetBinError( ibin + nbinshift));
     }

   for(int ibin=nbinshift; ibin<=hsum4->GetNbinsX(); ibin++)
     {
       hm->SetBinContent( ibin, hsum4->GetBinContent( ibin - nbinshift));
       hm->SetBinError( ibin, hsum4->GetBinError(ibin - nbinshift));
     }

   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*hsum1->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   double blife=exp(-(T_a)/(2*lifetime));
   double deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo_histogram_sum_ratio", "h_ratio", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));
       
      for(int ibin=hsum1->FindBin(0);ibin<=hsum1->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo histogram sum numerator", "h_num", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));

      h_num->Add(h1,1);
      h_num->Add(h2,1);
      h_num->Add(hp,-1);
      h_num->Add(hm,-1);

      h_deno->Add(h1,1);
      h_deno->Add(h2,1);
      h_deno->Add(hp,1);
      h_deno->Add(hm,1);
     
      h_num->Divide(h_deno);

     }
    //Assign error bars

   
 if(useerrorbars)
     {
       
     if(useFR)
       {
	 long double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,tdeno,ea,eap,eam,eb,ebp,ebm,ec,ecp,ecm,ed,edp,edm;
	 for(int ibin=hsum1->FindBin(0); ibin<=hsum1->GetNbinsX()-nbinshift; ibin++)
	   {
	     ea=hosum1->GetBinContent(2*ibin);
	     eam=hosum1->GetBinContent(2*ibin-1);
	     eap=hosum1->GetBinContent(2*ibin+1);
	     eb=hosum2->GetBinContent(2*ibin);
	     ebm=hosum2->GetBinContent(2*ibin-1);
	     ebp=hosum2->GetBinContent(2*ibin+1);
	     ec=hosum3->GetBinContent(2*ibin);
	     ecm=hosum3->GetBinContent(2*ibin-1);
	     ecp=hosum3->GetBinContent(2*ibin+1);
	     ed=hosum4->GetBinContent(2*ibin);
	     edm=hosum4->GetBinContent(2*ibin-1);
	     edp=hosum4->GetBinContent(2*ibin+1);

	     tdeno=pow((2*ea + eam + eap + 2*eb + ebm + ebp + 2*ec + ecm + ecp + 2*ed + edm + edp),2);

	     t1=((4*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum1->GetBinError(2*ibin);
	     t2=((2*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum1->GetBinError(2*ibin-1);
	     t3=((2*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum1->GetBinError(2*ibin+1);
	     t4=((4*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum2->GetBinError(2*ibin);
	     t5=((2*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum2->GetBinError(2*ibin-1);
	     t6=((2*(2*ec + ecm + ecp + 2*ed + edm + edp) )/tdeno)*hosum2->GetBinError(2*ibin+1);
	     t7=((-4*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum3->GetBinError(2*ibin);
	     t8=((-2*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum3->GetBinError(2*ibin-1);
	     t9=((-2*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum3->GetBinError(2*ibin+1);
	     t10=((-4*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum4->GetBinError(2*ibin);
	     t11=((-2*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum4->GetBinError(2*ibin-1);
	     t12=((-2*(2*ea + eam + eap + 2*eb + ebm + ebp) )/tdeno)*hosum4->GetBinError(2*ibin+1);

	     if(usesetbins)
             {
	      h_ratio->SetBinError( ibin, sqrt( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9) + (t10*t10) + (t11*t11) + (t12*t12) ) );
	      //cout<<"ibin "<<ibin<<"errorbar "<<sqrt( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9) + (t10*t10) + (t11*t11) + (t12*t12) )<<endl;
             }
             else
             {
	      h_num->SetBinError(ibin,  sqrt( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9) + (t10*t10) + (t11*t11) + (t12*t12) ) );
             }

	   }
       }
     else
      {
       long double t1,t2,t3, t4;
      for(int ibin=hsum1->FindBin(0); ibin<=hsum1->GetNbinsX()-nbinshift; ibin++)
      {
	t1=2*(h1->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

	t2=2*(h2->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t3=-2*(hp->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t4=-2*(hm->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( (t1*t2) + (t2*t2) + (t3*t3) + (t4*t4) ) );
       }
      }
     }
    }

 //send a common histogram name for fitting
 if(usesetbins)
      {
       hcalo=(TH1D*)h_ratio->Clone();
      }
   else
      {
       hcalo=(TH1D*)h_num->Clone();
      }



 /*  foutput=new TFile("ratio_4hist_randomized.root","new");
 h_ratio->Write();
 foutput->Close();
 */

 return hcalo;
  }

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(hcalo->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
   
    for(Int_t i=hcalo->FindBin(fit_start); i<=hcalo->FindBin(fit_stop); i++)
    {
      for(Int_t j=hcalo->FindBin(fit_start); j<=hcalo->FindBin(fit_stop); j++)
	{
	    if(j>=i-5 && j<=i+5)
	    {
	      ch=ch+((hcalo->GetBinContent(i))-(fuser->Eval(hcalo->GetBinCenter(i))))*cov[i][j]*((hcalo->GetBinContent(j))-(fuser->Eval(hcalo->GetBinCenter(j))));
	    }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}



void run3no_randomized_oneseed_starttime()
{

  
  _file[1]=TFile::Open(root_file_name);
  
  hsum1= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_0_0")));
  hsum2= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_1_0")));
  hsum3= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_2_0")));
  hsum4= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_3_0")));

  
  
  /*  hsum1->Sumw2(kTRUE);
  for(int i=1;i<=24;i++)
    {
      hsum1->Add(hr1[i],1);
    }
  hsum2->Sumw2(kTRUE);
   for(int i=1;i<=24;i++)
    {
      hsum2->Add(hr2[i],1);
    }
  hsum3->Sumw2(kTRUE);
   for(int i=1;i<=24;i++)
    {
      hsum3->Add(hr3[i],1);
    }
   hsum4->Sumw2(kTRUE);
 for(int i=1;i<=24;i++)
    {
      hsum4->Add(hr4[i],1);
    }
  */
 hsum2->SetLineColor(kRed);
 hsum3->SetLineColor(kBlack);
 hsum4->SetLineColor(kGreen);

 int start_bin=400;
 hsum2->Scale(hsum1->Integral(start_bin,16800)/hsum2->Integral(start_bin,16800));
 hsum3->Scale(hsum1->Integral(start_bin,16800)/hsum3->Integral(start_bin,16800));
 hsum4->Scale(hsum1->Integral(start_bin,16800)/hsum4->Integral(start_bin,16800));
 
  cout<<"Normalized from bin "<<start_bin<<endl;
 
   double binwidth=hsum1->GetBinWidth(1);
   int nbins=hsum1->GetNbinsX();
   double binlowedge=hsum1->GetBinLowEdge(1);
   double binhighedge=hsum1->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   
     hsum1->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum2->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum3->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum4->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     
   cout<<hsum1->GetBinWidth(1)<<" "<<hsum1->GetNbinsX()<<" "<<hsum1->GetBinLowEdge(1)<<" "<<(hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1))<<" "<<endl;

   hsum1->Rebin(4);
   hsum1->Scale(0.25);

   hsum2->Rebin(4);
   hsum2->Scale(0.25);

   hsum3->Rebin(4);
   hsum3->Scale(0.25);

   hsum4->Rebin(4);
   hsum4->Scale(0.25);


    hr_noshift1=(TH1D*)hsum1->Clone();
    hr_noshift2=(TH1D*)hsum2->Clone();
    hr_noshift3=(TH1D*)hsum3->Clone();
    hr_noshift4=(TH1D*)hsum4->Clone();
	
    hr_shift1=(TH1D*)hr_noshift1->Clone();
    hr_shift1->Reset();
    hr_shift2=(TH1D*)hr_noshift2->Clone();
    hr_shift2->Reset();
    hr_shift3=(TH1D*)hr_noshift3->Clone();
    hr_shift3->Reset();
    hr_shift4=(TH1D*)hr_noshift4->Clone();
    hr_shift4->Reset();


    for(int ibin=1;ibin<=hr_shift1->GetNbinsX();ibin++)
	{
	  hr_shift1->SetBinContent(ibin,hr_noshift1->GetBinContent(ibin+1));
	  hr_shift1->SetBinError(ibin,hr_noshift1->GetBinError(ibin+1));
	  hr_shift2->SetBinContent(ibin,hr_noshift2->GetBinContent(ibin+1));
	  hr_shift2->SetBinError(ibin,hr_noshift2->GetBinError(ibin+1));
          hr_shift3->SetBinContent(ibin,hr_noshift3->GetBinContent(ibin+1));
	  hr_shift3->SetBinError(ibin,hr_noshift3->GetBinError(ibin+1));
          hr_shift4->SetBinContent(ibin,hr_noshift4->GetBinContent(ibin+1));
	  hr_shift4->SetBinError(ibin,hr_noshift4->GetBinError(ibin+1));
	}

	hcomp_sum1=(TH1D*)hsum1->Clone();
	hcomp_sum1->Reset();
	hcomp_sum2=(TH1D*)hsum2->Clone();
	hcomp_sum2->Reset();
	hcomp_sum3=(TH1D*)hsum3->Clone();
	hcomp_sum3->Reset();
	hcomp_sum4=(TH1D*)hsum4->Clone();
	hcomp_sum4->Reset();


	for(int ibin=1;ibin<=hsum2->GetNbinsX();ibin++)
	  {
	    hcomp_sum1->SetBinContent(ibin,0.5*(hr_noshift1->GetBinContent(ibin)+hr_shift1->GetBinContent(ibin)));
	    hcomp_sum1->SetBinError(ibin,0.5*sqrt((hr_noshift1->GetBinError(ibin)*hr_noshift1->GetBinError(ibin))+(hr_shift1->GetBinError(ibin)*hr_shift1->GetBinError(ibin))));
	    hcomp_sum2->SetBinContent(ibin,0.5*(hr_noshift2->GetBinContent(ibin)+hr_shift2->GetBinContent(ibin)));
	    hcomp_sum2->SetBinError(ibin,0.5*sqrt((hr_noshift2->GetBinError(ibin)*hr_noshift2->GetBinError(ibin))+(hr_shift2->GetBinError(ibin)*hr_shift2->GetBinError(ibin))));
	    hcomp_sum3->SetBinContent(ibin,0.5*(hr_noshift3->GetBinContent(ibin)+hr_shift3->GetBinContent(ibin)));
	    hcomp_sum3->SetBinError(ibin,0.5*sqrt((hr_noshift3->GetBinError(ibin)*hr_noshift3->GetBinError(ibin))+(hr_shift3->GetBinError(ibin)*hr_shift3->GetBinError(ibin))));
	    hcomp_sum4->SetBinContent(ibin,0.5*(hr_noshift4->GetBinContent(ibin)+hr_shift4->GetBinContent(ibin)));
	    hcomp_sum4->SetBinError(ibin,0.5*sqrt((hr_noshift4->GetBinError(ibin)*hr_noshift4->GetBinError(ibin))+(hr_shift4->GetBinError(ibin)*hr_shift4->GetBinError(ibin))));
	  }

	
	hcomp_sum1->Rebin(2);
	hcomp_sum1->Scale(0.5);
	hcomp_sum2->Rebin(2);
	hcomp_sum2->Scale(0.5);
	hcomp_sum3->Rebin(2);
	hcomp_sum3->Scale(0.5);
	hcomp_sum4->Rebin(2);
	hcomp_sum4->Scale(0.5);

   
   
 
     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A", "Blind R", "Phi");
     fit_func3->SetNpx(1000000);

     fit_func7= new TF1("fprec7", fprec7,  30000,  309000, 7);
     fit_func7->SetParNames("A", "Blind R", "Phi");
     fit_func7->SetParName(3,"A_cbo");
     fit_func7->SetParName(4,"Tau_cbo");
     fit_func7->SetParName(5,"omega_cbo");
     fit_func7->SetParName(6,"phi_cbo");
     fit_func7->SetNpx(1000000);

     fit_func11= new TF1("fprec11", fprec11,  30000,  309000, 11);
     fit_func11->SetParNames("A", "Blind R", "Phi");
     fit_func11->SetParName(3,"A_cbo");
     fit_func11->SetParName(4,"Tau_cbo");
     fit_func11->SetParName(5,"omega_cbo");
     fit_func11->SetParName(6,"phi_cbo");
     fit_func11->SetParName(7,"A2_cbo");
     fit_func11->SetParName(8,"phi2_cbo");
     fit_func11->SetParName(9,"A3_cbo");
     fit_func11->SetParName(10,"phi3_cbo");
     fit_func11->SetNpx(1000000);

     fit_func15= new TF1("fprec15", fprec15,  30000,  309000, 15);
     fit_func15->SetParNames("A", "Blind R", "Phi");
     fit_func15->SetParName(3,"A_cbo");
     fit_func15->SetParName(4,"Tau_cbo");
     fit_func15->SetParName(5,"omega_cbo");
     fit_func15->SetParName(6,"phi_cbo");
     fit_func15->SetParName(7,"A2_cbo");
     fit_func15->SetParName(8,"phi2_cbo");
     fit_func15->SetParName(9,"A3_cbo");
     fit_func15->SetParName(10,"phi3_cbo");
     fit_func15->SetParName(11,"A_vw");
     fit_func15->SetParName(12,"Tau_vw");
     fit_func15->SetParName(13,"omega_vw");
     fit_func15->SetParName(14,"phi_vw");
     fit_func15->SetNpx(1000000);
    
     fit_func19= new TF1("fprec19", fprec19,  30000,  309000, 19);
     fit_func19->SetParNames("A", "Blind R", "Phi");
     fit_func19->SetParName(3,"A_cbo");
     fit_func19->SetParName(4,"Tau_cbo");
     fit_func19->SetParName(5,"omega_cbo");
     fit_func19->SetParName(6,"phi_cbo");
     fit_func19->SetParName(7,"A2_cbo");
     fit_func19->SetParName(8,"phi2_cbo");
     fit_func19->SetParName(9,"A3_cbo");
     fit_func19->SetParName(10,"phi3_cbo");
     fit_func19->SetParName(11,"A_vw");
     fit_func19->SetParName(12,"Tau_vw");
     fit_func19->SetParName(13,"omega_vw");
     fit_func19->SetParName(14,"phi_vw");
     fit_func19->SetParName(15,"A_y");
     fit_func19->SetParName(16,"Tau_y");
     fit_func19->SetParName(17,"omega_y");
     fit_func19->SetParName(18,"phi_y");
     fit_func19->SetNpx(1000000);


  
       fit_func23= new TF1("fprec23", fprec23,  30000,  309000, 23);
     fit_func23->SetParNames("A", "Blind R", "Phi");
     fit_func23->SetParName(3,"A_cbo_N");
     fit_func23->SetParName(4,"Tau_cbo");
     fit_func23->SetParName(5,"omega_cbo");
     fit_func23->SetParName(6,"phi_cbo_N");
     fit_func23->SetParName(7,"A_cbo_A");
     fit_func23->SetParName(8,"phi_cbo_A");
     fit_func23->SetParName(9,"A_cbo_phi");
     fit_func23->SetParName(10,"phi_cbo_phi");
     fit_func23->SetParName(11,"A_vw");
     fit_func23->SetParName(12,"Tau_vw");
     fit_func23->SetParName(13,"omega_vw");
     fit_func23->SetParName(14,"phi_vw");
     fit_func23->SetParName(15,"A_y");
     fit_func23->SetParName(16,"Tau_y");
     fit_func23->SetParName(17,"omega_y");
     fit_func23->SetParName(18,"phi_y");
     fit_func23->SetParName(19,"A_2cbo_N");
     fit_func23->SetParName(20,"Tau_2cbo");
     fit_func23->SetParName(21,"omega_2cbo");
     fit_func23->SetParName(22,"phi_2cbo_N");
     fit_func23->SetNpx(1000000);

	 
	fit_func3->SetParameters(0.23, 0.0, 4.01);


	hcalo=construct_rhist_rand(hcomp_sum1,hcomp_sum2,hcomp_sum3,hcomp_sum4,hsum1,hsum2,hsum3,hsum4);
  
	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");

	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");


        hcalo->GetXaxis()->SetRangeUser(fit_start,fit_stop);

        hcalo->Fit("fprec3","RE","",fit_start,fit_stop);
	//	hcalo->Draw();
	//	fit_func3->Draw();
	//	hcalo->Draw("same");
	
	
   
   
        fit_func7->SetParameter(0, fit_func3->GetParameter(0));
        fit_func7->SetParameter(1, fit_func3->GetParameter(1));
        fit_func7->SetParameter(2, fit_func3->GetParameter(2));
        fit_func7->SetParameter(3, 0.0025);
        fit_func7->SetParameter(4, 250000);
        fit_func7->SetParameter(5, 0.002331);
        fit_func7->SetParameter(6, 2.1);

   
      	hcalo->Fit("fprec7","RE","",fit_start,fit_stop);
      		
   
   
    
       fit_func11->SetParameter(0, fit_func7->GetParameter(0));
       fit_func11->SetParameter(1, fit_func7->GetParameter(1));
       fit_func11->SetParameter(2, fit_func7->GetParameter(2));
       fit_func11->SetParameter(3, fit_func7->GetParameter(3));
       fit_func11->SetParameter(4, fit_func7->GetParameter(4));
       fit_func11->SetParameter(5, fit_func7->GetParameter(5));
       fit_func11->SetParameter(6, fit_func7->GetParameter(6));
       fit_func11->SetParameter(7, 0.0);
       fit_func11->SetParameter(8, -6);
       fit_func11->SetParameter(9, 0.0);
       fit_func11->SetParameter(10, 2);


    
       hcalo->Fit("fprec11","RE","",fit_start,fit_stop);
      		
       
   
      
       fit_func15->SetParameter(0, fit_func11->GetParameter(0));
       fit_func15->SetParameter(1, fit_func11->GetParameter(1));
       fit_func15->SetParameter(2, fit_func11->GetParameter(2));
       fit_func15->SetParameter(3, fit_func11->GetParameter(3));
       fit_func15->SetParameter(4, fit_func11->GetParameter(4));
       fit_func15->SetParameter(5, fit_func11->GetParameter(5));
       fit_func15->SetParameter(6, fit_func11->GetParameter(6));
       fit_func15->SetParameter(7, fit_func11->GetParameter(7));
       fit_func15->SetParameter(8, fit_func11->GetParameter(8));
       fit_func15->SetParameter(9, fit_func11->GetParameter(9));
       fit_func15->SetParameter(10, fit_func11->GetParameter(10));
       fit_func15->SetParameter(11, 0.01);
       fit_func15->SetParameter(12, 29000);
       fit_func15->SetParameter(13, 0.01402);
       fit_func15->SetParameter(14, 3.5);

   
       hcalo->Fit("fprec15","RE","",fit_start,fit_stop);
      		

       //  return;

       fit_func19->SetParameter(0, fit_func15->GetParameter(0));
       fit_func19->SetParameter(1, fit_func15->GetParameter(1));
       fit_func19->SetParameter(2, fit_func15->GetParameter(2));
       fit_func19->SetParameter(3, fit_func15->GetParameter(3));
       fit_func19->SetParameter(4, fit_func15->GetParameter(4));
       fit_func19->SetParameter(5, fit_func15->GetParameter(5));
       fit_func19->SetParameter(6, fit_func15->GetParameter(6));
       fit_func19->SetParameter(7, fit_func15->GetParameter(7));
       fit_func19->SetParameter(8, fit_func15->GetParameter(8));
       fit_func19->SetParameter(9, fit_func15->GetParameter(9));
       fit_func19->SetParameter(10, fit_func15->GetParameter(10));
       fit_func19->SetParameter(11, fit_func15->GetParameter(11));
       fit_func19->SetParameter(12, fit_func15->GetParameter(12));
       fit_func19->SetParameter(13, fit_func15->GetParameter(13));
       fit_func19->SetParameter(14, fit_func15->GetParameter(14));
       fit_func19->SetParameter(15, 0.001);
       fit_func19->SetParameter(16, 300000);
       fit_func19->SetParameter(17, 0.01392);
       fit_func19->SetParameter(18, 3.1);

      
       hcalo->Fit("fprec19","RE","",fit_start,fit_stop);
       


       fit_func23->SetParameter(0, fit_func19->GetParameter(0));
       fit_func23->SetParameter(1, fit_func19->GetParameter(1));
       fit_func23->SetParameter(2, fit_func19->GetParameter(2));
       fit_func23->SetParameter(3, fit_func19->GetParameter(3));
       fit_func23->SetParameter(4, fit_func19->GetParameter(4));
       fit_func23->SetParameter(5, fit_func19->GetParameter(5));
       fit_func23->SetParameter(6, fit_func19->GetParameter(6));
       fit_func23->SetParameter(7, fit_func19->GetParameter(7));
       fit_func23->SetParameter(8, fit_func19->GetParameter(8));
       fit_func23->SetParameter(9, fit_func19->GetParameter(9));
       fit_func23->SetParameter(10, fit_func19->GetParameter(10));
       fit_func23->SetParameter(11, fit_func19->GetParameter(11));
       fit_func23->SetParameter(12, fit_func19->GetParameter(12));
       fit_func23->SetParameter(13, fit_func19->GetParameter(13));
       fit_func23->SetParameter(14, fit_func19->GetParameter(14));
       fit_func23->SetParameter(15, fit_func19->GetParameter(11));
       fit_func23->SetParameter(16, fit_func19->GetParameter(12));
       fit_func23->SetParameter(17, fit_func19->GetParameter(13));
       fit_func23->SetParameter(18, fit_func19->GetParameter(14));
       fit_func23->FixParameter(19, 0.000);
       fit_func23->FixParameter(20, fit_func19->GetParameter(4)/2);
       fit_func23->FixParameter(21, fit_func19->GetParameter(5)*2);
       fit_func23->FixParameter(22, 0);

      
       hcalo->Fit("fprec23","RE","",fit_start,fit_stop);

       
	
       gStyle->SetOptFit(1111);
            cauto=new TCanvas("auto-correlation", "auto-correlation");
       cauto->Divide(6,4);

        cfit=new TCanvas("23 parameter fit","23 parameter fit");


       
          for(int i=1; i<=18;i++)
       {

	cout<<"fit start time"<<fit_start<<" ns"<<endl; 

   
          

       fit_func23->SetParameter(0, fit_func23->GetParameter(0));
       fit_func23->SetParameter(1, fit_func23->GetParameter(1));
       fit_func23->SetParameter(2, fit_func23->GetParameter(2));
       fit_func23->SetParameter(3, fit_func23->GetParameter(3));
       fit_func23->SetParameter(4, fit_func23->GetParameter(4));
       fit_func23->SetParameter(5, fit_func23->GetParameter(5));
       fit_func23->SetParameter(6, fit_func23->GetParameter(6));
       fit_func23->SetParameter(7, fit_func23->GetParameter(7));
       fit_func23->SetParameter(8, fit_func23->GetParameter(8));
       fit_func23->SetParameter(9, fit_func23->GetParameter(9));
       fit_func23->SetParameter(10, fit_func23->GetParameter(10));
       fit_func23->SetParameter(11, fit_func23->GetParameter(11));
       fit_func23->SetParameter(12, fit_func23->GetParameter(12));
       fit_func23->SetParameter(13, fit_func23->GetParameter(13));
       fit_func23->SetParameter(14, fit_func23->GetParameter(14));
       fit_func23->SetParameter(15, fit_func23->GetParameter(15));
       fit_func23->SetParameter(16, fit_func23->GetParameter(16));
       fit_func23->SetParameter(17, fit_func23->GetParameter(17));
       fit_func23->SetParameter(18, fit_func23->GetParameter(18));
       fit_func23->FixParameter(19, fit_func23->GetParameter(19));
       fit_func23->FixParameter(20, fit_func23->GetParameter(20));
       fit_func23->FixParameter(21, fit_func23->GetParameter(21));
       fit_func23->FixParameter(22, fit_func23->GetParameter(22));

        hcalo->Fit("fprec23","RE","",fit_start,fit_stop);
	

       	h_res[i]= new TH1D("residual histogram ", "h_res", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalo->GetBinContent(ibin)- fit_func23->Eval( hcalo->GetBinCenter(ibin) ) );
          if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.05);
        h_res[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(h_res[i]->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.05);
       hfft[i]->GetYaxis()->SetLabelSize(0.05);

     
       hFFTsq[i] = (TH1D*)hfft[i]->Clone();
       hFFTsq[i]->Reset();
       for (int ib = 1; ib <= hFFTsq[i]->GetNbinsX(); ib++)
        {
         hFFTsq[i]->SetBinContent( ib, hfft[i]->GetBinContent(ib)*hfft[i]->GetBinContent(ib));
        }


       hAuto[i] = hFFTsq[i]->FFT(hAuto[i], "MAG");
       //hAuto[i]->SetTitle("auto-correlation");
       hAuto[i]->GetXaxis()->SetLabelSize(0.04);
       hAuto[i]->GetXaxis()->SetTitle("bins");
       hAuto[i]->GetXaxis()->SetRangeUser(0,h_res[i]->GetNbinsX()/2);
       hAuto[i]->SetLineColor(kRed);
       hAuto[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       cauto->cd(i);
       hAuto[i]->Draw("HIST");

       corr_coeff=hAuto[i]->GetBinContent(2)/hAuto[i]->GetBinContent(1);

       hfft[i]->Reset();
       h_res[i]->Reset();

       data2=ratio_cov_matrix(hcalo);

       TMatrixD cov2(mdim,mdim);

       cov2.SetMatrixArray(data2.GetArray());

       cov2.ResizeTo(hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), -1);

       cov2.SetTol(1.e-23);
       Double_t det1;
       cov2.Invert(&det1);

       cov.ResizeTo(hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), -1);

;
       for(int ix=hcalo->FindBin(fit_start);ix<=cov.GetRowUpb();ix++)
	 {
	   for(int jx=hcalo->FindBin(fit_start);jx<=cov.GetColUpb();jx++)
	     {
	       cov[ix][jx]=cov2[ix][jx];
	     }
	 }
   
      printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hcalo);
       gFitter->SetFCN(chi2);

      
       cfit->cd(1);

       hcalo->Fit("fprec23","RUE","",fit_start,fit_stop);

       countfcn=0;
       gStyle->SetOptFit(1111);
 
       gettimeofday(&t_end, NULL);
       printf("QRatio fit duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));


       A[m]=fit_func23->GetParameter(0);
       blindR[m]=fit_func23->GetParameter(1);
       phi_0[m]=fit_func23->GetParameter(2);
       A_cbo_N[m]=TMath::Abs(fit_func23->GetParameter(3));
       tau_cbo[m]=fit_func23->GetParameter(4);
       omega_cbo[m]=fit_func23->GetParameter(5);
       phi_cbo_N[m]=fit_func23->GetParameter(6);
       A_cbo_A[m]=TMath::Abs(fit_func23->GetParameter(7));
       phi_cbo_A[m]=fit_func23->GetParameter(8);
       A_cbo_phi[m]=TMath::Abs(fit_func23->GetParameter(9));
       phi_cbo_phi[m]=fit_func23->GetParameter(10);
       A_vw[m]=TMath::Abs(fit_func23->GetParameter(11));
       tau_vw[m]=fit_func23->GetParameter(12);
       omega_vw[m]=fit_func23->GetParameter(13);
       phi_vw[m]=fit_func23->GetParameter(14);
       A_y[m]=TMath::Abs(fit_func23->GetParameter(15));
       tau_y[m]=fit_func23->GetParameter(16);
       omega_y[m]=fit_func23->GetParameter(17);
       phi_y[m]=fit_func23->GetParameter(18);
       A_2cbo[m]=TMath::Abs(fit_func23->GetParameter(19));
       tau_2cbo[m]=fit_func23->GetParameter(20);
       omega_2cbo[m]=fit_func23->GetParameter(21);
       phi_2cbo[m]=fit_func23->GetParameter(22);

       dA[m]=fit_func23->GetParError(0);
       dblindR[m]=fit_func23->GetParError(1);
       dphi_0[m]=fit_func23->GetParError(2);
       dA_cbo_N[m]=fit_func23->GetParError(3);
       dtau_cbo[m]=fit_func23->GetParError(4);
       domega_cbo[m]=fit_func23->GetParError(5);
       dphi_cbo_N[m]=fit_func23->GetParError(6);
       dA_cbo_A[m]=fit_func23->GetParError(7);
       dphi_cbo_A[m]=fit_func23->GetParError(8);
       dA_cbo_phi[m]=fit_func23->GetParError(9);
       dphi_cbo_phi[m]=fit_func23->GetParError(10);
       dA_vw[m]=fit_func23->GetParError(11);
       dtau_vw[m]=fit_func23->GetParError(12);
       domega_vw[m]=fit_func23->GetParError(13);
       dphi_vw[m]=fit_func23->GetParError(14);
       dA_y[m]=fit_func23->GetParError(15);
       dtau_y[m]=fit_func23->GetParError(16);
       domega_y[m]=fit_func23->GetParError(17);
       dphi_y[m]=fit_func23->GetParError(18);
       dA_2cbo[m]=fit_func23->GetParError(19);
       dtau_2cbo[m]=fit_func23->GetParError(20);
       domega_2cbo[m]=fit_func23->GetParError(21);
       dphi_2cbo[m]=fit_func23->GetParError(22);

       chisq[m]=fit_func23->GetChisquare()/fit_func23->GetNDF();


       while(phi_cbo_N[m] < 0.0)
	 {
	   phi_cbo_N[m] = phi_cbo_N[m] + (2*(TMath::Pi()));
	 }

       while(phi_cbo_A[m] < 0.0)
	 {
	   phi_cbo_A[m] = phi_cbo_A[m] + (2*(TMath::Pi()));
	 }

       while(phi_cbo_phi[m] < 0.0)
	 {
	   phi_cbo_phi[m] = phi_cbo_phi[m] + (2*(TMath::Pi()));
	 }
       
       while(phi_vw[m] < 0.0)
	 {
	   phi_vw[m] = phi_vw[m] + (2*(TMath::Pi()));
	 }

       while(phi_y[m] < 0.0)
	 {
	   phi_y[m] = phi_y[m] + (2*(TMath::Pi()));
	 }

       while(phi_2cbo[m] < 0.0)
	 {
	   phi_2cbo[m] = phi_2cbo[m] + (2*(TMath::Pi()));
	 }

         while(phi_cbo_N[m] >  (2*(TMath::Pi())))
	 {
	   phi_cbo_N[m] = phi_cbo_N[m] - (2*(TMath::Pi()));
	 }

       while(phi_cbo_A[m] >  (2*(TMath::Pi())))
	 {
	   phi_cbo_A[m] = phi_cbo_A[m] - (2*(TMath::Pi()));
	 }

       while(phi_cbo_phi[m] >  (2*(TMath::Pi())))
	 {
	   phi_cbo_phi[m] = phi_cbo_phi[m] - (2*(TMath::Pi()));
	 }

       while(phi_vw[m] >  (2*(TMath::Pi())))
	 {
	   phi_vw[m] = phi_vw[m] - (2*(TMath::Pi()));
	 }

       while(phi_y[m] >  (2*(TMath::Pi())))
	 {
	   phi_y[m] = phi_y[m] - (2*(TMath::Pi()));
	 }

       while(phi_2cbo[m] >  (2*(TMath::Pi())))
	 {
	   phi_2cbo[m] = phi_2cbo[m] - (2*(TMath::Pi()));
	 }

       if(m==0)
	 {
	   cout<<"m = "<<m<<endl;
	   kband1p[m]=A[m];
	   kband1m[m]=A[m];
	   kband2p[m]=blindR[m];
	   kband2m[m]=blindR[m];
	   kband3p[m]=phi_0[m];
	   kband3m[m]=phi_0[m];
           kband4p[m]=A_cbo_N[m];
	   kband4m[m]=A_cbo_N[m];
           kband5p[m]=tau_cbo[m];
	   kband5m[m]=tau_cbo[m];
           kband6p[m]=omega_cbo[m];
	   kband6m[m]=omega_cbo[m];
           kband7p[m]=phi_cbo_N[m];
	   kband7m[m]=phi_cbo_N[m];
           kband8p[m]=A_cbo_A[m];
	   kband8m[m]=A_cbo_A[m];
           kband9p[m]=phi_cbo_A[m];
	   kband9m[m]=phi_cbo_A[m];
           kband10p[m]=A_cbo_phi[m];
	   kband10m[m]=A_cbo_phi[m];
           kband11p[m]=phi_cbo_phi[m];
	   kband11m[m]=phi_cbo_phi[m];
           kband12p[m]=A_vw[m];
	   kband12m[m]=A_vw[m];
           kband13p[m]=tau_vw[m];
	   kband13m[m]=tau_vw[m];
           kband14p[m]=omega_vw[m];
	   kband14m[m]=omega_vw[m];
           kband15p[m]=phi_vw[m];
	   kband15m[m]=phi_vw[m];
           kband16p[m]=A_y[m];
	   kband16m[m]=A_y[m];
           kband17p[m]=tau_y[m];
	   kband17m[m]=tau_y[m];
           kband18p[m]=omega_y[m];
	   kband18m[m]=omega_y[m];
           kband19p[m]=phi_y[m];
	   kband19m[m]=phi_y[m];
           kband20p[m]=A_2cbo[m];
	   kband20m[m]=A_2cbo[m];
           kband21p[m]=tau_2cbo[m];
	   kband21m[m]=tau_2cbo[m];
           kband22p[m]=omega_2cbo[m];
	   kband22m[m]=omega_2cbo[m];
           kband23p[m]=phi_2cbo[m];
	   kband23m[m]=phi_2cbo[m];

	 }

         else
	 {
	   cout<<"m = "<<m<<endl;
	   kband1p[m]=kband1p[0] + sqrt( TMath::Abs(dA[m]*dA[m] - dA[0]*dA[0]) );
	   kband1m[m]=kband1m[0] - sqrt( TMath::Abs(dA[m]*dA[m] - dA[0]*dA[0]) );
	   kband2p[m]=kband2p[0] + sqrt( TMath::Abs(dblindR[m]*dblindR[m] - dblindR[0]*dblindR[0]) );
					 kband2m[m]=kband2m[0] - sqrt( TMath::Abs(dblindR[m]*dblindR[m] - dblindR[0]*dblindR[0]) );
					 kband3p[m]=kband3p[0] + sqrt( TMath::Abs(dphi_0[m]*dphi_0[m] - dphi_0[0]*dphi_0[0]) );
					 kband3m[m]=kband3m[0] - sqrt( TMath::Abs(dphi_0[m]*dphi_0[m] - dphi_0[0]*dphi_0[0]) );
					 kband4p[m]=kband4p[0] + sqrt( TMath::Abs(dA_cbo_N[m]*dA_cbo_N[m] - dA_cbo_N[0]*dA_cbo_N[0]) );
					 kband4m[m]=kband4m[0] - sqrt( TMath::Abs(dA_cbo_N[m]*dA_cbo_N[m] - dA_cbo_N[0]*dA_cbo_N[0]) );
					 kband5p[m]=kband5p[0] + sqrt( TMath::Abs(dtau_cbo[m]*dtau_cbo[m] - dtau_cbo[0]*dtau_cbo[0]) );
					 kband5m[m]=kband5m[0] - sqrt( TMath::Abs(dtau_cbo[m]*dtau_cbo[m] - dtau_cbo[0]*dtau_cbo[0]) );
					 kband6p[m]=kband6p[0] + sqrt( TMath::Abs(domega_cbo[m]*domega_cbo[m] - domega_cbo[0]*domega_cbo[0]) );
					 kband6m[m]=kband6m[0] - sqrt( TMath::Abs(domega_cbo[m]*domega_cbo[m] - domega_cbo[0]*domega_cbo[0]) );
					 kband7p[m]=kband7p[0] + sqrt( TMath::Abs(dphi_cbo_N[m]*dphi_cbo_N[m] - dphi_cbo_N[0]*dphi_cbo_N[0]) );
					 kband7m[m]=kband7m[0] - sqrt( TMath::Abs(dphi_cbo_N[m]*dphi_cbo_N[m] - dphi_cbo_N[0]*dphi_cbo_N[0]) );
					 kband8p[m]=kband8p[0] + sqrt( TMath::Abs(dA_cbo_A[m]*dA_cbo_A[m] - dA_cbo_A[0]*dA_cbo_A[0]) );
					 kband8m[m]=kband8m[0] - sqrt( TMath::Abs(dA_cbo_A[m]*dA_cbo_A[m] - dA_cbo_A[0]*dA_cbo_A[0]) );
					 kband9p[m]=kband9p[0] + sqrt( TMath::Abs(dphi_cbo_A[m]*dphi_cbo_A[m] - dphi_cbo_A[0]*dphi_cbo_A[0]) );
					 kband9m[m]=kband9m[0] - sqrt( TMath::Abs(dphi_cbo_A[m]*dphi_cbo_A[m] - dphi_cbo_A[0]*dphi_cbo_A[0]) );
					 kband10p[m]=kband10p[0] + sqrt( TMath::Abs(dA_cbo_phi[m]*dA_cbo_phi[m] - dA_cbo_phi[0]*dA_cbo_phi[0]) );
					   kband10m[m]=kband10m[0] - sqrt( TMath::Abs(dA_cbo_phi[m]*dA_cbo_phi[m] - dA_cbo_phi[0]*dA_cbo_phi[0]) );
					   kband11p[m]=kband11p[0] + sqrt( TMath::Abs(dphi_cbo_phi[m]*dphi_cbo_phi[m] - dphi_cbo_phi[0]*dphi_cbo_phi[0]) );
					   kband11m[m]=kband11m[0] - sqrt( TMath::Abs(dphi_cbo_phi[m]*dphi_cbo_phi[m] - dphi_cbo_phi[0]*dphi_cbo_phi[0]) );
					   kband12p[m]=kband12p[0] + sqrt( TMath::Abs(dA_vw[m]*dA_vw[m] - dA_vw[0]*dA_vw[0]) );
					   kband12m[m]=kband12m[0] - sqrt( TMath::Abs(dA_vw[m]*dA_vw[m] - dA_vw[0]*dA_vw[0]) );
					   kband13p[m]=kband13p[0] + sqrt( TMath::Abs(dtau_vw[m]*dtau_vw[m] - dtau_vw[0]*dtau_vw[0]) );
					   kband13m[m]=kband13m[0] - sqrt( TMath::Abs(dtau_vw[m]*dtau_vw[m] - dtau_vw[0]*dtau_vw[0]) );
					   kband14p[m]=kband14p[0] + sqrt( TMath::Abs(domega_vw[m]*domega_vw[m] - domega_vw[0]*domega_vw[0]) );
					   kband14m[m]=kband14m[0] - sqrt( TMath::Abs(domega_vw[m]*domega_vw[m] - domega_vw[0]*domega_vw[0]) );
					   kband15p[m]=kband15p[0] + sqrt( TMath::Abs(dphi_vw[m]*dphi_vw[m] - dphi_vw[0]*dphi_vw[0]) );
					   kband15m[m]=kband15m[0] - sqrt( TMath::Abs(dphi_vw[m]*dphi_vw[m] - dphi_vw[0]*dphi_vw[0]) );
					   kband16p[m]=kband16p[0] + sqrt( TMath::Abs(dA_y[m]*dA_y[m] - dA_y[0]*dA_y[0]) );
					   kband16m[m]=kband16m[0] - sqrt( TMath::Abs(dA_y[m]*dA_y[m] - dA_y[0]*dA_y[0]) );
					   kband17p[m]=kband17p[0] + sqrt( TMath::Abs(dtau_y[m]*dtau_y[m] - dtau_y[0]*dtau_y[0]) );
					   kband17m[m]=kband17m[0] - sqrt( TMath::Abs(dtau_y[m]*dtau_y[m] - dtau_y[0]*dtau_y[0]) );
					   kband18p[m]=kband18p[0] + sqrt( TMath::Abs(domega_y[m]*domega_y[m] - domega_y[0]*domega_y[0]) );
					   kband18m[m]=kband18m[0] - sqrt( TMath::Abs(domega_y[m]*domega_y[m] - domega_y[0]*domega_y[0]) );
					   kband19p[m]=kband19p[0] + sqrt( TMath::Abs(dphi_y[m]*dphi_y[m] - dphi_y[0]*dphi_y[0]) );
					   kband19m[m]=kband19m[0] - sqrt( TMath::Abs(dphi_y[m]*dphi_y[m] - dphi_y[0]*dphi_y[0]) );
					   kband20p[m]=kband20p[0] + sqrt( TMath::Abs(dA_2cbo[m]*dA_2cbo[m] - dA_2cbo[0]*dA_2cbo[0]) );
					   kband20m[m]=kband20m[0] - sqrt( TMath::Abs(dA_2cbo[m]*dA_2cbo[m] - dA_2cbo[0]*dA_2cbo[0]) );
					   kband21p[m]=kband21p[0] + sqrt( TMath::Abs(dtau_2cbo[m]*dtau_2cbo[m] - dtau_2cbo[0]*dtau_2cbo[0]) );
					   kband21m[m]=kband21m[0] - sqrt( TMath::Abs(dtau_2cbo[m]*dtau_2cbo[m] - dtau_2cbo[0]*dtau_2cbo[0]) );
					   kband22p[m]=kband22p[0] + sqrt( TMath::Abs(domega_2cbo[m]*domega_2cbo[m] - domega_2cbo[0]*domega_2cbo[0]) );
					   kband22m[m]=kband22m[0] - sqrt( TMath::Abs(domega_2cbo[m]*domega_2cbo[m] - domega_2cbo[0]*domega_2cbo[0]) );
					   kband23p[m]=kband23p[0] + sqrt( TMath::Abs(dphi_2cbo[m]*dphi_2cbo[m] - dphi_2cbo[0]*dphi_2cbo[0]) );
					   kband23m[m]=kband23m[0] - sqrt( TMath::Abs(dphi_2cbo[m]*dphi_2cbo[m] - dphi_2cbo[0]*dphi_2cbo[0]) );

	 }

       n[m]=fit_start;
       m=m+1;
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalo->GetBinContent(ibin)- fit_func23->Eval( hcalo->GetBinCenter(ibin) ) );
          if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.05);
        h_res[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hcalo->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.05);
       hfft[i]->GetYaxis()->SetLabelSize(0.05);

       fit_start=fit_start+1000;
 
   }
     
  c2 = new TCanvas("c2","run3no start time residuals");
  c2->Divide(6,3);
  c3 = new TCanvas("c3","run3no start time fft");
  c3->Divide(6,3);
   
 for(int i=1;i<=18;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hfft[i]->Draw();
   }
    c1=new TCanvas("c1","run3no blind R vs fit start time");
    c1->Divide(2,2);
    c1->cd(1);
    gr1=new TGraphErrors(m,n,chisq,0,0);
    gr1->SetTitle("chisq vs fit start time");
    gr1->GetXaxis()->SetTitle("Start time [ns]");
    gr1->GetYaxis()->SetTitle("reduced chi-sq");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    

     c1->cd(2);
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle("A vs fit start time");
    gr2->GetXaxis()->SetTitle("Start time [ns]");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw();
    gkband1p=new TGraph(m,n,kband1p);
    gkband1m=new TGraph(m,n,kband1m);
    gkband1m->SetLineColor(kBlue);
    gkband1p->SetLineColor(kBlue);
    gkband1p->Draw("same");
    gkband1m->Draw("same");

     c1->cd(3);
    gr3=new TGraphErrors(m,n,blindR,0,dblindR);
    gr3->SetTitle("blindR vs fit start time");
    gr3->GetXaxis()->SetTitle("Start time [ns]");
    gr3->GetYaxis()->SetTitle("Blind R [ppm]");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw();
    gkband2p=new TGraph(m,n,kband2p);
    gkband2m=new TGraph(m,n,kband2m);
    gkband2m->SetLineColor(kBlue);
    gkband2p->SetLineColor(kBlue);
    gkband2p->Draw("same");
    gkband2m->Draw("same");


     c1->cd(4);
    gr4=new TGraphErrors(m,n,phi_0,0,dphi_0);
    gr4->SetTitle("phase vs fit start time");
    gr4->GetXaxis()->SetTitle("Start time [ns]");
    gr4->GetYaxis()->SetTitle("phi_0");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->Draw();
    gkband3p=new TGraph(m,n,kband3p);
    gkband3m=new TGraph(m,n,kband3m);
    gkband3m->SetLineColor(kBlue);
    gkband3p->SetLineColor(kBlue);
    gkband3p->Draw("same");
    gkband3m->Draw("same");


    c4=new TCanvas("c4","run3no cbo vs fit start time");
    c4->Divide(2,2);
    c4->cd(1);
    gr8=new TGraphErrors(m,n,A_cbo_N,0,dA_cbo_N);
    gr8->SetTitle("A_cbo_N vs fit start time");
    gr8->GetXaxis()->SetTitle("Start time [ns]");
    gr8->GetYaxis()->SetTitle("A_cbo_N");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kRed);
    gr8->Draw();
    gkband4p=new TGraph(m,n,kband4p);
    gkband4m=new TGraph(m,n,kband4m);
    gkband4m->SetLineColor(kBlue);
    gkband4p->SetLineColor(kBlue);
    gkband4p->Draw("same");
    gkband4m->Draw("same");


     c4->cd(2);
    gr5=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr5->SetTitle("tau_cbo_ vs fit start time");
    gr5->GetXaxis()->SetTitle("Start time [ns]");
    gr5->GetYaxis()->SetTitle("tau_cbo");
    gr5->SetMarkerStyle(20);
    gr5->SetLineColor(kRed);
    gr5->Draw();
    gkband5p=new TGraph(m,n,kband5p);
    gkband5m=new TGraph(m,n,kband5m);
    gkband5m->SetLineColor(kBlue);
    gkband5p->SetLineColor(kBlue);
    gkband5p->Draw("same");
    gkband5m->Draw("same");


     c4->cd(3);
    gr6=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr6->SetTitle("omega_cbo vs fit start time");
    gr6->GetXaxis()->SetTitle("Start time [ns]");
    gr6->GetYaxis()->SetTitle("omega_cbo");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kRed);
    gr6->Draw();
    gkband6p=new TGraph(m,n,kband6p);
    gkband6m=new TGraph(m,n,kband6m);
    gkband6m->SetLineColor(kBlue);
    gkband6p->SetLineColor(kBlue);
    gkband6p->Draw("same");
    gkband6m->Draw("same");


     c4->cd(4);
    gr7=new TGraphErrors(m,n,phi_cbo_N,0,dphi_cbo_N);
    gr7->SetTitle("phi_cbo_N vs fit start time");
    gr7->GetXaxis()->SetTitle("Start time [ns]");
    gr7->GetYaxis()->SetTitle("phi_cbo_N");
    gr7->SetMarkerStyle(20);
    gr7->SetLineColor(kRed);
    gr7->Draw();
    gkband7p=new TGraph(m,n,kband7p);
    gkband7m=new TGraph(m,n,kband7m);
    gkband7m->SetLineColor(kBlue);
    gkband7p->SetLineColor(kBlue);
    gkband7p->Draw("same");
    gkband7m->Draw("same");


    c6=new TCanvas("c6","run3no A,phi,cbo vs fit start time");
    c6->Divide(2,2);
    c6->cd(1);
    gr13=new TGraphErrors(m,n,A_cbo_A,0,dA_cbo_A);
    gr13->SetTitle("A_cbo_A vs fit start time");
    gr13->GetXaxis()->SetTitle("Start time [ns]");
    gr13->GetYaxis()->SetTitle("A_cbo_A");
    gr13->SetMarkerStyle(20);
    gr13->SetLineColor(kRed);
    gr13->Draw();
    gkband8p=new TGraph(m,n,kband8p);
    gkband8m=new TGraph(m,n,kband8m);
    gkband8m->SetLineColor(kBlue);
    gkband8p->SetLineColor(kBlue);
    gkband8p->Draw("same");
    gkband8m->Draw("same");

     c6->cd(2);
    gr14=new TGraphErrors(m,n,phi_cbo_A,0,dphi_cbo_A);
    gr14->SetTitle("phi_cbo_A vs fit start time");
    gr14->GetXaxis()->SetTitle("Start time [ns]");
    gr14->GetYaxis()->SetTitle("phi_cbo_A");
    gr14->SetMarkerStyle(20);
    gr14->SetLineColor(kRed);
    gr14->Draw();
    gkband9p=new TGraph(m,n,kband9p);
    gkband9m=new TGraph(m,n,kband9m);
    gkband9m->SetLineColor(kBlue);
    gkband9p->SetLineColor(kBlue);
    gkband9p->Draw("same");
    gkband9m->Draw("same");

     c6->cd(3);
    gr15=new TGraphErrors(m,n,A_cbo_phi,0,dA_cbo_phi);
    gr15->SetTitle("A_cbo_phi vs fit start time");
    gr15->GetXaxis()->SetTitle("Start time [ns]");
    gr15->GetYaxis()->SetTitle("A_cbo_phi");
    gr15->SetMarkerStyle(20);
    gr15->SetLineColor(kRed);
    gr15->Draw();
    gkband10p=new TGraph(m,n,kband10p);
    gkband10m=new TGraph(m,n,kband10m);
    gkband10m->SetLineColor(kBlue);
    gkband10p->SetLineColor(kBlue);
    gkband10p->Draw("same");
    gkband10m->Draw("same");

     c6->cd(4);
    gr16=new TGraphErrors(m,n,phi_cbo_phi,0,dphi_cbo_phi);
    gr16->SetTitle("phi_cbo_phi vs fit start time");
    gr16->GetXaxis()->SetTitle("Start time [ns]");
    gr16->GetYaxis()->SetTitle("phi_cbo_phi");
    gr16->SetMarkerStyle(20);
    gr16->SetLineColor(kRed);
    gr16->Draw();
    gkband11p=new TGraph(m,n,kband11p);
    gkband11m=new TGraph(m,n,kband11m);
    gkband11m->SetLineColor(kBlue);
    gkband11p->SetLineColor(kBlue);
    gkband11p->Draw("same");
    gkband11m->Draw("same");


     c7=new TCanvas("c7","run3no 2cbo vs fit start time");
     c7->Divide(2,2);
    c7->cd(1);
    gr17=new TGraphErrors(m,n,A_2cbo,0,dA_2cbo);
    gr17->SetTitle("A_2cbo vs fit start time");
    gr17->GetXaxis()->SetTitle("Start time [ns]");
    gr17->GetYaxis()->SetTitle("A_2cbo");
    gr17->SetMarkerStyle(20);
    gr17->SetLineColor(kRed);
    gr17->Draw();
    gkband20p=new TGraph(m,n,kband20p);
    gkband20m=new TGraph(m,n,kband20m);
    gkband20m->SetLineColor(kBlue);
    gkband20p->SetLineColor(kBlue);
    gkband20p->Draw("same");
    gkband20m->Draw("same");


     c7->cd(2);
    gr18=new TGraphErrors(m,n,tau_2cbo,0,dtau_2cbo);
    gr18->SetTitle("tau_2cbo_ vs fit start time");
    gr18->GetXaxis()->SetTitle("Start time [ns]");
    gr18->GetYaxis()->SetTitle("tau_2cbo");
    gr18->SetMarkerStyle(20);
    gr18->SetLineColor(kRed);
    gr18->Draw();
    gkband21p=new TGraph(m,n,kband21p);
    gkband21m=new TGraph(m,n,kband21m);
    gkband21m->SetLineColor(kBlue);
    gkband21p->SetLineColor(kBlue);
    gkband21p->Draw("same");
    gkband21m->Draw("same");


     c7->cd(3);
    gr19=new TGraphErrors(m,n,omega_2cbo,0,domega_2cbo);
    gr19->SetTitle("omega_2cbo vs fit start time");
    gr19->GetXaxis()->SetTitle("Start time [ns]");
    gr19->GetYaxis()->SetTitle("omega_2cbo");
    gr19->SetMarkerStyle(20);
    gr19->SetLineColor(kRed);
    gr19->Draw();
    gkband22p=new TGraph(m,n,kband22p);
    gkband22m=new TGraph(m,n,kband22m);
    gkband22m->SetLineColor(kBlue);
    gkband22p->SetLineColor(kBlue);
    gkband22p->Draw("same");
    gkband22m->Draw("same");


     c7->cd(4);
    gr20=new TGraphErrors(m,n,phi_2cbo,0,dphi_2cbo);
    gr20->SetTitle("phi_2cbo vs fit start time");
    gr20->GetXaxis()->SetTitle("Start time [ns]");
    gr20->GetYaxis()->SetTitle("phi_2cbo");
    gr20->SetMarkerStyle(20);
    gr20->SetLineColor(kRed);
    gr20->Draw();
    gkband23p=new TGraph(m,n,kband23p);
    gkband23m=new TGraph(m,n,kband23m);
    gkband23m->SetLineColor(kBlue);
    gkband23p->SetLineColor(kBlue);
    gkband23p->Draw("same");
    gkband23m->Draw("same");


    c8=new TCanvas("c8","run3no vw vs fit start time");
    c8->Divide(2,2);
    c8->cd(1);
    gr21=new TGraphErrors(m,n,A_vw,0,dA_vw);
    gr21->SetTitle("A_vw vs fit start time");
    gr21->GetXaxis()->SetTitle("Start time [ns]");
    gr21->GetYaxis()->SetTitle("A_vw");
    gr21->SetMarkerStyle(20);
    gr21->SetLineColor(kRed);
    gr21->Draw();
    gkband12p=new TGraph(m,n,kband12p);
    gkband12m=new TGraph(m,n,kband12m);
    gkband12m->SetLineColor(kBlue);
    gkband12p->SetLineColor(kBlue);
    gkband12p->Draw("same");
    gkband12m->Draw("same");


    c8->cd(2);
    gr22=new TGraphErrors(m,n,tau_vw,0,dtau_vw);
    gr22->SetTitle("tau_vw vs fit start time");
    gr22->GetXaxis()->SetTitle("Start time [ns]");
    gr22->GetYaxis()->SetTitle("tau_vw");
    gr22->SetMarkerStyle(20);
    gr22->SetLineColor(kRed);
    gr22->Draw();
    gkband13p=new TGraph(m,n,kband13p);
    gkband13m=new TGraph(m,n,kband13m);
    gkband13m->SetLineColor(kBlue);
    gkband13p->SetLineColor(kBlue);
    gkband13p->Draw("same");
    gkband13m->Draw("same");


    c8->cd(3);
    gr23=new TGraphErrors(m,n,omega_vw,0,domega_vw);
    gr23->SetTitle("omega_vw vs fit start time");
    gr23->GetXaxis()->SetTitle("Start time [ns]");
    gr23->GetYaxis()->SetTitle("omega_vw");
    gr23->SetMarkerStyle(20);
    gr23->SetLineColor(kRed);
    gr23->Draw();
    gkband14p=new TGraph(m,n,kband14p);
    gkband14m=new TGraph(m,n,kband14m);
    gkband14m->SetLineColor(kBlue);
    gkband14p->SetLineColor(kBlue);
    gkband14p->Draw("same");
    gkband14m->Draw("same");


    c8->cd(4);
    gr24=new TGraphErrors(m,n,phi_vw,0,dphi_vw);
    gr24->SetTitle("phi_vw vs fit start time");
    gr24->GetXaxis()->SetTitle("Start time [ns]");
    gr24->GetYaxis()->SetTitle("phi_vw");
    gr24->SetMarkerStyle(20);
    gr24->SetLineColor(kRed);
    gr24->Draw();
    gkband15p=new TGraph(m,n,kband15p);
    gkband15m=new TGraph(m,n,kband15m);
    gkband15m->SetLineColor(kBlue);
    gkband15p->SetLineColor(kBlue);
    gkband15p->Draw("same");
    gkband15m->Draw("same");


    c9=new TCanvas("c9","run3no y vs fit start time");
    c9->Divide(2,2);
    c9->cd(1);
    gr25=new TGraphErrors(m,n,A_y,0,dA_y);
    gr25->SetTitle("A_y vs fit start time");
    gr25->GetXaxis()->SetTitle("Start time [ns]");
    gr25->GetYaxis()->SetTitle("A_y");
    gr25->SetMarkerStyle(20);
    gr25->SetLineColor(kRed);
    gr25->Draw();
    gkband16p=new TGraph(m,n,kband16p);
    gkband16m=new TGraph(m,n,kband16m);
    gkband16m->SetLineColor(kBlue);
    gkband16p->SetLineColor(kBlue);
    gkband16p->Draw("same");
    gkband16m->Draw("same");


    c9->cd(2);
    gr26=new TGraphErrors(m,n,tau_y,0,dtau_y);
    gr26->SetTitle("tau_y vs fit start time");
    gr26->GetXaxis()->SetTitle("Start time [ns]");
    gr26->GetYaxis()->SetTitle("tau_y");
    gr26->SetMarkerStyle(20);
    gr26->SetLineColor(kRed);
    gr26->Draw();
    gkband17p=new TGraph(m,n,kband17p);
    gkband17m=new TGraph(m,n,kband17m);
    gkband17m->SetLineColor(kBlue);
    gkband17p->SetLineColor(kBlue);
    gkband17p->Draw("same");
    gkband17m->Draw("same");


    c9->cd(3);
    gr27=new TGraphErrors(m,n,omega_y,0,domega_y);
    gr27->SetTitle("omega_y vs fit start time");
    gr27->GetXaxis()->SetTitle("Start time [ns]");
    gr27->GetYaxis()->SetTitle("omega_y");
    gr27->SetMarkerStyle(20);
    gr27->SetLineColor(kRed);
    gr27->Draw();
    gkband18p=new TGraph(m,n,kband18p);
    gkband18m=new TGraph(m,n,kband18m);
    gkband18m->SetLineColor(kBlue);
    gkband18p->SetLineColor(kBlue);
    gkband18p->Draw("same");
    gkband18m->Draw("same");


    c9->cd(4);
    gr28=new TGraphErrors(m,n,phi_y,0,dphi_y);
    gr28->SetTitle("phi_y vs fit start time");
    gr28->GetXaxis()->SetTitle("Start time [ns]");
    gr28->GetYaxis()->SetTitle("phi_y");
    gr28->SetMarkerStyle(20);
    gr28->SetLineColor(kRed);
    gr28->Draw();
    gkband19p=new TGraph(m,n,kband19p);
    gkband19m=new TGraph(m,n,kband19m);
    gkband19m->SetLineColor(kBlue);
    gkband19p->SetLineColor(kBlue);
    gkband19p->Draw("same");
    gkband19m->Draw("same");


}
