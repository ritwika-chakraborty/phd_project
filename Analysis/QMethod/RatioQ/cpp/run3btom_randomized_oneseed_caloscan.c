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
char root_file_name[128] = "r3btom_w8_xtalreject_calos.root";
//char root_file_name[128] = "run2C_thresh300_calosum_ratio.root";
//char root_file_name[128] = "run2_thresh300_fbfDQC.root";
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*hr1[25],*hr2[25],*hr3[25],*hr4[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[50], *hr_sum, *hsum1,*hsum2,*hsum3,*hsum4,*hcomp_sum1[25],*hcomp_sum2[25],*hcomp_sum3[25],*hcomp_sum4[25],*hr_shift1,*hr_noshift1,*hr_shift2,*hr_noshift2,*hr_shift3,*hr_noshift3,*hr_shift4,*hr_noshift4;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*c11,*c12,*c13,*cauto,*cfit;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28,*gr29, *gr30, *gr31;
TGraph *gkband1p,*gkband1m,*gkband2p,*gkband2m,*gkband3p,*gkband3m,*gkband4p,*gkband4m,*gkband5p,*gkband5m,*gkband6p,*gkband6m,*gkband7p,*gkband7m,*gkband8p,*gkband8m,*gkband9p,*gkband9m,*gkband10p,*gkband10m,*gkband11p,*gkband11m,*gkband12p,*gkband12m,*gkband13p,*gkband13m,*gkband14p,*gkband14m,*gkband15p,*gkband15m,*gkband16p,*gkband16m,*gkband17p,*gkband17m,*gkband18p,*gkband18m,*gkband19p,*gkband19m,*gkband20p,*gkband20m,*gkband21p,*gkband21m,*gkband22p,*gkband22m,*gkband23p,*gkband23m;
TH1 *hfft[50],*hFFTsq[50],*hAuto[50];
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
bool useFR=false;
Int_t m=0;
Int_t corr_coeff;
Int_t countfcn=0;

Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi_0[50], dphi_0[50], A_cbo_N[50], tau_cbo[50], omega_cbo[50], phi_cbo_N[50], A_cbo_A[50], phi_cbo_A[50], A_cbo_phi[50], phi_cbo_phi[50], A_vw[50], tau_vw[50], omega_vw[50], phi_vw[50], A_y[50], tau_y[50], omega_y[50], phi_y[50], A_2cbo[50], tau_2cbo[50], omega_2cbo[50], phi_2cbo[50],  dA_cbo_N[50], dtau_cbo[50], domega_cbo[50], dphi_cbo_N[50], dA_cbo_A[50], dphi_cbo_A[50], dA_cbo_phi[50], dphi_cbo_phi[50], dA_vw[50], dtau_vw[50], domega_vw[50], dphi_vw[50], dA_y[50], dtau_y[50], domega_y[50], dphi_y[50], dA_2cbo[50], dtau_2cbo[50], domega_2cbo[50], dphi_2cbo[50], n[50],kband1p[50],kband1m[50],kband2p[50],kband2m[50],kband3p[50],kband3m[50],kband4p[50],kband4m[50],kband5p[50],kband5m[50],kband6p[50],kband6m[50],kband7p[50],kband7m[50],kband8p[50],kband8m[50],kband9p[50],kband9m[50],kband10p[50],kband10m[50],kband11p[50],kband11m[50],kband12p[50],kband12m[50],kband13p[50],kband13m[50],kband14p[50],kband14m[50],kband15p[50],kband15m[50],kband16p[50],kband16m[50],kband17p[50],kband17m[50],kband18p[50],kband18m[50],kband19p[50],kband19m[50],kband20p[50],kband20m[50],kband21p[50],kband21m[50],kband22p[50],kband22m[50],kband23p[50],kband23m[50], rval3[50], drval3[50], rval11[50], drval11[50], rval23[50], drval23[50];

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
  //  cov.SetMatrixArray(data.GetArray());


  // cov.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
   
  //   cov.Print();
  //  cov.SetTol(1.e-23);
  //  Double_t det1;
  //  cov.Invert(&det1);
  //   cout<<mdim<<endl;
    //cout<<"Matrix inverted "<<i0<<endl;
  // return cov;
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

    double asym_2cbo= par[11];

    double tau_2cbo = par[12];

    double omega_2cbo = par[13];

    double phi_2cbo = par[14];
    
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

    return (2*f*Ncbo - ff*Ncbof - fb*Ncbob)/(2*f*Ncbo + ff*Ncbof + fb*Ncbob);

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

    double asym_2cbo= par[11];

    double tau_2cbo = par[12];

    double omega_2cbo = par[13];

    double phi_2cbo = par[14];

    double asym_vw= par[15];

    double tau_vw = par[16];

    double omega_vw = par[17];

    double phi_vw = par[18];
    
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

    return (2*f*Ncbo*Nvw - ff*Ncbof*Nvwf - fb*Ncbob*Nvwb)/(2*f*Ncbo*Nvw + ff*Ncbof*Nvwf + fb*Ncbob*Nvwb);

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

    double asym_2cbo= par[11];

    double tau_2cbo = par[12];

    double omega_2cbo = par[13];

    double phi_2cbo = par[14];

    double asym_vw= par[15];

    double tau_vw = par[16];

    double omega_vw = par[17];

    double phi_vw = par[18];

    double asym_vbo= par[19];

    double tau_vbo = par[20];

    double omega_vbo = par[21];

    double phi_vbo = par[22];
    
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




TH1D* construct_rhist_rand(TH1D *hsum1, TH1D *hsum2, TH1D *hsum3, TH1D *hsum4, TH1D *hosum1, TH1D *hosum2, TH1D *hosum3, TH1D *hosum4)
{
 


   //create the component histograms of the ratio histograms   
   h1=(TH1D*)hsum1->Clone();
   h2=(TH1D*)hsum2->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

   //nbinshift=(0.5*T_a_true)/h1->GetBinWidth(1);
   nbinshift=15;
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
 
    cout<<hsum1->GetBinWidth(100)<<endl;
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
	      //cout<<m<<endl;
	    
	     }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}


void run3btom_randomized_oneseed_caloscan()
{

  
  _file[1]=TFile::Open(root_file_name);

    for(int i=1;i<=24;i++)
    {
     hr1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_0_2", i)));

     hr2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_1_2", i)));

     hr3[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_2_2", i)));

     hr4[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_3_2", i)));
    }
  
  /*  hsum1= new TH1D("calo histogram ratio sum 1", "hsum1", hr1[1]->GetNbinsX(), 100001, 352001);
  hsum2= new TH1D("calo histogram ratio sum 2", "hsum2", hr2[1]->GetNbinsX(), 100001, 352001);
  hsum3= new TH1D("calo histogram ratio sum 3", "hsum3", hr3[1]->GetNbinsX(), 100001, 352001);
  hsum4= new TH1D("calo histogram ratio sum 4", "hsum4", hr4[1]->GetNbinsX(), 100001, 352001);
  */

  
  
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
  // hsum2->SetLineColor(kRed);
  //hsum3->SetLineColor(kBlack);
  //hsum4->SetLineColor(kGreen);

  int start_bin=400;
  double binwidth=hr1[1]->GetBinWidth(1);
  int nbins=hr1[1]->GetNbinsX();
  double binlowedge=hr1[1]->GetBinLowEdge(1);
  double binhighedge=hr1[1]->GetBinLowEdge(nbins)+binwidth;
   
  cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   


 for(int i=1;i<=24;i++)
   {
   
    hr2[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr2[i]->Integral(start_bin,16800));
    hr3[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr3[i]->Integral(start_bin,16800));
    hr4[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr4[i]->Integral(start_bin,16800));
 
    // cout<<"Normalized from bin "<<start_bin<<endl;
 
    //h_sum= new TH1D("calo histogram sum", "h_sum", hsum1->GetNbinsX(), 100001, 352001);
  //hrsum= new TH1D("calo histogram sum ratio", "calosum ratio", hsum1->GetNbinsX(), 100001, 352001);

     hr1[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr2[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr3[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr4[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     
     //cout<<hr1[i]->GetBinWidth(1)<<" "<<hr1[i]->GetNbinsX()<<" "<<hr1[i]->GetBinLowEdge(1)<<" "<<(hr1[i]->GetBinLowEdge(hr1[i]->GetNbinsX())+hr1[i]->GetBinWidth(1))<<" "<<endl;

    hr1[i]->Rebin(4);
    hr1[i]->Scale(0.25);
    
    hr2[i]->Rebin(4);
    hr2[i]->Scale(0.25);

    hr3[i]->Rebin(4);
    hr3[i]->Scale(0.25);

    hr4[i]->Rebin(4);
    hr4[i]->Scale(0.25);
   
    if(useFR)
   {
    hr_noshift1=(TH1D*)hr1[i]->Clone();
    hr_noshift2=(TH1D*)hr2[i]->Clone();
    hr_noshift3=(TH1D*)hr3[i]->Clone();
    hr_noshift4=(TH1D*)hr4[i]->Clone();
	
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

	hcomp_sum1[i]=(TH1D*)hr1[i]->Clone();
	hcomp_sum1[i]->Reset();
	hcomp_sum2[i]=(TH1D*)hr2[i]->Clone();
	hcomp_sum2[i]->Reset();
	hcomp_sum3[i]=(TH1D*)hr3[i]->Clone();
	hcomp_sum3[i]->Reset();
	hcomp_sum4[i]=(TH1D*)hr4[i]->Clone();
	hcomp_sum4[i]->Reset();


	for(int ibin=1;ibin<=hr2[i]->GetNbinsX();ibin++)
	  {
	    hcomp_sum1[i]->SetBinContent(ibin,0.5*(hr_noshift1->GetBinContent(ibin)+hr_shift1->GetBinContent(ibin)));
	    hcomp_sum1[i]->SetBinError(ibin,0.5*sqrt((hr_noshift1->GetBinError(ibin)*hr_noshift1->GetBinError(ibin))+(hr_shift1->GetBinError(ibin)*hr_shift1->GetBinError(ibin))));
	    hcomp_sum2[i]->SetBinContent(ibin,0.5*(hr_noshift2->GetBinContent(ibin)+hr_shift2->GetBinContent(ibin)));
	    hcomp_sum2[i]->SetBinError(ibin,0.5*sqrt((hr_noshift2->GetBinError(ibin)*hr_noshift2->GetBinError(ibin))+(hr_shift2->GetBinError(ibin)*hr_shift2->GetBinError(ibin))));
	    hcomp_sum3[i]->SetBinContent(ibin,0.5*(hr_noshift3->GetBinContent(ibin)+hr_shift3->GetBinContent(ibin)));
	    hcomp_sum3[i]->SetBinError(ibin,0.5*sqrt((hr_noshift3->GetBinError(ibin)*hr_noshift3->GetBinError(ibin))+(hr_shift3->GetBinError(ibin)*hr_shift3->GetBinError(ibin))));
	    hcomp_sum4[i]->SetBinContent(ibin,0.5*(hr_noshift4->GetBinContent(ibin)+hr_shift4->GetBinContent(ibin)));
	    hcomp_sum4[i]->SetBinError(ibin,0.5*sqrt((hr_noshift4->GetBinError(ibin)*hr_noshift4->GetBinError(ibin))+(hr_shift4->GetBinError(ibin)*hr_shift4->GetBinError(ibin))));
	  }

	
	hcomp_sum1[i]->Rebin(2);
	hcomp_sum1[i]->Scale(0.5);
	hcomp_sum2[i]->Rebin(2);
	hcomp_sum2[i]->Scale(0.5);
	hcomp_sum3[i]->Rebin(2);
	hcomp_sum3[i]->Scale(0.5);
	hcomp_sum4[i]->Rebin(2);
	hcomp_sum4[i]->Scale(0.5);

	hr_shift1->Reset();
	hr_noshift1->Reset();
	hr_shift2->Reset();
	hr_noshift2->Reset();
	hr_shift3->Reset();
	hr_noshift3->Reset();
        hr_shift4->Reset();
	hr_noshift4->Reset();
   }
  else
   {
     hcomp_sum1[i]=(TH1D*)hr1[i]->Clone();
     hcomp_sum2[i]=(TH1D*)hr2[i]->Clone();
     hcomp_sum3[i]=(TH1D*)hr3[i]->Clone();
     hcomp_sum4[i]=(TH1D*)hr4[i]->Clone();

     hcomp_sum1[i]->Rebin(2);
     hcomp_sum1[i]->Scale(0.5);
     hcomp_sum2[i]->Rebin(2);
     hcomp_sum2[i]->Scale(0.5);
     hcomp_sum3[i]->Rebin(2);
     hcomp_sum3[i]->Scale(0.5);
     hcomp_sum4[i]->Rebin(2);
     hcomp_sum4[i]->Scale(0.5);

   }

  }


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
     fit_func15->SetParName(11,"A_2cbo");
     fit_func15->SetParName(12,"Tau_2cbo");
     fit_func15->SetParName(13,"omega_2cbo");
     fit_func15->SetParName(14,"phi_2cbo");
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
     fit_func19->SetParName(11,"A_2cbo");
     fit_func19->SetParName(12,"Tau_2cbo");
     fit_func19->SetParName(13,"omega_2cbo");
     fit_func19->SetParName(14,"phi_2cbo");
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
     fit_func23->SetParName(11,"A_2cbo");
     fit_func23->SetParName(12,"Tau_2cbo");
     fit_func23->SetParName(13,"omega_2cbo");
     fit_func23->SetParName(14,"phi_2cbo");
     fit_func23->SetParName(15,"A_y");
     fit_func23->SetParName(16,"Tau_y");
     fit_func23->SetParName(17,"omega_y");
     fit_func23->SetParName(18,"phi_y");
     fit_func23->SetParName(19,"A_vw");
     fit_func23->SetParName(20,"Tau_vw");
     fit_func23->SetParName(21,"omega_vw");
     fit_func23->SetParName(22,"phi_vw");
     fit_func23->SetNpx(1000000);


 


       cauto=new TCanvas("auto-correlation", "auto-correlation");
       cauto->Divide(6,4);

        cfit=new TCanvas("23 parameter fit","23 parameter fit");

	 

      		
     //hr_sum=(TH1D*)hcalo->Clone();
       
          for(int i=1; i<=24;i++)
       {

	cout<<"calo number "<<i<<" ns"<<endl; 
        
	 
	//	fit_func3->SetParameters(0.22, 0.0, 2.22);

		fit_func3->SetParameters(0.23, 0.0, 4.01);

		/*	hr1[i]->Rebin(2);
		hr1[i]->Scale(0.5);
		hr2[i]->Rebin(2);
		hr2[i]->Scale(0.5);
                hr3[i]->Rebin(2);
		hr3[i]->Scale(0.5);
                hr4[i]->Rebin(2);
		hr4[i]->Scale(0.5);
		hcalo=construct_rhist_rand(hr1[i],hr2[i],hr3[i],hr4[i]);*/


	     	hcalo=construct_rhist_rand(hcomp_sum1[i],hcomp_sum2[i],hcomp_sum3[i],hcomp_sum4[i], hr1[i], hr2[i], hr3[i], hr4[i]);
  
	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");


        hcalo->GetXaxis()->SetRangeUser(fit_start,fit_stop);



        hcalo->Fit("fprec3","RE","",fit_start,fit_stop);

	rval3[m]=fit_func3->GetParameter(1);
	drval3[m]=fit_func3->GetParError(1);
	
	/*	if(useFR)
	  {
           fit_func3->SetParameter(0, fit_func3->GetParameter(0));
           fit_func3->SetParameter(1, fit_func3->GetParameter(1));
           fit_func3->SetParameter(2, fit_func3->GetParameter(2));
	    
           hcalo->Fit("fprec3","RUE","",fit_start,fit_stop);
	  }
	*/	
         //hcalo->Draw();
	//	fit_func3->Draw();
	//	hcalo->Draw("same");
	
	
   
   
        fit_func7->SetParameter(0, fit_func3->GetParameter(0));
        fit_func7->SetParameter(1, fit_func3->GetParameter(1));
        fit_func7->SetParameter(2, fit_func3->GetParameter(2));
        fit_func7->SetParameter(3, 0.002);
	fit_func7->SetParameter(4, 240000);
        fit_func7->SetParameter(5, 0.0023287);
        fit_func7->SetParameter(6, 1.9+((6.28*(i-1))/24));

   
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

       rval11[m]=fit_func11->GetParameter(1);
       drval11[m]=fit_func11->GetParError(1);
      

       /*   if(useFR)
	  {
           fit_func11->SetParameter(0, fit_func7->GetParameter(0));
           fit_func11->SetParameter(1, fit_func7->GetParameter(1));
           fit_func11->SetParameter(2, fit_func7->GetParameter(2));
	   fit_func11->SetParameter(3, fit_func7->GetParameter(3));
           fit_func11->SetParameter(4, fit_func7->GetParameter(4));
           fit_func11->SetParameter(5, fit_func7->GetParameter(5));
	   fit_func11->SetParameter(6, fit_func7->GetParameter(6));
           fit_func11->SetParameter(7, fit_func7->GetParameter(7));
           fit_func11->SetParameter(8, fit_func7->GetParameter(8));
	   fit_func11->SetParameter(9, fit_func7->GetParameter(9));
           fit_func11->SetParameter(10, fit_func7->GetParameter(10));
	    
           hcalo->Fit("fprec11","RUE","",fit_start,fit_stop);
	  }
      		
       */  
   
      
       /*   fit_func15->SetParameter(0, fit_func11->GetParameter(0));
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
       fit_func15->SetParameter(11, 0.0005);
       fit_func15->FixParameter(12, fit_func11->GetParameter(4)/2);
       fit_func15->FixParameter(13, fit_func11->GetParameter(5)*2);
       fit_func15->SetParameter(14, 5.1);

   
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
       fit_func19->FixParameter(12, fit_func15->GetParameter(12));
       fit_func19->FixParameter(13, fit_func15->GetParameter(13));
       fit_func19->SetParameter(14, fit_func15->GetParameter(14));
       fit_func19->SetParameter(15, 0.001);
       fit_func19->FixParameter(16, 75360);
       fit_func19->FixParameter(17, 0.01389);
       fit_func19->SetParameter(18, 1.1);

      
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
       fit_func23->FixParameter(12, fit_func19->GetParameter(12));
       fit_func23->FixParameter(13, fit_func19->GetParameter(13));
       fit_func23->SetParameter(14, fit_func19->GetParameter(14));
       fit_func23->SetParameter(15, fit_func19->GetParameter(15));
       fit_func23->FixParameter(16, fit_func19->GetParameter(16));
       fit_func23->FixParameter(17, fit_func19->GetParameter(17));
       fit_func23->SetParameter(18, fit_func19->GetParameter(18));
       fit_func23->SetParameter(19, 0.001);
       fit_func23->FixParameter(20, 5500);
       fit_func23->FixParameter(21, 0.0140986);
       fit_func23->SetParameter(22, 3.7+((6.28*(i-1))/24));

       
      
       hcalo->Fit("fprec23","RE","",fit_start,fit_stop);

       rval23[m]=fit_func23->GetParameter(1);
       drval23[m]=fit_func23->GetParError(1);
	
       gStyle->SetOptFit(1111);

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

       //  cov2.ResizeTo(hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), -1);
   
  //   cov.Print();
       //  cov2.SetTol(1.e-23);


       // cov2=ratio_cov_matrix(hcalo);
       //cout<<"cov2[1000][1000] "<<cov2[1000][1000]<<endl;

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
       fit_func23->FixParameter(12, fit_func23->GetParameter(12));
       fit_func23->FixParameter(13, fit_func23->GetParameter(13));
       fit_func23->SetParameter(14, fit_func23->GetParameter(14));
       fit_func23->SetParameter(15, fit_func23->GetParameter(15));
       fit_func23->FixParameter(16, fit_func23->GetParameter(16));
       fit_func23->FixParameter(17, fit_func23->GetParameter(17));
       fit_func23->SetParameter(18, fit_func23->GetParameter(18));
       fit_func23->SetParameter(19, fit_func23->GetParameter(19));
       fit_func23->FixParameter(20, fit_func23->GetParameter(20));
       fit_func23->FixParameter(21, fit_func23->GetParameter(21));
       fit_func23->SetParameter(22, fit_func23->GetParameter(22));


  
       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hcalo);
       gFitter->SetFCN(chi2);

      
       cfit->cd(1);	
       hcalo->Fit("fprec11","RUE","",fit_start,fit_stop);
       countfcn=0;
       gStyle->SetOptFit(1111);

       gettimeofday(&t_end, NULL);
       printf("QRatio fit duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
     
       //cout<<" setting ch to zero "<<endl;
       //ch=0;
 
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
       A_2cbo[m]=TMath::Abs(fit_func23->GetParameter(11));
       tau_2cbo[m]=fit_func23->GetParameter(12);
       omega_2cbo[m]=fit_func23->GetParameter(13);
       phi_2cbo[m]=fit_func23->GetParameter(14);
       A_y[m]=TMath::Abs(fit_func23->GetParameter(15));
       tau_y[m]=fit_func23->GetParameter(16);
       omega_y[m]=fit_func23->GetParameter(17);
       phi_y[m]=fit_func23->GetParameter(18);
       A_vw[m]=TMath::Abs(fit_func23->GetParameter(19));
       tau_vw[m]=fit_func23->GetParameter(20);
       omega_vw[m]=fit_func23->GetParameter(21);
       phi_vw[m]=fit_func23->GetParameter(22);
       
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
       dA_2cbo[m]=fit_func23->GetParError(11);
       dtau_2cbo[m]=fit_func23->GetParError(12);
       domega_2cbo[m]=fit_func23->GetParError(13);
       dphi_2cbo[m]=fit_func23->GetParError(14);
       dA_y[m]=fit_func23->GetParError(15);
       dtau_y[m]=fit_func23->GetParError(16);
       domega_y[m]=fit_func23->GetParError(17);
       dphi_y[m]=fit_func23->GetParError(18);
       dA_vw[m]=fit_func23->GetParError(19);
       dtau_vw[m]=fit_func23->GetParError(20);
       domega_vw[m]=fit_func23->GetParError(21);
       dphi_vw[m]=fit_func23->GetParError(22);
       
       // chisq[m]=fit_func23->GetChisquare()/fit_func23->GetNDF();


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
       */
            n[m]=i;
       m=m+1;
      		
      h_res[i]= new TH1D("residual histogram ", "h_res", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalo->GetBinContent(ibin)- fit_func11->Eval( hcalo->GetBinCenter(ibin) ) );
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

       
 
   }
	  	   
  c2 = new TCanvas("c2","run 3btom calo scan residuals");
  c2->Divide(6,4);
  c3 = new TCanvas("c3","run 3btom calo scan fft");
  c3->Divide(6,4);
   
 for(int i=1;i<=24;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     //  hfft[i]->GetXaxis()->SetRangeUser(2.0,2.5);
     //hfft[i]->GetYaxis()->SetRangeUser(-4,100);
     hfft[i]->Draw();
   }
 
    c1=new TCanvas("c1","run 3btom blind R vs calorimeter");
    c1->Divide(2,2);
    c1->cd(1);
    gr1=new TGraphErrors(m,n,chisq,0,0);
    gr1->SetTitle("chisq vs caloriemeter");
    gr1->GetXaxis()->SetTitle("calo index");
    gr1->GetYaxis()->SetTitle("reduced chi-sq");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kBlue);
    gr1->Draw("ap");
    

     c1->cd(2);
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle("A vs caloriemeter");
    gr2->GetXaxis()->SetTitle("calo index");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kBlue);
    gr2->Draw("ap");
  
     c1->cd(3);
    gr3=new TGraphErrors(m,n,blindR,0,dblindR);
    gr3->SetTitle("blindR vs caloriemeter");
    gr3->GetXaxis()->SetTitle("calo index");
    gr3->GetYaxis()->SetTitle("Blind R [ppm]");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kBlue);
    gr3->Draw("ap");
 
     c1->cd(4);
    gr4=new TGraphErrors(m,n,phi_0,0,dphi_0);
    gr4->SetTitle("phase vs caloriemeter");
    gr4->GetXaxis()->SetTitle("calo index");
    gr4->GetYaxis()->SetTitle("phi_0");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kBlue);
    gr4->Draw("ap");

    c4=new TCanvas("c4","run 3btom cbo vs caloriemeter");
    c4->Divide(2,2);
    c4->cd(1);
    gr8=new TGraphErrors(m,n,A_cbo_N,0,dA_cbo_N);
    gr8->SetTitle("A_cbo_N vs caloriemeter");
    gr8->GetXaxis()->SetTitle("calo index");
    gr8->GetYaxis()->SetTitle("A_cbo_N");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kBlue);
    gr8->Draw("ap");
  
     c4->cd(2);
    gr5=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr5->SetTitle("tau_cbo_ vs caloriemeter");
    gr5->GetXaxis()->SetTitle("calo index");
    gr5->GetYaxis()->SetTitle("tau_cbo");
    gr5->SetMarkerStyle(20);
    gr5->SetLineColor(kBlue);
    gr5->Draw("ap");
 
     c4->cd(4);
    gr6=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr6->SetTitle("omega_cbo vs caloriemeter");
    gr6->GetXaxis()->SetTitle("calo index");
    gr6->GetYaxis()->SetTitle("omega_cbo");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kBlue);
    gr6->Draw("ap");
 
     c4->cd(3);
    gr7=new TGraphErrors(m,n,phi_cbo_N,0,dphi_cbo_N);
    gr7->SetTitle("phi_cbo_N vs caloriemeter");
    gr7->GetXaxis()->SetTitle("calo index");
    gr7->GetYaxis()->SetTitle("phi_cbo_N");
    gr7->SetMarkerStyle(20);
    gr7->SetLineColor(kBlue);
    gr7->Draw("ap");
 
    c6=new TCanvas("c6","run3 btom A,phi,cbo vs caloriemeter");
    c6->Divide(2,2);
    c6->cd(1);
    gr13=new TGraphErrors(m,n,A_cbo_A,0,dA_cbo_A);
    gr13->SetTitle("A_cbo_A vs caloriemeter");
    gr13->GetXaxis()->SetTitle("calo index");
    gr13->GetYaxis()->SetTitle("A_cbo_A");
    gr13->SetMarkerStyle(20);
    gr13->SetLineColor(kBlue);
    gr13->Draw("ap");

    c6->cd(3);
    gr14=new TGraphErrors(m,n,phi_cbo_A,0,dphi_cbo_A);
    gr14->SetTitle("phi_cbo_A vs caloriemeter");
    gr14->GetXaxis()->SetTitle("calo index");
    gr14->GetYaxis()->SetTitle("phi_cbo_A");
    gr14->SetMarkerStyle(20);
    gr14->SetLineColor(kBlue);
    gr14->Draw("ap");

    c6->cd(2);
    gr15=new TGraphErrors(m,n,A_cbo_phi,0,dA_cbo_phi);
    gr15->SetTitle("A_cbo_phi vs caloriemeter");
    gr15->GetXaxis()->SetTitle("calo index");
    gr15->GetYaxis()->SetTitle("A_cbo_phi");
    gr15->SetMarkerStyle(20);
    gr15->SetLineColor(kBlue);
    gr15->Draw("ap");

    c6->cd(4);
    gr16=new TGraphErrors(m,n,phi_cbo_phi,0,dphi_cbo_phi);
    gr16->SetTitle("phi_cbo_phi vs caloriemeter");
    gr16->GetXaxis()->SetTitle("calo index");
    gr16->GetYaxis()->SetTitle("phi_cbo_phi");
    gr16->SetMarkerStyle(20);
    gr16->SetLineColor(kBlue);
    gr16->Draw("ap");
 

     c7=new TCanvas("c7","run 3btom 2cbo vs caloriemeter");
     c7->Divide(2,2);
    c7->cd(1);
    gr17=new TGraphErrors(m,n,A_2cbo,0,dA_2cbo);
    gr17->SetTitle("A_2cbo vs caloriemeter");
    gr17->GetXaxis()->SetTitle("calo index");
    gr17->GetYaxis()->SetTitle("A_2cbo");
    gr17->SetMarkerStyle(20);
    gr17->SetLineColor(kBlue);
    gr17->Draw("ap");
 

     c7->cd(2);
    gr18=new TGraphErrors(m,n,tau_2cbo,0,dtau_2cbo);
    gr18->SetTitle("tau_2cbo_ vs caloriemeter");
    gr18->GetXaxis()->SetTitle("calo index");
    gr18->GetYaxis()->SetTitle("tau_2cbo");
    gr18->SetMarkerStyle(20);
    gr18->SetLineColor(kBlue);
    gr18->Draw("ap");
 

     c7->cd(4);
    gr19=new TGraphErrors(m,n,omega_2cbo,0,domega_2cbo);
    gr19->SetTitle("omega_2cbo vs caloriemeter");
    gr19->GetXaxis()->SetTitle("calo index");
    gr19->GetYaxis()->SetTitle("omega_2cbo");
    gr19->SetMarkerStyle(20);
    gr19->SetLineColor(kBlue);
    gr19->Draw("ap");
 
     c7->cd(3);
    gr20=new TGraphErrors(m,n,phi_2cbo,0,dphi_2cbo);
    gr20->SetTitle("phi_2cbo vs caloriemeter");
    gr20->GetXaxis()->SetTitle("calo index");
    gr20->GetYaxis()->SetTitle("phi_2cbo");
    gr20->SetMarkerStyle(20);
    gr20->SetLineColor(kBlue);
    gr20->Draw("ap");
 

    c8=new TCanvas("c8","run 3btom vw vs caloriemeter");
    c8->Divide(2,2);
    c8->cd(1);
    gr21=new TGraphErrors(m,n,A_vw,0,dA_vw);
    gr21->SetTitle("A_vw vs calo index");
    gr21->GetXaxis()->SetTitle("calo index");
    gr21->GetYaxis()->SetTitle("A_vw");
    gr21->SetMarkerStyle(20);
    gr21->SetLineColor(kBlue);
    gr21->Draw("ap");
  

    c8->cd(2);
    gr22=new TGraphErrors(m,n,tau_vw,0,dtau_vw);
    gr22->SetTitle("tau_vw vs calo index");
    gr22->GetXaxis()->SetTitle("calo index");
    gr22->GetYaxis()->SetTitle("tau_vw");
    gr22->SetMarkerStyle(20);
    gr22->SetLineColor(kBlue);
    gr22->Draw("ap");
  
    c8->cd(4);
    gr23=new TGraphErrors(m,n,omega_vw,0,domega_vw);
    gr23->SetTitle("omega_vw vs calo index");
    gr23->GetXaxis()->SetTitle("calo index");
    gr23->GetYaxis()->SetTitle("omega_vw");
    gr23->SetMarkerStyle(20);
    gr23->SetLineColor(kBlue);
    gr23->Draw("ap");
  

    c8->cd(3);
    gr24=new TGraphErrors(m,n,phi_vw,0,dphi_vw);
    gr24->SetTitle("phi_vw vs calo index");
    gr24->GetXaxis()->SetTitle("calo index");
    gr24->GetYaxis()->SetTitle("phi_vw");
    gr24->SetMarkerStyle(20);
    gr24->SetLineColor(kBlue);
    gr24->Draw("ap");
 
    
    c9=new TCanvas("c9","run 3btom y vs caloriemeter");
    c9->Divide(2,2);
    c9->cd(1);
    gr25=new TGraphErrors(m,n,A_y,0,dA_y);
    gr25->SetTitle("A_y vs calo index");
    gr25->GetXaxis()->SetTitle("calo index");
    gr25->GetYaxis()->SetTitle("A_y");
    gr25->SetMarkerStyle(20);
    gr25->SetLineColor(kBlue);
    gr25->Draw("ap");
 

    c9->cd(2);
    gr26=new TGraphErrors(m,n,tau_y,0,dtau_y);
    gr26->SetTitle("tau_y vs calo index");
    gr26->GetXaxis()->SetTitle("calo index");
    gr26->GetYaxis()->SetTitle("tau_y");
    gr26->SetMarkerStyle(20);
    gr26->SetLineColor(kBlue);
    gr26->Draw("ap");
 

    c9->cd(4);
    gr27=new TGraphErrors(m,n,omega_y,0,domega_y);
    gr27->SetTitle("omega_y vs calo index");
    gr27->GetXaxis()->SetTitle("calo index");
    gr27->GetYaxis()->SetTitle("omega_y");
    gr27->SetMarkerStyle(20);
    gr27->SetLineColor(kBlue);
    gr27->Draw("ap");
 
    c9->cd(3);
    gr28=new TGraphErrors(m,n,phi_y,0,dphi_y);
    gr28->SetTitle("phi_y vs calo index");
    gr28->GetXaxis()->SetTitle("calo index");
    gr28->GetYaxis()->SetTitle("phi_y");
    gr28->SetMarkerStyle(20);
    gr28->SetLineColor(kBlue);
    gr28->Draw("ap");

    c10=new TCanvas("c10","run 3btom R vs caloriemeter");
    gr29=new TGraphErrors(m,n,rval3,0,drval3);
    gr29->SetTitle("R(3 par) vs calo index");
    gr29->GetXaxis()->SetTitle("calo index");
    gr29->GetYaxis()->SetTitle("R");
    gr29->SetMarkerStyle(20);
    gr29->SetLineColor(kBlue);
    gr29->Draw("ap");

    c11=new TCanvas("c11","run 3btom R vs caloriemeter");
    gr30=new TGraphErrors(m,n,rval11,0,drval11);
    gr30->SetTitle("R(11 par) vs calo index");
    gr30->GetXaxis()->SetTitle("calo index");
    gr30->GetYaxis()->SetTitle("R");
    gr30->SetMarkerStyle(20);
    gr30->SetLineColor(kBlue);
    gr30->Draw("ap");


    c12=new TCanvas("c12","run 3btom R vs caloriemeter");
    gr31=new TGraphErrors(m,n,rval23,0,drval23);
    gr31->SetTitle("R(23 par) vs calo index");
    gr31->GetXaxis()->SetTitle("calo index");
    gr31->GetYaxis()->SetTitle("R");
    gr31->SetMarkerStyle(20);
    gr31->SetLineColor(kBlue);
    gr31->Draw("ap");


    

    
}
