#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include <sys/time.h>
#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );



float toddiff(struct timeval *tod2, struct timeval *tod1)
{
  float fdt, fmudt;
  long long t1, t2, mut1, mut2;
  long long dt, mudt;
  t1 = tod1->tv_sec;
  mut1 = tod1->tv_usec;
  t2 = tod2->tv_sec;
  mut2 = tod2->tv_usec;
  dt = (t2 - t1);
  mudt = (mut2 - mut1);
  fdt = (float)dt;
  fmudt = (float)mudt;
  return 1.0e6 * fdt + fmudt;
}


TVirtualFitter *gFitter;
TCanvas *cfit3,*cfit7,*cfit11,*cfit15,*cfit19,*cfit23,*cfit24,*cauto;
TH1D *h_res3,*h_res7,*h_res11,*h_res15,*h_res19,*h_res20,*h_res23, *hcomp;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24, *hFFTsq, *hAuto;
TH1F *hlm;
Double_t fit_start, fit_stop;
Int_t countfcn=0;
Int_t m;

Double_t fprec3(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

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

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

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

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

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

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

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

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));

    double Nvbob=(1+ asym_vbo*exp(-(time - T_a/2)/tau_vbo)*cos(omega_vbo*(time - T_a/2) + phi_vbo));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

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

    double asym_pr= par[15];

    double tau_pr = par[16];

    double omega_pr = par[17];

    double phi_pr = par[18];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));

    double Nvbob=(1+ asym_vbo*exp(-(time - T_a/2)/tau_vbo)*cos(omega_vbo*(time - T_a/2) + phi_vbo));

    double Npr=(1+ asym_pr*exp(-time/tau_pr)*cos(omega_pr*time + phi_pr));

    double Nprf=(1+ asym_pr*exp(-(time + T_a/2)/tau_pr)*cos(omega_pr*(time + T_a/2) + phi_pr));

    double Nprb=(1+ asym_pr*exp(-(time - T_a/2)/tau_pr)*cos(omega_pr*(time - T_a/2) + phi_pr));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw*Nvbo*Npr - ff*Ncbof*Nvwf*Nvbof*Nprf - fb*Ncbob*Nvwb*Nvbob*Nprb)/(2*f*Ncbo*Nvw*Nvbo*Npr + ff*Ncbof*Nvwf*Nvbof*Nprf + fb*Ncbob*Nvwb*Nvbob*Nprb);

}



Double_t fprec2(Double_t *x, Double_t *par)
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

    double lost_muon_amp = par[15];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double Nlm= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double Nlmf= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time + T_a/2) * 1.0e-3)));

    double Nlmb= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time - T_a/2) * 1.0e-3)));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));
    
     return (2*f*Ncbo*Nvw*Nlm - ff*Ncbof*Nvwf*Nlmf - fb*Ncbob*Nvwb*Nlmb)/(2*f*Ncbo*Nvw*Nlm + ff*Ncbof*Nvwf*Nlmf + fb*Ncbob*Nvwb*Nlmb);
    //  return 2*(1+asym*cos(omega_a*time + phi)) - (1+asym*cos(omega_a*(time+ (T_a+20)/2) + phi)) - (1+asym*cos(omega_a*(time - (T_a+20)/2) + phi));
    // return (2*f - ff - fb)/(2*f + ff + fb);
      // /(2*f + ff + fb);
}

 void automated_ratio_fit()

 {
 

 // cout<<r3->Rndm()<<" "<<endl;
   // T_a = T_a + 20;
   //gStyle->SetOptFit(1111);
   
 cfit3=new TCanvas("cfit","3 parameter ratio fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");



     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A", "Blind R", "Phi");
     fit_func3->SetNpx(1000000);


     fit_start=30000;
     fit_stop=300000;

     fit_func3->SetParameters(0.2, 0.0, 3.14/2);

     cfit3->Divide(1,3);
     cfit3->cd(1);


	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");
        hcalo->GetXaxis()->SetRangeUser(fit_start,fit_stop);


    printf("start Q-method analyzer\n");
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
     hcalo->Fit("fprec3","RE","",fit_start,fit_stop);
      gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));

      
		
      h_res3= new TH1D("residual histogram 3", "h_res3", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit3->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func3->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res3->SetBinContent(ibin, (res)  );
     
      }

	h_res3->GetYaxis()->SetTitle("Energy [MeV]");
	h_res3->GetXaxis()->SetTitle("time [ns]");
	h_res3->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res3->GetXaxis()->SetLabelSize(0.04);
	h_res3->GetYaxis()->SetLabelSize(0.04);
        h_res3->Draw();

   cfit3->cd(3);
   hfft3 = h_res3->FFT(hfft3, "MAG");
   hfft3->SetLineColor(kBlack);
   hfft3->SetBins(hcalo->GetNbinsX(),0,1000/hcalo->GetBinWidth(1));
   hfft3->GetXaxis()->SetRangeUser(0,1000/(2*hcalo->GetBinWidth(1)));
   hfft3->GetXaxis()->SetTitle("Freq [MHz]");
   hfft3->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft3->GetXaxis()->SetLabelSize(0.04);
   hfft3->GetYaxis()->SetLabelSize(0.04);
   hfft3->Draw();
   
   
     fit_func7= new TF1("fprec7", fprec7,  30000,  309000, 7);
     fit_func7->SetParNames("A", "Blind R", "Phi");
     fit_func7->SetParName(3,"A_cbo");
     fit_func7->SetParName(4,"Tau_cbo");
     fit_func7->SetParName(5,"omega_cbo");
     fit_func7->SetParName(6,"phi_cbo");
     fit_func7->SetNpx(1000000);

     fit_func7->SetParameter(0, fit_func3->GetParameter(0));
     fit_func7->SetParameter(1, fit_func3->GetParameter(1));
     fit_func7->SetParameter(2, fit_func3->GetParameter(2));
     fit_func7->SetParameter(3, 0.0025);
     fit_func7->SetParameter(4, 235000);
     fit_func7->SetParameter(5, 0.00234023);
     fit_func7->SetParameter(6, 0.48);

     cfit7=new TCanvas("cfit 7","7 parameter ratio fit"); 
     cfit7->Divide(1,3);
     cfit7->cd(1);

     hcalo->Fit("fprec7","RE","",fit_start,fit_stop);
      		
      h_res7= new TH1D("residual histogram 7", "h_res7", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit7->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func7->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res7->SetBinContent(ibin, (res)  );
     
      }

	h_res7->GetYaxis()->SetTitle("ADC counts");
	h_res7->GetXaxis()->SetTitle("time [ns]");	
        h_res7->Draw();

   cfit7->cd(3);
   hfft7 = h_res7->FFT(hfft7, "MAG");
   hfft7->SetLineColor(kBlack);
   hfft7->SetBins(hcalo->GetNbinsX(),0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1))));
   hfft7->GetXaxis()->SetRangeUser(0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)))/2);
   hfft7->GetXaxis()->SetTitle("Freq [GHz]");
   hfft7->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft7->Draw();
   
   
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

     fit_func11->SetParameter(0, fit_func7->GetParameter(0));
     fit_func11->SetParameter(1, fit_func7->GetParameter(1));
     fit_func11->SetParameter(2, fit_func7->GetParameter(2));
     fit_func11->SetParameter(3, fit_func7->GetParameter(3));
     fit_func11->SetParameter(4, fit_func7->GetParameter(4));
     fit_func11->SetParameter(5, fit_func7->GetParameter(5));
     fit_func11->SetParameter(6, fit_func7->GetParameter(6));
     /*    fit_func11->SetParameter(7, 0.0008);
     fit_func11->SetParameter(8, 0.97);
     fit_func11->SetParameter(9, 0.0002);
     fit_func11->SetParameter(10, -2.6);
     */
     fit_func11->SetParameter(7, 0.0001);
     fit_func11->SetParameter(8, 0.1);
     fit_func11->SetParameter(9, 0.0002);
     fit_func11->SetParameter(10, 0.1);


     cfit11=new TCanvas("cfit 11","11 parameter ratio fit"); 
     cfit11->Divide(1,3);
     cfit11->cd(1);

     hcalo->Fit("fprec11","RE","",fit_start,fit_stop);
      		
      h_res11= new TH1D("residual histogram 11", "h_res11", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit11->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func11->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res11->SetBinContent(ibin, (res)  );
     
      }

	h_res11->GetYaxis()->SetTitle("ADC counts");
	h_res11->GetXaxis()->SetTitle("time [ns]");	
        h_res11->Draw();

   cfit11->cd(3);
   hfft11 = h_res11->FFT(hfft11, "MAG");
   hfft11->SetLineColor(kBlack);
   hfft11->SetBins(hcalo->GetNbinsX(),0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1))));
   hfft11->GetXaxis()->SetRangeUser(0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)))/2);
   hfft11->GetXaxis()->SetTitle("Freq [GHz]");
   hfft11->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft11->Draw();
   
   
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
     fit_func15->SetParameter(11, 0.009);
     fit_func15->SetParameter(12, 16200);
     fit_func15->SetParameter(13, 0.01419594);
     fit_func15->SetParameter(14, 9.6);

     cfit15=new TCanvas("cfit 15","15 parameter ratio fit"); 
     cfit15->Divide(1,3);
     cfit15->cd(1);

     hcalo->Fit("fprec15","RE","",fit_start,fit_stop);
      		
      h_res15= new TH1D("residual histogram 15", "h_res15", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit15->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func15->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res15->SetBinContent(ibin, (res)  );
     
      }

	h_res15->GetYaxis()->SetTitle("ADC counts");
	h_res15->GetXaxis()->SetTitle("time [ns]");	
        h_res15->Draw();

   cfit15->cd(3);
   hfft15 = h_res15->FFT(hfft15, "MAG");
   hfft15->SetLineColor(kBlack);
   hfft15->SetBins(hcalo->GetNbinsX(),0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1))));
   hfft15->GetXaxis()->SetRangeUser(0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)))/2);
   hfft15->GetXaxis()->SetTitle("Freq [GHz]");
   hfft15->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft15->Draw();
   
   
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
     fit_func19->SetParName(15,"A_vbo");
     fit_func19->SetParName(16,"Tau_vbo");
     fit_func19->SetParName(17,"omega_vbo");
     fit_func19->SetParName(18,"phi_vbo");

     fit_func19->SetNpx(1000000);

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
     fit_func19->SetParameter(15, 0.000545);
     fit_func19->SetParameter(16, 85003.6);
     fit_func19->SetParameter(17, 0.01392276);
     fit_func19->SetParameter(18, 0.5907);

     cfit19=new TCanvas("cfit 19","19 parameter ratio fit"); 
     cfit19->Divide(1,3);
     cfit19->cd(1);
     gStyle->SetOptFit(1111);
     //gPrintViaErrorHandler = kTRUE;
     //gErrorIgnoreLevel = kWarning;
     hcalo->Fit("fprec19","RE","",fit_start,fit_stop);
      		
     h_res19= new TH1D("residual histogram 19", "h_res19", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit19->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func19->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res19->SetBinContent(ibin, (res)  );
     
      }

	h_res19->GetYaxis()->SetTitle("Energy [MeV]");
	h_res19->GetXaxis()->SetTitle("time [ns]");
	h_res19->GetXaxis()->SetRangeUser(fit_start,fit_stop);
	h_res19->GetXaxis()->SetLabelSize(0.04);
	h_res19->GetYaxis()->SetLabelSize(0.04);
        h_res19->Draw();

   cfit19->cd(3);
   hfft19 = h_res19->FFT(hfft19, "MAG");
   hfft19->SetLineColor(kBlack);
   hfft19->SetBins(hcalo->GetNbinsX(),0,1000/hcalo->GetBinWidth(1));
   hfft19->GetXaxis()->SetRangeUser(0,1000/(2*hcalo->GetBinWidth(1)));
   hfft19->GetXaxis()->SetTitle("Freq [MHz]");
   hfft19->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft19->GetXaxis()->SetLabelSize(0.04);
   hfft19->GetYaxis()->SetLabelSize(0.04);
   hfft19->Draw();

   cauto=new TCanvas("auto-correlation", "auto-correlation");
   hFFTsq = (TH1D*)hfft19->Clone();
   hFFTsq->Reset();
   for (int ib = 1; ib <= hFFTsq->GetNbinsX(); ib++)
      {
       hFFTsq->SetBinContent( ib, hfft19->GetBinContent(ib)*hfft19->GetBinContent(ib));
      }


       hAuto = hFFTsq->FFT(hAuto, "MAG");
       hAuto->SetTitle("auto-correlation");
       hAuto->GetXaxis()->SetLabelSize(0.04);
       hAuto->GetXaxis()->SetTitle("bins");
       hAuto->GetXaxis()->SetRangeUser(0,h_res19->GetNbinsX()/2);
       hAuto->SetLineColor(kRed);
       hAuto->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hAuto->Draw("HIST");
 

    return;
     fit_func23= new TF1("fprec23", fprec23,  30000,  309000, 23);
     fit_func23->SetParNames("A", "Blind R", "Phi");
     fit_func23->SetParName(3,"A_cbo");
     fit_func23->SetParName(4,"Tau_cbo");
     fit_func23->SetParName(5,"omega_cbo");
     fit_func23->SetParName(6,"phi_cbo");
     fit_func23->SetParName(7,"A2_cbo");
     fit_func23->SetParName(8,"phi2_cbo");
     fit_func23->SetParName(9,"A3_cbo");
     fit_func23->SetParName(10,"phi3_cbo");
     fit_func23->SetParName(11,"A_vw");
     fit_func23->SetParName(12,"Tau_vw");
     fit_func23->SetParName(13,"omega_vw");
     fit_func23->SetParName(14,"phi_vw");
     fit_func23->SetParName(15,"A_vbo");
     fit_func23->SetParName(16,"Tau_vbo");
     fit_func23->SetParName(17,"omega_vbo");
     fit_func23->SetParName(18,"phi_vbo");
     fit_func23->SetParName(19,"A_pr");
     fit_func23->SetParName(20,"Tau_pr");
     fit_func23->SetParName(21,"omega_pr");
     fit_func23->SetParName(22,"phi_pr");
     

     fit_func23->SetNpx(1000000);

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
     fit_func23->SetParameter(15, fit_func19->GetParameter(15));
     fit_func23->SetParameter(16, fit_func19->GetParameter(16));
     fit_func23->SetParameter(17, fit_func19->GetParameter(17));
     fit_func23->SetParameter(18, fit_func19->GetParameter(18));
     fit_func23->SetParameter(19,0.01);
     fit_func23->SetParameter(20,16000);
     fit_func23->SetParameter(21,0.01085);
     fit_func23->SetParameter(22,1.5);

     cfit23=new TCanvas("cfit 23","23 parameter ratio fit"); 
     cfit23->Divide(1,3);
     cfit23->cd(1);

     hcalo->Fit("fprec23","RE","",fit_start,fit_stop);
      		
     h_res23= new TH1D("residual histogram 23", "h_res23", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit23->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func23->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res23->SetBinContent(ibin, (res)  );
     
      }

	h_res23->GetYaxis()->SetTitle("ADC counts");
	h_res23->GetXaxis()->SetTitle("time [ns]");	
        h_res23->Draw();

   cfit23->cd(3);
   hfft23 = h_res23->FFT(hfft23, "MAG");
   hfft23->SetLineColor(kBlack);
   hfft23->SetBins(hcalo->GetNbinsX(),0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1))));
   hfft23->GetXaxis()->SetRangeUser(0,hcalo->GetNbinsX()/((hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)))/2);
   hfft23->GetXaxis()->SetTitle("Freq [GHz]");
   hfft23->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft23->Draw();


 
   
 }
