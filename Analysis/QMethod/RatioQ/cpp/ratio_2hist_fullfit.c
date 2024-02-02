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



//char root_file_name[128] = "highstat_wigsim_output_statfluc_binintegral_e16.root";
//char root_file_name[128] = "run2_thresh300_fbfDQC.root";
char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
TVirtualFitter *gFitter;
TFile *_file[25];
TFile *foutput;
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res, *hr_sum, *hcov_14_p, *hcov_14_m, *hcov_28_p, *hcov_28_m;
TCanvas *c,*c1, *c2,*c3, *cnew, *cauto;
TGraphErrors *gr1;
TH1 *hfft, *hFFTsq, *hAuto;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=29926.250;
Double_t  fit_stop=299926.25;
//Double_t fit_start=181300;
//Double_t fit_stop=196000;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
Int_t m=0;
Double_t blindR[25],dblindR[25],n[25];
Int_t countfcn=0;


TH2D *hcov;
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
TMatrixD cov(mdim,mdim),cov2(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim));

TMatrixD ratio_cov_matrix(TH1D *hr_sum, TH1D *h_sum){

  /* for (Int_t i = 0; i < mdim*mdim; i++) {
    const Int_t ir = i/mdim;
    const Int_t ic = i%mdim;
    data[i] = 0.0;
    if ( ir == ic ) data[i] = 1;
    if ( ir == ic-1 || ir == ic+1 ) data[i] = 7;
  }
  h.SetMatrixArray(data.GetArray());

  h.Print();

  Double_t det1;
  h.Invert(&det1);

  h.Print();
  cout<<h[1][2]<<" "<<endl;*/

  double covar, dprod1, dprod2, dprod3, dprod4, dprod5, dI11, dI12, dI13, dI21, dI22, dI23;

  double u,v,w,x,z,s,t,y;

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


	//	cout<<i<<" "<<ir<<" "<<ic<<endl;
	
        w=h_sum->GetBinContent(ir);
        x=h_sum->GetBinContent(ir+nbinshift);
        y=h_sum->GetBinContent(ir-nbinshift);
        z=h_sum->GetBinContent(ir+(2*nbinshift));

        
	
        data[i] = 0.0;
 
        if(ir==ic)
	  {
	    //  data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir))/(hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	    //	    if(data[i]==0)
	    //  {
	    //	cout<<"diag ele is zero at "<<ir<<" "<<ic<<endl;
	    //	}
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }
       
       	if(ir==ic+nbinshift)
	  {
	    //  cout<<i<<" "<<ir<<" "<<ic<<endl;
	    
	    dprod1=(-((4*x*(x - z))/((w + x)*(w + x)*(w + x)*(x + z))));

	    dprod2=(-((4*(w*w*w*z - x*x*x*z + w*w*z*(3*x + 2*z) + w*(-x*x*x + 3*x*z*z + z*z*z)))/((w + x)*(w + x)*(w + x)*(x + z)*(x + z)*(x + z))));
	      
	    dprod3=((4*(w - x)*x)/((w + x)*(x + z)*(x + z)*(x + z)));

	    dI11=(-(4*x)/((w + x)*(w + x)*(w + x)));

	    dI12=(4*w)/((w + x)*(w + x)*(w + x));

	    dI21=(-(4*z)/((x + z)*(x + z)*(x + z)));

	    dI22=(4*x)/((x + z)*(x + z)*(x + z));
	      
	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir+nbinshift) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift)) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift))*(hr_sum->GetBinContent(ir+nbinshift) + 0.5*dI21*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI22*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift)));

	    //  data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    data[i]=covar;
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

	  if(ir==ic-nbinshift)
	  {

	    dprod1=-((4*(x + y)*(-w*w*w + 3*w*x*y + x*y*(x + y)))/((w + x)*(w + x)*(w + x)*(w + y)*(w + y)*(w + y)));

	    dprod2=-((4*w*(w - y))/((w + x)*(w + x)*(w + x)*(w + y)));
	      
	    dprod3=-((4*w*(w - x))/((w + x)*(w + y)*(w + y)*(w + y)));
	      
	    dI11=-(4*x)/((w + x)*(w + x)*(w + x));
	      
	    dI12=(4*w)/((w + x)*(w + x)*(w + x));

	    dI21=(4*y)/((w + y)*(w + y)*(w + y));
	      
	    dI22=-(4*w)/((w + y)*(w + y)*(w + y));

	      
	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir-nbinshift) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift)) * (hr_sum->GetBinContent(ir-nbinshift) + 0.5*dI21*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI22*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift));

	    // data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	    data[i]=covar;
	  }

		

      }
  cout<<"flgs "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());

  cov.ResizeTo(hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
   
  //  cov.Print();
    cov.SetTol(1.e-23);
    Double_t det1;
    cov.Invert(&det1);
    //  cov.Print();
    /*     hcov = new TH2D("hcov", "cov matrix hist", 2100, 0.0, 2100.0, 2100, 0.0, 2100.0);
  for(int irow=44; irow<=2060; irow++)
    {
      for(int icol=44; icol<=2060; icol++)
	{
	  hcov->SetBinContent(irow,icol,cov(irow,icol));
	}
    }
    */
    cout<<mdim<<endl;
    //cout<<"Matrix inverted "<<i0<<endl;
  return cov;

}


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

    return (f - ff)/(f + ff);
    

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

    return (f*Ncbo - ff*Ncbof)/(f*Ncbo + ff*Ncbof);

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

    return (f*Ncbo - ff*Ncbof)/(f*Ncbo + ff*Ncbof);

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

    return (f*Ncbo*Nvw - ff*Ncbof*Nvwf)/(f*Ncbo*Nvw + ff*Ncbof*Nvwf);

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

    return (f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof)/(f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof);

}



TH1D* construct_rhist_copy(TH1D *h_sum)
{
 
   h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
   h_sum->Scale(0.125);
   

   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

     nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
      /*   if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   // nbinshift=20;
   cout<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=h_sum->GetNbinsX(); ibin++)             
     {
       hp->SetBinContent( ibin, h_sum->GetBinContent( ibin + nbinshift));
       hp->SetBinError(ibin, h_sum->GetBinError( ibin + nbinshift));
     }

   for(int ibin=nbinshift; ibin<=h_sum->GetNbinsX(); ibin++)
     {
       hm->SetBinContent( ibin, h_sum->GetBinContent( ibin - nbinshift));
       hm->SetBinError( ibin, h_sum->GetBinError(ibin - nbinshift));
     }
 

   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*h_sum->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   //double blife=exp(-(T_a)/(2*lifetime));
   double deno=1+flife;

   h1->Scale(1/deno);
   //h2->Scale(1/deno);
   hp->Scale(flife/deno);
   // hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)-hp->GetBinContent(ibin))/(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo histogram sum numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

      h_num->Add(h1,1);
      // h_num->Add(h2,1);
      h_num->Add(hp,-1);
      // h_num->Add(hm,-1);

      h_deno->Add(h1,1);
      //h_deno->Add(h2,1);
      h_deno->Add(hp,1);
      //h_deno->Add(hm,1);
     
      h_num->Divide(h_deno);

     }
    //Assign error bars
   
 if(useerrorbars)
     {
       
       long double t1,t2;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
       	t1=2*(h1->GetBinError(ibin))*(hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+hp->GetBinContent(ibin))*(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));

	t2=-2*(hp->GetBinError(ibin))*(hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+hp->GetBinContent(ibin))*(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));
	


 	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt((t1*t1) + (t2*t2) ));
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( (t1*t1) + (t2*t2) ) );
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

 //  hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //  hcalo->Scale(0.125);

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

   if(hr_sum->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
   
   for(Int_t i=hr_sum->FindBin(fit_start); i<=hr_sum->FindBin(fit_stop); i++)
    {
      for(Int_t j=hr_sum->FindBin(fit_start); j<=hr_sum->FindBin(fit_stop); j++)
	{
	  //  for(m=0; m<=105; m++)
	  // {
	  //  if(i==j || i==j+15 || i==j-15 || i==j+30 || i==j-30 || i==j+45 || i==j-45 || i==j+60 || i==j-60 || i==j+75 || i==j-75 || i==j+90 || i==j-90 || i==j+105 || i== j-105 || i==j+120 || i==j-120 || i==j+135 || i==j-135 || i==j+150 || i==j+165 || i==j-165 || i==j+180 || i==j-180)
	   if(cov2[i][j]!=0)
	   {
	     ch=ch+((hr_sum->GetBinContent(i))-(fuser->Eval(hr_sum->GetBinCenter(i))))*cov2[i][j]*((hr_sum->GetBinContent(j))-(fuser->Eval(hr_sum->GetBinCenter(j))));
	   }
	      //cout<<m<<endl;
	      // ch=ch+1;
	      // cout<<i<<" "<<j<<endl;
	   // }
	     // }
	  // cout<<(hr_sum->GetBinCenter(i+hr_sum->FindBin(fit_start)))<<endl;
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}


void ratio_2hist_fullfit()
{

    _file[1]=TFile::Open(root_file_name);
  //  _file[1]->GetObject("QRatioFillByFillAnalyzerDB",dir[1]);
  // dir[1]->GetObject("qHist_1D_sig_sum_1",qHist_1D[1]);
  // _file[1]->GetObject("hwiggle",qHist_1D[1]);
  // hsum1=(TH1D*)qHist_1D[1]->Clone();

     for(int i=1;i<=24;i++)
     {
      h[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
      }



      h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
      h_sum->Sumw2(kTRUE);
   for(Int_t i=1; i<=24; i++)
    {
     h_sum->Add(h[i],1); 
    }

    //  h_sum=(TH1D*)(_file[1]->FindObjectAny("hwiggle"));
    
   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   for(int i=1;i<=24;i++)
    {
    h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
    }

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

  


    







     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A","R","phi");
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
     fit_func19->SetParName(15,"A_vbo");
     fit_func19->SetParName(16,"Tau_vbo");
     fit_func19->SetParName(17,"omega_vbo");
     fit_func19->SetParName(18,"phi_vbo");
     fit_func19->SetNpx(1000000);


       




     
 

        hr_sum = construct_rhist_copy(h_sum);

	//  	h_sum->Rebin(8);

       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hcalo);
       gFitter->SetFCN(chi2);  

	 
       fit_func3->SetParameters(0.229,0.0, 2.23);

       
	hr_sum->GetYaxis()->SetTitle("ADC counts");
	hr_sum->GetXaxis()->SetTitle("time [ns]");
        hr_sum->GetXaxis()->SetRangeUser(fit_start,fit_stop);

	//  	cov2.ResizeTo(h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), -1);
	 cov2.ResizeTo(h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), -1);
   
  //   cov.Print();
        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr_sum,h_sum);
	//	cov2.Print();

      	hcov = new TH2D("hcov", "cov matrix hist", 2100, 0.0, 2100.0, 2100, 0.0, 2100.0);

	 hcov_14_p = new TH1D("hcov_14_p", "cov matrix hist 14p", 2100, 0.0, 2100.0);
	 hcov_14_m = new TH1D("hcov_14_m", "cov matrix hist 14m", 2100, 0.0, 2100.0);
	 //hcov_28_p = new TH1D("hcov_28_p", "cov matrix hist 28p", 2100, 0.0, 2100.0);
	 //hcov_28_m = new TH1D("hcov_28_m", "cov matrix hist 28m", 2100, 0.0, 2100.0);

	    
	for(int irow=hr_sum->FindBin(fit_start); irow<=hr_sum->FindBin(fit_stop); irow++)
        {
	  for(int icol=hr_sum->FindBin(fit_start); icol<=hr_sum->FindBin(fit_stop); icol++)
	 {
	  hcov->SetBinContent(irow,icol,cov2(irow,icol));
	 }
        }
	    


      for(int irow=hr_sum->FindBin(fit_start); irow<=hr_sum->FindBin(fit_stop); irow++)
      {
       for(int icol=hr_sum->FindBin(fit_start); icol<=hr_sum->FindBin(fit_stop); icol++)
	{
	  if(irow==icol)
	    {
	      hcov_14_p->SetBinContent(irow,hcov->GetBinContent(irow,irow + nbinshift));
	      hcov_14_m->SetBinContent(irow,hcov->GetBinContent(irow,irow - nbinshift));
	      //  hcov_28_p->SetBinContent(irow,hcov->GetBinContent(irow,irow + (2*nbinshift)));
	      //hcov_28_m->SetBinContent(irow,hcov->GetBinContent(irow,irow - (2*nbinshift)));
	    }
	}
    }

   hcov_14_m->SetLineColor(kRed);
   //  hcov_28_p->SetLineColor(kBlack);
   //hcov_28_m->SetLineColor(kMagenta);

   /*   foutput=new TFile("correlation_matrix_nbinshift_14.root","new");
   hcov_14_p->Write();
   hcov_14_m->Write();
  // hcov_28_p->Write();
  // hcov_28_m->Write();
   foutput->Close();
   */
   hcov_14_p->Draw();
   hcov_14_m->Draw("same");
   // hcov_28_p->Draw("same");
   // hcov_28_m->Draw("same");
   
   // return;
   //  hr_sum->Fit("fprec3","RE","",fit_start,fit_stop);
 
      
	/*	fit_func7->SetParameter(0, 0.229);
        fit_func7->SetParameter(1, 0.0);
        fit_func7->SetParameter(2, 2.224);
        fit_func7->SetParameter(3, 0.00256);
        fit_func7->SetParameter(4, 221100);
        fit_func7->SetParameter(5, 0.00234);
        fit_func7->SetParameter(6, 0.464);
	*/
	fit_func7->SetParameter(0, 0.229);
        fit_func7->SetParameter(1, 0.0);
        fit_func7->SetParameter(2, 2.224);
        fit_func7->SetParameter(3, 0.00238);
        fit_func7->SetParameter(4, 256000);
        fit_func7->SetParameter(5, 0.00234);
        fit_func7->SetParameter(6, 0.412);

   
	//   hr_sum->Fit("fprec7","RE","",fit_start,fit_stop);
      		
   
   
    
	/*   fit_func11->SetParameter(0, 0.2291);
       fit_func11->SetParameter(1, 0.0);
       fit_func11->SetParameter(2, 2.225);
       fit_func11->SetParameter(3, 0.00226);
       fit_func11->SetParameter(4, 290570);
       fit_func11->SetParameter(5, 0.0023404);
       fit_func11->SetParameter(6, 0.454);
       fit_func11->SetParameter(7, 0.0006);
       fit_func11->SetParameter(8, 6.7);
       fit_func11->SetParameter(9, 0.00014);
       fit_func11->SetParameter(10, -3.8);
	*/
	fit_func11->SetParameter(0, 0.229);
        fit_func11->SetParameter(1, 0.0);
        fit_func11->SetParameter(2, 2.224);
        fit_func11->SetParameter(3, 0.00238);
        fit_func11->SetParameter(4, 256000);
        fit_func11->SetParameter(5, 0.00234);
        fit_func11->SetParameter(6, 0.412);
	fit_func11->SetParameter(7, 0.0006);
        fit_func11->SetParameter(8, 6.7);
        fit_func11->SetParameter(9, 0.00014);
        fit_func11->SetParameter(10, -3.8);


   

    
	//   hr_sum->Fit("fprec11","RE","",fit_start,fit_stop);
      		
       
   
    
	/*   fit_func15->SetParameter(0, 0.2291);
       fit_func15->SetParameter(1, 0.0);
       fit_func15->SetParameter(2, 2.224);
       fit_func15->SetParameter(3, 0.00254);
       fit_func15->SetParameter(4, 224300);
       fit_func15->SetParameter(5, 0.0023404);
       fit_func15->SetParameter(6, 0.459);
       fit_func15->SetParameter(7, 0.0006);
       fit_func15->SetParameter(8, 6.7);
       fit_func15->SetParameter(9, -0.00015);
       fit_func15->SetParameter(10, -6.9);
       fit_func15->SetParameter(11, 0.008);
       fit_func15->SetParameter(12, 17400);
       fit_func15->SetParameter(13, 0.01403);
       fit_func15->SetParameter(14, 3);
	*/

	fit_func15->SetParameter(0, 0.229);
        fit_func15->SetParameter(1, 0.0);
        fit_func15->SetParameter(2, 2.224);
        fit_func15->SetParameter(3, 0.00238);
        fit_func15->SetParameter(4, 256000);
        fit_func15->SetParameter(5, 0.00234);
        fit_func15->SetParameter(6, 0.412);
	fit_func15->SetParameter(7, 0.0007);
        fit_func15->SetParameter(8, 6.6);
        fit_func15->SetParameter(9, -0.0001);
        fit_func15->SetParameter(10, 0.05);
	fit_func15->SetParameter(11, 0.008);
        fit_func15->SetParameter(12, 25400);
        fit_func15->SetParameter(13, 0.01403);
        fit_func15->SetParameter(14, 3);


   
	//    hr_sum->Fit("fprec15","RE","",fit_start,fit_stop);
      		
       
   
          

       fit_func19->SetParameter(0, 0.2291);
       fit_func19->SetParameter(1, 0.0);
       fit_func19->SetParameter(2, 2.224);
       fit_func19->SetParameter(3, 0.00255);
       fit_func19->SetParameter(4, 224460);
       fit_func19->SetParameter(5, 0.0023404);
       fit_func19->SetParameter(6, 0.46);
       fit_func19->SetParameter(7, 0.00066);
       fit_func19->SetParameter(8, 0.39);
       fit_func19->SetParameter(9, -0.000143);
       fit_func19->SetParameter(10, -19.5);
       fit_func19->SetParameter(11, 0.0036);
       fit_func19->SetParameter(12, 25144);
       fit_func19->SetParameter(13, 0.01403);
       fit_func19->SetParameter(14, 2.79);
       fit_func19->SetParameter(15, 0.00022);
       fit_func19->SetParameter(16, 335000);
       fit_func19->SetParameter(17, 0.01393);
       fit_func19->SetParameter(18, 6.27);
	

       /*   	fit_func19->SetParameter(0, 0.229);
        fit_func19->SetParameter(1, 0.0);
        fit_func19->SetParameter(2, 2.224);
        fit_func19->SetParameter(3, 0.00238);
        fit_func19->SetParameter(4, 256000);
        fit_func19->SetParameter(5, 0.00234);
        fit_func19->SetParameter(6, 0.412);
	fit_func19->SetParameter(7, 0.0007);
        fit_func19->SetParameter(8, 6.6);
        fit_func19->SetParameter(9, -0.0001);
        fit_func19->SetParameter(10, 0.05);
	fit_func19->SetParameter(11, 0.006);
        fit_func19->SetParameter(12, 26900);
        fit_func19->SetParameter(13, 0.01403);
        fit_func19->SetParameter(14, 1.67);
	fit_func19->SetParameter(15, 0.00022);
        fit_func19->SetParameter(16, 335000);
        fit_func19->SetParameter(17, 0.01393);
        fit_func19->SetParameter(18, 6.27);
       */

      
        hr_sum->Fit("fprec19","RUE","",fit_start,fit_stop);
       
	
    
   
   
     	
       gettimeofday(&t_end, NULL);
       printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
       gStyle->SetOptFit(1111);
      		
       h_res= new TH1D("residual histogram ", "h_res", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_sum->GetBinContent(ibin)- fit_func19->Eval( hr_sum->GetBinCenter(ibin) ) );
          if(hr_sum->GetBinError(ibin)!=0){res=(res/hr_sum->GetBinError(ibin));}
          h_res->SetBinContent(ibin, (res)  );
     
        }

        h_res->GetYaxis()->SetTitle("Energy [MeV]");
	h_res->GetXaxis()->SetTitle("time [ns]");
	h_res->GetXaxis()->SetRangeUser(fit_start, fit_stop);
        h_res->GetXaxis()->SetLabelSize(0.05);
        h_res->GetYaxis()->SetLabelSize(0.05);


 
       hfft = h_res->FFT(hfft, "MAG");
       hfft->SetLineColor(kBlack);
       hfft->SetBins(hr_sum->GetNbinsX(),0,1000/h_res->GetBinWidth(1));
       hfft->GetXaxis()->SetRangeUser(0,1000/(2*h_res->GetBinWidth(1)));
       hfft->GetXaxis()->SetTitle("Freq [MHz]");
       hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft->GetXaxis()->SetLabelSize(0.05);
       hfft->GetYaxis()->SetLabelSize(0.05);
 
       cnew=new TCanvas("cnew","ratio fit");
       cnew->Divide(1,3);
       cnew->cd(1);
       hr_sum->Draw();
       cnew->cd(2);
       h_res->Draw();
       cnew->cd(3);
       hfft->Draw();

        cauto=new TCanvas("auto-correlation", "auto-correlation");
	   hFFTsq = (TH1D*)hfft->Clone();
           hFFTsq->Reset();
           for (int ib = 1; ib <= hFFTsq->GetNbinsX(); ib++)
	     {
	       hFFTsq->SetBinContent( ib, hfft->GetBinContent(ib)*hfft->GetBinContent(ib));
	     }


       hAuto = hFFTsq->FFT(hAuto, "MAG");
       hAuto->SetTitle("auto-correlation");
       hAuto->GetXaxis()->SetLabelSize(0.04);
       hAuto->GetXaxis()->SetTitle("bins");
       hAuto->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/2);
       hAuto->SetLineColor(kRed);
       hAuto->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hAuto->Draw("HIST");




   
}
