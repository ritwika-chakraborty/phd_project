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
char root_file_name[128] = "run2_thresh300_fbfDQC.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed.root";
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[50], *hr_sum, *hr_calosum, *h_res_calosum;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28;
TGraph *gkband1p,*gkband1m,*gkband2p,*gkband2m,*gkband3p,*gkband3m,*gkband4p,*gkband4m,*gkband5p,*gkband5m,*gkband6p,*gkband6m,*gkband7p,*gkband7m,*gkband8p,*gkband8m,*gkband9p,*gkband9m,*gkband10p,*gkband10m,*gkband11p,*gkband11m,*gkband12p,*gkband12m,*gkband13p,*gkband13m,*gkband14p,*gkband14m,*gkband15p,*gkband15m,*gkband16p,*gkband16m,*gkband17p,*gkband17m,*gkband18p,*gkband18m,*gkband19p,*gkband19m,*gkband20p,*gkband20m,*gkband21p,*gkband21m,*gkband22p,*gkband22m,*gkband23p,*gkband23m;
TH1 *hfft[50], *hfft_calosum;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
//Double_t fit_start=29926.250;
//Double_t fit_stop=299926.25;
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
Int_t m=0;
Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi_0[50], dphi_0[50], A_cbo_N[50], tau_cbo[50], omega_cbo[50], phi_cbo_N[50], A_cbo_A[50], phi_cbo_A[50], A_cbo_phi[50], phi_cbo_phi[50], A_vw[50], tau_vw[50], omega_vw[50], phi_vw[50], A_y[50], tau_y[50], omega_y[50], phi_y[50], A_2cbo[50], tau_2cbo[50], omega_2cbo[50], phi_2cbo[50],  dA_cbo_N[50], dtau_cbo[50], domega_cbo[50], dphi_cbo_N[50], dA_cbo_A[50], dphi_cbo_A[50], dA_cbo_phi[50], dphi_cbo_phi[50], dA_vw[50], dtau_vw[50], domega_vw[50], dphi_vw[50], dA_y[50], dtau_y[50], domega_y[50], dphi_y[50], dA_2cbo[50], dtau_2cbo[50], domega_2cbo[50], dphi_2cbo[50], n[50],kband1p[50],kband1m[50],kband2p[50],kband2m[50],kband3p[50],kband3m[50],kband4p[50],kband4m[50],kband5p[50],kband5m[50],kband6p[50],kband6m[50],kband7p[50],kband7m[50],kband8p[50],kband8m[50],kband9p[50],kband9m[50],kband10p[50],kband10m[50],kband11p[50],kband11m[50],kband12p[50],kband12m[50],kband13p[50],kband13m[50],kband14p[50],kband14m[50],kband15p[50],kband15m[50],kband16p[50],kband16m[50],kband17p[50],kband17m[50],kband18p[50],kband18m[50],kband19p[50],kband19m[50],kband20p[50],kband20m[50],kband21p[50],kband21m[50],kband22p[50],kband22m[50],kband23p[50],kband23m[50];
Int_t countfcn=0;
Double_t ch=0;

TFile *foutput;
TH2D *hcov;
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
TMatrixD cov2(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim));

TMatrixD ratio_cov_matrix(TH1D *hr_sum, TH1D *h_sum, Double_t fit_start){
  TMatrixD cov(mdim,mdim);
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

  double u,v,w,x,z,s,t;

  for(Int_t k=hr_sum->FindBin(10000); k<=hr_sum->GetNbinsX(); k++)
    {
      if(hr_sum->GetBinError(k)!=0)
	{ i0=k;
	  break;}
    }
  //cout<<"i0 is"<<i0<<" "<<endl;
  
   for (Int_t i = 0; i < (mdim)*(mdim); i++)
      {
	const Int_t ir = (i)/(mdim);
	const Int_t ic = (i)%(mdim);
        
          u=h_sum->GetBinContent(ir);
          v=h_sum->GetBinContent(ir+nbinshift);
          w=h_sum->GetBinContent(ir-nbinshift);
          x=h_sum->GetBinContent(ir+(2*nbinshift));
          z=h_sum->GetBinContent(ir-(2*nbinshift));
          s=h_sum->GetBinContent(ir+(3*nbinshift));
          t=h_sum->GetBinContent(ir-(3*nbinshift));
        
	
        data[i] = 0.0;
 
        if(ir==ic)
	  {
	    //   data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir))/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	     data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	     /*	    if(data[i]==0)
	      {
		cout<<"diag ele is zero at "<<ir<<" "<<ic<<endl;
		}*/
	     
	     //if(ir==247){cout<<i<<" "<<ir<<" "<<ic<<" "<<data[i]<<endl;}
	  }

	    	if(ir==ic+nbinshift)
	  {

	    dprod1=(8*(-25*v*v*v*v + 2*u*u*u*(5*v + w) + 6*u*u*(v + w)*x + 2*w*x*x*x - v*v*v*(35*w + 12*x) + v*v*(-11*w*w - 16*w*x + 4*x*x) - v*(w*w*w + 4*w*w*x - 4*w*x*x - 2*x*x*x) - 6*u*(5*v*v*v + 6*v*v*w - w*x*x + v*(w*w - x*x))))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x));

	    dprod2=-((8*(25*u*u*u*u - 2*(v + w)*(v + w)*(v + w)*x + u*u*u*(30*v + 12*w + 35*x) + u*(-10*v*v*v - 6*v*v*w - 6*v*w*w - 2*w*w*w - 4*w*w*x + 6*v*x*x + 4*w*x*x + x*x*x) + u*u*(-4*w*w + 16*w*x + x*(36*v + 11*x))))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x)));
	      
	    dprod3=-((8*u*(u - 2*v + x))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*v + x)));

	    dprod4=-((8*v*(-2*u + v + w))/((2*u + v + w)*(u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x)));

	    dI11=-((16*(v + w))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)));

	    dI12=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

	    dI13=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

	    dI21=(8*v)/((u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x));

	    dI22=-((16*(u + x))/((u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x)));

	    dI23=(8*v)/((u + 2*v + x)*(u + 2*v + x)*(u + 2*v + x));

	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir+nbinshift) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dprod4*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift))-(hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI13*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift))*(hr_sum->GetBinContent(ir+nbinshift) + 0.5*dI21*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI22*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI23*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift)));

	    //  data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	       data[i]=covar;
	    // data[i]=(-2.0/3.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //    data[i]=(-2.0/3.0);
	       //     cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

		if(ir==ic-nbinshift)
	  {

	    dprod1=(8*(-v*v*v*w + 2*u*u*u*(v + 5*w) + 6*u*u*(v + w)*z - v*v*w*(11*w + 4*z) - 6*u*(v*v*w + 6*v*w*w + 5*w*w*w - v*z*z - w*z*z) + v*(-35*w*w*w - 16*w*w*z + 4*w*z*z + 2*z*z*z) + w*(-25*w*w*w - 12*w*w*z + 4*w*z*z + 2*z*z*z)))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z));

	    dprod2=-((8*u*(u - 2*w + z))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*w + z)));
	      
	    dprod3=-((8*(25*u*u*u*u - 2*(v + w)*(v + w)*(v + w)*z + u*u*u*(12*v + 30*w + 35*z) + u*u*(-4*v*v + 16*v*z + z*(36*w + 11*z)) + u*(-2*v*v*v - 10*w*w*w + 6*w*z*z + z*z*z - 2*v*v*(3*w + 2*z) + v*(-6*w*w + 4*z*z))))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z)));

	    dprod4=(8*(2*u - v - w)*w)/((2*u + v + w)*(u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z));
	      
	    dI11=-((16*(v + w))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)));
	      
	    dI12=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

	    dI13=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

	    dI21=(8*w)/((u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z));
	      
	    dI22=-((16*(u + z))/((u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z)));

	    dI23=(8*w)/((u + 2*w + z)*(u + 2*w + z)*(u + 2*w + z));
	      
	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir-nbinshift) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dprod4*h_sum->GetBinError(ir-(2*nbinshift))*h_sum->GetBinError(ir-(2*nbinshift))-(hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI13*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift))*(hr_sum->GetBinContent(ir-nbinshift) + 0.5*dI21*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI22*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dI23*h_sum->GetBinError(ir-(2*nbinshift))*h_sum->GetBinError(ir-(2*nbinshift)));

	    //  data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    
	     data[i]=covar;
	     //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	    //  data[i]=(-2.0/3.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //   data[i]=(-2.0/3.0);
	  }

		if(ir==ic+(2*nbinshift))
	      {

		dprod1=(16*(v + w)*(s + v - 2*x))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(s + v + 2*x));

		dprod2=-((8*(s*s*s*u - 8*u*u*u*x + (v + w)*(v + w)*(v + w)*x + s*s*u*(3*v + 2*x) - 4*u*u*x*(3*v + w + 4*x) + s*u*(3*v*v - 4*x*(2*u + w + x)) + u*(v*v*v - 12*v*x*x + 2*x*(w*w - 4*w*x - 4*x*x))))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x)));

		dprod3=-((8*u*(s + v - 2*x))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(s + v + 2*x)));

		dprod4=(16*(s + v)*(-2*u + v + w))/((2*u + v + w)*(s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x));

		dprod5=(8*(2*u - v - w)*x)/((2*u + v + w)*(s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x));

	        dI11=-((16*(v + w))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)));

		dI12=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

		dI13=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

		dI21=(8*x)/((s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x));

		dI22=(8*x)/((s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x));

		dI23=-((16*(s + v))/((s + v + 2*x)*(s + v + 2*x)*(s + v + 2*x)));

		covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir+(2*nbinshift)) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dprod4*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift)) +  0.5*dprod5*h_sum->GetBinError(ir+(3*nbinshift))*h_sum->GetBinError(ir+(3*nbinshift)) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI13*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift))*(hr_sum->GetBinContent(ir+(2*nbinshift)) + 0.5*dI21*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI22*h_sum->GetBinError(ir+(3*nbinshift))*h_sum->GetBinError(ir+(3*nbinshift)) + 0.5*dI23*h_sum->GetBinError(ir+(2*nbinshift))*h_sum->GetBinError(ir+(2*nbinshift)));

		// 	data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	        
	       	  data[i]=covar;
		  //	cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		//	data[i]=(1.0/6.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	data[i]=(1.0/6.0);
	      }

		if(ir==ic-(2*nbinshift))
	      {

		dprod1=(16*(v + w)*(t + w - 2*z))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(t + w + 2*z));

		dprod2=-((8*u*(t + w - 2*z))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(t + w + 2*z)));

		dprod3=-((8*(t*t*t*u - 8*u*u*u*z + (v + w)*(v + w)*(v + w)*z + t*t*u*(3*w + 2*z) - 4*u*u*z*(v + 3*w + 4*z) + t*u*(3*w*w - 4*z*(2*u + v + z)) + u*(w*w*w - 12*w*z*z + 2*z*(v*v - 4*v*z - 4*z*z))))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)*(t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z)));

		dprod4=-((16*(2*u - v - w)*(t + w))/((2*u + v + w)*(t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z)));

		dprod5=(8*(2*u - v - w)*z)/((2*u + v + w)*(t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z));

		dI11=-((16*(v + w))/((2*u + v + w)*(2*u + v + w)*(2*u + v + w)));

		dI12=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

		dI13=(8*u)/((2*u + v + w)*(2*u + v + w)*(2*u + v + w));

		dI21=(8*z)/((t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z));

		dI22=(8*z)/((t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z));

		dI23=-((16*(t + w))/((t + w + 2*z)*(t + w + 2*z)*(t + w + 2*z)));

		covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir-(2*nbinshift)) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dprod3*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dprod4*h_sum->GetBinError(ir-(2*nbinshift))*h_sum->GetBinError(ir-(2*nbinshift)) +  0.5*dprod5*h_sum->GetBinError(ir-(3*nbinshift))*h_sum->GetBinError(ir-(3*nbinshift)) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+nbinshift)*h_sum->GetBinError(ir+nbinshift) + 0.5*dI13*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift))*(hr_sum->GetBinContent(ir-(2*nbinshift)) + 0.5*dI21*h_sum->GetBinError(ir-nbinshift)*h_sum->GetBinError(ir-nbinshift) + 0.5*dI22*h_sum->GetBinError(ir-(3*nbinshift))*h_sum->GetBinError(ir-(3*nbinshift)) + 0.5*dI23*h_sum->GetBinError(ir-(2*nbinshift))*h_sum->GetBinError(ir-(2*nbinshift)));

		// 	data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
			
			data[i]=covar;
			//		cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		// 	data[i]=(1.0/6.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	data[i]=(1.0/6.0);
			
	      }
	 
		//	cout<<data[i]<<endl;	
      }
  cout<<"flgs "<<flg<<" "<<endl;
  // cout<<data[518947]<<endl;
  cov.SetMatrixArray(data.GetArray());
  //cout<<data[518947]<<endl;
  // cov.Print();
  //cout<<"matrix element "<<cov[hr_sum->FindBin(fit_start)][hr_sum->FindBin(fit_start)]<<endl;
  //  cov.ResizeTo(h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), -1);
  //cout<<"start bin "<<hr_sum->FindBin(fit_start)<<endl;
  cov.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
  //cout<<"matrix element "<<cov[hr_sum->FindBin(fit_start)][hr_sum->FindBin(fit_start)]<<endl;
  //   cov.Print();
    cov.SetTol(1.e-23);
    Double_t det1;
    cov.Invert(&det1);
    //   cov.Print();
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

    double asym_2cbo= par[11];

    double tau_2cbo = par[12];

    double omega_2cbo = par[13];

    double phi_2cbo = par[14];
    
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

    double N2cbo=(asym_2cbo*exp(-time/tau_2cbo)*cos(omega_2cbo*time + phi_2cbo));

    double N2cbof=(asym_2cbo*exp(-(time + T_a/2)/tau_2cbo)*cos(omega_2cbo*(time + T_a/2) + phi_2cbo));

    double N2cbob=(asym_2cbo*exp(-(time - T_a/2)/tau_2cbo)*cos(omega_2cbo*(time - T_a/2) + phi_2cbo));

    Ncbo=Ncbo+N2cbo;
    Ncbof=Ncbof+N2cbof;
    Ncbob=Ncbob+N2cbob;


    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

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

    double N2cbo=(asym_2cbo*exp(-time/tau_2cbo)*cos(omega_2cbo*time + phi_2cbo));

    double N2cbof=(asym_2cbo*exp(-(time + T_a/2)/tau_2cbo)*cos(omega_2cbo*(time + T_a/2) + phi_2cbo));

    double N2cbob=(asym_2cbo*exp(-(time - T_a/2)/tau_2cbo)*cos(omega_2cbo*(time - T_a/2) + phi_2cbo));

    Ncbo=Ncbo+N2cbo;
    Ncbof=Ncbof+N2cbof;
    Ncbob=Ncbob+N2cbob;

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

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

    double Ncbo=(1 + asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

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

    double N2cbo=(asym_2cbo*exp(-time/tau_2cbo)*cos(omega_2cbo*time + phi_2cbo));

    double N2cbof=(asym_2cbo*exp(-(time + T_a/2)/tau_2cbo)*cos(omega_2cbo*(time + T_a/2) + phi_2cbo));

    double N2cbob=(asym_2cbo*exp(-(time - T_a/2)/tau_2cbo)*cos(omega_2cbo*(time - T_a/2) + phi_2cbo));

    Ncbo=Ncbo+N2cbo;
    Ncbof=Ncbof+N2cbof;
    Ncbob=Ncbob+N2cbob;

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof - fb*Ncbob*Nvwb*Nvbob)/(2*f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof + fb*Ncbob*Nvwb*Nvbob);

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
   //  nbinshift=15;
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
   double blife=exp(-(T_a)/(2*lifetime));
   double deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo histogram sum numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

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
       
       long double t1,t2,t3,t4;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
       	t1=2*(h1->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

	t2=2*(h2->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t3=-2*(hp->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t4=-2*(hm->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

 	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt( ((t1 + t2)*(t1 + t2)) + (t3*t3) + (t4*t4) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( ((t1 + t2)*(t1 + t2)) + (t3*t3) + (t4*t4) ) );
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

 //hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //hcalo->Scale(0.125);

 return hcalo;
  }

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
 
   if(hr_sum->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }

  f = 0;
  ch=0;
  //  countfcn=0;

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
	      //cout<<m<<endl;
	    
	     }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}


void ratio_calo_scan()
{
  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzerDB",dir[1]);
   dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  // _file[1]->GetObject("hwiggle",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

    _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzerDB",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzerDB",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzerDB",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzerDB",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

    _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzerDB",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzerDB",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzerDB",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzerDB",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzerDB",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

  _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzerDB",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzerDB",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzerDB",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzerDB",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzerDB",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzerDB",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzerDB",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzerDB",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzerDB",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzerDB",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzerDB",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzerDB",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzerDB",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
  _file[24]->GetObject("QFillByFillAnalyzerDB",dir[24]);
  dir[24]->GetObject("qHist1D_sig_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();


  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  hr_sum= new TH1D("calo histogram sum ratio", "h_ratio", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
  for(Int_t i=1; i<=24; i++)
    {
      h_sum->Add(h[i],1); 
    }

   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   for(int i=1;i<=24;i++)
     {
      h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

  
   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));





    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");



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
     fit_func19->SetParName(15,"A_vw");
     fit_func19->SetParName(16,"Tau_vw");
     fit_func19->SetParName(17,"omega_vw");
     fit_func19->SetParName(18,"phi_vw");

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
     fit_func23->SetParName(15,"A_vw");
     fit_func23->SetParName(16,"Tau_vw");
     fit_func23->SetParName(17,"omega_vw");
     fit_func23->SetParName(18,"phi_vw");
     fit_func23->SetParName(19,"A_y");
     fit_func23->SetParName(20,"Tau_y");
     fit_func23->SetParName(21,"omega_y");
     fit_func23->SetParName(22,"phi_y");
     fit_func23->SetNpx(1000000);


     /*   fit_func19->SetParameter(0, 0.22);
       fit_func19->SetParameter(1, 0.0);
       fit_func19->SetParameter(2, 2.214);
       fit_func19->SetParameter(3, 0.02028);
       fit_func19->SetParameter(4, 239700);
       fit_func19->SetParameter(5, 0.00234071);
       fit_func19->SetParameter(6, 4.16);
       fit_func19->SetParameter(7, 0.0217);
       fit_func19->SetParameter(8, 2.1);
       fit_func19->SetParameter(9, -0.003);
       fit_func19->SetParameter(10, -23.5);
       fit_func19->SetParameter(11, 0.000036);
       fit_func19->SetParameter(12, 55144);
       fit_func19->SetParameter(13, 0.01403);
       fit_func19->SetParameter(14, 2.79);
       fit_func19->SetParameter(15, 0.0000022);
       fit_func19->SetParameter(16, 335000);
       fit_func19->SetParameter(17, 0.01393);
       fit_func19->SetParameter(18, 6.27);


       fit_func23->SetParameter(0, 0.22189);
       fit_func23->SetParameter(1, 0.0);
       fit_func23->SetParameter(2, 2.2143);
       fit_func23->SetParameter(3, 0.02027);
       fit_func23->SetParameter(4, 240400);
       fit_func23->SetParameter(5, 0.0023407);
       fit_func23->SetParameter(6, 4.16-0.262);
       fit_func23->SetParameter(7, 0.0217);
       fit_func23->SetParameter(8, 2.17-0.262);
       fit_func23->SetParameter(9, -0.0029);
       fit_func23->SetParameter(10, -23.4-0.262);
       fit_func23->SetParameter(11, 0.001779);
       fit_func23->SetParameter(12, 32500);
       fit_func23->FixParameter(13, 0.01404);
       fit_func23->SetParameter(14, 3.26);
       fit_func23->SetParameter(15, 0.0002);
       fit_func23->SetParameter(16, 82400);
       fit_func23->FixParameter(17, 0.01392);
       fit_func23->SetParameter(18, 7.77);
       fit_func23->SetParameter(19, 0.00048);
       fit_func23->SetParameter(20, 128000);
       fit_func23->SetParameter(21, 0.00468);
       fit_func23->SetParameter(22, 0.67);
     

       fit_func7->SetParameter(0, 0.2291);
      fit_func7->SetParameter(1, 0.0);
      fit_func7->SetParameter(2, 2.22);
      fit_func7->SetParameter(3, 0.025);
      fit_func7->SetParameter(4, 257000);
      fit_func7->SetParameter(5, 0.0023406);
      fit_func7->SetParameter(6, 0.45);
     
          fit_func7->SetParameter(0, 0.2219);
        fit_func7->SetParameter(1, 0.0);
        fit_func7->SetParameter(2, 2.214);
        fit_func7->SetParameter(3, 0.0205);
        fit_func7->SetParameter(4, 239600);
        fit_func7->SetParameter(5, 0.00234059);
        fit_func7->SetParameter(6, 4.19);
     */   
       
     fit_func19->SetNpx(1000000);

     


   
     
     for(int i=1; i<=24;i++)
       {
	hr_sum->Reset();
	cout<<"calo "<<i<<endl; 
        
	 
	//	fit_func3->SetParameters(0.22, 0.0, 2.22);

	hr_sum = construct_rhist_copy(h[i]);
  
	hr_sum->GetYaxis()->SetTitle("ADC counts");
	hr_sum->GetXaxis()->SetTitle("time [ns]");
        hr_sum->GetXaxis()->SetRangeUser(fit_start,fit_stop);




		
 	cov2.ResizeTo(h[i]->FindBin(fit_start), h[i]->FindBin(fit_stop), h[i]->FindBin(fit_start), h[i]->FindBin(fit_stop), -1);
   

        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr_sum,h[i], fit_start);

   
          

	/*   fit_func23->SetParameter(0, fit_func23->GetParameter(0));
       fit_func23->SetParameter(1, fit_func23->GetParameter(1));
       fit_func23->SetParameter(2, fit_func23->GetParameter(2));
       fit_func23->SetParameter(3, fit_func23->GetParameter(3));
       fit_func23->SetParameter(4, fit_func23->GetParameter(4));
       fit_func23->SetParameter(5, fit_func23->GetParameter(5));
       fit_func23->SetParameter(6, fit_func23->GetParameter(6)+0.262);
       fit_func23->SetParameter(7, fit_func23->GetParameter(7));
       fit_func23->SetParameter(8, fit_func23->GetParameter(8)+0.262);
       fit_func23->SetParameter(9, fit_func23->GetParameter(9));
       fit_func23->SetParameter(10, fit_func23->GetParameter(10)+0.262);
       fit_func23->SetParameter(11, fit_func23->GetParameter(11));
       fit_func23->SetParameter(12, fit_func23->GetParameter(12));
       fit_func23->SetParameter(13, fit_func23->GetParameter(13));
       fit_func23->SetParameter(14, fit_func23->GetParameter(14));
       fit_func23->SetParameter(15, fit_func23->GetParameter(15));
       fit_func23->SetParameter(16, fit_func23->GetParameter(16));
       fit_func23->SetParameter(17, fit_func23->GetParameter(17));
       fit_func23->SetParameter(18, fit_func23->GetParameter(18));
       fit_func23->SetParameter(19, fit_func23->GetParameter(19));
       fit_func23->SetParameter(20, fit_func23->GetParameter(20));
       fit_func23->SetParameter(21, fit_func23->GetParameter(21));
       fit_func23->SetParameter(22, fit_func23->GetParameter(22));
       	

       	  fit_func19->SetParameter(0, fit_func19->GetParameter(0));
       fit_func19->SetParameter(1, fit_func19->GetParameter(1));
       fit_func19->SetParameter(2, fit_func19->GetParameter(2));
       fit_func19->SetParameter(3, fit_func19->GetParameter(3));
       fit_func19->SetParameter(4, fit_func19->GetParameter(4));
       fit_func19->SetParameter(5, fit_func19->GetParameter(5));
       fit_func19->SetParameter(6, fit_func19->GetParameter(6));
       fit_func19->SetParameter(7, fit_func19->GetParameter(7));
       fit_func19->SetParameter(8, fit_func19->GetParameter(8));
       fit_func19->SetParameter(9, fit_func19->GetParameter(9));
       fit_func19->SetParameter(10, fit_func19->GetParameter(10));
       fit_func19->SetParameter(11, fit_func19->GetParameter(11));
       fit_func19->SetParameter(12, fit_func19->GetParameter(12));
       fit_func19->SetParameter(13, fit_func19->GetParameter(13));
       fit_func19->SetParameter(14, fit_func19->GetParameter(14));
       fit_func19->SetParameter(15, fit_func19->GetParameter(15));
       fit_func19->SetParameter(16, fit_func19->GetParameter(16));
       fit_func19->SetParameter(17, fit_func19->GetParameter(17));
       fit_func19->SetParameter(18, fit_func19->GetParameter(18));
	*/
	

       /*   fit_func7->SetParameter(0, fit_func7->GetParameter(0));
       fit_func7->SetParameter(1, fit_func7->GetParameter(1));
       fit_func7->SetParameter(2, fit_func7->GetParameter(2));
       fit_func7->SetParameter(3, fit_func7->GetParameter(3));
       fit_func7->SetParameter(4, fit_func7->GetParameter(4));
       fit_func7->SetParameter(5, fit_func7->GetParameter(5));
       fit_func7->SetParameter(6, fit_func7->GetParameter(6));
       */
          	fit_func3->SetParameters(0.22, 0.0, 2.22);


        hr_sum->Fit("fprec3","RE","",fit_start,fit_stop);
 
      
		
    
   
   
        fit_func7->SetParameter(0, fit_func3->GetParameter(0));
        fit_func7->SetParameter(1, fit_func3->GetParameter(1));
        fit_func7->SetParameter(2, fit_func3->GetParameter(2));
        fit_func7->SetParameter(3, 0.0204);
        fit_func7->FixParameter(4, 243936.1);
        fit_func7->FixParameter(5, 0.0023403791);
        fit_func7->SetParameter(6, 2.1);

   
   	hr_sum->Fit("fprec7","RE","",fit_start,fit_stop);
      		
   
   
    
       fit_func11->SetParameter(0, fit_func7->GetParameter(0));
       fit_func11->SetParameter(1, fit_func7->GetParameter(1));
       fit_func11->SetParameter(2, fit_func7->GetParameter(2));
       fit_func11->SetParameter(3, fit_func7->GetParameter(3));
       fit_func11->FixParameter(4, fit_func7->GetParameter(4));
       fit_func11->FixParameter(5, fit_func7->GetParameter(5));
       fit_func11->SetParameter(6, fit_func7->GetParameter(6));
       fit_func11->SetParameter(7, 0.02);
       fit_func11->SetParameter(8, 2.1);
       fit_func11->SetParameter(9, 0.003);
       fit_func11->SetParameter(10, -1.4);


    
       hr_sum->Fit("fprec11","RE","",fit_start,fit_stop);
      		
       
   
    
       fit_func15->SetParameter(0, fit_func11->GetParameter(0));
       fit_func15->SetParameter(1, fit_func11->GetParameter(1));
       fit_func15->SetParameter(2, fit_func11->GetParameter(2));
       fit_func15->SetParameter(3, fit_func11->GetParameter(3));
       fit_func15->FixParameter(4, fit_func11->GetParameter(4));
       fit_func15->FixParameter(5, fit_func11->GetParameter(5));
       fit_func15->SetParameter(6, fit_func11->GetParameter(6));
       fit_func15->SetParameter(7, fit_func11->GetParameter(7));
       fit_func15->SetParameter(8, fit_func11->GetParameter(8));
       fit_func15->SetParameter(9, fit_func11->GetParameter(9));
       fit_func15->SetParameter(10, fit_func11->GetParameter(10));
       fit_func15->SetParameter(11, 0.0005);
       fit_func15->FixParameter(12, fit_func11->GetParameter(4)/2);
       fit_func15->FixParameter(13, 2*fit_func11->GetParameter(5));
       fit_func15->SetParameter(14, 0.67);

   
        hr_sum->Fit("fprec15","RE","",fit_start,fit_stop);
      		

        
  

       fit_func19->SetParameter(0, fit_func15->GetParameter(0));
       fit_func19->SetParameter(1, fit_func15->GetParameter(1));
       fit_func19->SetParameter(2, fit_func15->GetParameter(2));
       fit_func19->SetParameter(3, fit_func15->GetParameter(3));
       fit_func19->FixParameter(4, fit_func15->GetParameter(4));
       fit_func19->FixParameter(5, fit_func15->GetParameter(5));
       fit_func19->SetParameter(6, fit_func15->GetParameter(6));
       fit_func19->SetParameter(7, fit_func15->GetParameter(7));
       fit_func19->SetParameter(8, fit_func15->GetParameter(8));
       fit_func19->SetParameter(9, fit_func15->GetParameter(9));
       fit_func19->SetParameter(10, fit_func15->GetParameter(10));
       fit_func19->SetParameter(11, fit_func15->GetParameter(11));
       fit_func19->FixParameter(12, fit_func15->GetParameter(12));
       fit_func19->FixParameter(13, fit_func15->GetParameter(13));
       fit_func19->SetParameter(14, fit_func15->GetParameter(14));
       fit_func19->SetParameter(15, 0.002);
       fit_func19->FixParameter(16, 28629);
       fit_func19->FixParameter(17, 0.014037648);
       fit_func19->SetParameter(18, -2.4);


       hr_sum->Fit(fit_func19,"RE","",fit_start,fit_stop);

       fit_func23->SetParameter(0, fit_func19->GetParameter(0));
       fit_func23->SetParameter(1, fit_func19->GetParameter(1));
       fit_func23->SetParameter(2, fit_func19->GetParameter(2));
       fit_func23->SetParameter(3, fit_func19->GetParameter(3));
       fit_func23->FixParameter(4, fit_func19->GetParameter(4));
       fit_func23->FixParameter(5, fit_func19->GetParameter(5));
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
       fit_func23->SetParameter(19, 0.01);
       fit_func23->FixParameter(20, 111541.29);
       fit_func23->FixParameter(21, 0.013931315);
       fit_func23->SetParameter(22, 2.1);



       hr_sum->Fit(fit_func23,"RE","",fit_start,fit_stop);

       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hr_sum);
       gFitter->SetFCN(chi2);  
	
       //   hr_sum->Fit("fprec23","RE","",fit_start,fit_stop);
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
       A_vw[m]=TMath::Abs(fit_func23->GetParameter(15));
       tau_vw[m]=fit_func23->GetParameter(16);
       omega_vw[m]=fit_func23->GetParameter(17);
       phi_vw[m]=fit_func23->GetParameter(18);
       A_y[m]=TMath::Abs(fit_func23->GetParameter(19));
       tau_y[m]=fit_func23->GetParameter(20);
       omega_y[m]=fit_func23->GetParameter(21);
       phi_y[m]=fit_func23->GetParameter(22);

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
       dA_vw[m]=fit_func23->GetParError(15);
       dtau_vw[m]=fit_func23->GetParError(16);
       domega_vw[m]=fit_func23->GetParError(17);
       dphi_vw[m]=fit_func23->GetParError(18);
       dA_y[m]=fit_func23->GetParError(19);
       dtau_y[m]=fit_func23->GetParError(20);
       domega_y[m]=fit_func23->GetParError(21);
       dphi_y[m]=fit_func23->GetParError(22);

       chisq[m]=fit_func23->GetChisquare()/fit_func23->GetNDF();


       /*    while(phi_cbo_N[m] < 0.0)
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
       /*     if(m==0)
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
*/
       n[m]=i;
       m=m+1;
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_sum->GetBinContent(ibin)- fit_func23->Eval( hr_sum->GetBinCenter(ibin) ) );
          if(hr_sum->GetBinError(ibin)!=0){res=(res/hr_sum->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.05);
        h_res[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hr_sum->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.05);
       hfft[i]->GetYaxis()->SetLabelSize(0.05);


       // fit_start=fit_start+1000;
 
   }

  	hr_calosum = construct_rhist_copy(h_sum);
  
	hr_calosum->GetYaxis()->SetTitle("ADC counts");
	hr_calosum->GetXaxis()->SetTitle("time [ns]");
        hr_calosum->GetXaxis()->SetRangeUser(fit_start,fit_stop);




		
 	cov2.ResizeTo(hr_calosum->FindBin(fit_start), hr_calosum->FindBin(fit_stop), hr_calosum->FindBin(fit_start), hr_calosum->FindBin(fit_stop), -1);
   

        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr_calosum,h_sum, fit_start);


          	fit_func3->SetParameters(0.228, 0.0, 2.23);


        hr_calosum->Fit("fprec3","RE","",fit_start,fit_stop);
 
  h_res_calosum= new TH1D("residual histogram calosum", "h_res calosum", hr_calosum->GetNbinsX(), hr_calosum->GetBinLowEdge(1), (hr_calosum->GetBinLowEdge(hr_calosum->GetNbinsX())+hr_calosum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_calosum->GetBinLowEdge(1))/hr_calosum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_calosum->GetBinLowEdge(1))/hr_calosum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_calosum->GetBinContent(ibin)- fit_func3->Eval( hr_calosum->GetBinCenter(ibin) ) );
          if(hr_calosum->GetBinError(ibin)!=0){res=(res/hr_calosum->GetBinError(ibin));}
          h_res_calosum->SetBinContent(ibin, (res)  );
     
        }

	h_res_calosum->GetYaxis()->SetTitle("Energy [MeV]");
	h_res_calosum->GetXaxis()->SetTitle("time [ns]");
	h_res_calosum->GetXaxis()->SetLabelSize(0.05);
        h_res_calosum->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft_calosum = h_res_calosum->FFT(hfft_calosum, "MAG");
       hfft_calosum->SetLineColor(kBlack);
       hfft_calosum->SetBins(hr_calosum->GetNbinsX(),0,1000/h_res_calosum->GetBinWidth(1));
       hfft_calosum->GetXaxis()->SetRangeUser(0,1000/(2*h_res_calosum->GetBinWidth(1)));
       hfft_calosum->GetXaxis()->SetTitle("Freq [MHz]");
       hfft_calosum->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft_calosum->GetXaxis()->SetLabelSize(0.05);
       hfft_calosum->GetYaxis()->SetLabelSize(0.05);
       hfft_calosum->SetLineColor(kRed);


     
     c2 = new TCanvas("c2","start time residuals");
     c2->Divide(4,6);
     c3 = new TCanvas("c3","start time fft");
     c3->Divide(4,6);
  
 
     hfft_calosum->GetXaxis()->SetRangeUser(2.1,2.4);
     hfft_calosum->GetYaxis()->SetRangeUser(0,1000);

 for(int i=1;i<=24;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hfft[i]->GetXaxis()->SetRangeUser(2.1,2.4);
     hfft[i]->GetYaxis()->SetRangeUser(0,1000);
     hfft[i]->Draw();
     hfft_calosum->Draw("same");
   }
    c1=new TCanvas("c1","blind R vs fit start time");
    c1->Divide(2,2);
    c1->cd(1);
    gr1=new TGraphErrors(m,n,chisq,0,0);
    gr1->SetTitle("chisq vs caloriemeter index");
    gr1->GetXaxis()->SetTitle("calo index");
    gr1->GetYaxis()->SetTitle("reduced chi-sq");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    

     c1->cd(2);
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle("A vs caloriemeter index");
    gr2->GetXaxis()->SetTitle("calo index");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw();
    /*  gkband1p=new TGraph(m,n,kband1p);
    gkband1m=new TGraph(m,n,kband1m);
    gkband1m->SetLineColor(kBlue);
    gkband1p->SetLineColor(kBlue);
    gkband1p->Draw("same");
    gkband1m->Draw("same");
    */
     c1->cd(3);
    gr3=new TGraphErrors(m,n,blindR,0,dblindR);
    gr3->SetTitle("blindR vs caloriemeter index");
    gr3->GetXaxis()->SetTitle("calo index");
    gr3->GetYaxis()->SetTitle("Blind R [ppm]");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw();
    /* gkband2p=new TGraph(m,n,kband2p);
    gkband2m=new TGraph(m,n,kband2m);
    gkband2m->SetLineColor(kBlue);
    gkband2p->SetLineColor(kBlue);
    gkband2p->Draw("same");
    gkband2m->Draw("same");
    */

     c1->cd(4);
    gr4=new TGraphErrors(m,n,phi_0,0,dphi_0);
    gr4->SetTitle("phase vs caloriemeter index");
    gr4->GetXaxis()->SetTitle("calo index");
    gr4->GetYaxis()->SetTitle("phi_0");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->Draw();
    /* gkband3p=new TGraph(m,n,kband3p);
    gkband3m=new TGraph(m,n,kband3m);
    gkband3m->SetLineColor(kBlue);
    gkband3p->SetLineColor(kBlue);
    gkband3p->Draw("same");
    gkband3m->Draw("same");
    */

    c4=new TCanvas("c4","cbo vs fit start time");
    c4->Divide(2,2);
    c4->cd(1);
    gr8=new TGraphErrors(m,n,A_cbo_N,0,dA_cbo_N);
    gr8->SetTitle("A_cbo_N vs caloriemeter index");
    gr8->GetXaxis()->SetTitle("calo index");
    gr8->GetYaxis()->SetTitle("A_cbo_N");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kRed);
    gr8->Draw();
    /* gkband4p=new TGraph(m,n,kband4p);
    gkband4m=new TGraph(m,n,kband4m);
    gkband4m->SetLineColor(kBlue);
    gkband4p->SetLineColor(kBlue);
    gkband4p->Draw("same");
    gkband4m->Draw("same");
    */

     c4->cd(2);
    gr5=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr5->SetTitle("tau_cbo_ vs caloriemeter index");
    gr5->GetXaxis()->SetTitle("calo index");
    gr5->GetYaxis()->SetTitle("tau_cbo");
    gr5->SetMarkerStyle(20);
    gr5->SetLineColor(kRed);
    gr5->Draw();
    /* gkband5p=new TGraph(m,n,kband5p);
    gkband5m=new TGraph(m,n,kband5m);
    gkband5m->SetLineColor(kBlue);
    gkband5p->SetLineColor(kBlue);
    gkband5p->Draw("same");
    gkband5m->Draw("same");
    */

     c4->cd(3);
    gr6=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr6->SetTitle("omega_cbo vs caloriemeter index");
    gr6->GetXaxis()->SetTitle("calo index");
    gr6->GetYaxis()->SetTitle("omega_cbo");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kRed);
    gr6->Draw();
    /*  gkband6p=new TGraph(m,n,kband6p);
    gkband6m=new TGraph(m,n,kband6m);
    gkband6m->SetLineColor(kBlue);
    gkband6p->SetLineColor(kBlue);
    gkband6p->Draw("same");
    gkband6m->Draw("same");
    */

     c4->cd(4);
    gr7=new TGraphErrors(m,n,phi_cbo_N,0,dphi_cbo_N);
    gr7->SetTitle("phi_cbo_N vs caloriemeter index");
    gr7->GetXaxis()->SetTitle("calo index");
    gr7->GetYaxis()->SetTitle("phi_cbo_N");
    gr7->SetMarkerStyle(20);
    gr7->SetLineColor(kRed);
    gr7->Draw();
    /*   gkband7p=new TGraph(m,n,kband7p);
    gkband7m=new TGraph(m,n,kband7m);
    gkband7m->SetLineColor(kBlue);
    gkband7p->SetLineColor(kBlue);
    gkband7p->Draw("same");
    gkband7m->Draw("same");
    */
    /*   c5=new TCanvas("c4","cbo vs fit start time");
    c5->Divide(2,2);
    c5->cd(1);
    gr9=new TGraphErrors(m,n,A_cbo_N,0,dA_cbo_N);
    gr9->SetTitle("A_cbo_N vs fit start time");
    gr9->GetXaxis()->SetTitle("Start time [ns]");
    gr9->GetYaxis()->SetTitle("A_cbo_N");
    gr9->SetMarkerStyle(20);
    gr9->SetLineColor(kRed);
    gr9->Draw();
    gkband8p=new TGraph(m,n,kband8p);
    gkband8m=new TGraph(m,n,kband8m);
    gkband8m->SetLineColor(kBlue);
    gkband8p->SetLineColor(kBlue);
    gkband8p->Draw("same");
    gkband8m->Draw("same");


     c5->cd(2);
    gr10=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr10->SetTitle("tau_cbo vs fit start time");
    gr10->GetXaxis()->SetTitle("Start time [ns]");
    gr10->GetYaxis()->SetTitle("tau_cbo");
    gr10->SetMarkerStyle(20);
    gr10->SetLineColor(kRed);
    gr10->Draw();

     c5->cd(3);
    gr11=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr11->SetTitle("omega_cbo vs fit start time");
    gr11->GetXaxis()->SetTitle("Start time [ns]");
    gr11->GetYaxis()->SetTitle("omega_cbo");
    gr11->SetMarkerStyle(20);
    gr11->SetLineColor(kRed);
    gr11->Draw();

     c5->cd(4);
    gr12=new TGraphErrors(m,n,phi_cbo_N,0,dphi_cbo_N);
    gr12->SetTitle("phi_cbo_N vs fit start time");
    gr12->GetXaxis()->SetTitle("Start time [ns]");
    gr12->GetYaxis()->SetTitle("phi_cbo_N");
    gr12->SetMarkerStyle(20);
    gr12->SetLineColor(kRed);
    gr12->Draw();
    */
    c6=new TCanvas("c6","A,phi,cbo vs caloriemeter index");
    c6->Divide(2,2);
    c6->cd(1);
    gr13=new TGraphErrors(m,n,A_cbo_A,0,dA_cbo_A);
    gr13->SetTitle("A_cbo_A vs fit caloriemeter index");
    gr13->GetXaxis()->SetTitle("calo index");
    gr13->GetYaxis()->SetTitle("A_cbo_A");
    gr13->SetMarkerStyle(20);
    gr13->SetLineColor(kRed);
    gr13->Draw();
    /*  gkband8p=new TGraph(m,n,kband8p);
    gkband8m=new TGraph(m,n,kband8m);
    gkband8m->SetLineColor(kBlue);
    gkband8p->SetLineColor(kBlue);
    gkband8p->Draw("same");
    gkband8m->Draw("same");
    */
     c6->cd(2);
    gr14=new TGraphErrors(m,n,phi_cbo_A,0,dphi_cbo_A);
    gr14->SetTitle("phi_cbo_A vs caloriemeter index");
    gr14->GetXaxis()->SetTitle("calo index");
    gr14->GetYaxis()->SetTitle("phi_cbo_A");
    gr14->SetMarkerStyle(20);
    gr14->SetLineColor(kRed);
    gr14->Draw();
    /*  gkband9p=new TGraph(m,n,kband9p);
    gkband9m=new TGraph(m,n,kband9m);
    gkband9m->SetLineColor(kBlue);
    gkband9p->SetLineColor(kBlue);
    gkband9p->Draw("same");
    gkband9m->Draw("same");
    */
     c6->cd(3);
    gr15=new TGraphErrors(m,n,A_cbo_phi,0,dA_cbo_phi);
    gr15->SetTitle("A_cbo_phi vs caloriemeter index");
    gr15->GetXaxis()->SetTitle("calo index");
    gr15->GetYaxis()->SetTitle("A_cbo_phi");
    gr15->SetMarkerStyle(20);
    gr15->SetLineColor(kRed);
    gr15->Draw();
    /*  gkband10p=new TGraph(m,n,kband10p);
    gkband10m=new TGraph(m,n,kband10m);
    gkband10m->SetLineColor(kBlue);
    gkband10p->SetLineColor(kBlue);
    gkband10p->Draw("same");
    gkband10m->Draw("same");
    */
     c6->cd(4);
    gr16=new TGraphErrors(m,n,phi_cbo_phi,0,dphi_cbo_phi);
    gr16->SetTitle("phi_cbo_phi vs caloriemeter index");
    gr16->GetXaxis()->SetTitle("calo index");
    gr16->GetYaxis()->SetTitle("phi_cbo_phi");
    gr16->SetMarkerStyle(20);
    gr16->SetLineColor(kRed);
    gr16->Draw();
    /*  gkband11p=new TGraph(m,n,kband11p);
    gkband11m=new TGraph(m,n,kband11m);
    gkband11m->SetLineColor(kBlue);
    gkband11p->SetLineColor(kBlue);
    gkband11p->Draw("same");
    gkband11m->Draw("same");
    */

     c7=new TCanvas("c7","2cbo vs fit start time");
     c7->Divide(2,2);
    c7->cd(1);
    gr17=new TGraphErrors(m,n,A_2cbo,0,dA_2cbo);
    gr17->SetTitle("A_2cbo vs caloriemeter index");
    gr17->GetXaxis()->SetTitle("calo index");
    gr17->GetYaxis()->SetTitle("A_2cbo");
    gr17->SetMarkerStyle(20);
    gr17->SetLineColor(kRed);
    gr17->Draw();
    /*  gkband20p=new TGraph(m,n,kband20p);
    gkband20m=new TGraph(m,n,kband20m);
    gkband20m->SetLineColor(kBlue);
    gkband20p->SetLineColor(kBlue);
    gkband20p->Draw("same");
    gkband20m->Draw("same");
    */

     c7->cd(2);
    gr18=new TGraphErrors(m,n,tau_2cbo,0,dtau_2cbo);
    gr18->SetTitle("tau_2cbo_ vs caloriemeter index");
    gr18->GetXaxis()->SetTitle("calo index");
    gr18->GetYaxis()->SetTitle("tau_2cbo");
    gr18->SetMarkerStyle(20);
    gr18->SetLineColor(kRed);
    gr18->Draw();
    /*  gkband21p=new TGraph(m,n,kband21p);
    gkband21m=new TGraph(m,n,kband21m);
    gkband21m->SetLineColor(kBlue);
    gkband21p->SetLineColor(kBlue);
    gkband21p->Draw("same");
    gkband21m->Draw("same");
    */

     c7->cd(3);
    gr19=new TGraphErrors(m,n,omega_2cbo,0,domega_2cbo);
    gr19->SetTitle("omega_2cbo vs caloriemeter index");
    gr19->GetXaxis()->SetTitle("calo index");
    gr19->GetYaxis()->SetTitle("omega_2cbo");
    gr19->SetMarkerStyle(20);
    gr19->SetLineColor(kRed);
    gr19->Draw();
    /*  gkband22p=new TGraph(m,n,kband22p);
    gkband22m=new TGraph(m,n,kband22m);
    gkband22m->SetLineColor(kBlue);
    gkband22p->SetLineColor(kBlue);
    gkband22p->Draw("same");
    gkband22m->Draw("same");
    */

     c7->cd(4);
    gr20=new TGraphErrors(m,n,phi_2cbo,0,dphi_2cbo);
    gr20->SetTitle("phi_2cbo vs caloriemeter index");
    gr20->GetXaxis()->SetTitle("calo index");
    gr20->GetYaxis()->SetTitle("phi_2cbo");
    gr20->SetMarkerStyle(20);
    gr20->SetLineColor(kRed);
    gr20->Draw();
    /*  gkband23p=new TGraph(m,n,kband23p);
    gkband23m=new TGraph(m,n,kband23m);
    gkband23m->SetLineColor(kBlue);
    gkband23p->SetLineColor(kBlue);
    gkband23p->Draw("same");
    gkband23m->Draw("same");
    */

    c8=new TCanvas("c8","vw vs fit start time");
    c8->Divide(2,2);
    c8->cd(1);
    gr21=new TGraphErrors(m,n,A_vw,0,dA_vw);
    gr21->SetTitle("A_vw vs caloriemeter index");
    gr21->GetXaxis()->SetTitle("calo index");
    gr21->GetYaxis()->SetTitle("A_vw");
    gr21->SetMarkerStyle(20);
    gr21->SetLineColor(kRed);
    gr21->Draw();
    /*  gkband12p=new TGraph(m,n,kband12p);
    gkband12m=new TGraph(m,n,kband12m);
    gkband12m->SetLineColor(kBlue);
    gkband12p->SetLineColor(kBlue);
    gkband12p->Draw("same");
    gkband12m->Draw("same");
    */

    c8->cd(2);
    gr22=new TGraphErrors(m,n,tau_vw,0,dtau_vw);
    gr22->SetTitle("tau_vw vs caloriemeter index");
    gr22->GetXaxis()->SetTitle("calo index");
    gr22->GetYaxis()->SetTitle("tau_vw");
    gr22->SetMarkerStyle(20);
    gr22->SetLineColor(kRed);
    gr22->Draw();
    /* gkband13p=new TGraph(m,n,kband13p);
    gkband13m=new TGraph(m,n,kband13m);
    gkband13m->SetLineColor(kBlue);
    gkband13p->SetLineColor(kBlue);
    gkband13p->Draw("same");
    gkband13m->Draw("same");
    */

    c8->cd(3);
    gr23=new TGraphErrors(m,n,omega_vw,0,domega_vw);
    gr23->SetTitle("omega_vw vs caloriemeter index");
    gr23->GetXaxis()->SetTitle("calo index");
    gr23->GetYaxis()->SetTitle("omega_vw");
    gr23->SetMarkerStyle(20);
    gr23->SetLineColor(kRed);
    gr23->Draw();
    /*  gkband14p=new TGraph(m,n,kband14p);
    gkband14m=new TGraph(m,n,kband14m);
    gkband14m->SetLineColor(kBlue);
    gkband14p->SetLineColor(kBlue);
    gkband14p->Draw("same");
    gkband14m->Draw("same");
    */

    c8->cd(4);
    gr24=new TGraphErrors(m,n,phi_vw,0,dphi_vw);
    gr24->SetTitle("phi_vw vs caloriemeter index");
    gr24->GetXaxis()->SetTitle("calo index");
    gr24->GetYaxis()->SetTitle("phi_vw");
    gr24->SetMarkerStyle(20);
    gr24->SetLineColor(kRed);
    gr24->Draw();
    /*  gkband15p=new TGraph(m,n,kband15p);
    gkband15m=new TGraph(m,n,kband15m);
    gkband15m->SetLineColor(kBlue);
    gkband15p->SetLineColor(kBlue);
    gkband15p->Draw("same");
    gkband15m->Draw("same");
    */

    c9=new TCanvas("c9","y vs fit start time");
    c9->Divide(2,2);
    c9->cd(1);
    gr25=new TGraphErrors(m,n,A_y,0,dA_y);
    gr25->SetTitle("A_y vs caloriemeter index");
    gr25->GetXaxis()->SetTitle("calo index");
    gr25->GetYaxis()->SetTitle("A_y");
    gr25->SetMarkerStyle(20);
    gr25->SetLineColor(kRed);
    gr25->Draw();
    /*  gkband16p=new TGraph(m,n,kband16p);
    gkband16m=new TGraph(m,n,kband16m);
    gkband16m->SetLineColor(kBlue);
    gkband16p->SetLineColor(kBlue);
    gkband16p->Draw("same");
    gkband16m->Draw("same");
    */

    c9->cd(2);
    gr26=new TGraphErrors(m,n,tau_y,0,dtau_y);
    gr26->SetTitle("tau_y vs caloriemeter index");
    gr26->GetXaxis()->SetTitle("calo index");
    gr26->GetYaxis()->SetTitle("tau_y");
    gr26->SetMarkerStyle(20);
    gr26->SetLineColor(kRed);
    gr26->Draw();
    /*  gkband17p=new TGraph(m,n,kband17p);
    gkband17m=new TGraph(m,n,kband17m);
    gkband17m->SetLineColor(kBlue);
    gkband17p->SetLineColor(kBlue);
    gkband17p->Draw("same");
    gkband17m->Draw("same");
    */

    c9->cd(3);
    gr27=new TGraphErrors(m,n,omega_y,0,domega_y);
    gr27->SetTitle("omega_y vs caloriemeter index");
    gr27->GetXaxis()->SetTitle("calo index");
    gr27->GetYaxis()->SetTitle("omega_y");
    gr27->SetMarkerStyle(20);
    gr27->SetLineColor(kRed);
    gr27->Draw();
    /* gkband18p=new TGraph(m,n,kband18p);
    gkband18m=new TGraph(m,n,kband18m);
    gkband18m->SetLineColor(kBlue);
    gkband18p->SetLineColor(kBlue);
    gkband18p->Draw("same");
    gkband18m->Draw("same");
    */

    c9->cd(4);
    gr28=new TGraphErrors(m,n,phi_y,0,dphi_y);
    gr28->SetTitle("phi_y vs caloriemeter index");
    gr28->GetXaxis()->SetTitle("calo index");
    gr28->GetYaxis()->SetTitle("phi_y");
    gr28->SetMarkerStyle(20);
    gr28->SetLineColor(kRed);
    gr28->Draw();
    /*  gkband19p=new TGraph(m,n,kband19p);
    gkband19m=new TGraph(m,n,kband19m);
    gkband19m->SetLineColor(kBlue);
    gkband19p->SetLineColor(kBlue);
    gkband19p->Draw("same");
    gkband19m->Draw("same");
    */

   
}
