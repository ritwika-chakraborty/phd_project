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
char root_file_name[128] = "r2_w8_calos.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
//char root_file_name[128] = "run3BtoM_thresh300_nofbfDQC_wndw_8_2.root";
char *histname = new char[10];
char *histname2 = new char[10];
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res, *hr_sum, *hcov_1_p, *hcov_1_m,*hcov_13_p,*hcov_13_m,*hcov_14_p, *hcov_14_m, *hcov_15_p, *hcov_15_m, *hcov_27_p, *hcov_27_m,*hcov_28_p, *hcov_28_m, *hcov_29_p, *hcov_29_m, *hr_noshift, *hr_shift, *hcomp[25],*hcomp_sum,*htemp, *hrm1[50],*hrm2[50],*hrm3[50],*hrm4[50];
TCanvas *c,*c1, *c2,*c3, *cauto;
TGraphErrors *gr1;
TH1 *hfft, *hFFTsq, *hAuto;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=30000;
Double_t fit_stop=300000;
//Double_t fit_stop=306076.25;
Double_t rawBinToNs = 1.25;
Double_t inject_time;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
Int_t m=0;
Double_t blindR[25],dblindR[25],n[25];
Int_t countfcn=0;

TFile *foutput;
TH2D *hcov[75];
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
//TMatrixD cov2(mdim,mdim);
TMatrixD cov(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim)),data2((mdim)*(mdim));

TArrayD ratio_cov_matrix(TH1D *hr_sum, TH1D *h_sum, Double_t fit_start){

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
  cout<<"i0 is"<<i0<<" "<<endl;
  cout<<"check "<<h_sum->GetBinError(250)<<endl;
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
	    data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	    //data[i]=1;
	    //    data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	    /*	    if(data[i]==0)
	      {
		cout<<"diag ele is zero at "<<ir<<" "<<ic;
		}*/
	    // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

	    	if(ir==ic-nbinshift)
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
	    //  data[i]=covar;
	     data[i]=(-0.675)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //    data[i]=(-2.0/3.0);
	    //data[i]=covar;
	       //   cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

		if(ir==ic+nbinshift)
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


	     data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	    //   data[i]=covar;
	    data[i]=(-0.675)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    // data[i]=(-2.0/3.0);
	    //data[i]=covar;
	  }

		if(ir==ic-(2*nbinshift))
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
        
		//data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	        //cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		//	  data[i]=covar;
			data[i]=(0.175)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	data[i]=(1.0/6.0);
		//data[i]=covar;
	      }

		if(ir==ic+(2*nbinshift))
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
	        
		//	data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		//	data[i]=covar;
		 	data[i]=(0.175)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	data[i]=(1.0/6.0);
		//data[i]=covar;	
	      }
	 
	
      }

    cout<<mdim<<endl;
    //cout<<"Matrix inverted "<<i0<<endl;
  return data;

}


TH1D* fast_rotation_smoothing(TH1D *h_sum)
{
    hr_noshift=(TH1D*)h_sum->Clone();
	
    hr_shift=(TH1D*)hr_noshift->Clone();
    hr_shift->Reset();

    for(int ibin=1;ibin<=hr_shift->GetNbinsX();ibin++)
       {
	hr_shift->SetBinContent(ibin,hr_noshift->GetBinContent(ibin+1));
	hr_shift->SetBinError(ibin,hr_noshift->GetBinError(ibin+1));
       }

    hcomp_sum=(TH1D*)h_sum->Clone();
    hcomp_sum->Reset();

    for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
       {
	hcomp_sum->SetBinContent(ibin,0.5*(hr_noshift->GetBinContent(ibin)+hr_shift->GetBinContent(ibin)));
	hcomp_sum->SetBinError(ibin,0.5*sqrt((hr_noshift->GetBinError(ibin)*hr_noshift->GetBinError(ibin))+(hr_shift->GetBinError(ibin)*hr_shift->GetBinError(ibin))));
       }

    return hcomp_sum;

}


TH1D* construct_rhist_copy(TH1D *h_sum, TH1D *ho_sum)
{
 
  //h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
  //h_sum->Scale(0.125);
   
  cout<<"check "<<h_sum->GetBinError(500)<<endl;
   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

   nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
      /*     if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   //nbinshift=15;
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
      h_ratio=new TH1D("calo_histogram_sum_ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo_histogram_sum_numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
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
       long double et1,det1,et2,det2,et3,det3,et4,det4,et5,det5,et6,det6,et7,det7,et8,det8,et9,det9;
       long double t1,t2,t3,t4,t5,t6,t7,t8,t9;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
	et1=ho_sum->GetBinContent(2*ibin);
        det1=ho_sum->GetBinError(2*ibin);
	et2=ho_sum->GetBinContent((2*ibin)+1);
	det2=ho_sum->GetBinError((2*ibin)+1);
	et3=ho_sum->GetBinContent((2*ibin)-1);
	det3=ho_sum->GetBinError((2*ibin)-1);
	et4=ho_sum->GetBinContent((2*ibin)+27);
	det4=ho_sum->GetBinError((2*ibin))+27;
	et5=ho_sum->GetBinContent((2*ibin)-27);
	det5=ho_sum->GetBinError((2*ibin)-27);
	et6=ho_sum->GetBinContent((2*ibin)+28);
	det6=ho_sum->GetBinError((2*ibin)+28);
	et7=ho_sum->GetBinContent((2*ibin)-28);
	det7=ho_sum->GetBinError((2*ibin)-28);
	et8=ho_sum->GetBinContent((2*ibin)+29);
	det8=ho_sum->GetBinError((2*ibin)+29);
	et9=ho_sum->GetBinContent((2*ibin)-29);
	det9=ho_sum->GetBinError((2*ibin)-29);

	
	t1=((8*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det1;

	t2=((4*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det2;
	
	t3=((4*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det3;
	
	t4=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det4;

	t5=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det5;

        t6=-((8*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det6;

	t7=-((8*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det7;

	t8=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det8;

	t9=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det9;
 	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt( ( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9)) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( ( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9))  ) );
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

 // hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //hcalo->Scale(0.125);

 /* foutput=new TFile("ratio_4hist_copy.root","new");
 h_ratio->Write();
 foutput->Close();
 */
 
 return hcalo;
  }


void ratio_noFR_cov_construct_percalo_FRcorr()
{

    _file[1]=TFile::Open(root_file_name);

    // ############# My root file 

    /*   for(int i=1;i<=24;i++)
    {
     h[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }
    */
    // ############### Tim's root file

      for(int i=1;i<=24;i++)
    {
     hrm1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_0_0", i)));
     hrm2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_1_0", i)));
     hrm3[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_2_0", i)));
     hrm4[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_3_0", i)));
    }
    

    /* hrm1[1]= (TH1D*)(_file[1]->FindObjectAny("qHist1D_sig_rm_sum_0_0"));
     hrm2[1]= (TH1D*)(_file[1]->FindObjectAny("qHist1D_sig_rm_sum_1_0"));
     hrm3[1]= (TH1D*)(_file[1]->FindObjectAny("qHist1D_sig_rm_sum_2_0"));
     hrm4[1]= (TH1D*)(_file[1]->FindObjectAny("qHist1D_sig_rm_sum_3_0"));
    */


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
    
    
  
  /*h[24]->Fit("gaus","","",104750,104850);
  TF1 *fit = h[24]->GetFunction("gaus");
  inject_time=fit->GetParameter(1);
  */
  inject_time=104800;
  
  h_sum= new TH1D("calo_histogram_sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
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
   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

    
   h_sum->Rebin(4);
   h_sum->Scale(0.25);

	for(int i=1;i<=24;i++)
	  {
	   h[i]->Rebin(4);
	   h[i]->Scale(0.25);

     
	   hcomp[i]=fast_rotation_smoothing(h[i]);

	    cout<<"binwidth check after smoothing "<<hcomp[i]->GetBinWidth(500)<<" "<<h[i]->GetBinWidth(500)<<endl;
            //hcomp[i]=(TH1D*)h[i]->Clone();
	    
            hcomp[i]->Rebin(2);
            hcomp[i]->Scale(0.5);

	    	  
            sprintf(histname2,"hcalo_ratio_%d",i);
	
           hr[i]= new TH1D(histname2, histname2, hcomp[i]->GetNbinsX(),hcomp[i]->GetBinLowEdge(1),hcomp[i]->GetBinLowEdge(hcomp[i]->GetNbinsX())+hcomp[i]->GetBinWidth(1));

	   cout<<"binwidth check before ratio construction "<<hcomp[i]->GetBinWidth(500)<<" "<<h[i]->GetBinWidth(500)<<endl;

	   htemp = construct_rhist_copy(hcomp[i],h[i]);

	   //h[i]->Rebin(8);
	   //h[i]->Scale(0.125);

	   for(int ibin=1;ibin<=htemp->GetNbinsX();ibin++)
	     {
	       hr[i]->SetBinContent(ibin,htemp->GetBinContent(ibin));
	       hr[i]->SetBinError(ibin,htemp->GetBinError(ibin));
	     }
           htemp->Reset();
           cout<<"check "<<h[i]->GetBinError(500)<<endl;

	  
	   hr[i]->GetYaxis()->SetTitle("ADC counts");
	   hr[i]->GetXaxis()->SetTitle("time [ns]");
           hr[i]->GetXaxis()->SetRangeUser(fit_start,fit_stop);

  
 
	   h[i]->Rebin(2);
	   h[i]->Scale(0.5);

	   //	   cov2.ResizeTo(hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), -1);
   

	   // cov2.SetTol(1.e-23);

	   cout<<"binwidth check before covariance calculation "<<hr[i]->GetBinWidth(500)<<" "<<h[i]->GetBinWidth(500)<<endl;

           data2=ratio_cov_matrix(hr[i],h[i],fit_start);

	   

	   TMatrixD cov2(mdim,mdim);
		
           cout<<"before cov[1000][1000] "<<cov2[1000][1000]<<endl;
           cov2.SetMatrixArray(data2.GetArray());
           cout<<"after cov[1000][1000] "<<cov2[1000][1000]<<endl;

           //cov2.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
	   cov2.ResizeTo(hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), -1);
   
	   // cov2.Print();
           cov2.SetTol(1.e-50);

	   Double_t det1; 
	   cov2.Invert(&det1);
	   cout<<"covariance matrix inverted "<<endl;


	   cout<<"calo "<<i<<endl; 


	   sprintf(histname,"hcov_%d",i);
	
           hcov[i]= new TH2D(histname, histname, 2100, 0.0, 2100.0, 2100, 0.0, 2100.0);


	    	      
           for(int irow=hr[i]->FindBin(fit_start); irow<=hr[i]->FindBin(fit_stop); irow++)
           {
	    for(int icol=hr[i]->FindBin(fit_start); icol<=hr[i]->FindBin(fit_stop); icol++)
	    {
	     hcov[i]->SetBinContent(irow,icol,cov2(irow,icol));
	    }
           }

	    
	   cov2.Clear();
	// cov.Clear("C");
        //cov.UnitMatrix();
 
   }
  
       
    foutput=new TFile("run2_4hcopy_ratio_noFR_percalo_cov_mat_FRcorr_tim_inverse.root","new");
   for(int i=1;i<=24;i++)
     {
      hr[i]->Write();
      hcov[i]->Write();
     }
   foutput->Close();
     
   
}
