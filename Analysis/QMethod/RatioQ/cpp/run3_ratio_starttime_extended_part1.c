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
char root_file_name[128] = "run3_thresh300_nofbfDQC_wndw_8.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed.root";
TVirtualFitter *gFitter;
TFile *_file[25];

TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res_1[50], *hr_sum;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28;

TH1 *hfft_1[50];
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
TH1D *chisq_1, *A_1, *blindR_1, *phi_0_1, *A_cbo_N_1, *tau_cbo_1, *omega_cbo_1, *phi_cbo_N_1, *A_cbo_A_1, *phi_cbo_A_1, *A_cbo_phi_1, *phi_cbo_phi_1, *A_vw_1, *tau_vw_1, *omega_vw_1, *phi_vw_1, *A_y_1, *tau_y_1, *omega_y_1, *phi_y_1, *A_2cbo_1, *tau_2cbo_1, *omega_2cbo_1, *phi_2cbo_1,*n_1;
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
	       data[i]=covar;
	    // data[i]=(-2.0/3.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //    data[i]=(-2.0/3.0);
	       //     cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
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

	    //  data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    
	     data[i]=covar;
	     //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	    //  data[i]=(-2.0/3.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	    //   data[i]=(-2.0/3.0);
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

		// 	data[i]=covar/sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	        
	       	  data[i]=covar;
		  //	cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		//	data[i]=(1.0/6.0)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
		//	data[i]=(1.0/6.0);
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

    double asym_2cbo= par[19];

    double tau_2cbo = par[20];

    double omega_2cbo = par[21];

    double phi_2cbo = par[22];
    
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


void run3_ratio_starttime_extended_part1()
{

  _file[1]=TFile::Open(root_file_name);

  for(int i=1;i<=24;i++)
    {
     h[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }



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



   A_1=new TH1D("A_1","A_1",3,0,3);
   blindR_1=new TH1D("blindR_1","blindR_1",3,0,3);
   phi_0_1=new TH1D("phi_0_1","phi_0_1",3,0,3);
   A_cbo_N_1=new TH1D("A_cbo_N_1","A_cbo_N_1",3,0,3);
   tau_cbo_1=new TH1D("tau_cbo_1","tau_cbo_1",3,0,3);
   omega_cbo_1=new TH1D("omega_cbo_1","omega_cbo_1",3,0,3);
   phi_cbo_N_1=new TH1D("phi_cbo_N_1","phi_cbo_N_1",3,0,3);
   A_cbo_A_1=new TH1D("A_cbo_A_1","A_cbo_A_1",3,0,3);
   phi_cbo_A_1=new TH1D("phi_cbo_A_1","phi_cbo_A_1",3,0,3);
   A_cbo_phi_1=new TH1D("A_cbo_phi_1","A_cbo_phi_1",3,0,3);
   phi_cbo_phi_1=new TH1D("phi_cbo_phi_1","phi_cbo_phi_1",3,0,3);
   A_vw_1=new TH1D("A_vw_1","A_vw_1",3,0,3);
   tau_vw_1=new TH1D("tau_vw_1","tau_vw_1",3,0,3);
   omega_vw_1=new TH1D("omega_vw_1","omega_vw_1",3,0,3);
   phi_vw_1=new TH1D("phi_vw_1","phi_vw_1",3,0,3);
   A_y_1=new TH1D("A_y_1","A_y_1",3,0,3);
   tau_y_1=new TH1D("tau_y_1","tau_y_1",3,0,3);
   omega_y_1=new TH1D("omega_y_1","omega_y_1",3,0,3);
   phi_y_1=new TH1D("phi_y_1","phi_y_1",3,0,3);
   A_2cbo_1=new TH1D("A_2cbo_1","A_2cbo_1",3,0,3);
   tau_2cbo_1=new TH1D("tau_2cbo_1","tau_2cbo_1",3,0,3);
   omega_2cbo_1=new TH1D("omega_2cbo_1","omega_2cbo_1",3,0,3);
   phi_2cbo_1=new TH1D("phi_2cbo_1","phi_2cbo_1",3,0,3);
   n_1=new TH1D("n_1","n_1",3,0,3);
   
   chisq_1=new TH1D("chisq_1","chisq_1",3,0,3);




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


       fit_func23->SetParameter(0, 0.2285);
       fit_func23->SetParameter(1, 0.0);
       fit_func23->SetParameter(2, 2.23);
       fit_func23->SetParameter(3, -0.00127);
       //fit_func23->SetParameter(4, 191200);
       //fit_func23->SetParameter(5, 0.002331);
       fit_func23->FixParameter(4, 191153);
       fit_func23->FixParameter(5, 0.00233128);
       fit_func23->SetParameter(6, -2.71);
       fit_func23->SetParameter(7, 0.00033);
       fit_func23->SetParameter(8, 38.3);
       fit_func23->SetParameter(9, -0.0000635);
       fit_func23->SetParameter(10, -29.3);
       fit_func23->SetParameter(11, 0.00487);
       //fit_func23->SetParameter(12, 22690);
       //fit_func23->SetParameter(13, 0.014063);
       fit_func23->FixParameter(12, 22688.5);
       fit_func23->FixParameter(13,  0.0140631);
       fit_func23->SetParameter(14, 4.7);
       fit_func23->SetParameter(15, -0.00064);
       //fit_func23->SetParameter(16, 74670);
       //fit_func23->SetParameter(17, 0.013895);
       fit_func23->FixParameter(16, 74647.8);
       fit_func23->FixParameter(17, 0.0138950);
       fit_func23->SetParameter(18, -2.888);
       fit_func23->SetParameter(19, -0.000031);
       fit_func23->FixParameter(20, 95639);
       fit_func23->FixParameter(21, 0.004662);
       fit_func23->SetParameter(22, -0.45);

     

     /*   fit_func7->SetParameter(0, 0.2291);
      fit_func7->SetParameter(1, 0.0);
      fit_func7->SetParameter(2, 2.22);
      fit_func7->SetParameter(3, 0.025);
      fit_func7->SetParameter(4, 257000);
      fit_func7->SetParameter(5, 0.0023406);
      fit_func7->SetParameter(6, 0.45);
     */
       /*    fit_func7->SetParameter(0, 0.2291);
        fit_func7->SetParameter(1, 0.0);
        fit_func7->SetParameter(2, 2.224);
        fit_func7->SetParameter(3, 0.00256);
        fit_func7->SetParameter(4, 221100);
        fit_func7->SetParameter(5, 0.00234032);
        fit_func7->SetParameter(6, 0.464);
       */
       
     fit_func23->SetNpx(1000000);

     hr_sum = construct_rhist_copy(h_sum);


   
     
     for(int i=1; i<=3;i++)
       {

	cout<<"fit start time"<<fit_start<<" ns"<<endl; 
        
	 
	//	fit_func3->SetParameters(0.22, 0.0, 2.22);

  
	hr_sum->GetYaxis()->SetTitle("ADC counts");
	hr_sum->GetXaxis()->SetTitle("time [ns]");
        hr_sum->GetXaxis()->SetRangeUser(fit_start,fit_stop);




		
 	cov2.ResizeTo(h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), h_sum->FindBin(fit_start), h_sum->FindBin(fit_stop), -1);
   

        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr_sum,h_sum, fit_start);

   
          

       fit_func23->SetParameter(0, fit_func23->GetParameter(0));
       fit_func23->SetParameter(1, fit_func23->GetParameter(1));
       fit_func23->SetParameter(2, fit_func23->GetParameter(2));
       fit_func23->SetParameter(3, fit_func23->GetParameter(3));
       fit_func23->FixParameter(4, fit_func23->GetParameter(4));
       fit_func23->FixParameter(5, fit_func23->GetParameter(5));
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
	


	

	/*   fit_func7->SetParameter(0, fit_func7->GetParameter(0));
       fit_func7->SetParameter(1, fit_func7->GetParameter(1));
       fit_func7->SetParameter(2, fit_func7->GetParameter(2));
       fit_func7->SetParameter(3, fit_func7->GetParameter(3));
       fit_func7->SetParameter(4, fit_func7->GetParameter(4));
       fit_func7->SetParameter(5, fit_func7->GetParameter(5));
       fit_func7->SetParameter(6, fit_func7->GetParameter(6));
	*/ 
       
       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hr_sum);
       gFitter->SetFCN(chi2);  
	
       hr_sum->Fit("fprec23","RUE","",fit_start,fit_stop);
       countfcn=0;
       gStyle->SetOptFit(1111);

          gettimeofday(&t_end, NULL);
       printf("QRatio fit duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
     
       //cout<<" setting ch to zero "<<endl;
       //ch=0;
 
       A_1->SetBinContent(i,fit_func23->GetParameter(0));
       blindR_1->SetBinContent(i,fit_func23->GetParameter(1));
       phi_0_1->SetBinContent(i,fit_func23->GetParameter(2));
       A_cbo_N_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(3)));
       tau_cbo_1->SetBinContent(i,fit_func23->GetParameter(4));
       omega_cbo_1->SetBinContent(i,fit_func23->GetParameter(5));
       phi_cbo_N_1->SetBinContent(i,fit_func23->GetParameter(6));
       A_cbo_A_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(7)));
       phi_cbo_A_1->SetBinContent(i,fit_func23->GetParameter(8));
       A_cbo_phi_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(9)));
       phi_cbo_phi_1->SetBinContent(i,fit_func23->GetParameter(10));
       A_vw_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(11)));
       tau_vw_1->SetBinContent(i,fit_func23->GetParameter(12));
       omega_vw_1->SetBinContent(i,fit_func23->GetParameter(13));
       phi_vw_1->SetBinContent(i,fit_func23->GetParameter(14));
       A_y_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(15)));
       tau_y_1->SetBinContent(i,fit_func23->GetParameter(16));
       omega_y_1->SetBinContent(i,fit_func23->GetParameter(17));
       phi_y_1->SetBinContent(i,fit_func23->GetParameter(18));
       A_2cbo_1->SetBinContent(i,TMath::Abs(fit_func23->GetParameter(19)));
       tau_2cbo_1->SetBinContent(i,fit_func23->GetParameter(20));
       omega_2cbo_1->SetBinContent(i,fit_func23->GetParameter(21));
       phi_2cbo_1->SetBinContent(i,fit_func23->GetParameter(22));

     
       chisq_1->SetBinContent(i,fit_func23->GetChisquare()/fit_func23->GetNDF());


       while(phi_cbo_N_1->GetBinContent(i) < 0.0)
	 {
	   phi_cbo_N_1->SetBinContent(i, phi_cbo_N_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }

       while(phi_cbo_A_1->GetBinContent(i) < 0.0)
	 {
	   phi_cbo_A_1->SetBinContent(i, phi_cbo_A_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }

       while(phi_cbo_phi_1->GetBinContent(i) < 0.0)
	 {
	   phi_cbo_phi_1->SetBinContent(i, phi_cbo_phi_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }
       
       while(phi_vw_1->GetBinContent(i) < 0.0)
	 {
	   phi_vw_1->SetBinContent(i, phi_vw_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }

       while(phi_y_1->GetBinContent(i) < 0.0)
	 {
	   phi_y_1->SetBinContent(i, phi_y_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }

       while(phi_2cbo_1->GetBinContent(i) < 0.0)
	 {
	   phi_2cbo_1->SetBinContent(i, phi_2cbo_1->GetBinContent(i) + (2*(TMath::Pi())));
	 }

         while(phi_cbo_N_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_cbo_N_1->SetBinContent(i, phi_cbo_N_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }

       while(phi_cbo_A_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_cbo_A_1->SetBinContent(i, phi_cbo_A_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }

       while(phi_cbo_phi_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_cbo_phi_1->SetBinContent(i, phi_cbo_phi_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }

       while(phi_vw_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_vw_1->SetBinContent(i, phi_vw_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }

       while(phi_y_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_y_1->SetBinContent(i, phi_y_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }

       while(phi_2cbo_1->GetBinContent(i) >  (2*(TMath::Pi())))
	 {
	   phi_2cbo_1->SetBinContent(i, phi_2cbo_1->GetBinContent(i) - (2*(TMath::Pi())));
	 }
       
       A_1->SetBinError(i,fit_func23->GetParError(0));
       blindR_1->SetBinError(i,fit_func23->GetParError(1));
       phi_0_1->SetBinError(i,fit_func23->GetParError(2));
       A_cbo_N_1->SetBinError(i,fit_func23->GetParError(3));
       tau_cbo_1->SetBinError(i,fit_func23->GetParError(4));
       omega_cbo_1->SetBinError(i,fit_func23->GetParError(5));
       phi_cbo_N_1->SetBinError(i,fit_func23->GetParError(6));
       A_cbo_A_1->SetBinError(i,fit_func23->GetParError(7));
       phi_cbo_A_1->SetBinError(i,fit_func23->GetParError(8));
       A_cbo_phi_1->SetBinError(i,fit_func23->GetParError(9));
       phi_cbo_phi_1->SetBinError(i,fit_func23->GetParError(10));
       A_vw_1->SetBinError(i,fit_func23->GetParError(11));
       tau_vw_1->SetBinError(i,fit_func23->GetParError(12));
       omega_vw_1->SetBinError(i,fit_func23->GetParError(13));
       phi_vw_1->SetBinError(i,fit_func23->GetParError(14));
       A_y_1->SetBinError(i,fit_func23->GetParError(15));
       tau_y_1->SetBinError(i,fit_func23->GetParError(16));
       omega_y_1->SetBinError(i,fit_func23->GetParError(17));
       phi_y_1->SetBinError(i,fit_func23->GetParError(18));
       A_2cbo_1->SetBinError(i,fit_func23->GetParError(19));
       tau_2cbo_1->SetBinError(i,fit_func23->GetParError(20));
       omega_2cbo_1->SetBinError(i,fit_func23->GetParError(21));
       phi_2cbo_1->SetBinError(i,fit_func23->GetParError(22));

     
       n_1->SetBinContent(i,fit_start);
       m=m+1;
      		
       h_res_1[i]= new TH1D("residual histogram 1 ", "h_res_1", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_sum->GetBinContent(ibin)- fit_func23->Eval( hr_sum->GetBinCenter(ibin) ) );
          if(hr_sum->GetBinError(ibin)!=0){res=(res/hr_sum->GetBinError(ibin));}
          h_res_1[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res_1[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res_1[i]->GetXaxis()->SetTitle("time [ns]");
	h_res_1[i]->GetXaxis()->SetLabelSize(0.05);
        h_res_1[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft_1[i] = h_res_1[i]->FFT(hfft_1[i], "MAG");
       hfft_1[i]->SetLineColor(kBlack);
       hfft_1[i]->SetBins(hr_sum->GetNbinsX(),0,1000/h_res_1[i]->GetBinWidth(1));
       hfft_1[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res_1[i]->GetBinWidth(1)));
       hfft_1[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft_1[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft_1[i]->GetXaxis()->SetLabelSize(0.05);
       hfft_1[i]->GetYaxis()->SetLabelSize(0.05);


       fit_start=fit_start+1000;
 
   }
     
  c2 = new TCanvas("c2","start time residuals");
  c2->Divide(3,1);
  c3 = new TCanvas("c3","start time fft");
  c3->Divide(3,1);
   
 for(int i=1;i<=3;i++)
   {
     c2->cd(i);
     h_res_1[i]->Draw();
     c3->cd(i);
     hfft_1[i]->Draw();
   }
 
  foutput=new TFile("run3_starttime_output_corr_fixed_extended_part1.root","new");
  
   A_1->Write();
   blindR_1->Write();
   phi_0_1->Write();
   A_cbo_N_1->Write();
   tau_cbo_1->Write();
   omega_cbo_1->Write();
   phi_cbo_N_1->Write();
   A_cbo_A_1->Write();
   phi_cbo_A_1->Write();
   A_cbo_phi_1->Write();
   phi_cbo_phi_1->Write();
   A_vw_1->Write();
   tau_vw_1->Write();
   omega_vw_1->Write();
   phi_vw_1->Write();
   A_y_1->Write();
   tau_y_1->Write();
   omega_y_1->Write();
   phi_y_1->Write();
   A_2cbo_1->Write();
   tau_2cbo_1->Write();
   omega_2cbo_1->Write();
   phi_2cbo_1->Write();
   n_1->Write();
   foutput->Close();
 


   
}
