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
//char root_file_name[128] = "run2g_corrOOF_thresh400.root";
char root_file_name[128] = "run2_thresh300_fbfDQC.root";
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[25], *hr_noshift[25], *hr_shift[25];
TCanvas *c,*c1, *c2,*c3,*c4;
TGraphErrors *gr1,*gr2;
TH1 *hfft[25];
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
Int_t m=0;
Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25];
Int_t countfcn=0;
Double_t ch;

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
  cout<<"i0 is"<<i0<<" "<<endl;
  
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
	     
	     if(ir==247){cout<<i<<" "<<ir<<" "<<ic<<" "<<data[i]<<endl;}
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
    // cov.Invert(&det1);
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
 
   //   h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
   //  h_sum->Scale(0.125);
   

   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

      nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
     if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
   
   //   nbinshift=20;
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
       
      long double t1,t2,t3;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
	t1=4*(h1->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	t2=-4*(hp->GetBinError(ibin))*(h1->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	t3=-4*(hm->GetBinError(ibin))*(h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
       if(usesetbins)
       {
	h_ratio->SetBinError( ibin, sqrt( (t1*t1) + (t2*t2) + (t3*t3) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	h_num->SetBinError(ibin,  sqrt( (t1*t2) + (t2*t2) + (t3*t3) ) );
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

 hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 hcalo->Scale(0.125);

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

void ratio_percalo_fit_FR()
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
     fit_func23->SetParName(19,"A_2cbo");
     fit_func23->SetParName(20,"Tau_2cbo");
     fit_func23->SetParName(21,"omega_2cbo");
     fit_func23->SetParName(22,"phi_2cbo");
     fit_func23->SetNpx(1000000);

        fit_func23->SetParameter(0, 0.2288);
       fit_func23->SetParameter(1, 0.0);
       fit_func23->SetParameter(2, 2.268);
       fit_func23->SetParameter(3, 0.002426);
       fit_func23->SetParameter(4, 254739);
       fit_func23->SetParameter(5, 0.00234065);
       fit_func23->SetParameter(6, 0.552);
       fit_func23->SetParameter(7, 0.0007);
       fit_func23->SetParameter(8, 6.7);
       fit_func23->SetParameter(9, -0.00007);
       fit_func23->SetParameter(10, -19.5);
       fit_func23->SetParameter(11, 0.008);
       fit_func23->SetParameter(12, 25090);
       fit_func23->SetParameter(13, 0.01403);
       fit_func23->SetParameter(14, 2.55);
       fit_func23->SetParameter(15, 0.00022);
       fit_func23->SetParameter(16, 335000);
       fit_func23->SetParameter(17, 0.01393);
       fit_func23->SetParameter(18, 6.27);
       fit_func23->SetParameter(19, 0.00012);
       fit_func23->SetParameter(20, 158000);
       fit_func23->SetParameter(21, 0.00469);
       fit_func23->SetParameter(22, 2.9);



     

     
     for(int i=1; i<=1;i++)
       {

	cout<<"Calo "<<i<<endl; 
        hr_noshift[i] = construct_rhist_copy(h[i]);

	hr_shift[i]=(TH1D*)hr_noshift[i]->Clone();
	hr_shift[i]->Reset();

	for(int ibin=1;ibin<=hr_shift[i]->GetNbinsX();ibin++)
	  {
	    hr_shift[i]->SetBinContent(ibin,hr_noshift[i]->GetBinContent(ibin+1));
	    hr_shift[i]->SetBinError(ibin,hr_noshift[i]->GetBinError(ibin+1));
	  }

	hr[i]=(TH1D*)hr_shift[i]->Clone();
	hr[i]->Reset();

	for(int ibin=1;ibin<=hr_shift[i]->GetNbinsX();ibin++)
	  {
	    hr[i]->SetBinContent(ibin,0.5*(hr_noshift[i]->GetBinContent(ibin)+hr_shift[i]->GetBinContent(ibin)));
	    hr[i]->SetBinError(ibin,0.5*(hr_noshift[i]->GetBinError(ibin)+hr_shift[i]->GetBinError(ibin)));
	  }

	
	hr[i]->Rebin(2);
	hr[i]->Scale(0.5);
	h[i]->Rebin(2);
	h[i]->Scale(0.5);

	 
	//	fit_func3->SetParameters(0.22, 0.0, 2.22);

  
	hr[i]->GetYaxis()->SetTitle("ADC counts");
	hr[i]->GetXaxis()->SetTitle("time [ns]");
        hr[i]->GetXaxis()->SetRangeUser(fit_start,fit_stop);



        printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hcalo);
       gFitter->SetFCN(chi2);  
		
 	cov2.ResizeTo(h[i]->FindBin(fit_start), h[i]->FindBin(fit_stop), h[i]->FindBin(fit_start), h[i]->FindBin(fit_stop), -1);
   

        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr[i],h[i], fit_start);

   
          

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
       fit_func23->SetParameter(19, fit_func23->GetParameter(19));
       fit_func23->SetParameter(20, fit_func23->GetParameter(20));
       fit_func23->SetParameter(21, fit_func23->GetParameter(21));
       fit_func23->SetParameter(22, fit_func23->GetParameter(22));


      
       hr[i]->Fit("fprec23","RE","",fit_start,fit_stop);
       gStyle->SetOptFit(1111);

       blindR[m]=fit_func23->GetParameter(1);
       phase[m]=fit_func23->GetParameter(2);
       //blindR[m]=fit_func19->GetChisquare()/fit_func19->GetNDF();
       dblindR[m]=fit_func23->GetParError(1);
       dphase[m]=fit_func23->GetParError(2);
       n[m]=i;
       m=m+1;
   
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hr[i]->GetNbinsX(), hr[i]->GetBinLowEdge(1), (hr[i]->GetBinLowEdge(hr[i]->GetNbinsX())+hr[i]->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr[i]->GetBinContent(ibin)- fit_func19->Eval( hr[i]->GetBinCenter(ibin) ) );
          if(hr[i]->GetBinError(ibin)!=0){res=(res/hr[i]->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("ADC counts");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");	
      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hr[i]->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.04);
       hfft[i]->GetYaxis()->SetLabelSize(0.04);

 
   }

  c2 = new TCanvas("c2","per calo residuals");
  c2->Divide(4,6);
  c3 = new TCanvas("c3","per calo fft");
  c3->Divide(4,6);
   
 for(int i=1;i<=24;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hfft[i]->Draw();
   }
    c1=new TCanvas("c1","Ratio fit blind R vs caloriemeter #");
    gr1=new TGraphErrors(m,n,blindR,0,dblindR);
    gr1->SetTitle("Ratio blind R vs calo");
    gr1->GetXaxis()->SetTitle("calo #");
    gr1->GetYaxis()->SetTitle("R [ppm]");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->GetYaxis()->SetRangeUser(-75,-15);
    //  gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr1->Draw();

    c4=new TCanvas("c4","Ratio fit phase vs caloriemeter #");
    gr2=new TGraphErrors(m,n,phase,0,dphase);
    gr2->SetTitle("g-2 phase vs calo");
    gr2->GetXaxis()->SetTitle("calo #");
    gr2->GetYaxis()->SetTitle("phi");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    //  gr1->GetYaxis()->SetRangeUser(-75,-15);
    //  gr1->GetYaxis()->SetRangeUser(0.8,1.3);
    gr2->Draw();

   
}
