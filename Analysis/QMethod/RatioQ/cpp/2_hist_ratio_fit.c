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
//char root_file_name[128] = "run2g_corrOOF_thresh400.root";
char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res, *hr_sum;
TCanvas *c,*c1, *c2,*c3;
TGraphErrors *gr1;
TH1 *hfft;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=181351.25;
Double_t  fit_stop=196051.25;
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
int mdim = 50;
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
	const Int_t ir = (i)/(mdim) + hr_sum->FindBin(fit_start);
	const Int_t ic = (i)%(mdim) + hr_sum->FindBin(fit_start);


	//	cout<<i<<" "<<ir<<" "<<ic<<endl;
	
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
	    // data[i]=(hcalo->GetBinError(ir)*hcalo->GetBinError(ir))/(hcalo->GetBinError(ir)*hcalo->GetBinError(ir));
	    data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	    /*	    if(data[i]==0)
	      {
		cout<<"diag ele is zero at "<<ir<<" "<<ic;
		}*/
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }
       
       	if(ir==ic+nbinshift)
	  {
	    cout<<i<<" "<<ir<<" "<<ic<<endl;
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

	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir+15) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dprod3*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dprod4*h_sum->GetBinError(ir+30)*h_sum->GetBinError(ir+30)-(hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI13*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15))*(hr_sum->GetBinContent(ir+15) + 0.5*dI21*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI22*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI23*h_sum->GetBinError(ir+30)*h_sum->GetBinError(ir+30));

	    //  data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
	    data[i]=covar;
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
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
	      
	    covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir-15) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dprod3*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dprod4*h_sum->GetBinError(ir-30)*h_sum->GetBinError(ir-30)-(hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI13*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15))*(hr_sum->GetBinContent(ir-15) + 0.5*dI21*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI22*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dI23*h_sum->GetBinError(ir-30)*h_sum->GetBinError(ir-30));

	    // data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
	    // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	     data[i]=covar;
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

		covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir+30) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dprod3*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dprod4*h_sum->GetBinError(ir+30)*h_sum->GetBinError(ir+30) +  0.5*dprod5*h_sum->GetBinError(ir+45)*h_sum->GetBinError(ir+45) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI13*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15))*(hr_sum->GetBinContent(ir+30) + 0.5*dI21*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI22*h_sum->GetBinError(ir+45)*h_sum->GetBinError(ir+45) + 0.5*dI23*h_sum->GetBinError(ir+30)*h_sum->GetBinError(ir+30));

		// data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
		 // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	       	data[i]=covar;
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

		covar= hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ir-30) + 0.5*dprod1*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dprod2*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dprod3*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dprod4*h_sum->GetBinError(ir-30)*h_sum->GetBinError(ir-30) +  0.5*dprod5*h_sum->GetBinError(ir-45)*h_sum->GetBinError(ir-45) - (hr_sum->GetBinContent(ir) + 0.5*dI11*h_sum->GetBinError(ir)*h_sum->GetBinError(ir) + 0.5*dI12*h_sum->GetBinError(ir+15)*h_sum->GetBinError(ir+15) + 0.5*dI13*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15))*(hr_sum->GetBinContent(ir-30) + 0.5*dI21*h_sum->GetBinError(ir-15)*h_sum->GetBinError(ir-15) + 0.5*dI22*h_sum->GetBinError(ir-45)*h_sum->GetBinError(ir-45) + 0.5*dI23*h_sum->GetBinError(ir-30)*h_sum->GetBinError(ir-30));

		//	data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
		//cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	       	data[i]=covar;
		      }
	

      }
  cout<<"flgs "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());

  // cov.ResizeTo(h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), -1);
   
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

    return (2*f - ff - fb)/(2*f + ff + fb);
    

}




TH1D* construct_rhist_copy(TH1D *h_sum)
{
 
   h_sum->Rebin(16);//This brings the bin width of the calo summed histograms to 150 ns
   h_sum->Scale(0.0625);
   

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

 // hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //hcalo->Scale(0.125);

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
   
   for(Int_t i=0; i<50; i++)
    {
      for(Int_t j=0; j<50; j++)
	{
	  //  for(m=0; m<=105; m++)
	  // {
	  //  if(i==j || i==j+15 || i==j-15 || i==j+30 || i==j-30 || i==j+45 || i==j-45 || i==j+60 || i==j-60 || i==j+75 || i==j-75 || i==j+90 || i==j-90 || i==j+105 || i== j-105 || i==j+120 || i==j-120 || i==j+135 || i==j-135 || i==j+150 || i==j+165 || i==j-165 || i==j+180 || i==j-180)
	  //  if(cov2[i][j]!=0)
	  // {
	   ch=ch+((hr_sum->GetBinContent(i+hr_sum->FindBin(fit_start)))-(fuser->Eval(hr_sum->GetBinCenter(i+hr_sum->FindBin(fit_start)))))*cov2[i][j]*((hr_sum->GetBinContent(j+hr_sum->FindBin(fit_start)))-(fuser->Eval(hr_sum->GetBinCenter(j+hr_sum->FindBin(fit_start)))));
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


void hist_2_ratio_fit()
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



       




     
 

        hr_sum = construct_rhist_copy(h_sum);

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
   
  //   cov.Print();
        cov2.SetTol(1.e-23);

        cov2=ratio_cov_matrix(hr_sum,h_sum);
	cov2.Print();

      	hcov = new TH2D("hcov", "cov matrix hist", 50, 0.0, 50.0, 50, 0.0, 50.0);
	    
	for(int irow=0; irow<50; irow++)
        {
         for(int icol=0; icol<50; icol++)
	 {
	  hcov->SetBinContent(irow,icol,cov2(irow,icol));
	 }
        }
	    
	//	return;   		
        hr_sum->Fit("fprec3","RE","",fit_start,fit_stop);
 
      
	cout<<"fit over!!"<<endl;		
    
   
   
     	
       gettimeofday(&t_end, NULL);
       printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
       gStyle->SetOptFit(1111);
      		
       h_res= new TH1D("residual histogram ", "h_res", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_sum->GetBinContent(ibin)- fit_func3->Eval( hr_sum->GetBinCenter(ibin) ) );
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
 
   



   
}
