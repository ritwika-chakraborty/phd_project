#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"


/*int mdim = 3;
TMatrixD h(3,3);
TArrayD data(100);
*/
TH2D *hcov;
int i0, iend;
int flg = 0;
int mdim = h1->GetNbinsX();
// int mdim = 10;
  // TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
  TMatrixD cov(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
  TArrayD data((mdim)*(mdim));

int ratio_cov_matrix(){

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
  
  for(Int_t k=hcalo->FindBin(10000); k<=hcalo->GetNbinsX(); k++)
    {
      if(hcalo->GetBinError(k)!=0)
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
	    // data[i]=(hcalo->GetBinError(ir)*hcalo->GetBinError(ir))/(hcalo->GetBinError(ir)*hcalo->GetBinError(ir));
	    data[i]=(hcalo->GetBinError(ir)*hcalo->GetBinError(ir));
	    /*	    if(data[i]==0)
	      {
		cout<<"diag ele is zero at "<<ir<<" "<<ic;
		}*/
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

	if(ir==ic+15)
	  {
            u=h1->GetBinContent(ir);
	    v=h1->GetBinContent(ir+15);
	    w=h1->GetBinContent(ir-15);
	    x=h1->GetBinContent(ir+30);
	    z=h1->GetBinContent(ir-30);
	    s=h1->GetBinContent(ir+45);
	    t=h1->GetBinContent(ir-45);
	    
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

	    covar= hcalo->GetBinContent(ir)*hcalo->GetBinContent(ir+15) + 0.5*dprod1*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dprod2*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dprod3*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dprod4*h1->GetBinError(ir+30)*h1->GetBinError(ir+30)-(hcalo->GetBinContent(ir) + 0.5*dI11*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI12*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI13*h1->GetBinError(ir-15)*h1->GetBinError(ir-15))*(hcalo->GetBinContent(ir+15) + 0.5*dI21*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI22*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI23*h1->GetBinError(ir+30)*h1->GetBinError(ir+30));

	    //  data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
	    data[i]=covar;
	    //  cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	  }

		if(ir==ic-15)
	  {
            u=h1->GetBinContent(ir);
	    v=h1->GetBinContent(ir+15);
	    w=h1->GetBinContent(ir-15);
	    x=h1->GetBinContent(ir+30);
	    z=h1->GetBinContent(ir-30);
	    s=h1->GetBinContent(ir+45);
	    t=h1->GetBinContent(ir-45);
	    
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
	      
	    covar= hcalo->GetBinContent(ir)*hcalo->GetBinContent(ir-15) + 0.5*dprod1*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dprod2*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dprod3*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dprod4*h1->GetBinError(ir-30)*h1->GetBinError(ir-30)-(hcalo->GetBinContent(ir) + 0.5*dI11*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI12*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI13*h1->GetBinError(ir-15)*h1->GetBinError(ir-15))*(hcalo->GetBinContent(ir-15) + 0.5*dI21*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI22*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dI23*h1->GetBinError(ir-30)*h1->GetBinError(ir-30));

	    // data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
	    // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
	    data[i]=covar;
	  }

		if(ir==ic+30)
	      {
                u=h1->GetBinContent(ir);
	        v=h1->GetBinContent(ir+15);
	        w=h1->GetBinContent(ir-15);
	        x=h1->GetBinContent(ir+30);
	        z=h1->GetBinContent(ir-30);
	        s=h1->GetBinContent(ir+45);
	        t=h1->GetBinContent(ir-45);

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

		covar= hcalo->GetBinContent(ir)*hcalo->GetBinContent(ir+30) + 0.5*dprod1*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dprod2*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dprod3*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dprod4*h1->GetBinError(ir+30)*h1->GetBinError(ir+30) +  0.5*dprod5*h1->GetBinError(ir+45)*h1->GetBinError(ir+45) - (hcalo->GetBinContent(ir) + 0.5*dI11*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI12*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI13*h1->GetBinError(ir-15)*h1->GetBinError(ir-15))*(hcalo->GetBinContent(ir+30) + 0.5*dI21*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI22*h1->GetBinError(ir+45)*h1->GetBinError(ir+45) + 0.5*dI23*h1->GetBinError(ir+30)*h1->GetBinError(ir+30));

		// data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
		 // cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		data[i]=covar;
	      }

	     if(ir==ic-30)
	      {
                u=h1->GetBinContent(ir);
	        v=h1->GetBinContent(ir+15);
	        w=h1->GetBinContent(ir-15);
	        x=h1->GetBinContent(ir+30);
	        z=h1->GetBinContent(ir-30);
	        s=h1->GetBinContent(ir+45);
	        t=h1->GetBinContent(ir-45);

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

		covar= hcalo->GetBinContent(ir)*hcalo->GetBinContent(ir-30) + 0.5*dprod1*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dprod2*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dprod3*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dprod4*h1->GetBinError(ir-30)*h1->GetBinError(ir-30) +  0.5*dprod5*h1->GetBinError(ir-45)*h1->GetBinError(ir-45) - (hcalo->GetBinContent(ir) + 0.5*dI11*h1->GetBinError(ir)*h1->GetBinError(ir) + 0.5*dI12*h1->GetBinError(ir+15)*h1->GetBinError(ir+15) + 0.5*dI13*h1->GetBinError(ir-15)*h1->GetBinError(ir-15))*(hcalo->GetBinContent(ir-30) + 0.5*dI21*h1->GetBinError(ir-15)*h1->GetBinError(ir-15) + 0.5*dI22*h1->GetBinError(ir-45)*h1->GetBinError(ir-45) + 0.5*dI23*h1->GetBinError(ir-30)*h1->GetBinError(ir-30));

		//	data[i]=covar/sqrt(hcalo->GetBinError(ir)*hcalo->GetBinError(ir)*hcalo->GetBinError(ic)*hcalo->GetBinError(ic));
		//cout<<ir<<" "<<ic<<" "<<data[i]<<endl;
		data[i]=covar;
	      }
	  

      }
  cout<<"flgs "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());

  cov.ResizeTo(h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), h_sum->FindBin(0)+4, mdim-(h_sum->FindBin(0)), -1);
   
  //   cov.Print();
    cov.SetTol(1.e-23);
    Double_t det1;
    cov.Invert(&det1);
    //   cov.Print();
     hcov = new TH2D("hcov", "cov matrix hist", 2100, 0.0, 2100.0, 2100, 0.0, 2100.0);
  for(int irow=44; irow<=2060; irow++)
    {
      for(int icol=44; icol<=2060; icol++)
	{
	  hcov->SetBinContent(irow,icol,cov(irow,icol));
	}
    }

    cout<<mdim<<endl;
    //cout<<"Matrix inverted "<<i0<<endl;
  return 0;

}
