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
int i0, iend;
int flg = 0;
int mdim = (hcalo->GetNbinsX()/2)-1;
// int mdim = 10;
  // TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
  TMatrixD cov(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
  TArrayD data((mdim)*(mdim));

int matrix(){

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
 
    
    if ( ir == ic && ic>=i0)
      {
	if((2*ir-1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = ((0.25)*(hcalo->GetBinError(2*ir-1))*(hcalo->GetBinError(2*ir-1))+(1)*(hcalo->GetBinError(2*ir))*(hcalo->GetBinError(2*ir))+(0.25)*(hcalo->GetBinError(2*ir+1))*(hcalo->GetBinError(2*ir+1)));
	 cout<<"diagonal "<<ir<<" "<<ic<<endl;
	 // cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<h[1]->GetBinError(2*ir+1)<<" "<<endl;
	  if(data[ir*mdim+ic]==0)
	    {
	      cout<<ir<<" "<<ic<<"breaking"<<endl;
      	      iend=ir;
	      break;
	    }
	  else
	    {
	      iend=mdim;
	    }
       }
      }   
    if ( ic == ir+1 && ic>=i0 && ir>=i0)
      {
        if((2*ir+1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = (0.25)*(hcalo->GetBinError(2*ir+1))*(hcalo->GetBinError(2*ir+1));
	 //	 	cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<((0.25)*(h[1]->GetBinError(2*ir+1))*(h[1]->GetBinError(2*ir+1)))<<" "<<endl;
	 cout<<"diagonal+1  "<<ir<<" "<<ic<<endl;
       }
      }
     if ( ic == ir-1 && ic>=i0 && ir>=i0)
      {
        if((2*ir-1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = (0.25)*(hcalo->GetBinError(2*ir-1))*(hcalo->GetBinError(2*ir-1));
	 //	 	cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<((0.25)*(h[1]->GetBinError(2*ir+1))*(h[1]->GetBinError(2*ir+1)))<<" "<<endl;
	 cout<<"diagonal+1  "<<ir<<" "<<ic<<endl;
       }
      }

    flg=i;
     }
  cout<<"flg "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());

  cov.ResizeTo(i0, iend-1, i0, iend-1, -1);
  
  //   cov.Print();

   Double_t det1;
   cov.Invert(&det1);
   //   cov.Print();
   cout<<"Matrix inverted "<<i0<<endl;
  return 0;
}
