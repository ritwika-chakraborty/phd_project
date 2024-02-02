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

//blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

//blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );



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
TCanvas *c11, *c12;
TF1  *fit_func1;
TH1D *h_res1, *h1[2], *hcomp;
TH1D *hsigma;
TH1 *hm1;

Int_t countfcn=0;

void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(hcomp->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
  
   for(Int_t i=hcomp->FindBin(fit_start); i<=hcomp->FindBin(fit_stop); i++)
    {
      for(Int_t j=hcomp->FindBin(fit_start); j<=hcomp->FindBin(fit_stop); j++)
	{
	  // if(countfcn<=1000)
	  // {
	   if(j<=i+5 && j>=i-5)
	   {
	     ch=ch+((hcomp->GetBinContent(i))-(fuser->Eval(hcomp->GetBinCenter(i))))*cov[i][j]*((hcomp->GetBinContent(j))-(fuser->Eval(hcomp->GetBinCenter(j))));
	    }
	   //  }
	   /*  else
	    {
	      if(j!=i)
		{ ch=ch+((hcomp->GetBinContent(i))-(fuser->Eval(hcomp->GetBinCenter(i))))*cov[i][j]*((hcomp->GetBinContent(j))-(fuser->Eval(hcomp->GetBinCenter(j))));}
		}*/
		}
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}

 void run2_calo_fr()

 {
 

  
   


 // cout<<r3->Rndm()<<" "<<endl;
  
 c11=new TCanvas("c11","9 parameter wiggle_fit");  
 TFile *_file[2];
 TDirectoryFile *dir[2];
 //  TFile *_file=TFile::Open("rc_hist_sum_60hr.root");

 fit_func1= new TF1("fprec", fprec,  30000,  309000, 18);
  fit_func1->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo", "A2_cbo","phi2_cbo");
  //, "A_vbo", "Tau_vbo", "omega_vbo", "phi_vbo");
  // );
     fit_func1->SetParName(11, "A3_cbo");
     fit_func1->SetParName(12, "phi3_cbo");
     fit_func1->SetParName(13, "lm_amp");
     fit_func1->SetParName(14, "A_vbo");
     fit_func1->SetParName(15, "Tau_vbo");
     fit_func1->SetParName(16, "omega_vbo");
     fit_func1->SetParName(17, "phi_vbo");

   
 fit_func1->SetNpx(1000000);



 //gStyle->SetOptFit(1111);
 //  c2->Divide(1,3);
 //  c2->cd(1);
 // h_sum->Rebin(8);




      fit_start=30000;
      fit_stop=300000;

      //   h_sum->Rebin(8);
   fit_func1->SetParameters(961159000, 64440.8, 0.2248, 0, 2.284, 0.02025, 240550, 0.00234054, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func1->SetParameter(9,0.02219);
     fit_func1->SetParameter(10,0);
     fit_func1->SetParameter(11, 0.00272);
     fit_func1->SetParameter(12,0);
     fit_func1->SetParameter(13,0);

     fit_func1->SetParameter(14, -0.00226);
     fit_func1->SetParameter(15, 25387);
     fit_func1->SetParameter(16, 0.1531);
     fit_func1->SetParameter(17, 0);
  
     
     //fit_func->FixParameter(5,0);
     //  fit_func->FixParameter(9,0);
     //fit_func->FixParameter(11,0);
     //fit_func->FixParameter(13,0);



     // gStyle->SetOptFit(1111);
  c11->Divide(1,3);
  c11->cd(1);
  //  h[1]->Rebin(4);
   
   h1[1]= new TH1D("180 degree shifted histogram", "h1[1]", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   h1[0]= new TH1D("-180 degree shifted histogram", "h1[0]", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   hcomp= new TH1D("superposed histogram", "h_calo1", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   for(Int_t k=1; k<=h1[1]->GetNbinsX(); k++)
     {
       h1[1]->SetBinContent(k,hcalo->GetBinContent(k+1));
       h1[0]->SetBinContent(k,hcalo->GetBinContent(k-1));
     }
   for(Int_t l=1; l<=hcomp->GetNbinsX(); l++)
     {
       hcomp->SetBinContent(l,(hcalo->GetBinContent(l)+h1[1]->GetBinContent(l))/2);
     }
   h1[1]->SetLineColor(kRed);
   hcomp->SetLineColor(kBlack);

   hcalo->GetYaxis()->SetRangeUser(-100., 60000000.);//rebin 0

  
   
  
   hcomp->Rebin(2);
   // hcomp->Fit("fprec","R","",fit_start,fit_stop);
   c11->cd(1);
       for(Int_t p=2; p<=hcomp->GetNbinsX(); p++)
     {
       hcomp->SetBinError(p,sqrt((1/4)*(hcalo->GetBinError(2*p-1))*(hcalo->GetBinError(2*p-1))+(1)*(hcalo->GetBinError(2*p))*(hcalo->GetBinError(2*p))+(1/4)*(hcalo->GetBinError(2*p+1))*(hcalo->GetBinError(2*p+1))));

       }

    hcomp->GetYaxis()->SetRangeUser(-100.,  30000000000./24);
	hcomp->GetYaxis()->SetTitle("ADC counts");
	hcomp->GetXaxis()->SetTitle("time [ns]");

   //cout<<(1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1))<<" "<<endl;
    printf("start Q-method analyzer\n");
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
     gFitter = TVirtualFitter::Fitter(hcomp);
     gFitter->SetFCN(chi2);  
     hcomp->Fit("fprec","RUE","",fit_start,fit_stop);
    // hcomp->Draw();
     // fit_func1->Draw("same");
      gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));

      


 
 c11->cd(2);
 h_res1= new TH1D("residual histogram fr", "h_res1", hcomp->GetNbinsX(),  hcomp->GetBinLowEdge(1), (hcomp->GetBinLowEdge(hcomp->GetNbinsX())+hcomp->GetBinWidth(1)));
  for (int ibin = ((fit_start + 1 -  hcomp->GetBinLowEdge(1))/hcomp->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 -  hcomp->GetBinLowEdge(1))/hcomp->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcomp->GetBinContent(ibin)- fit_func1->Eval( hcomp->GetBinCenter(ibin) ) );
      if(hcomp->GetBinError(ibin)!=0){res=(res/hcomp->GetBinError(ibin));}
      h_res1->SetBinContent(ibin, (res)  );
     
      }
    h_res1->GetYaxis()->SetTitle("ADC counts");
    h_res1->GetXaxis()->SetTitle("time [ns]");	
    h_res1->Draw();
    c11->cd(3);
    hm1 = h_res1->FFT(hm1, "MAG");
    hm1->SetLineColor(kBlack);
     hm1->SetBins(h_res1->GetNbinsX(),0,h_res1->GetNbinsX()/((h_res1->GetBinLowEdge(h_res1->GetNbinsX())+h_res1->GetBinWidth(1))));
   hm1->GetXaxis()->SetRangeUser(0,h_res1->GetNbinsX()/((h_res1->GetBinLowEdge(h_res1->GetNbinsX())+h_res1->GetBinWidth(1)))/2);
   hm1->GetXaxis()->SetTitle("Freq [GHz]");
   hm1->GetYaxis()->SetTitle("FFT Mag [arb. units]");

    hm1->Draw();
 
   c12=new TCanvas("c12","sigma"); 
   hsigma= new TH1D("sigma histogram", "hsigma", 10000, -5, 5);
   for(Int_t u=((fit_start + 1 - 100001)/hcomp->GetBinWidth(1))+1; u<=((fit_stop +1 - 100001)/hcomp->GetBinWidth(1))+1; u++)
     {
       for(Int_t v=1; v<=hsigma->GetNbinsX(); v++)
	 {
           if((int (h_res1->GetBinContent(u)*1000))==(int (hsigma->GetBinCenter(v)*1000)))
	     {
	       //	       cout<<"voila "<<v<<" "<<hsigma->GetBinCenter(v)<<endl;
	       hsigma->AddBinContent(v,1);
	     }
	 }
     }
   hsigma->Rebin(100);
   hsigma->Draw();
   hsigma->Fit("gaus");
   // hsigma->GetStartTime();
   
  
 }









