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
TCanvas *c4, *c5;
TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *h_res, *h1[2], *hcomp;
TH1D *hsigma;
TH1 *hm;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Int_t countfcn=0;

Double_t wiggle(Double_t *x, Double_t *par)
{
  Double_t f_x=par[0]*exp(-x[0]/par[1])*(1+par[2]*cos(par[3]*x[0]+par[4]));
  return f_x;
}

Double_t cbo(Double_t *x, Double_t *par)
{
  Double_t f_x=(1-par[0]*exp(-x[0]/par[1])*(cos(par[2]*x[0]+par[3])));
  return f_x;
}

double paramToFreq(double blindedValue)
{

  double unblindedR = blindedValue - deltaR;
  return 2 * TMath::Pi() * refFreq * (1 + (unblindedR * precisionR));
}


Double_t fprec(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A1_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi1_cbo = par[8];

    double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double time = x[0];
    
    double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym* (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo)))*cos(omega_a*time + phi)) * (1-A1_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi1_cbo)));

    }

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

 void run2_fr_thresh_850()

 {
   Double_t chi[40],cal[40],norm[40],dcal[40];
   Double_t n[40];
   Int_t m=0;
   Double_t N;

  
   
  TString blindingString("Ritwika's blinding!"); // your blinding string
  TRandom3 *r3 = new TRandom3(Hash(blindingString));
  deltaR = 24.*2*(r3->Rndm() - 0.5); // random number +- 24 ppm

 // cout<<r3->Rndm()<<" "<<endl;
  
 c4=new TCanvas("c4","9 parameter wiggle_fit");  
 TFile *_file[2];
 TDirectoryFile *dir[2];
 //  TFile *_file=TFile::Open("rc_hist_sum_60hr.root");



 fit_func= new TF1("fprec", fprec,  114000,  350000, 11);
 fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A1_cbo", "Tau_cbo", "omega_cbo", "phi1_cbo", "A2_cbo", "phi2_cbo"); 
 fit_func->SetNpx(1000000);

 // fit_func->SetParameters(15800000000/1.1, 51568, 0.2, 135.54, 3.14/2, 0.01, 128325, 0.0029, 1.57, 0.001, 1.5);
 //  fit_func->SetParameters(1580000000/0.85, 51568, 0.2, 135.54, 3.14/2, 0.01, 128325, 0.0029, 1.57, 0.0, 1.5);//rebin 0
 // fit_func->SetParameters(15800000000/4.25, 51568, 0.2, 135.54, 3.14/2, 0.005, 148325, 0.0029, 1.57, 0.0, 3.14/2);//rebin 2
 // fit_func->SetParameters(15800000000/2.1, 51568, 0.2, 135.54, 3.14/2, 0.001, 128325, 0.0029, 1.57, 0.0, 3.14/2);//rebin 4
 //  fit_func->SetParameters(1580000000/2.5, 51568, 0.2, 135.54, 3.14/2, 0.01, 128325, 0.0029, 3.14/2);
   fit_func->SetParameters(1580000000/3.1, 51568, 0.2, 135.54, 3.14/2, 0.01, 128325, 0.0029, 1.57, 0.0, 1.5);//calo 1
  


  gStyle->SetOptFit(1111);
  c4->Divide(1,3);
  c4->cd(1);
  //  h[1]->Rebin(4);
   
   h1[1]= new TH1D("180 degree shifted histogram", "h1[1]", h[1]->GetNbinsX(), 100001, 352001);
   h1[0]= new TH1D("-180 degree shifted histogram", "h1[0]", h[1]->GetNbinsX(), 100001, 352001);
   hcomp= new TH1D("superposed histogram", "hcomp", h[1]->GetNbinsX(), 100001, 352001);   
   for(Int_t k=1; k<=h1[1]->GetNbinsX(); k++)
     {
       h1[1]->SetBinContent(k,h[1]->GetBinContent(k+1));
       h1[0]->SetBinContent(k,h[1]->GetBinContent(k-1));
     }
   for(Int_t l=1; l<=hcomp->GetNbinsX(); l++)
     {
       hcomp->SetBinContent(l,(h[1]->GetBinContent(l)+h1[1]->GetBinContent(l))/2);
     }
   h1[1]->SetLineColor(kRed);
   hcomp->SetLineColor(kBlack);

   h[1]->GetYaxis()->SetRangeUser(-100., 60000000.);//rebin 0

  
   
  
   hcomp->Rebin(2);
   // hcomp->Fit("fprec","R","",fit_start,fit_stop);
   c4->cd(1);
       for(Int_t p=2; p<=hcomp->GetNbinsX(); p++)
     {
       hcomp->SetBinError(p,sqrt((1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1))));
       // cout<<(1/2)*sqrt((1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1)))<<" "<<endl;
       // cout<<(1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1))<<" "<<endl;
       }

    hcomp->GetYaxis()->SetRangeUser(-100., 120000000.);

   //cout<<(1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1))<<" "<<endl;
    printf("start Q-method analyzer\n");
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
     gFitter = TVirtualFitter::Fitter(hcomp);
     gFitter->SetFCN(chi2);  
     hcomp->Fit("fprec","RUE","",fit_start,fit_stop);
    //   hcomp->Draw();
    //  fit_func->Draw("same");
     gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));

      

 cbo_func= new TF1("cbo", cbo,  fit_start,  fit_stop, 4);
 cbo_func->SetParNames("A_cbo", "Tau_cbo", "omega_cbo", "Phi"); 
 cbo_func->SetNpx(1000000);
 cbo_func->SetParameters(699000000, 22060, 0.0029, 3.14/2);
 // h_res->Fit("cbo","R","",130000,330000);
 // cbo_func->Draw("same");
 
 c4->cd(2);
 h_res= new TH1D("residual histogram", "h_res", hcomp->GetNbinsX(), 100001, 352001);
  for (int ibin = ((fit_start + 1 - 100001)/hcomp->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/hcomp->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcomp->GetBinContent(ibin)- fit_func->Eval( hcomp->GetBinCenter(ibin) ) );
      if(hcomp->GetBinError(ibin)!=0){res=(res/hcomp->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();
   c4->cd(3);
    hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   hm->Draw();
 
   c5=new TCanvas("c5","sigma"); 
   hsigma= new TH1D("sigma histogram", "hsigma", 10000, -5, 5);
   for(Int_t u=((fit_start + 1 - 100001)/hcomp->GetBinWidth(1))+1; u<=((fit_stop +1 - 100001)/hcomp->GetBinWidth(1))+1; u++)
     {
       for(Int_t v=1; v<=hsigma->GetNbinsX(); v++)
	 {
           if((int (h_res->GetBinContent(u)*1000))==(int (hsigma->GetBinCenter(v)*1000)))
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









