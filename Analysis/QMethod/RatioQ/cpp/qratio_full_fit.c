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


TVirtualFitter *gFitter;
TCanvas *cfit;
TF1  *fit_func1;
TH1D *h_res, *hcomp;
TF1 *fit_func, *fit_func2;
TH1 *hfft, *hfft2;
TH1F *hlm;
Double_t fit_start, fit_stop;
Int_t countfcn=0;
Int_t m;



Double_t fprec(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    return asym*cos(omega_a*time + phi);
}

Double_t fprec2(Double_t *x, Double_t *par)
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

    double lost_muon_amp = par[15];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

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

    double Nlm= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double Nlmf= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time + T_a/2) * 1.0e-3)));

    double Nlmb= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time - T_a/2) * 1.0e-3)));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));
    
     return (2*f*Ncbo*Nvw*Nlm - ff*Ncbof*Nvwf*Nlmf - fb*Ncbob*Nvwb*Nlmb)/(2*f*Ncbo*Nvw*Nlm + ff*Ncbof*Nvwf*Nlmf + fb*Ncbob*Nvwb*Nlmb);
    //  return 2*(1+asym*cos(omega_a*time + phi)) - (1+asym*cos(omega_a*(time+ (T_a+20)/2) + phi)) - (1+asym*cos(omega_a*(time - (T_a+20)/2) + phi));
    // return (2*f - ff - fb)/(2*f + ff + fb);
      // /(2*f + ff + fb);
}


void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(hcalo->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
   
   for(Int_t i=hcalo->FindBin(fit_start); i<=hcalo->FindBin(fit_stop); i++)
    {
      for(Int_t j=hcalo->FindBin(fit_start); j<=hcalo->FindBin(fit_stop); j++)
	{
	  //  for(m=0; m<=105; m++)
	  // {
	  //  if(i==j || i==j+15 || i==j-15 || i==j+30 || i==j-30 || i==j+45 || i==j-45 || i==j+60 || i==j-60 || i==j+75 || i==j-75 || i==j+90 || i==j-90 || i==j+105 || i== j-105 || i==j+120 || i==j-120 || i==j+135 || i==j-135 || i==j+150 || i==j+165 || i==j-165 || i==j+180 || i==j-180)
	  if(cov[i][j]!=0)
	     {
	      ch=ch+((hcalo->GetBinContent(i))-(fuser->Eval(hcalo->GetBinCenter(i))))*cov[i][j]*((hcalo->GetBinContent(j))-(fuser->Eval(hcalo->GetBinCenter(j))));
	      //cout<<m<<endl;
	    
	     }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}

 void qratio_full_fit()

 {
 

 // cout<<r3->Rndm()<<" "<<endl;
   // T_a = T_a + 20;
   
 cfit=new TCanvas("cfit","full ratio fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");



      fit_func= new TF1("fprec", fprec,  30000,  309000, 3);
     fit_func->SetParNames("A", "Blind R", "Phi");
     fit_func->SetNpx(1000000);

        fit_func2= new TF1("fprec2", fprec2,  30000,  309000, 16);
	fit_func2->SetParNames("A", "Blind R", "Phi", "A_cbo", "tau_cbo", "omega_cbo", "phi_cbo","asym_cbo_A","asym_cbo_phi","phi_cbo_A","phi_cbo_phi");
	fit_func2->SetParName(11,"A_vw");
	fit_func2->SetParName(12,"tau_vw");
	fit_func2->SetParName(13,"omega_vw");
	fit_func2->SetParName(14,"phi_vw");
	fit_func2->SetParName(15,"lm_amp");

     fit_func2->SetNpx(1000000);
     fit_func2->SetParameters(0.23, 0.0, 2.18, 0.02,237000,0.00234,-2.2);
     fit_func2->SetParameter(7,-0.02);
     fit_func2->SetParameter(8,4.9);
     fit_func2->SetParameter(9,0.0);
     fit_func2->SetParameter(10,0.0);
     fit_func2->SetParameter(11,0.002);
     fit_func2->SetParameter(12,1000);
     fit_func2->SetParameter(13,0.1535);
     fit_func2->SetParameter(14,2.24);
     fit_func2->SetParameter(15,0.0);

     //  fit_func2->FixParameter(3,0);
     // fit_func2->FixParameter(4,0);
     // fit_func2->FixParameter(5,0);
     // fit_func2->FixParameter(6,0);
     // fit_func2->FixParameter(7,0);
     // fit_func2->FixParameter(8,0);
     //  fit_func2->FixParameter(9,0);
     // fit_func2->FixParameter(10,0);
     //fit_func2->FixParameter(11,0);
     // fit_func2->FixParameter(12,0);
      fit_func2->FixParameter(13,0.1535);
     // fit_func2->FixParameter(14,0);
      fit_func2->FixParameter(15,0);

     //  fit_func2->Draw();
     //return;


     //  gStyle->SetOptFit(1111);

     fit_start=30000;
     fit_stop=300000;

     fit_func->SetParameters(0.2, 0.0, 3.14/2);
     cfit->Divide(1,3);
     cfit->cd(1);

     //h_num->GetYaxis()->SetRangeUser(-100., 30000000000.);
	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h_num->Draw();
	//  fit_func->Draw("same");


	//	return;


    printf("start Q-method analyzer\n");
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
     gFitter = TVirtualFitter::Fitter(hcalo);
     gFitter->SetFCN(chi2);  
     hcalo->Fit("fprec2","RUE","",fit_start,fit_stop);
     //hcalo->Draw();
     //fit_func2->Draw("same");
      gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));

      
		
		h_res= new TH1D("residual histogram", "h_res", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


        cfit->cd(2);
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcalo->GetBinContent(ibin)- fit_func2->Eval( hcalo->GetBinCenter(ibin) ) );
         if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");	
        h_res->Draw();

   cfit->cd(3);
   hfft = h_res->FFT(hfft, "MAG");
   hfft->SetLineColor(kBlack);
   hfft->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hfft->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hfft->GetXaxis()->SetTitle("Freq [GHz]");
   hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft->Draw();


 
   
 }









