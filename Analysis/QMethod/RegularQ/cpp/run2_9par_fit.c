#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );


TCanvas *c2;
TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *h_res;
TH1 *hm;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;

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

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

    double time = x[0];
    
    double omega_a = rawBinToNs * 1.e-3 * getBlinded.paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi)) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo)));

    }

 void run2_9par_fit()

 {
   Double_t chi[40],cal[40],norm[40],dcal[40];
   Double_t n[40];
   Int_t m=0;
   Double_t N;

   h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);
   
   /*  TString blindingString("Ritwika's blinding!"); // your blinding string
  TRandom3 *r3 = new TRandom3(Hash(blindingString));
  deltaR = 24.*2*(r3->Rndm() - 0.5); // random number +- 24 ppm
   */
  // cout<<r3->Rndm()<<" "<<endl;
  
 c2=new TCanvas("c2","5 parameter wiggle_fit");  
 TFile *_file[2];
  TDirectoryFile *dir[2];
  //  TFile *_file=TFile::Open("rc_hist_sum_60hr.root");



 fit_func= new TF1("fprec", fprec,  114000,  330000, 9);
 fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo"); 
 fit_func->SetNpx(1000000);

  fit_func->SetParameters(15800000000/1.1, 51568, 0.2, 0, 3.14/2, 0.01, 128325, 0.0029, 3.14/2);
 //  fit_func->SetParameters(1580000000/2.5, 51568, 0.2, 135.54, 3.14/2, 0.01, 128325, 0.0029, 3.14/2);


 gStyle->SetOptFit(1111);
  c2->Divide(1,3);
  c2->cd(1);
 h_sum->Rebin(8);

  // fit_func->SetParLimits(4, -6.28, 6.28);
   h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
    h_sum->Fit("fprec","R","",130000,330000);
   
    //    h_sum->Draw();
    //  fit_func->Draw("same");

  for (int ibin = 251; ibin <= 2100; ++ibin)
     {
       // if(ibin>=120 && ibin<=1900)
	 // {
      double res =  (h_sum->GetBinContent(ibin)- fit_func->Eval( h_sum->GetBinCenter(ibin) ) );
      h_res->SetBinContent(ibin, res  );
     
     }
  
   c2->cd(2);
   h_res->GetYaxis()->SetRangeUser(-6000000., 6000000.);
   //  h_res->GetXaxis()->SetRangeUser(125000, 350000);
   h_res->Draw();


 cbo_func= new TF1("cbo", cbo,  130000,  330000, 4);
 cbo_func->SetParNames("A_cbo", "Tau_cbo", "omega_cbo", "Phi"); 
 cbo_func->SetNpx(1000000);
 cbo_func->SetParameters(699000000, 22060, 0.0029, 3.14/2);
 // h_res->Fit("cbo","R","",130000,330000);
 // cbo_func->Draw("same");
 
 c2->cd(3);
 // hm= new TH1("FFT histogram", "hm", 2100, 100001, 352001);
 hm = h_res->FFT(hm, "MAG");
 hm->SetLineColor(kBlack);
  hm->Draw();
 //  cbo_func->Draw("same"); 
 //  h_res->Fit("cbo","R","",170000,330000); 
  /*   c1->cd(2);
 // h[1]->Fit("wiggle","R","",500,2500);
  h[1]->Rebin(8);
 fit_func->SetParameters(15800000000, 51568, 0.2, 135.54, 0);
 h[1]->Fit("fprec","R","",114000,330000);
 h[1]->GetYaxis()->SetRangeUser(-100., 3000000000.);
 h[1]->Draw();*/ 
 }
