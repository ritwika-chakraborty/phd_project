#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TCanvas *c1;

TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5;

void ratiofit_omega_a_vs_binshift()
{
  Double_t omega_cbo[40], domega_cbo[40], tau_cbo[40], dtau_cbo[40], A_cbo[40], dA_cbo[40], phi_cbo[40], dphi_cbo[40], n[40];
  Int_t m=0;

   n[m]=10;
   A_cbo[m]=2.29070e-01;
   dA_cbo[m]=2.08640e-05;
   omega_cbo[m]=-4.84037e+01;
   domega_cbo[m]=1.16244e+00;
   phi_cbo[m]=2.22505e+00;
   dphi_cbo[m]=1.78602e-04;
   m=m+1;

  n[m]=11;
  A_cbo[m]=2.29070e-01;
   dA_cbo[m]=1.88504e-05;
  omega_cbo[m]=-4.78090e+01;
  domega_cbo[m]=1.05066e+00;
    phi_cbo[m]=2.22493e+00;
   dphi_cbo[m]=1.61459e-04;
  m=m+1;

  n[m]=12;
 A_cbo[m]=2.29071e-01;
   dA_cbo[m]=1.74664e-05;
  omega_cbo[m]=-4.71457e+01;
  domega_cbo[m]=9.75135e-01;
    phi_cbo[m]=2.22480e+00;
   dphi_cbo[m]=1.49875e-04;
   m=m+1;

  n[m]=13;
 A_cbo[m]=2.29073e-01;
   dA_cbo[m]=1.65910e-05;
  omega_cbo[m]=-4.66220e+01;
  domega_cbo[m]=9.28022e-01;
    phi_cbo[m]=2.22469e+00;
   dphi_cbo[m]=1.42649e-04;
   m=m+1;

  n[m]=14;
 A_cbo[m]=2.29075e-01;
   dA_cbo[m]=1.61582e-05;
  omega_cbo[m]=-4.61337e+01;
  domega_cbo[m]=9.04926e-01;
    phi_cbo[m]=2.22459e+00;
   dphi_cbo[m]=1.39106e-04;
   m=m+1;

  n[m]=15;
A_cbo[m]=2.29078e-01;
   dA_cbo[m]=1.61382e-05;
  omega_cbo[m]=-4.57073e+01;
  domega_cbo[m]=9.03828e-01;
    phi_cbo[m]=2.22450e+00;
   dphi_cbo[m]=1.38937e-04;
   m=m+1;

   n[m]=16;
A_cbo[m]=2.29080e-01;
   dA_cbo[m]=1.65298e-05;
  omega_cbo[m]=-4.54458e+01;
  domega_cbo[m]=9.24636e-01;
    phi_cbo[m]=2.22444e+00;
   dphi_cbo[m]=1.42129e-04;
   m=m+1;

  n[m]=17;
A_cbo[m]=2.29080e-01;
   dA_cbo[m]=1.73600e-05;
  omega_cbo[m]=-4.52233e+01;
  domega_cbo[m]=9.69167e-01;
    phi_cbo[m]=2.22439e+00;
   dphi_cbo[m]=1.48959e-04;
   m=m+1;

    n[m]=18;
A_cbo[m]=2.29080e-01;
   dA_cbo[m]=1.86912e-05;
  omega_cbo[m]=-4.51330e+01;
  domega_cbo[m]=1.04155e+00;
    phi_cbo[m]=2.22437e+00;
   dphi_cbo[m]=1.60061e-04;
   m=m+1;
  
  n[m]=19;
A_cbo[m]=2.29078e-01;
   dA_cbo[m]=2.06386e-05;
  omega_cbo[m]=-4.52248e+01;
  domega_cbo[m]=1.14921e+00;
    phi_cbo[m]=2.22439e+00;
   dphi_cbo[m]=1.76573e-04;
   m=m+1;
  
  n[m]=20;
A_cbo[m]=2.29075e-01;
   dA_cbo[m]=2.34041e-05;
  omega_cbo[m]=-4.54970e+01;
  domega_cbo[m]=1.30485e+00;
    phi_cbo[m]=2.22445e+00;
   dphi_cbo[m]=2.00441e-04;
   m=m+1;
  
  
   
    c1=new TCanvas("c1","Phase vs nbinshift");
    gr1=new TGraphErrors(m,n,phi_cbo,0,dphi_cbo);
    gr1->SetTitle("Phase vs nbinshift");
    gr1->GetXaxis()->SetTitle("nbinshift");
    gr1->GetYaxis()->SetTitle("phi");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    
   
    
}
