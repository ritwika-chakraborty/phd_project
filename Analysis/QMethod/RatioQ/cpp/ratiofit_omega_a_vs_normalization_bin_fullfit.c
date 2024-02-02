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

void ratiofit_omega_a_vs_normalization_bin_fullfit()
{
  Double_t omega_cbo[40], domega_cbo[40], tau_cbo[40], dtau_cbo[40], A_cbo[40], dA_cbo[40], phi_cbo[40], dphi_cbo[40], n[40];
  Int_t m=0;

   n[m]=300;
   omega_cbo[m]=-4.52131e+01;
   domega_cbo[m]=1.05058e+00;
   m=m+1;

    n[m]=310;
   omega_cbo[m]=-4.52131e+01;
   domega_cbo[m]=1.05058e+00;
   m=m+1;

    n[m]=320;
   omega_cbo[m]=-4.51634e+01;
   domega_cbo[m]=1.05051e+00;
   m=m+1;

    n[m]=330;
   omega_cbo[m]=-4.51203e+01;
   domega_cbo[m]=1.05049e+00;
   m=m+1;

    n[m]=340;
   omega_cbo[m]=-4.51681e+01;
   domega_cbo[m]=1.05052e+00;
   m=m+1;

    n[m]=350;
   omega_cbo[m]=-4.51794e+01;
   domega_cbo[m]=1.05052e+00;
   m=m+1;

    n[m]=360;
   omega_cbo[m]=-4.51843e+01;
   domega_cbo[m]=1.05054e+00;
   m=m+1;

    n[m]=370;
   omega_cbo[m]=-4.51869e+01;
   domega_cbo[m]=1.05054e+00;
   m=m+1;

    n[m]=380;
   omega_cbo[m]=-4.51852e+01;
   domega_cbo[m]=1.05054e+00;
   m=m+1;

    n[m]=390;
   omega_cbo[m]=-4.51863e+01;
   domega_cbo[m]=1.05055e+00;
   m=m+1;

    n[m]=400;
   omega_cbo[m]=-4.51858e+01;
   domega_cbo[m]=1.05054e+00;
   m=m+1;

 
  
   
    c1=new TCanvas("c1","R vs normalization start bin");
    gr1=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr1->SetTitle("R vs Nomalization");
    gr1->GetXaxis()->SetTitle("norm start bin");
    gr1->GetYaxis()->SetTitle("R");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    
   
    
}
