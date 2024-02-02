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

void ratiofit_omega_a_vs_binshift_fullfit()
{
  Double_t omega_cbo[40], domega_cbo[40], tau_cbo[40], dtau_cbo[40], A_cbo[40], dA_cbo[40], phi_cbo[40], dphi_cbo[40], n[40];
  Int_t m=0;

   n[m]=10;
   A_cbo[m]=2.29077e-01;
   dA_cbo[m]=2.08877e-05;
   omega_cbo[m]=-4.51128e+01;
   domega_cbo[m]=1.16405e+00;
   phi_cbo[m]=2.22438e+00;
   dphi_cbo[m]=1.79037e-04;
   m=m+1;

  n[m]=11;
  A_cbo[m]=2.29078e-01;
   dA_cbo[m]=1.88708e-05;
  omega_cbo[m]=-4.51110e+01;
  domega_cbo[m]=1.05212e+00;
    phi_cbo[m]=2.22438e+00;
   dphi_cbo[m]=1.61860e-04;
  m=m+1;

    n[m]=12;
 A_cbo[m]=2.29079e-01;
   dA_cbo[m]=1.74852e-05;
  omega_cbo[m]=-4.50937e+01;
  domega_cbo[m]=9.76412e-01;
    phi_cbo[m]=2.22438e+00;
   dphi_cbo[m]=1.50223e-04;
   m=m+1;
  
  n[m]=13;
 A_cbo[m]=2.29080e-01;
   dA_cbo[m]=1.66092e-05;
  omega_cbo[m]=-4.51182e+01;
  domega_cbo[m]=9.29337e-01;
    phi_cbo[m]=2.22438e+00;
   dphi_cbo[m]=1.43002e-04;
   m=m+1;

  n[m]=14;
 A_cbo[m]=2.29080e-01;
   dA_cbo[m]=1.61755e-05;
  omega_cbo[m]=-4.51255e+01;
  domega_cbo[m]=9.06240e-01;
    phi_cbo[m]=2.22438e+00;
   dphi_cbo[m]=1.39458e-04;
   m=m+1;

     n[m]=15;
A_cbo[m]=2.29081e-01;
   dA_cbo[m]=1.61546e-05;
  omega_cbo[m]=-4.51606e+01;
  domega_cbo[m]=9.05113e-01;
    phi_cbo[m]=2.22439e+00;
   dphi_cbo[m]=1.39271e-04;
   m=m+1;
   
   n[m]=16;
A_cbo[m]=2.29082e-01;
   dA_cbo[m]=1.65452e-05;
  omega_cbo[m]=-4.52363e+01;
  domega_cbo[m]=9.25754e-01;
    phi_cbo[m]=2.22440e+00;
   dphi_cbo[m]=1.42431e-04;
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
A_cbo[m]=2.29081e-01;
   dA_cbo[m]=2.34267e-05;
  omega_cbo[m]=-4.52789e+01;
  domega_cbo[m]=1.30618e+00;
    phi_cbo[m]=2.22440e+00;
   dphi_cbo[m]=2.00801e-04;
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
