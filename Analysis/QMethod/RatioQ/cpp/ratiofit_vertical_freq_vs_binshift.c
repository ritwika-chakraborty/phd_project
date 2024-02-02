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

void ratiofit_vertical_freq_vs_binshift()
{
  Double_t omega_cbo[40], domega_cbo[40], tau_cbo[40], dtau_cbo[40], A_cbo[40], dA_cbo[40], phi_cbo[40], dphi_cbo[40], n[40];
  Int_t m=0;

  n[m]=10;
  A_cbo[m]=7.94027e-03;
  dA_cbo[m]=1.81655e-03;
  tau_cbo[m]=1.72234e+04;
  dtau_cbo[m]=1.75897e+03;
  omega_cbo[m]=1.40329e-02;
  domega_cbo[m]=6.14943e-06;
  phi_cbo[m]=1.54767e+01;
  dphi_cbo[m]=2.34893e-01;
  m=m+1;

  n[m]=11;
  A_cbo[m]=9.04335e-03;
  dA_cbo[m]=2.32428e-03;
  tau_cbo[m]=1.62954e+04;
  dtau_cbo[m]=1.81494e+03;
  omega_cbo[m]=1.40231e-02;
  domega_cbo[m]=8.00132e-06;
  phi_cbo[m]=1.58929e+01;
  dphi_cbo[m]=3.04604e-01;
  m=m+1;

  /*  n[m]=12;
  A_cbo[m]=1.02004e-04;
  dA_cbo[m]=2.28298e-05;
  tau_cbo[m]=3.46544e+05;
  dtau_cbo[m]=2.95517e+05;
  omega_cbo[m]=1.52427e-02;
  domega_cbo[m]=2.74382e-06;
  phi_cbo[m]=-2.22529e+01;
  dphi_cbo[m]=2.47422e-01;
  m=m+1;
  */
  n[m]=13;
  A_cbo[m]=7.99907e-03;
  dA_cbo[m]=1.79627e-03;
  tau_cbo[m]=1.72217e+04;
  dtau_cbo[m]=1.72785e+03;
  omega_cbo[m]=1.40318e-02;
  domega_cbo[m]=6.00534e-06;
  phi_cbo[m]=1.55125e+01;
  dphi_cbo[m]=2.28892e-01;
  m=m+1;

  n[m]=14;
  A_cbo[m]=9.04143e-03;
  dA_cbo[m]=2.32645e-03;
  tau_cbo[m]=1.63285e+04;
  dtau_cbo[m]=1.83309e+03;
  omega_cbo[m]=1.40213e-02;
  domega_cbo[m]=8.22597e-06;
  phi_cbo[m]=1.59645e+01;
  dphi_cbo[m]=3.13215e-01;
  m=m+1;

  /*  n[m]=15;
   A_cbo[m]=-1.73491e-01;
  dA_cbo[m]=4.55017e-01;
  tau_cbo[m]=2.51225e+08;
  dtau_cbo[m]=2.58505e+07;
  omega_cbo[m]=1.39567e-02;
  domega_cbo[m]=1.38208e-05;
  phi_cbo[m]=1.99650e+01;
  dphi_cbo[m]=1.47555e+00;
  m=m+1;
  */
   n[m]=16;
  A_cbo[m]=8.00205e-03;
  dA_cbo[m]=1.71537e-03;
  tau_cbo[m]=1.71915e+04;
  dtau_cbo[m]=1.64592e+03;
  omega_cbo[m]=1.40314e-02;
  domega_cbo[m]=5.90592e-06;
  phi_cbo[m]=1.55254e+01;
  dphi_cbo[m]=2.24659e-01;
  m=m+1;

   n[m]=17;
  A_cbo[m]=9.22316e-03;
  dA_cbo[m]=2.41653e-03;
  tau_cbo[m]=1.61174e+04;
  dtau_cbo[m]=1.83696e+03;
  omega_cbo[m]=1.40188e-02;
  domega_cbo[m]=8.90373e-06;
  phi_cbo[m]=1.60763e+01;
  dphi_cbo[m]=3.37945e-01;
  m=m+1;

  /*    n[m]=18;
  A_cbo[m]=3.95225e-04;
  dA_cbo[m]=2.94316e-04;
  tau_cbo[m]=3.78993e+09;
  dtau_cbo[m]=6.80684e+08;
  omega_cbo[m]=1.41315e-02;
  domega_cbo[m]=4.86360e-06;
  phi_cbo[m]=1.24225e+01;
  dphi_cbo[m]=6.87409e-01;
  m=m+1;
  */
  n[m]=19;
  A_cbo[m]=7.89273e-03;
  dA_cbo[m]=1.75039e-03;
  tau_cbo[m]=1.71913e+04;
  dtau_cbo[m]=1.70496e+03;
  omega_cbo[m]=1.40325e-02;
  domega_cbo[m]=5.85686e-06;
  phi_cbo[m]=1.54882e+01;
  dphi_cbo[m]=2.22331e-01;
  m=m+1;
  
  n[m]=20;
  A_cbo[m]=9.43695e-03;
  dA_cbo[m]=2.40735e-03;
  tau_cbo[m]=1.58895e+04;
  dtau_cbo[m]=1.75451e+03;
  omega_cbo[m]=1.40194e-02;
  domega_cbo[m]=9.27998e-06;
  phi_cbo[m]=1.60628e+01;
  dphi_cbo[m]=3.50716e-01;
  m=m+1;
  
  
   
    c1=new TCanvas("c1","horizontal cbo vs nbinshift");
    c1->Divide(2,2);
    c1->cd(1);
    gr1=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr1->SetTitle("vw freq vs nbinshift");
    gr1->GetXaxis()->SetTitle("nbinshift");
    gr1->GetYaxis()->SetTitle("vw_cbo");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    
    c1->cd(2);
    gr2=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr2->SetTitle("vw lifetime vs nbinshift");
    gr2->GetXaxis()->SetTitle("nbinshift");
    gr2->GetYaxis()->SetTitle("tau_vw");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw();

    c1->cd(3);
    gr3=new TGraphErrors(m,n,A_cbo,0,dA_cbo);
    gr3->SetTitle("vw amplitude vs nbinshift");
    gr3->GetXaxis()->SetTitle("nbinshift");
    gr3->GetYaxis()->SetTitle("A_vw");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw();

    c1->cd(4);
    gr4=new TGraphErrors(m,n,phi_cbo,0,dphi_cbo);
    gr4->SetTitle("vw phase vs nbinshift");
    gr4->GetXaxis()->SetTitle("nbinshift");
    gr4->GetYaxis()->SetTitle("phi_vw");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->Draw();

    
}
