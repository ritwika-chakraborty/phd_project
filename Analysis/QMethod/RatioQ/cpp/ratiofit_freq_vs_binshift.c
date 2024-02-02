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

void ratiofit_freq_vs_binshift()
{
  Double_t omega_cbo[40], domega_cbo[40], tau_cbo[40], dtau_cbo[40], A_cbo[40], dA_cbo[40], phi_cbo[40], dphi_cbo[40], n[40];
  Int_t m=0;

  n[m]=10;
  A_cbo[m]=2.44064e-03;
  dA_cbo[m]=4.94593e-05;
  tau_cbo[m]=2.45889e+05;
  dtau_cbo[m]=1.47170e+04;
  omega_cbo[m]=2.34017e-03;
  domega_cbo[m]=2.34242e-07;
  phi_cbo[m]=4.95288e-01;
  dphi_cbo[m]=1.96600e-02;
  m=m+1;

  n[m]=11;
  A_cbo[m]=2.44565e-03;
  dA_cbo[m]=5.47659e-05;
  tau_cbo[m]=2.44228e+05;
  dtau_cbo[m]=1.60710e+04;
  omega_cbo[m]=2.34018e-03;
  domega_cbo[m]=2.58537e-07;
  phi_cbo[m]=4.95022e-01;
  dphi_cbo[m]=2.16684e-02;
  m=m+1;

  n[m]=12;
  A_cbo[m]=2.45440e-03;
  dA_cbo[m]=6.51678e-05;
  tau_cbo[m]=2.41943e+05;
  dtau_cbo[m]=1.87439e+04;
  omega_cbo[m]=2.34022e-03;
  domega_cbo[m]=3.06319e-07;
  phi_cbo[m]=4.91351e-01;
  dphi_cbo[m]=2.56251e-02;
  m=m+1;

  n[m]=13;
  A_cbo[m]=2.46005e-03;
  dA_cbo[m]=8.42445e-05;
  tau_cbo[m]=2.40702e+05;
  dtau_cbo[m]=2.39593e+04;
  omega_cbo[m]=2.34025e-03;
  domega_cbo[m]=3.94681e-07;
  phi_cbo[m]=4.89750e-01;
  dphi_cbo[m]=3.29883e-02;
  m=m+1;

  n[m]=14;
  A_cbo[m]=2.47257e-03;
  dA_cbo[m]=1.21707e-04;
  tau_cbo[m]=2.38652e+05;
  dtau_cbo[m]=3.39264e+04;
  omega_cbo[m]=2.34027e-03;
  domega_cbo[m]=5.67011e-07;
  phi_cbo[m]=4.87168e-01;
  dphi_cbo[m]=4.73270e-02;
  m=m+1;

  n[m]=15;
   A_cbo[m]=2.50282e-03;
  dA_cbo[m]=2.06367e-04;
  tau_cbo[m]=2.34629e+05;
  dtau_cbo[m]=5.51229e+04;
  omega_cbo[m]=2.34022e-03;
  domega_cbo[m]=9.52584e-07;
  phi_cbo[m]=4.86548e-01;
  dphi_cbo[m]=7.92699e-02;
  m=m+1;

   n[m]=16;
  A_cbo[m]=2.55291e-03;
  dA_cbo[m]=4.47975e-04;
  tau_cbo[m]=2.31306e+05;
  dtau_cbo[m]=1.13936e+05;
  omega_cbo[m]=2.33933e-03;
  domega_cbo[m]=2.08120e-06;
  phi_cbo[m]=5.40146e-01;
  dphi_cbo[m]=1.72535e-01;
  m=m+1;

  /*  n[m]=17;
  A_cbo[m]=2.02076e-03;
  dA_cbo[m]=5.35134e-04;
  tau_cbo[m]=2.81260e+09;
  dtau_cbo[m]=1.40630e+09;
  omega_cbo[m]=2.33009e-03;
  domega_cbo[m]=4.67164e-06;
  phi_cbo[m]=1.08669e+00;
  dphi_cbo[m]=5.14226e-01;
  m=m+1;

    n[m]=18;
  A_cbo[m]=1.25143e-01;
  dA_cbo[m]=1.87989e-01;
  tau_cbo[m]=3.21457e+10;
  dtau_cbo[m]=5.43369e+09;
  omega_cbo[m]=2.32913e-03;
  domega_cbo[m]=3.07161e-06;
  phi_cbo[m]=3.92841e+00;
  dphi_cbo[m]=7.07597e-01;
  m=m+1;
  
  n[m]=19;
  A_cbo[m]=1.87965e-03;
  dA_cbo[m]=3.60523e-04;
  tau_cbo[m]=5.82923e+09;
  dtau_cbo[m]=1.41421e+00;
  omega_cbo[m]=2.34507e-03;
  domega_cbo[m]=3.24516e-06;
  phi_cbo[m]=2.78852e-01;
  dphi_cbo[m]=3.52955e-01;
  m=m+1;
  */ 
  n[m]=20;
  A_cbo[m]=2.22462e-03;
  dA_cbo[m]=3.00657e-04;
  tau_cbo[m]=3.65043e+05;
  dtau_cbo[m]=2.00793e+05;
  omega_cbo[m]=2.34186e-03;
  domega_cbo[m]=1.45614e-06;
  phi_cbo[m]=4.09004e-01;
  dphi_cbo[m]=1.31924e-01;
  m=m+1;
  
  
   
    c1=new TCanvas("c1","horizontal cbo vs nbinshift");
    c1->Divide(2,2);
    c1->cd(1);
    gr1=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr1->SetTitle("cbo freq vs nbinshift");
    gr1->GetXaxis()->SetTitle("nbinshift");
    gr1->GetYaxis()->SetTitle("omega_cbo");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    
    c1->cd(2);
    gr2=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr2->SetTitle("cbo lifetime vs nbinshift");
    gr2->GetXaxis()->SetTitle("nbinshift");
    gr2->GetYaxis()->SetTitle("tau_cbo");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw();

    c1->cd(3);
    gr3=new TGraphErrors(m,n,A_cbo,0,dA_cbo);
    gr3->SetTitle("cbo amplitude vs nbinshift");
    gr3->GetXaxis()->SetTitle("nbinshift");
    gr3->GetYaxis()->SetTitle("A_cbo");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw();

    c1->cd(4);
    gr4=new TGraphErrors(m,n,phi_cbo,0,dphi_cbo);
    gr4->SetTitle("cbo phase vs nbinshift");
    gr4->GetXaxis()->SetTitle("nbinshift");
    gr4->GetYaxis()->SetTitle("phi_cbo");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->Draw();

    
}
