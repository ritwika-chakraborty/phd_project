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

void run2_cov_graph()
{
  Double_t chi[40], time[40], n[40];
  Int_t m=0;

  n[m]=0 ;
  time[m]=107.27 ;
  chi[m]=1.181 ;
  m=m+1;

   n[m]=1 ;
  time[m]=300.5 ;
  chi[m]=1.105 ;
  m=m+1;

   n[m]=2 ;
  time[m]=144.97 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=3 ;
  time[m]=148.46 ;
  chi[m]=1.106;
  m=m+1;

   n[m]=4 ;
  time[m]=159.17 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=5 ;
  time[m]=154.09 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=10 ;
  time[m]=188.88 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=15 ;
  time[m]=214.05 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=20 ;
  time[m]=244.24 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=40 ;
  time[m]=334.49 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=60 ;
  time[m]=437.08 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=80 ;
  time[m]=532.57 ;
  chi[m]=1.106 ;
  m=m+1;

   n[m]=100 ;
  time[m]=618.22 ;
  chi[m]=1.106 ;
  m=m+1;
  
    c1=new TCanvas("c1","chi-square, time to fit, number of rows on either side of diagonals");
    c1->Divide(1,2);
      c1->cd(1);
    gr1=new TGraphErrors(m,n,time,0,0);
    gr1->SetTitle("time vs diag+n");
    gr1->GetXaxis()->SetTitle("diag+n");
    gr1->GetYaxis()->SetTitle("time (s)");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    
    c1->cd(2);
     gr2=new TGraphErrors(m,n,chi,0,0);
    gr2->SetTitle("chi-square vs diag+n");
    gr2->GetXaxis()->SetTitle("diag+n");
    gr2->GetYaxis()->SetTitle("chi-square");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kBlue);
    gr2->Draw();
}
