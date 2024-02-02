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

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_odd[25], *h_even[25], *h_diff[25], *h_average[25], *stat_hist[25], *hstat[25];
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res;
TH1 *hm;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t fit_start, fit_stop;

char root_file_name1[128] = "run2c_dqc_thresh_300.root";
char root_file_name2[128] = "run2c_dqc_24727_24739.root";
char root_file_name3[128] = "run2d_dqc_26013_100sr.root";
char root_file_name4[128] = "run2d_dqc_26013_subrun136.root";
char root_file_name5[128] = "run2c_1fill.root";


void diff_nonlin_statnorm()

 {
   Double_t chi[50],life[50],dlife[50],freq[50],dfreq[50],norm[50],dnorm[50],phi[50],dphi[50], p0[50], dp0[50], average_trace[50], diff_trace[50], p1[50], p2[50], p3[50], p4[50], p5[50];
   Double_t n[50], n1[50], n2[50], n3[50], n4[50], n5[50];
   Int_t m=0;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[1]=TFile::Open(root_file_name1);
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
  dir[1]->GetObject("Hseqnumstats",stat_hist[1]);
  hstat[1]=(TH1D*)stat_hist[1]->Clone();

  _file[2]=TFile::Open(root_file_name2);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_1_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();
   dir[2]->GetObject("Hseqnumstats",stat_hist[2]);
  hstat[2]=(TH1D*)stat_hist[2]->Clone();


   _file[3]=TFile::Open(root_file_name3);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_1_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();
   dir[3]->GetObject("Hseqnumstats",stat_hist[3]);
  hstat[3]=(TH1D*)stat_hist[3]->Clone();


   _file[4]=TFile::Open(root_file_name4);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_1_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();
   dir[4]->GetObject("Hseqnumstats",stat_hist[4]);
  hstat[4]=(TH1D*)stat_hist[4]->Clone();


   _file[5]=TFile::Open(root_file_name5);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_1_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();
   dir[5]->GetObject("Hseqnumstats",stat_hist[5]);
  hstat[5]=(TH1D*)stat_hist[5]->Clone();


   
  gStyle->SetOptFit(1111);
 
  // c1->Divide(4,6);


   
  for(Int_t j=1; j<=5; j++)
    {
  
	 h_odd[j]= new TH1D("raw trace h_odd", "raw trace h_odd", h[j]->GetNbinsX()/2, 100001, 352001);
    h_even[j]= new TH1D("raw trace h_even", "raw trace h_even", h[j]->GetNbinsX()/2, 100001, 352001);
    h_diff[j]= new TH1D("raw trace odd-even diff", "h_diff", h[j]->GetNbinsX()/2, 100001, 352001);
    h_average[j]= new TH1D("raw trace odd-even average", "h_average", h[j]->GetNbinsX()/2, 100001, 352001);
     for(Int_t i=1; i<=h[j]->GetNbinsX(); i++)
       {
        if(i%2==0)
	 {
	   h_even[j]->SetBinContent(i/2, h[j]->GetBinContent(i));
	 }
        else
	 {
           h_odd[j]->SetBinContent((i+1)/2, h[j]->GetBinContent(i));
	 }
        }
    } 
  //   cout<<"j ="<<j<<endl;
     //return;
   for(Int_t j=1; j<=5; j++)
      {	
     for(Int_t i=1; i<=h[j]->GetNbinsX()/2; i++)
       {

	 h_diff[j]->SetBinContent(i,(h_odd[j]->GetBinContent(i)-h_even[j]->GetBinContent(i)));
	 h_average[j]->SetBinContent(i,(h_even[j]->GetBinContent(i)+h_odd[j]->GetBinContent(i))/2);
     }

      }

          
   for(Int_t k=1; k<=5;k++)
     {
       p1[m]=h_diff[k]->Integral(6667,8334)/(2*hstat[k]->Integral());
       
	 n1[m]=k;
         m=m+1;
     }
  
    c2=new TCanvas("c2","odd even difference parameter");   
      gr1=new TGraphErrors(m,n1,p1,0,0);
    gr1->SetTitle("normalized difference vs dataset size ");
    gr1->GetXaxis()->SetTitle("dataset size");
    gr1->GetYaxis()->SetTitle("Inetgral(odd-even)/no. of fills");
    gr1->SetLineWidth(7);
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kBlue);
    gr1->Draw();
    
    
       
 }
