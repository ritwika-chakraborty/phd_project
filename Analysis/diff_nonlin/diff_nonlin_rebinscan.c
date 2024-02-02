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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_odd[25], *h_even[25], *h_diff[25], *h_average[25], *stat_hist[25], *hstat[25], *hr[5];
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

//char root_file_name1[128] = "run2d_dqc_thresh_300.root";
char root_file_name1[128] = "run2d_1000sr_thresh_300.root";
char root_file_name2[128] = "run2d_1000sr_thresh_225.root";
char root_file_name3[128] = "run2d_1000sr_thresh_150.root";
char root_file_name4[128] = "run2d_1000sr_thresh_80.root";
char root_file_name5[128] = "run2d_1000sr_thresh_10.root";


void diff_nonlin_rebinscan()

 {
   Double_t chi[50],life[50],dlife[50],freq[50],dfreq[50],norm[50],dnorm[50],phi[50],dphi[50], p0[50], dp0[50], average_trace[50], diff_trace[50], p1[50], p2[50], p3[50], p4[50], p5[50], dp1[50];
   Double_t n[50], n1[50], n2[50], n3[50], n4[50], n5[50];
   Int_t m=0;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[1]=TFile::Open(root_file_name1);
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
  dir[1]->GetObject("Hseqnumstats",stat_hist[1]);
  hstat[1]=(TH1D*)stat_hist[1]->Clone();

  _file[2]=TFile::Open(root_file_name2);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_1_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();
   dir[2]->GetObject("Hseqnumstats",stat_hist[2]);
  hstat[2]=(TH1D*)stat_hist[2]->Clone();


   _file[3]=TFile::Open(root_file_name3);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_1_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();
   dir[3]->GetObject("Hseqnumstats",stat_hist[3]);
  hstat[3]=(TH1D*)stat_hist[3]->Clone();


   _file[4]=TFile::Open(root_file_name4);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_1_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();
   dir[4]->GetObject("Hseqnumstats",stat_hist[4]);
  hstat[4]=(TH1D*)stat_hist[4]->Clone();


   _file[5]=TFile::Open(root_file_name5);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_1_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();
   dir[5]->GetObject("Hseqnumstats",stat_hist[5]);
  hstat[5]=(TH1D*)stat_hist[5]->Clone();


   
  gStyle->SetOptFit(1111);
 
  // c1->Divide(4,6);

  hr[1]= new TH1D("wiggle rebin 1", "wiggle rebin 1", h[1]->GetNbinsX(), 100001, 352001);
  for(Int_t i=1;i<=hr[1]->GetNbinsX();i++)
    {
      hr[1]->SetBinContent(i,h[1]->GetBinContent(i));
      hr[1]->SetBinError(i,h[1]->GetBinError(i));
    }
  h[1]->Rebin(2); //30
  h[1]->Scale(0.5);
  hr[2]= new TH1D("wiggle rebin 2", "wiggle rebin 2", h[1]->GetNbinsX(), 100001, 352001);
   for(Int_t i=1;i<=hr[2]->GetNbinsX();i++)
    {
      hr[2]->SetBinContent(i,h[1]->GetBinContent(i));
      hr[2]->SetBinError(i,h[1]->GetBinError(i));
    }

  h[1]->Rebin(2); //60
  h[1]->Scale(0.5);
  hr[3]= new TH1D("wiggle rebin 4", "wiggle rebin 4", h[1]->GetNbinsX(), 100001, 352001);
   for(Int_t i=1;i<=hr[3]->GetNbinsX();i++)
    {
      hr[3]->SetBinContent(i,h[1]->GetBinContent(i));
      hr[3]->SetBinError(i,h[1]->GetBinError(i));
    }

  h[1]->Rebin(2); //120
  h[1]->Scale(0.5);
  hr[4]= new TH1D("wiggle rebin 8", "wiggle rebin 8", h[1]->GetNbinsX(), 100001, 352001);
   for(Int_t i=1;i<=hr[4]->GetNbinsX();i++)
    {
      hr[4]->SetBinContent(i,h[1]->GetBinContent(i));
      hr[4]->SetBinError(i,h[1]->GetBinError(i));
    }



  int k=1;
   
  for(Int_t j=1; j<=4; j++)
    {
  
      h_odd[j]= new TH1D("raw trace h_odd", "raw trace h_odd", hr[j]->GetNbinsX()/2, 100001, 352001);
      h_even[j]= new TH1D("raw trace h_even", "raw trace h_even", hr[j]->GetNbinsX()/2, 100001, 352001);
      h_diff[j]= new TH1D("raw trace odd-even diff", "h_diff", hr[j]->GetNbinsX()/2, 100001, 352001);
      h_average[j]= new TH1D("raw trace odd-even average", "h_average", hr[j]->GetNbinsX()/2, 100001, 352001);
      for(Int_t i=1; i<=hr[j]->GetNbinsX(); i++)
       {
        if(i%2==0)
	 {
	   h_even[j]->SetBinContent(i/2, hr[j]->GetBinContent(i));
	   h_even[j]->SetBinError(i/2, hr[j]->GetBinError(i));
	 }
        else
	 {
           h_odd[j]->SetBinContent((i+1)/2, hr[j]->GetBinContent(i));
	   h_odd[j]->SetBinError((i+1)/2, hr[j]->GetBinError(i));
	 }
        }
     k=k*2;
    } 
  //   cout<<"j ="<<j<<endl;
     //return;
  for(Int_t j=1; j<=4; j++)
    {
     for(Int_t i=1; i<=hr[j]->GetNbinsX()/2; i++)
       {

	 h_diff[j]->SetBinContent(i,(h_even[j]->GetBinContent(i)-h_odd[j]->GetBinContent(i)));
	 h_diff[j]->SetBinError(i,sqrt(h_even[j]->GetBinError(i)*h_even[j]->GetBinError(i)+h_odd[j]->GetBinError(i)*h_odd[j]->GetBinError(i)));
	 h_average[j]->SetBinContent(i,(h_even[j]->GetBinContent(i)+h_odd[j]->GetBinContent(i))/2);
     }
    }
 
  double err;
          
       for(Int_t k=1; k<=4;k++)
     {
       //  p1[m]=h_diff[k]->Integral(6667,8334)/(2*54*1668*hstat[k]->Integral());
       p1[m]=h_diff[k]->Integral(h_diff[k]->FindBin(300000),h_diff[k]->FindBin(350000))/ (1-h_diff[k]->FindBin(300000)+h_diff[k]->FindBin(350000));
       err=0;
       for(int l=h_diff[k]->FindBin(300000);l<=h_diff[k]->FindBin(350000);l++)
	 {
	   err=err+(h_diff[k]->GetBinError(l)*h_diff[k]->GetBinError(l));
	 }
       dp1[m]=sqrt(err)/ (1-h_diff[k]->FindBin(300000)+h_diff[k]->FindBin(350000));
       if(m==0){n1[m]=1;}
       if(m==1){n1[m]=2;}
       if(m==2){n1[m]=4;}
       if(m==3){n1[m]=8;}
         m=m+1;
     }
    
    c2=new TCanvas("c2","odd even difference parameter");   
      gr1=new TGraphErrors(m,n1,p1,0,dp1);
    gr1->SetTitle("Picket fence amplitude vs rebin factor");
    gr1->GetXaxis()->SetTitle("Rebin factor");
    gr1->GetYaxis()->SetTitle("Average amplitude [MeV]");
    gr1->SetLineWidth(7);
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kBlue);
    gr1->Draw();
    
    
     
 }
