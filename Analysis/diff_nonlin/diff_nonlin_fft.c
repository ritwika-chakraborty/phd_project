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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_odd[25], *h_even[25], *h_diff[25], *h_average[25], *stat_hist[25], *hstat[25], *h_diff_fft[25];
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res;
TH1 *hm[6];
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;

char root_file_name1[128] = "run2c_dqc_thresh_300.root";
char root_file_name2[128] = "run2c_dqc_24727_24739.root";
char root_file_name3[128] = "run2d_dqc_26013_100sr.root";
char root_file_name4[128] = "run2d_dqc_26013_subrun136.root";
char root_file_name5[128] = "run2c_1fill.root";
char root_file_name6[128] = "run2c_dqc_thesh_300_oldcongif.root";


void diff_nonlin_fft()

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
  dir[1]->GetObject("qHist1D_3_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
  dir[1]->GetObject("Hseqnumstats",stat_hist[1]);
  hstat[1]=(TH1D*)stat_hist[1]->Clone();

  _file[2]=TFile::Open(root_file_name2);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_3_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();
   dir[2]->GetObject("Hseqnumstats",stat_hist[2]);
  hstat[2]=(TH1D*)stat_hist[2]->Clone();


   _file[3]=TFile::Open(root_file_name3);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();
   dir[3]->GetObject("Hseqnumstats",stat_hist[3]);
  hstat[3]=(TH1D*)stat_hist[3]->Clone();


   _file[4]=TFile::Open(root_file_name4);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_3_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();
   dir[4]->GetObject("Hseqnumstats",stat_hist[4]);
  hstat[4]=(TH1D*)stat_hist[4]->Clone();


   _file[5]=TFile::Open(root_file_name5);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_3_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();
   dir[5]->GetObject("Hseqnumstats",stat_hist[5]);
  hstat[5]=(TH1D*)stat_hist[5]->Clone();


   
  gStyle->SetOptFit(1111);

  double binwidth=h[1]->GetBinWidth(1);
  int nbins=h[1]->GetNbinsX();
  double binlowedge=h[1]->GetBinLowEdge(1);
  double binhighedge=h[1]->GetBinLowEdge(nbins)+binwidth;

      
  // c1->Divide(4,6);
 for(Int_t j=1; j<=5; j++)
    {
      
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

      h[j]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

    }

   
  for(Int_t j=1; j<=5; j++)
    {
  
	 h_odd[j]= new TH1D("raw trace h_odd", "raw trace h_odd", h[j]->GetNbinsX()/2, h[j]->GetBinLowEdge(1), h[j]->GetBinLowEdge(h[j]->GetNbinsX())+h[j]->GetBinWidth(1));
    h_even[j]= new TH1D("raw trace h_even", "raw trace h_even", h[j]->GetNbinsX()/2, h[j]->GetBinLowEdge(1),  h[j]->GetBinLowEdge(h[j]->GetNbinsX())+h[j]->GetBinWidth(1));
    h_diff[j]= new TH1D("raw trace odd-even diff", "h_diff", h[j]->GetNbinsX()/2, h[j]->GetBinLowEdge(1),  h[j]->GetBinLowEdge(h[j]->GetNbinsX())+h[j]->GetBinWidth(1));
    h_average[j]= new TH1D("raw trace odd-even average", "h_average", h[j]->GetNbinsX()/2, h[j]->GetBinLowEdge(1),  h[j]->GetBinLowEdge(h[j]->GetNbinsX())+h[j]->GetBinWidth(1));
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
   
   //fourier transform of the raw histogram
   /*  for(Int_t j=1; j<=5; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334*2, 255000, 305025);
       for(Int_t i=13920; i<=16587; i++)
       {
         h_diff_fft[j]->SetBinContent(i-13919,h[j]->GetBinContent(i));
	 
       }

      }
   */
   //fourier transform of the difference histogram
   /*    for(Int_t j=1; j<=5; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334, 255000, 305025);
       for(Int_t i=6960; i<=8293; i++)
       {
         h_diff_fft[j]->SetBinContent(i-6959,h_diff[j]->GetBinContent(i));
	 
       }

       }*/
     //fourier transform of the rebiined difference
        for(Int_t j=1; j<=5; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334/2, 255000, 305025);
       h_diff[j]->Rebin(2);
       for(Int_t i=3480; i<=4146; i++)
       {
         h_diff_fft[j]->SetBinContent(i-3479,h_diff[j]->GetBinContent(i));
	 
       }

       }
   

     //fourier transform of the rebiined difference
   /*    for(Int_t j=1; j<=5; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334/8, 255000, 305025);
       h_diff[j]->Rebin(4);
       for(Int_t i=1740; i<=2073; i++)
       {
         h_diff_fft[j]->SetBinContent(i-1739,h_diff[j]->GetBinContent(i));
	 
       }

       }
   */
          
   /* for(Int_t k=1; k<=5;k++)
     {
       p1[m]=h_diff[k]->Integral(6667,8334)/(hstat[k]->Integral());
       
	 n1[m]=k;
         m=m+1;
     }
   */
    c2=new TCanvas("c2","odd even difference parameter");
    c2->Divide(1,5);

     for(Int_t j=1; j<=5; j++)
      {

         hm[j]=h_diff_fft[j]->FFT(hm[j],"MAG");
	 hm[j]->SetBins(h_diff_fft[j]->GetNbinsX(),0,1/h_diff_fft[j]->GetBinWidth(1));
	 //	 hm[j]->SetBins(h_diff_fft[j]->GetNbinsX(),0,h_diff_fft[j]->GetNbinsX()/((h_diff_fft[j]->GetBinLowEdge(h_diff_fft[j]->GetNbinsX())+h_diff_fft[j]->GetBinWidth(1))));
	 //   hm[j]->GetXaxis()->SetRangeUser(0,1/(1.99*h_diff_fft[j]->GetBinWidth(1)));
         hm[j]->GetXaxis()->SetTitle("Freq [GHz]");
         hm[j]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
 	 hm[j]->SetLineColor(kBlack);
      	 hm[j]->GetYaxis()->SetRangeUser(0,350000000);
	 c2->cd(j);
	 hm[j]->Draw();


      }
     /*      hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm->GetXaxis()->SetTitle("Freq [GHz]");
   hm->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm->Draw();     
     */
    /*  gr1=new TGraphErrors(m,n1,p1,0,0);
    gr1->SetTitle("normalized difference vs dataset size ");
    gr1->GetXaxis()->SetTitle("dataset size");
    gr1->GetYaxis()->SetTitle("Inetgral(odd-even)/no. of fills");
    gr1->SetLineWidth(7);
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kBlue);
    gr1->Draw();
    
    */
       
 }
