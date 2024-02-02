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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_odd[25], *h_even[25], *h_diff[25], *h_average[25], *hstat, *h_diff_fft[25];
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res;
TH1 *hm[25];
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;

char root_file_name[128] = "run2c_dqc_thresh_300.root";


void diff_nonlin_fft_hslice()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40], p0[40], dp0[40], average_trace[40], diff_trace[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[0]=TFile::Open(root_file_name);
  _file[0]->GetObject("QFillByFillAnalyzer",dir[0]);
  dir[0]->GetObject("Hseqnumstats",qHist_1D[0]);
  hstat=(TH1D*)qHist_1D[0]->Clone();

 
  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_hslice_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_hslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_hslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_hslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_hslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_hslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

 
  Double_t sum=0;
  
 for(Int_t k=1; k<=6; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
  c1->Divide(1,6);
  
 double binwidth=h[1]->GetBinWidth(1);
  int nbins=h[1]->GetNbinsX();
  double binlowedge=h[1]->GetBinLowEdge(1);
  double binhighedge=h[1]->GetBinLowEdge(nbins)+binwidth;


 for(Int_t j=1; j<=6; j++)
    {
      
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

      h[j]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

    }

   
  for(Int_t j=1; j<=6; j++)
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
   for(Int_t j=1; j<=6; j++)
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
       for(Int_t j=1; j<=6; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334, 255000, 305025);
       for(Int_t i=6960; i<=8293; i++)
       {
         h_diff_fft[j]->SetBinContent(i-6959,h_diff[j]->GetBinContent(i));
	 
       }

       }

     for(Int_t j=1; j<=6; j++)
      {

         hm[j]=h_diff_fft[j]->FFT(hm[j],"MAG");
	 hm[j]->SetBins(h_diff_fft[j]->GetNbinsX(),0,1/h_diff_fft[j]->GetBinWidth(1));
	 //	 hm[j]->SetBins(h_diff_fft[j]->GetNbinsX(),0,h_diff_fft[j]->GetNbinsX()/((h_diff_fft[j]->GetBinLowEdge(h_diff_fft[j]->GetNbinsX())+h_diff_fft[j]->GetBinWidth(1))));
	 //   hm[j]->GetXaxis()->SetRangeUser(0,1/(1.99*h_diff_fft[j]->GetBinWidth(1)));
         hm[j]->GetXaxis()->SetTitle("Freq [GHz]");
         hm[j]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
 	 hm[j]->SetLineColor(kBlack);
      	 hm[j]->GetYaxis()->SetRangeUser(0,300000000);
	 c1->cd(j);
	 hm[j]->Draw();


      }
  
 }
