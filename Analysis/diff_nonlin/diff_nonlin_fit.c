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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *h_odd[25], *h_even[25], *h_diff[25], *h_average[25], *stat_hist[25], *h_diff_fft[25], *hstat;
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

char root_file_name[128] = "run2c_dqc_thresh_300.root";
char root_file_name2[128] = "run2c_dqc_24727_24739.root";
char root_file_name3[128] = "run2d_dqc_26013_100sr.root";
char root_file_name4[128] = "run2d_dqc_26013_subrun136.root";
char root_file_name5[128] = "run2c_1fill.root";

Double_t fpicket(Double_t *x, Double_t *par)
{
  Double_t f_x=par[0]*(1+par[1]*sin(par[2]*x[0]+par[3]))*(1+par[4]*sin(par[5]*x[0]+par[6]));
  return f_x;
}

void diff_nonlin_fit()

 {
   Double_t chi[50],life[50],dlife[50],freq[50],dfreq[50],norm[50],dnorm[50],phi[50],dphi[50], p0[50], dp0[50], average_trace[50], diff_trace[50], p1[50], p2[50], p3[50], p4[50], p5[50];
   Double_t n[50], n1[50], n2[50], n3[50], n4[50], n5[50], A1[50], dA1[50], A2[50], dA2[50], phi1[50], dphi1[50], phi2[50], dphi2[50];
   Int_t m=0;
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
  dir[1]->GetObject("qHist1D_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();


   
  gStyle->SetOptFit(1111);

  double binwidth=h[1]->GetBinWidth(1);
  int nbins=h[1]->GetNbinsX();
  double binlowedge=h[1]->GetBinLowEdge(1);
  double binhighedge=h[1]->GetBinLowEdge(nbins)+binwidth;

      
  // c1->Divide(4,6);
 for(Int_t j=1; j<=24; j++)
    {
      
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

      h[j]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
      cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

    }

   
  for(Int_t j=1; j<=24; j++)
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
   for(Int_t j=1; j<=24; j++)
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
      for(Int_t j=1; j<=24; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334, 255000, 305025);
       for(Int_t i=6960; i<=8293; i++)
       {
         h_diff_fft[j]->SetBinContent(i-6959,h_diff[j]->GetBinContent(i));
	 
       }

       }
     //fourier transform of the rebiined difference
   /*    for(Int_t j=1; j<=5; j++)
      {
       h_diff_fft[j]= new TH1D("hdiff_fft", "h_diff_fft", 1334/2, 255000, 305025);
       h_diff[j]->Rebin(2);
       for(Int_t i=3480; i<=4146; i++)
       {
         h_diff_fft[j]->SetBinContent(i-3479,h_diff[j]->GetBinContent(i));
	 
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
    c2->Divide(1,3);

    c2->cd(1);
     fit_func= new TF1("fpicket", fpicket,  255000,  305025, 7);
     fit_func->SetParNames("N","A1", "omega_1", "phi1", "A2", "omega_2", "phi2");
 
     // fit_func->FixParameter(2,0.083733333);
     //fit_func->FixParameter(3,11/7);
     //fit_func->SetParLimits(3,-44/7,44/7);
     //fit_func->SetParLimits(6,-44/7,44/7);
     fit_func->SetNpx(10000);
     fit_start=255000;
     fit_stop=305025;
     //   h_diff_fft[1]->Draw();
     //fit_func->Draw("same");
     int j=1;
    fit_func->SetParameters(-62365000, 0.043, 0.083733333, 11/7, 0.01, 0.041904, 0);
    cout<<"fitting calo "<<j<<" "<<endl;
    h_diff_fft[j]->Fit("fpicket","R","",255000,305025);
    A1[j]=fit_func->GetParameter(1);
    phi1[j]=fit_func->GetParameter(3);
    A2[j]=fit_func->GetParameter(4);
    phi2[j]=fit_func->GetParameter(6);

     dA1[j]=fit_func->GetParError(1);
    dphi1[j]=fit_func->GetParError(3);
    dA2[j]=fit_func->GetParError(4);
    dphi2[j]=fit_func->GetParError(6);
    n[j]=j;
    j=j+1;

        fit_func->SetParameters(-11805000, 0.03, 0.083733333, 0, 0.04, 0.041904, 0);
    cout<<"fitting calo "<<j<<" "<<endl;
    //    fit_func->FixParameter(2,0.083733333);
    //fit_func->FixParameter(5,0.041904);
    h_diff_fft[j]->Fit("fpicket","R","",255000,305025);
    //h_diff_fft[j]->Draw();
    // fit_func->Draw("same");
    A1[j]=fit_func->GetParameter(1);
    phi1[j]=fit_func->GetParameter(3);
    A2[j]=fit_func->GetParameter(4);
    phi2[j]=fit_func->GetParameter(6);

     dA1[j]=fit_func->GetParError(1);
    dphi1[j]=fit_func->GetParError(3);
    dA2[j]=fit_func->GetParError(4);
    dphi2[j]=fit_func->GetParError(6);
    n[j]=j;
    j=j+1;


    	h_res= new TH1D("residual histogram", "h_res", h_diff_fft[1]->GetNbinsX(), h_diff_fft[1]->GetBinLowEdge(1), (h_diff_fft[1]->GetBinLowEdge(h_diff_fft[1]->GetNbinsX())+h_diff_fft[1]->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- h_diff_fft[1]->GetBinLowEdge(1))/h_diff_fft[1]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_diff_fft[1]->GetBinLowEdge(1))/h_diff_fft[1]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_diff_fft[1]->GetBinContent(ibin)- fit_func->Eval( h_diff_fft[1]->GetBinCenter(ibin) ) );
      if(h_diff_fft[1]->GetBinError(ibin)!=0){res=(res/h_diff_fft[1]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	
        h_res->Draw();

  c2->cd(3);
     

         hm[1]=h_res->FFT(hm[1],"MAG");
	 	 hm[1]->SetBins(h_res->GetNbinsX(),0,1/h_res->GetBinWidth(1));
		 //    hm[1]->GetXaxis()->SetRangeUser(0,h_diff_fft[1]->GetNbinsX()/((h_diff_fft[1]->GetBinLowEdge(h_diff_fft[1]->GetNbinsX())+h_diff_fft[1]->GetBinWidth(1)))/2);

	 hm[1]->SetLineColor(kBlack);
	 hm[1]->GetYaxis()->SetRangeUser(0,10000);
	 hm[1]->Draw();


     
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
