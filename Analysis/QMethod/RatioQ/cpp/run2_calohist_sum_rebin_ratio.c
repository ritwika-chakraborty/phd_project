#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"


TCanvas *c2, *c3;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[50], *h[50], *h_sum, *h_sum2;
Double_t fit_start, fit_stop;
//TH1 *hm;

void run2_calohist_sum_rebin_ratio()

 {
   /* Double_t chi[40],cal[40],norm[40],dcal[40];
   Double_t n[40];
   Int_t m=0;
   Double_t N;*/

 c2=new TCanvas("c2","5 parameter wiggle_fit");
 
 TFile *_file[50];
 TDirectoryFile *dir[50];

  _file[1]=TFile::Open("run2_982sr_corr.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2_982sr_corr.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2_982sr_corr.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2_982sr_corr.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2_982sr_corr.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2_982sr_corr.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open("run2_982sr_corr.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open("run2_982sr_corr.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open("run2_982sr_corr.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open("run2_982sr_corr.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open("run2_982sr_corr.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open("run2_982sr_corr.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open("run2_982sr_corr.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open("run2_982sr_corr.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open("run2_982sr_corr.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open("run2_982sr_corr.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open("run2_982sr_corr.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open("run2_982sr_corr.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open("run2_982sr_corr.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open("run2_982sr_corr.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open("run2_982sr_corr.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open("run2_982sr_corr.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open("run2_982sr_corr.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open("run2_982sr_corr.root");
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_sig_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=24; i++)
    {
      h_sum->Add(h[i],1); 
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;
   
  _file[25]=TFile::Open("run2_982sr_corr.root");
  _file[25]->GetObject("QFillByFillAnalyzer",dir[25]);
  dir[25]->GetObject("qHist1D_sig_1_2",qHist_1D[25]);
  h[25]=(TH1D*)qHist_1D[25]->Clone();

  _file[26]=TFile::Open("run2_982sr_corr.root");
  _file[26]->GetObject("QFillByFillAnalyzer",dir[26]);
  dir[26]->GetObject("qHist1D_sig_2_2",qHist_1D[26]);
  h[26]=(TH1D*)qHist_1D[26]->Clone();

   _file[27]=TFile::Open("run2_982sr_corr.root");
  _file[27]->GetObject("QFillByFillAnalyzer",dir[27]);
  dir[27]->GetObject("qHist1D_sig_3_2",qHist_1D[27]);
  h[27]=(TH1D*)qHist_1D[27]->Clone();

   _file[28]=TFile::Open("run2_982sr_corr.root");
  _file[28]->GetObject("QFillByFillAnalyzer",dir[28]);
  dir[28]->GetObject("qHist1D_sig_4_2",qHist_1D[28]);
  h[28]=(TH1D*)qHist_1D[28]->Clone();

   _file[29]=TFile::Open("run2_982sr_corr.root");
  _file[29]->GetObject("QFillByFillAnalyzer",dir[29]);
  dir[29]->GetObject("qHist1D_sig_5_2",qHist_1D[29]);
  h[29]=(TH1D*)qHist_1D[29]->Clone();

   _file[30]=TFile::Open("run2_982sr_corr.root");
  _file[30]->GetObject("QFillByFillAnalyzer",dir[30]);
  dir[30]->GetObject("qHist1D_sig_6_2",qHist_1D[30]);
  h[30]=(TH1D*)qHist_1D[30]->Clone();

   _file[31]=TFile::Open("run2_982sr_corr.root");
  _file[31]->GetObject("QFillByFillAnalyzer",dir[31]);
  dir[31]->GetObject("qHist1D_sig_7_2",qHist_1D[31]);
  h[31]=(TH1D*)qHist_1D[31]->Clone();

   _file[32]=TFile::Open("run2_982sr_corr.root");
  _file[32]->GetObject("QFillByFillAnalyzer",dir[32]);
  dir[32]->GetObject("qHist1D_sig_8_2",qHist_1D[32]);
  h[32]=(TH1D*)qHist_1D[32]->Clone();

   _file[33]=TFile::Open("run2_982sr_corr.root");
  _file[33]->GetObject("QFillByFillAnalyzer",dir[33]);
  dir[33]->GetObject("qHist1D_sig_9_2",qHist_1D[33]);
  h[33]=(TH1D*)qHist_1D[33]->Clone();

   _file[34]=TFile::Open("run2_982sr_corr.root");
  _file[34]->GetObject("QFillByFillAnalyzer",dir[34]);
  dir[34]->GetObject("qHist1D_sig_10_2",qHist_1D[34]);
  h[34]=(TH1D*)qHist_1D[34]->Clone();

   _file[35]=TFile::Open("run2_982sr_corr.root");
  _file[35]->GetObject("QFillByFillAnalyzer",dir[35]);
  dir[35]->GetObject("qHist1D_sig_11_2",qHist_1D[35]);
  h[35]=(TH1D*)qHist_1D[35]->Clone();

   _file[36]=TFile::Open("run2_982sr_corr.root");
  _file[36]->GetObject("QFillByFillAnalyzer",dir[36]);
  dir[36]->GetObject("qHist1D_sig_12_2",qHist_1D[36]);
  h[36]=(TH1D*)qHist_1D[36]->Clone();

   _file[37]=TFile::Open("run2_982sr_corr.root");
  _file[37]->GetObject("QFillByFillAnalyzer",dir[37]);
  dir[37]->GetObject("qHist1D_sig_13_2",qHist_1D[37]);
  h[37]=(TH1D*)qHist_1D[37]->Clone();

   _file[38]=TFile::Open("run2_982sr_corr.root");
  _file[38]->GetObject("QFillByFillAnalyzer",dir[38]);
  dir[38]->GetObject("qHist1D_sig_14_2",qHist_1D[38]);
  h[38]=(TH1D*)qHist_1D[38]->Clone();

   _file[39]=TFile::Open("run2_982sr_corr.root");
  _file[39]->GetObject("QFillByFillAnalyzer",dir[39]);
  dir[39]->GetObject("qHist1D_sig_15_2",qHist_1D[39]);
  h[39]=(TH1D*)qHist_1D[39]->Clone();

   _file[40]=TFile::Open("run2_982sr_corr.root");
  _file[40]->GetObject("QFillByFillAnalyzer",dir[40]);
  dir[40]->GetObject("qHist1D_sig_16_2",qHist_1D[40]);
  h[40]=(TH1D*)qHist_1D[40]->Clone();

   _file[41]=TFile::Open("run2_982sr_corr.root");
  _file[41]->GetObject("QFillByFillAnalyzer",dir[41]);
  dir[41]->GetObject("qHist1D_sig_17_2",qHist_1D[41]);
  h[41]=(TH1D*)qHist_1D[41]->Clone();

   _file[42]=TFile::Open("run2_982sr_corr.root");
  _file[42]->GetObject("QFillByFillAnalyzer",dir[42]);
  dir[42]->GetObject("qHist1D_sig_18_2",qHist_1D[42]);
  h[42]=(TH1D*)qHist_1D[42]->Clone();

   _file[43]=TFile::Open("run2_982sr_corr.root");
  _file[43]->GetObject("QFillByFillAnalyzer",dir[43]);
  dir[43]->GetObject("qHist1D_sig_19_2",qHist_1D[43]);
  h[43]=(TH1D*)qHist_1D[43]->Clone();

   _file[44]=TFile::Open("run2_982sr_corr.root");
  _file[44]->GetObject("QFillByFillAnalyzer",dir[44]);
  dir[44]->GetObject("qHist1D_sig_20_2",qHist_1D[44]);
  h[44]=(TH1D*)qHist_1D[44]->Clone();

   _file[45]=TFile::Open("run2_982sr_corr.root");
  _file[45]->GetObject("QFillByFillAnalyzer",dir[45]);
  dir[45]->GetObject("qHist1D_sig_21_2",qHist_1D[45]);
  h[45]=(TH1D*)qHist_1D[45]->Clone();

   _file[46]=TFile::Open("run2_982sr_corr.root");
  _file[46]->GetObject("QFillByFillAnalyzer",dir[46]);
  dir[46]->GetObject("qHist1D_sig_22_2",qHist_1D[46]);
  h[46]=(TH1D*)qHist_1D[46]->Clone();

   _file[47]=TFile::Open("run2_982sr_corr.root");
  _file[47]->GetObject("QFillByFillAnalyzer",dir[47]);
  dir[47]->GetObject("qHist1D_sig_23_2",qHist_1D[47]);
  h[47]=(TH1D*)qHist_1D[47]->Clone();

   _file[48]=TFile::Open("run2_982sr_corr.root");
  _file[48]->GetObject("QFillByFillAnalyzer",dir[48]);
  dir[48]->GetObject("qHist1D_sig_24_2",qHist_1D[48]);
  h[48]=(TH1D*)qHist_1D[48]->Clone();

  Double_t sum2=0;
  
 for(Int_t k=25; k<=48; k++)
   {
    sum2=sum2+h[k]->GetBinContent(500);
   }

 cout<<sum2<<" "<<endl;
   c3= new TCanvas("c3", "difference");  
  gStyle->SetOptFit(1111);
 
   
  h_sum2= new TH1D("calo histogram sum", "h_sum", h[25]->GetNbinsX(), 100001, 352001);
  h_sum2->Sumw2(kTRUE);
    
   for(Int_t i=25; i<=48; i++)
    {
      h_sum2->Add(h[i],1); 
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   //   h_sum2->Draw();
   h_sum->Rebin(4);
   h_sum->Divide(h_sum2);
   h_sum->Draw();
   

   cout<<h_sum2->GetBinContent(500)<<" "<<endl;
 
 
  
 }
