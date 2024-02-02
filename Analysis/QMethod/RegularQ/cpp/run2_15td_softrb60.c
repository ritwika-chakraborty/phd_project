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
TH1D *qHist_1D[25], *h[25], *h_sum, *h_diff;
Double_t fit_start, fit_stop;
//TH1 *hm;

void run2_15td_softrb60()

 {
   /* Double_t chi[40],cal[40],norm[40],dcal[40];
   Double_t n[40];
   Int_t m=0;
   Double_t N;*/

 c2=new TCanvas("c2","5 parameter wiggle_fit");
 
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[1]=TFile::Open("run2_982sr_corr.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_2",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2_982sr_corr.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_2",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2_982sr_corr.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2_982sr_corr.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_2",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2_982sr_corr.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_2",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2_982sr_corr.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_2",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open("run2_982sr_corr.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_2",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open("run2_982sr_corr.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_2",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open("run2_982sr_corr.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_2",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open("run2_982sr_corr.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_2",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open("run2_982sr_corr.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_2",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open("run2_982sr_corr.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_2",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open("run2_982sr_corr.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_2",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open("run2_982sr_corr.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_2",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open("run2_982sr_corr.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_2",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open("run2_982sr_corr.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_2",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open("run2_982sr_corr.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_2",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open("run2_982sr_corr.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_2",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open("run2_982sr_corr.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_2",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open("run2_982sr_corr.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_2",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open("run2_982sr_corr.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_2",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open("run2_982sr_corr.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_2",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open("run2_982sr_corr.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_2",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open("run2_982sr_corr.root");
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_sig_24_2",qHist_1D[24]);
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
    fit_start=140000;
   fit_stop=300000;
 }
