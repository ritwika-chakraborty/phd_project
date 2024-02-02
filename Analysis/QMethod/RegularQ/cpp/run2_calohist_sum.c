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
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D0[25], *h0[25], *h_sum0;
//TH1 *hm;

void run2_calohist_sum()

 {
   /*Double_t chi[40],cal[40],norm[40],dcal[40];
   Double_t n[40];
   Int_t m=0;
   Double_t N;*/

 c1=new TCanvas("c1","5 parameter wiggle_fit");  
 TFile *_file0[25];
 TDirectoryFile *dir0[25];

  _file0[1]=TFile::Open("run2_982sr_prelim.root");
  _file0[1]->GetObject("QFillByFillAnalyzer",dir0[1]);
  dir0[1]->GetObject("qHist1D_sig_1_0",qHist_1D0[1]);
  h0[1]=(TH1D*)qHist_1D0[1]->Clone();

  _file0[2]=TFile::Open("run2_982sr_prelim.root");
  _file0[2]->GetObject("QFillByFillAnalyzer",dir0[2]);
  dir0[2]->GetObject("qHist1D_sig_2_0",qHist_1D0[2]);
  h0[2]=(TH1D*)qHist_1D0[2]->Clone();

   _file0[3]=TFile::Open("run2_982sr_prelim.root");
  _file0[3]->GetObject("QFillByFillAnalyzer",dir0[3]);
  dir0[3]->GetObject("qHist1D_sig_3_0",qHist_1D0[3]);
  h0[3]=(TH1D*)qHist_1D0[3]->Clone();

   _file0[4]=TFile::Open("run2_982sr_prelim.root");
  _file0[4]->GetObject("QFillByFillAnalyzer",dir0[4]);
  dir0[4]->GetObject("qHist1D_sig_4_0",qHist_1D0[4]);
  h0[4]=(TH1D*)qHist_1D0[4]->Clone();

   _file0[5]=TFile::Open("run2_982sr_prelim.root");
  _file0[5]->GetObject("QFillByFillAnalyzer",dir0[5]);
  dir0[5]->GetObject("qHist1D_sig_5_0",qHist_1D0[5]);
  h0[5]=(TH1D*)qHist_1D0[5]->Clone();

   _file0[6]=TFile::Open("run2_982sr_prelim.root");
  _file0[6]->GetObject("QFillByFillAnalyzer",dir0[6]);
  dir0[6]->GetObject("qHist1D_sig_6_0",qHist_1D0[6]);
  h0[6]=(TH1D*)qHist_1D0[6]->Clone();

   _file0[7]=TFile::Open("run2_982sr_prelim.root");
  _file0[7]->GetObject("QFillByFillAnalyzer",dir0[7]);
  dir0[7]->GetObject("qHist1D_sig_7_0",qHist_1D0[7]);
  h0[7]=(TH1D*)qHist_1D0[7]->Clone();

   _file0[8]=TFile::Open("run2_982sr_prelim.root");
  _file0[8]->GetObject("QFillByFillAnalyzer",dir0[8]);
  dir0[8]->GetObject("qHist1D_sig_8_0",qHist_1D0[8]);
  h0[8]=(TH1D*)qHist_1D0[8]->Clone();

   _file0[9]=TFile::Open("run2_982sr_prelim.root");
  _file0[9]->GetObject("QFillByFillAnalyzer",dir0[9]);
  dir0[9]->GetObject("qHist1D_sig_9_0",qHist_1D0[9]);
  h0[9]=(TH1D*)qHist_1D0[9]->Clone();

   _file0[10]=TFile::Open("run2_982sr_prelim.root");
  _file0[10]->GetObject("QFillByFillAnalyzer",dir0[10]);
  dir0[10]->GetObject("qHist1D_sig_10_0",qHist_1D0[10]);
  h0[10]=(TH1D*)qHist_1D0[10]->Clone();

   _file0[11]=TFile::Open("run2_982sr_prelim.root");
  _file0[11]->GetObject("QFillByFillAnalyzer",dir0[11]);
  dir0[11]->GetObject("qHist1D_sig_11_0",qHist_1D0[11]);
  h0[11]=(TH1D*)qHist_1D0[11]->Clone();

   _file0[12]=TFile::Open("run2_982sr_prelim.root");
  _file0[12]->GetObject("QFillByFillAnalyzer",dir0[12]);
  dir0[12]->GetObject("qHist1D_sig_12_0",qHist_1D0[12]);
  h0[12]=(TH1D*)qHist_1D0[12]->Clone();

   _file0[13]=TFile::Open("run2_982sr_prelim.root");
  _file0[13]->GetObject("QFillByFillAnalyzer",dir0[13]);
  dir0[13]->GetObject("qHist1D_sig_13_0",qHist_1D0[13]);
  h0[13]=(TH1D*)qHist_1D0[13]->Clone();

   _file0[14]=TFile::Open("run2_982sr_prelim.root");
  _file0[14]->GetObject("QFillByFillAnalyzer",dir0[14]);
  dir0[14]->GetObject("qHist1D_sig_14_0",qHist_1D0[14]);
  h0[14]=(TH1D*)qHist_1D0[14]->Clone();

   _file0[15]=TFile::Open("run2_982sr_prelim.root");
  _file0[15]->GetObject("QFillByFillAnalyzer",dir0[15]);
  dir0[15]->GetObject("qHist1D_sig_15_0",qHist_1D0[15]);
  h0[15]=(TH1D*)qHist_1D0[15]->Clone();

   _file0[16]=TFile::Open("run2_982sr_prelim.root");
  _file0[16]->GetObject("QFillByFillAnalyzer",dir0[16]);
  dir0[16]->GetObject("qHist1D_sig_16_0",qHist_1D0[16]);
  h0[16]=(TH1D*)qHist_1D0[16]->Clone();

   _file0[17]=TFile::Open("run2_982sr_prelim.root");
  _file0[17]->GetObject("QFillByFillAnalyzer",dir0[17]);
  dir0[17]->GetObject("qHist1D_sig_17_0",qHist_1D0[17]);
  h0[17]=(TH1D*)qHist_1D0[17]->Clone();

   _file0[18]=TFile::Open("run2_982sr_prelim.root");
  _file0[18]->GetObject("QFillByFillAnalyzer",dir0[18]);
  dir0[18]->GetObject("qHist1D_sig_18_0",qHist_1D0[18]);
  h0[18]=(TH1D*)qHist_1D0[18]->Clone();

   _file0[19]=TFile::Open("run2_982sr_prelim.root");
  _file0[19]->GetObject("QFillByFillAnalyzer",dir0[19]);
  dir0[19]->GetObject("qHist1D_sig_19_0",qHist_1D0[19]);
  h0[19]=(TH1D*)qHist_1D0[19]->Clone();

   _file0[20]=TFile::Open("run2_982sr_prelim.root");
  _file0[20]->GetObject("QFillByFillAnalyzer",dir0[20]);
  dir0[20]->GetObject("qHist1D_sig_20_0",qHist_1D0[20]);
  h0[20]=(TH1D*)qHist_1D0[20]->Clone();

   _file0[21]=TFile::Open("run2_982sr_prelim.root");
  _file0[21]->GetObject("QFillByFillAnalyzer",dir0[21]);
  dir0[21]->GetObject("qHist1D_sig_21_0",qHist_1D0[21]);
  h0[21]=(TH1D*)qHist_1D0[21]->Clone();

   _file0[22]=TFile::Open("run2_982sr_prelim.root");
  _file0[22]->GetObject("QFillByFillAnalyzer",dir0[22]);
  dir0[22]->GetObject("qHist1D_sig_22_0",qHist_1D0[22]);
  h0[22]=(TH1D*)qHist_1D0[22]->Clone();

   _file0[23]=TFile::Open("run2_982sr_prelim.root");
  _file0[23]->GetObject("QFillByFillAnalyzer",dir0[23]);
  dir0[23]->GetObject("qHist1D_sig_23_0",qHist_1D0[23]);
  h0[23]=(TH1D*)qHist_1D0[23]->Clone();

   _file0[24]=TFile::Open("run2_982sr_prelim.root");
  _file0[24]->GetObject("QFillByFillAnalyzer",dir0[24]);
  dir0[24]->GetObject("qHist1D_sig_24_0",qHist_1D0[24]);
  h0[24]=(TH1D*)qHist_1D0[24]->Clone();

  Double_t sum0=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum0=sum0+h0[k]->GetBinContent(500);
   }

 cout<<sum0<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum0= new TH1D("calo histogram sum", "h_sum0", h0[1]->GetNbinsX(), 100001, 352001);
  h_sum0->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=24; i++)
    {
      h_sum0->Add(h0[i],1); 
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_sum0->Draw();

   cout<<h_sum0->GetBinContent(500)<<" "<<endl;
   h0[1]->Rebin(4);
 }
