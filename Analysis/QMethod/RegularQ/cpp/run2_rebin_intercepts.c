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
TH1D *qHist_1D[25], *qHist_1D1[25], *h[25], *h1[25], *h_sum;
//TH1 *hm;

TCanvas *c2, *c3, *c4;
TF1 *f1;
//TH1D *h_res;
//TH1 *hm;
TGraphErrors *gr_int, *gr6;




void run2_rebin_intercepts()

 {
   Double_t intercept[40],dintercept[40];
   long double slope[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

   f1 = new TF1("f1","[0]+[1]*x",100001,352001);
   f1->SetParameters(1,0);

  _file[1]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
  dir[1]->GetObject("qHist1D_sig_1_2",qHist_1D1[1]);
  h1[1]=(TH1D*)qHist_1D1[1]->Clone();
  h[1]->Rebin(4);
  h[1]->Divide(h1[1]);
  h[1]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

  _file[2]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();
  dir[2]->GetObject("qHist1D_sig_2_2",qHist_1D1[2]);
  h1[2]=(TH1D*)qHist_1D1[2]->Clone();
  h[2]->Rebin(4);
  h[2]->Divide(h1[2]);
  h[2]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[3]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();
  dir[3]->GetObject("qHist1D_sig_3_2",qHist_1D1[3]);
  h1[3]=(TH1D*)qHist_1D1[3]->Clone();
  h[3]->Rebin(4);
  h[3]->Divide(h1[3]);
  h[3]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[4]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();
  dir[4]->GetObject("qHist1D_sig_4_2",qHist_1D1[4]);
  h1[4]=(TH1D*)qHist_1D1[4]->Clone();
  h[4]->Rebin(4);
  h[4]->Divide(h1[4]);
  h[4]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[5]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();
   dir[5]->GetObject("qHist1D_sig_5_2",qHist_1D1[5]);
   h1[5]=(TH1D*)qHist_1D1[5]->Clone();
  h[5]->Rebin(4);
  h[5]->Divide(h1[5]);
  h[5]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[6]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();
  dir[6]->GetObject("qHist1D_sig_6_2",qHist_1D1[6]);
  h1[6]=(TH1D*)qHist_1D1[6]->Clone();
  h[6]->Rebin(4);
  h[6]->Divide(h1[6]);
  h[6]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[7]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();
  dir[7]->GetObject("qHist1D_sig_7_2",qHist_1D1[7]);
  h1[7]=(TH1D*)qHist_1D1[7]->Clone();
  h[7]->Rebin(4);
  h[7]->Divide(h1[7]);
  h[7]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[8]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();
  dir[8]->GetObject("qHist1D_sig_8_2",qHist_1D1[8]);
   h1[8]=(TH1D*)qHist_1D1[8]->Clone();
  h[8]->Rebin(4);
  h[8]->Divide(h1[8]);
  h[8]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[9]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();
  dir[9]->GetObject("qHist1D_sig_9_2",qHist_1D1[9]);
   h1[9]=(TH1D*)qHist_1D1[9]->Clone();
  h[9]->Rebin(4);
  h[9]->Divide(h1[9]);
  h[9]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[10]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();
  dir[10]->GetObject("qHist1D_sig_10_2",qHist_1D1[10]);
  h1[10]=(TH1D*)qHist_1D1[10]->Clone();
  h[10]->Rebin(4);
  h[10]->Divide(h1[10]);
  h[10]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[11]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();
  dir[11]->GetObject("qHist1D_sig_11_2",qHist_1D1[11]);
  h1[11]=(TH1D*)qHist_1D1[11]->Clone();
  h[11]->Rebin(4);
  h[11]->Divide(h1[11]);
  h[11]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[12]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();
  dir[12]->GetObject("qHist1D_sig_12_2",qHist_1D1[12]);
  h1[12]=(TH1D*)qHist_1D1[12]->Clone();
  h[12]->Rebin(4);
  h[12]->Divide(h1[12]);
  h[12]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[13]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();
  dir[13]->GetObject("qHist1D_sig_13_2",qHist_1D1[13]);
  h1[13]=(TH1D*)qHist_1D1[13]->Clone();
  h[13]->Rebin(4);
  h[13]->Divide(h1[13]);
  h[13]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[14]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();
  dir[14]->GetObject("qHist1D_sig_14_2",qHist_1D1[14]);
  h1[14]=(TH1D*)qHist_1D1[14]->Clone();
  h[14]->Rebin(4);
  h[14]->Divide(h1[14]);
  h[14]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[15]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();
  dir[15]->GetObject("qHist1D_sig_15_2",qHist_1D1[15]);
  h1[15]=(TH1D*)qHist_1D1[15]->Clone();
  h[15]->Rebin(4);
  h[15]->Divide(h1[15]);
  h[15]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[16]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();
  dir[16]->GetObject("qHist1D_sig_16_2",qHist_1D1[16]);
  h1[16]=(TH1D*)qHist_1D1[16]->Clone();
  h[16]->Rebin(4);
  h[16]->Divide(h1[16]);
  h[16]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[17]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();
  dir[17]->GetObject("qHist1D_sig_17_2",qHist_1D1[17]);
  h1[17]=(TH1D*)qHist_1D1[17]->Clone();
  h[17]->Rebin(4);
  h[17]->Divide(h1[17]);
  h[17]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[18]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();
  dir[18]->GetObject("qHist1D_sig_18_2",qHist_1D1[18]);
  h1[18]=(TH1D*)qHist_1D1[18]->Clone();
  h[18]->Rebin(4);
  h[18]->Divide(h1[18]);
  h[18]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[19]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();
  dir[19]->GetObject("qHist1D_sig_19_2",qHist_1D1[19]);
  h1[19]=(TH1D*)qHist_1D1[19]->Clone();
  h[19]->Rebin(4);
  h[19]->Divide(h1[19]);
  h[19]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[20]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();
  dir[20]->GetObject("qHist1D_sig_20_2",qHist_1D1[20]);
  h1[20]=(TH1D*)qHist_1D1[20]->Clone();
  h[20]->Rebin(4);
  h[20]->Divide(h1[20]);
  h[20]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;
  
   _file[21]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();
  dir[21]->GetObject("qHist1D_sig_21_2",qHist_1D1[21]);
  h1[21]=(TH1D*)qHist_1D1[21]->Clone();
  h[21]->Rebin(4);
  h[21]->Divide(h1[21]);
  h[21]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[22]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();
  dir[22]->GetObject("qHist1D_sig_22_2",qHist_1D1[22]);
  h1[22]=(TH1D*)qHist_1D1[22]->Clone();
  h[22]->Rebin(4);
  h[22]->Divide(h1[22]);
  h[22]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[23]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();
  dir[23]->GetObject("qHist1D_sig_23_2",qHist_1D1[23]);
  h1[23]=(TH1D*)qHist_1D1[23]->Clone();
  h[23]->Rebin(4);
  h[23]->Divide(h1[23]);
  h[23]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  m=m+1;

   _file[24]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_sig_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();
  dir[24]->GetObject("qHist1D_sig_24_2",qHist_1D1[24]);
  h1[24]=(TH1D*)qHist_1D1[24]->Clone();
  h[24]->Rebin(4);
  h[24]->Divide(h1[24]);
  h[24]->Fit("f1","","",140000, 350000);
  n[m]=m;
  intercept[m]=f1->GetParameter(0);
  dintercept[m]=f1->GetParError(0);
  slope[m]=f1->GetParameter(1);
  //  m=m+1;

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  /*  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=24; i++)
    {
      h_sum->Add(h[i],1); 
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;
  */

   

   
      
        c4=new TCanvas("c4","chi-squared vs calo#");
	//	c4->Divide(1,2);
	//	c4->cd(1);
     gr_int=new TGraphErrors(m,n,intercept,0,dintercept);
    gr_int->SetTitle("intercept of rebin ratio fit vs calo#");
    gr_int->GetXaxis()->SetTitle("calo #");
    gr_int->GetYaxis()->SetTitle("intercept");
     gr_int->SetMarkerStyle(20);
    gr_int->Draw();
    /* c4->cd(2);
      gr6=new TGraphErrors(m,n,slope,0,0);
    gr6->SetTitle("slope of rebin ratio fit vs calo#");
    gr6->GetXaxis()->SetTitle("calo #");
    gr6->GetYaxis()->SetTitle("slope");
     gr6->SetMarkerStyle(20);
    gr6->Draw();
    */	
     	            
     
      
 }
