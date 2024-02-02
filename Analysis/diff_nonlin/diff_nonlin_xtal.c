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
TH1D *qHist_1D[55], *h[55], *h_sum, *hcalo, *h_odd[55], *h_even[55], *h_diff[55], *h_average[55], *hstat;
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

//char root_file_name[128] = "run2d_dqc_26013_100sr_calo6_xtals.root";
//char root_file_name[128] = "run2d_dqc_26013_100sr_calo6_xtals.root";
//char root_file_name[128] = "run2d_dqc_26013_100sr_calo12_xtals.root";
//char root_file_name[128] = "run2d_dqc_26013_100sr_calo18_xtals.root";
char root_file_name[128] = "run2d_dqc_26013_100sr_calo23_xtals.root";


void diff_nonlin_xtal()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40], p0[56], dp0[40], average_trace[40], diff_trace[40];
   Double_t n[56];
   Int_t m=0;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[56];
 TDirectoryFile *dir[56];

  _file[0]=TFile::Open(root_file_name);
  _file[0]->GetObject("QFillByFillAnalyzer",dir[0]);
  dir[0]->GetObject("qHist1D_xtal_0",qHist_1D[0]);
  h[0]=(TH1D*)qHist_1D[0]->Clone();


  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_xtal_1",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_xtal_2",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_xtal_3",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_xtal_4",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_xtal_5",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_xtal_6",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_xtal_7",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_xtal_8",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_xtal_9",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_xtal_10",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_xtal_11",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_xtal_12",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_xtal_13",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_xtal_14",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_xtal_15",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_xtal_16",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_xtal_17",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_xtal_18",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_xtal_19",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_xtal_20",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_xtal_21",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_xtal_22",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_xtal_23",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_xtal_24",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();

    _file[25]=TFile::Open(root_file_name);
  _file[25]->GetObject("QFillByFillAnalyzer",dir[25]);
  dir[25]->GetObject("qHist1D_xtal_25",qHist_1D[25]);
  h[25]=(TH1D*)qHist_1D[25]->Clone();

     _file[26]=TFile::Open(root_file_name);
  _file[26]->GetObject("QFillByFillAnalyzer",dir[26]);
  dir[26]->GetObject("qHist1D_xtal_26",qHist_1D[26]);
  h[26]=(TH1D*)qHist_1D[26]->Clone();

     _file[27]=TFile::Open(root_file_name);
  _file[27]->GetObject("QFillByFillAnalyzer",dir[27]);
  dir[27]->GetObject("qHist1D_xtal_27",qHist_1D[27]);
  h[27]=(TH1D*)qHist_1D[27]->Clone();

     _file[28]=TFile::Open(root_file_name);
  _file[28]->GetObject("QFillByFillAnalyzer",dir[28]);
  dir[28]->GetObject("qHist1D_xtal_28",qHist_1D[28]);
  h[28]=(TH1D*)qHist_1D[28]->Clone();

     _file[29]=TFile::Open(root_file_name);
  _file[29]->GetObject("QFillByFillAnalyzer",dir[29]);
  dir[29]->GetObject("qHist1D_xtal_29",qHist_1D[29]);
  h[29]=(TH1D*)qHist_1D[29]->Clone();

     _file[30]=TFile::Open(root_file_name);
  _file[30]->GetObject("QFillByFillAnalyzer",dir[30]);
  dir[30]->GetObject("qHist1D_xtal_30",qHist_1D[30]);
  h[30]=(TH1D*)qHist_1D[30]->Clone();

     _file[31]=TFile::Open(root_file_name);
  _file[31]->GetObject("QFillByFillAnalyzer",dir[31]);
  dir[31]->GetObject("qHist1D_xtal_31",qHist_1D[31]);
  h[31]=(TH1D*)qHist_1D[31]->Clone();

     _file[32]=TFile::Open(root_file_name);
  _file[32]->GetObject("QFillByFillAnalyzer",dir[32]);
  dir[32]->GetObject("qHist1D_xtal_32",qHist_1D[32]);
  h[32]=(TH1D*)qHist_1D[32]->Clone();

     _file[33]=TFile::Open(root_file_name);
  _file[33]->GetObject("QFillByFillAnalyzer",dir[33]);
  dir[33]->GetObject("qHist1D_xtal_33",qHist_1D[33]);
  h[33]=(TH1D*)qHist_1D[33]->Clone();

     _file[34]=TFile::Open(root_file_name);
  _file[34]->GetObject("QFillByFillAnalyzer",dir[34]);
  dir[34]->GetObject("qHist1D_xtal_34",qHist_1D[34]);
  h[34]=(TH1D*)qHist_1D[34]->Clone();

     _file[35]=TFile::Open(root_file_name);
  _file[35]->GetObject("QFillByFillAnalyzer",dir[35]);
  dir[35]->GetObject("qHist1D_xtal_35",qHist_1D[35]);
  h[35]=(TH1D*)qHist_1D[35]->Clone();

     _file[36]=TFile::Open(root_file_name);
  _file[36]->GetObject("QFillByFillAnalyzer",dir[36]);
  dir[36]->GetObject("qHist1D_xtal_36",qHist_1D[36]);
  h[36]=(TH1D*)qHist_1D[36]->Clone();

     _file[37]=TFile::Open(root_file_name);
  _file[37]->GetObject("QFillByFillAnalyzer",dir[37]);
  dir[37]->GetObject("qHist1D_xtal_37",qHist_1D[37]);
  h[37]=(TH1D*)qHist_1D[37]->Clone();

     _file[38]=TFile::Open(root_file_name);
  _file[38]->GetObject("QFillByFillAnalyzer",dir[38]);
  dir[38]->GetObject("qHist1D_xtal_38",qHist_1D[38]);
  h[38]=(TH1D*)qHist_1D[38]->Clone();

      _file[39]=TFile::Open(root_file_name);
  _file[39]->GetObject("QFillByFillAnalyzer",dir[39]);
  dir[39]->GetObject("qHist1D_xtal_39",qHist_1D[39]);
  h[39]=(TH1D*)qHist_1D[39]->Clone();

     _file[40]=TFile::Open(root_file_name);
  _file[40]->GetObject("QFillByFillAnalyzer",dir[40]);
  dir[40]->GetObject("qHist1D_xtal_40",qHist_1D[40]);
  h[40]=(TH1D*)qHist_1D[40]->Clone();

     _file[41]=TFile::Open(root_file_name);
  _file[41]->GetObject("QFillByFillAnalyzer",dir[41]);
  dir[41]->GetObject("qHist1D_xtal_41",qHist_1D[41]);
  h[41]=(TH1D*)qHist_1D[41]->Clone();

     _file[42]=TFile::Open(root_file_name);
  _file[42]->GetObject("QFillByFillAnalyzer",dir[42]);
  dir[42]->GetObject("qHist1D_xtal_42",qHist_1D[42]);
  h[42]=(TH1D*)qHist_1D[42]->Clone();

     _file[43]=TFile::Open(root_file_name);
  _file[43]->GetObject("QFillByFillAnalyzer",dir[43]);
  dir[43]->GetObject("qHist1D_xtal_43",qHist_1D[43]);
  h[43]=(TH1D*)qHist_1D[43]->Clone();

     _file[44]=TFile::Open(root_file_name);
  _file[44]->GetObject("QFillByFillAnalyzer",dir[44]);
  dir[44]->GetObject("qHist1D_xtal_44",qHist_1D[44]);
  h[44]=(TH1D*)qHist_1D[44]->Clone();

     _file[45]=TFile::Open(root_file_name);
  _file[45]->GetObject("QFillByFillAnalyzer",dir[45]);
  dir[45]->GetObject("qHist1D_xtal_45",qHist_1D[45]);
  h[45]=(TH1D*)qHist_1D[45]->Clone();

     _file[46]=TFile::Open(root_file_name);
  _file[46]->GetObject("QFillByFillAnalyzer",dir[46]);
  dir[46]->GetObject("qHist1D_xtal_46",qHist_1D[46]);
  h[46]=(TH1D*)qHist_1D[46]->Clone();

     _file[47]=TFile::Open(root_file_name);
  _file[47]->GetObject("QFillByFillAnalyzer",dir[47]);
  dir[47]->GetObject("qHist1D_xtal_47",qHist_1D[47]);
  h[47]=(TH1D*)qHist_1D[47]->Clone();

     _file[48]=TFile::Open(root_file_name);
  _file[48]->GetObject("QFillByFillAnalyzer",dir[48]);
  dir[48]->GetObject("qHist1D_xtal_48",qHist_1D[48]);
  h[48]=(TH1D*)qHist_1D[48]->Clone();

     _file[49]=TFile::Open(root_file_name);
  _file[49]->GetObject("QFillByFillAnalyzer",dir[49]);
  dir[49]->GetObject("qHist1D_xtal_49",qHist_1D[49]);
  h[49]=(TH1D*)qHist_1D[49]->Clone();

     _file[50]=TFile::Open(root_file_name);
  _file[50]->GetObject("QFillByFillAnalyzer",dir[50]);
  dir[50]->GetObject("qHist1D_xtal_50",qHist_1D[50]);
  h[50]=(TH1D*)qHist_1D[50]->Clone();

     _file[51]=TFile::Open(root_file_name);
  _file[51]->GetObject("QFillByFillAnalyzer",dir[51]);
  dir[51]->GetObject("qHist1D_xtal_51",qHist_1D[51]);
  h[51]=(TH1D*)qHist_1D[51]->Clone();

     _file[52]=TFile::Open(root_file_name);
  _file[52]->GetObject("QFillByFillAnalyzer",dir[52]);
  dir[52]->GetObject("qHist1D_xtal_52",qHist_1D[52]);
  h[52]=(TH1D*)qHist_1D[52]->Clone();

     _file[53]=TFile::Open(root_file_name);
  _file[53]->GetObject("QFillByFillAnalyzer",dir[53]);
  dir[53]->GetObject("qHist1D_xtal_53",qHist_1D[53]);
  h[53]=(TH1D*)qHist_1D[53]->Clone();

    _file[54]=TFile::Open(root_file_name);
  _file[54]->GetObject("QFillByFillAnalyzer",dir[54]);
  dir[54]->GetObject("Hseqnumstats",qHist_1D[54]);
  hstat=(TH1D*)qHist_1D[54]->Clone();




  Double_t sum=0;
  
 for(Int_t k=0; k<=53; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
  c1->Divide(9,6);


   
    for(Int_t j=0; j<=53; j++)
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

   for(Int_t j=0; j<=53; j++)
      {	
     for(Int_t i=1; i<=h[j]->GetNbinsX()/2; i++)
       {
	 // diff_trace[j]=0;
	 //average_trace[j]=0;
	 h_diff[j]->SetBinContent(i,(h_odd[j]->GetBinContent(i)-h_even[j]->GetBinContent(i)));
	 h_average[j]->SetBinContent(i,(h_even[j]->GetBinContent(i)+h_odd[j]->GetBinContent(i))/2);
	 //diff_trace[j]=diff_trace[j]/average_trace[j];
	 //h_diff[j]->SetBinContent(i, diff_trace[j]);
       }
     //h_diff[j]->Divide(h_average[j]);
      c1->cd(j+1);
      //  h_diff[j]->Fit("pol1","R","",150000, 350000);
      n[m]=j;
      // TF1 *myfit = (TF1*) h_diff[j]->GetFunction("pol1");
      p0[m]=h_diff[j]->Integral(6667,8334)/(1664*hstat->Integral());
     //dp0[m]=myfit->GetParError(0);
      //h_diff[j]->GetYaxis()->SetRangeUser(-4, 4);
      h_diff[j]->GetYaxis()->SetRangeUser(-1500,1000);
     h_diff[j]->Draw();
     m=m+1;
      }

   c2=new TCanvas("c2","odd even difference parameter");   
    gr1=new TGraphErrors(m,n,p0,0,0);
    gr1->SetTitle("Amplitude vs xtal #; Calo 24");
    gr1->GetXaxis()->SetTitle("xtal #");
    gr1->GetYaxis()->SetTitle("Amplitude");
     gr1->SetMarkerStyle(20);
     gr1->SetLineColor(kRed);
     gr1->SetLineWidth(7);
     gr1->GetYaxis()->SetRangeUser(-2,2);
    gr1->Draw();
   
       
 }
