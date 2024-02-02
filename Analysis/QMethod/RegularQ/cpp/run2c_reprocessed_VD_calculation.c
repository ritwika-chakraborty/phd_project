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

char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *qHist_hslice_1D[25], *h[25], *hs[25], *h_sum, *hcalo;
//TH1 *hm;

TH1F *hlm, h0;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func, *VD_func;
TH1D *h_res[7],*hdiff;
TH1 *hm;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;

Double_t fVD(Double_t *x, Double_t *par)
{
  double A = par[0];

  double tau_A = par[1];

  //  double C = par[2];

  // double R = par[2];

  //double phi= par[3];

  //double asym= par[4];

  //  double omega_a = rawBinToNs * 1.e-3 * getBlinded.paramToFreq(R);

  double time = x[0];

  //return (1+A * exp(-time/tau_A) * (1+asym*cos(omega_a*time+phi))) ;

  return 1+(A * exp(-time/tau_A));

  //return (1+asym*cos(omega_a*time+phi));
}



Double_t fprec(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];

    double A_cbo = par[5];

    double Tau_cbo = par[6];

    double omega_cbo = par[7];

    double phi_cbo = par[8];

     double A2_cbo = par[9];

    double phi2_cbo = par[10];

    double A3_cbo = par[11];

    double phi3_cbo = par[12];

    double lost_muon_amp = par[13];

    double A_vbo = par[14];

    double Tau_vbo = par[15];

    double omega_vbo = par[16];

    double phi_vbo = par[17];

    double A_vbo2 = par[18];

    double Tau_vbo2 = par[19];

    double omega_vbo2 = par[20];

    double phi_vbo2 = par[21];
    
    
    double time = x[0];
    
    double omega_a = rawBinToNs * 1.e-3 * getBlinded.paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3*1.25)))  * (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) * (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));


   
    }

void run2c_reprocessed_VD_calculation()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40],norm1[40],dnorm1[40],norm2[40],dnorm2[40],norm3[40],dnorm3[40],norm4[40],dnorm4[40],norm5[40],dnorm5[40], norm6[40], norm7[40], norm8[40], norm9[40], norm10[40], dnorm6[40], dnorm7[40], dnorm8[40], dnorm9[40], dnorm10[40], diff[40],ratio[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

   Double_t A[10], tau_A[10];

 
 TFile *_fileh[7];
 TDirectoryFile *dirh[7];

 TFile *_file[25];
 TDirectoryFile *dir[25];

  _fileh[0]=TFile::Open("r2c.root");
  _fileh[0]->GetObject("hI",hlm);
    hlm->SetName("hlm");
 
  _fileh[1]=TFile::Open(root_file_name);
  _fileh[1]->GetObject("QFillByFillAnalyzerDB",dirh[1]);
  dirh[1]->GetObject("qHist1D_sig_hslice_0",qHist_hslice_1D[1]);
  hs[1]=(TH1D*)qHist_hslice_1D[1]->Clone();

  _fileh[2]=TFile::Open(root_file_name);
  _fileh[2]->GetObject("QFillByFillAnalyzerDB",dirh[2]);
  dirh[2]->GetObject("qHist1D_sig_hslice_1",qHist_hslice_1D[2]);
  hs[2]=(TH1D*)qHist_hslice_1D[2]->Clone();

   _fileh[3]=TFile::Open(root_file_name);
  _fileh[3]->GetObject("QFillByFillAnalyzerDB",dirh[3]);
  dirh[3]->GetObject("qHist1D_sig_hslice_2",qHist_hslice_1D[3]);
  hs[3]=(TH1D*)qHist_hslice_1D[3]->Clone();

   _fileh[4]=TFile::Open(root_file_name);
  _fileh[4]->GetObject("QFillByFillAnalyzerDB",dirh[4]);
  dirh[4]->GetObject("qHist1D_sig_hslice_3",qHist_hslice_1D[4]);
  hs[4]=(TH1D*)qHist_hslice_1D[4]->Clone();

   _fileh[5]=TFile::Open(root_file_name);
  _fileh[5]->GetObject("QFillByFillAnalyzerDB",dirh[5]);
  dirh[5]->GetObject("qHist1D_sig_hslice_4",qHist_hslice_1D[5]);
  hs[5]=(TH1D*)qHist_hslice_1D[5]->Clone();

   _fileh[6]=TFile::Open(root_file_name);
  _fileh[6]->GetObject("QFillByFillAnalyzerDB",dirh[6]);
  dirh[6]->GetObject("qHist1D_sig_hslice_5",qHist_hslice_1D[6]);
  hs[6]=(TH1D*)qHist_hslice_1D[6]->Clone();

    _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzerDB",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzerDB",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzerDB",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzerDB",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzerDB",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzerDB",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzerDB",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzerDB",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzerDB",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzerDB",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzerDB",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzerDB",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzerDB",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzerDB",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzerDB",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzerDB",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzerDB",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzerDB",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzerDB",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzerDB",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzerDB",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzerDB",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzerDB",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
  _file[24]->GetObject("QFillByFillAnalyzerDB",dir[24]);
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
   //  h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;

    double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   
      for(Int_t i=1; i<=6; i++)
     {
      hs[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }   

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   
   h_sum->Rebin(8);
   h_sum->Draw();
   h_sum->Rebin(30);
   h_sum->Scale(0.03333);
   
   c2= new TCanvas("c2","hslice/hsum ratios"); 
   c2->Divide(1,6);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0000);

    VD_func= new TF1("fVD", fVD,  1000,  350000,2);
    VD_func->SetParNames("A", "tau_A", "constant_bg");
	VD_func->SetNpx(1000000);
	VD_func->SetParameters(0.5, 500000);
   fit_start=25000;
     fit_stop=300000;
   //   fit_stop=200000;
   
   for(Int_t i=1; i<=6; i++)
     {
       hs[i]->Rebin(8);
       hs[i]->Rebin(30);
       hs[i]->Scale(0.03333);
       
       
       hs[i]->Scale(h_sum->Integral(67,16800)/hs[i]->Integral(67,16800));
       hs[i]->Divide(h_sum);
       c2->cd(i);
      
       //  hs[i]->Fit("fVD","R","",fit_start,fit_stop);
        hs[i]->Draw("hist");

     }
   
   
           c3= new TCanvas("c3","hslice/hsum ratios fits");
	   c3->Divide(1,6);
	   c3->cd(1);
	   hs[1]->Fit("fVD","R","",fit_start,fit_stop);
	   hs[1]->GetYaxis()->SetRangeUser(0.8,1.2);
	   hs[1]->Draw();
	   A[1]=VD_func->GetParameter(0);
	   tau_A[1]=VD_func->GetParameter(1);
	   //  VD_func->Draw("same");

	   
       	VD_func->SetParameters(0.5, 50000);
     
   
   c3->cd(2);
     hs[2]->Fit("fVD","R","",fit_start,fit_stop);
     hs[2]->GetYaxis()->SetRangeUser(0.8,1.2);
     hs[2]->Draw();
     A[2]=VD_func->GetParameter(0);
     tau_A[2]=VD_func->GetParameter(1);

	   // VD_func->Draw("same");
     

               	VD_func->SetParameters(0.1, 50000);
	       
   
   c3->cd(3);
   	   hs[3]->Fit("fVD","R","",fit_start,fit_stop);
   	   hs[3]->GetYaxis()->SetRangeUser(0.8,1.2);
           hs[3]->Draw();
	   A[3]=VD_func->GetParameter(0);
	   tau_A[3]=VD_func->GetParameter(1);

	   // VD_func->Draw("same");
  
   	        	VD_func->SetParameters(-0.1, 50000);
              
		
   c3->cd(4);
   	   hs[4]->Fit("fVD","R","",fit_start,fit_stop);
   	   hs[4]->GetYaxis()->SetRangeUser(0.8,1.2);
	   hs[4]->Draw();
	   A[4]=VD_func->GetParameter(0);
	   tau_A[4]=VD_func->GetParameter(1);

	   // VD_func->Draw("same");

       	VD_func->SetParameters(-0.5, 50000);
		        
   
   c3->cd(5);
   	   hs[5]->Fit("fVD","R","",fit_start,fit_stop);
   	   hs[5]->GetYaxis()->SetRangeUser(0.8,1.2);
	   hs[5]->Draw();
	   A[5]=VD_func->GetParameter(0);
	   tau_A[5]=VD_func->GetParameter(1);

	   // VD_func->Draw("same");




               	VD_func->SetParameters(-0.5, 50000);

   
   c3->cd(6);
    hs[6]->Fit("fVD","R","",fit_start,fit_stop);
    hs[6]->GetYaxis()->SetRangeUser(0.8,1.2);
	   hs[6]->Draw();
	   A[6]=VD_func->GetParameter(0);
	   tau_A[6]=VD_func->GetParameter(1);

	   // VD_func->Draw("same");
 
	   double amp_tau_sum=0;
	   double amp_sum=0;
	   for(Int_t k=1;k<=6; k++)
	     {
	       amp_tau_sum=amp_tau_sum+(A[k]*tau_A[k]);
	       amp_sum=amp_sum+A[k];
	     }

	   cout<<"weighted mean lifetime is "<<amp_tau_sum/amp_sum<<endl;
 }
