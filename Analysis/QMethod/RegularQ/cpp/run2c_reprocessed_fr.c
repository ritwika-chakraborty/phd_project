#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TCanvas.h"
#include <sys/time.h>
#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );

char root_file_name[128] = "run2C_thresh400_reprocessed.root";

float toddiff(struct timeval *tod2, struct timeval *tod1)
{
  float fdt, fmudt;
  long long t1, t2, mut1, mut2;
  long long dt, mudt;
  t1 = tod1->tv_sec;
  mut1 = tod1->tv_usec;
  t2 = tod2->tv_sec;
  mut2 = tod2->tv_usec;
  dt = (t2 - t1);
  mudt = (mut2 - mut1);
  fdt = (float)dt;
  fmudt = (float)mudt;
  return 1.0e6 * fdt + fmudt;
}


TVirtualFitter *gFitter;

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo, *hcomp, *h1[2];
//TH1 *hm;
TH1F *hlm, h0;
TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res;
TH1 *hm;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t fit_start, fit_stop;
Int_t countfcn=0;
Int_t i0;
TMatrixD cov(2099,2099);



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

    double A_pr = par[22];

    double Tau_pr = par[23];

    double omega_pr = par[24];

    double phi_pr = par[25];

    double A_vd = par[26];

    double tau_vd = par[27];

    double C_vd = par[28];
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    double Ncbox = (1 + A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi_cbo)));

    double Acbox =  (1 + A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi2_cbo)));

    double phicbox = A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time + phi3_cbo));

    double Nvoy1 =  (1 + A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time + phi_vbo)));

    double Nvoy2 = (1 + A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time + phi_vbo2)));

    double muloss =  ( 1 - lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    // double pedring =  (1 + A_pr * exp(-time/Tau_pr) * (cos(omega_pr*time+phi_pr)));

    //double vertdrift =  (1 + (A_vd * exp(-time/tau_vd)) + C_vd);
    
    //  return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)))  * (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) * (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2))) * (1-A_pr * exp(-time/Tau_pr) * (cos(omega_pr*time+phi_pr))) * (1+(A_vd * exp(-time/tau_vd)) + C_vd);

    return norm  * Ncbox * Nvoy1 * Nvoy2 * muloss * exp(-time/life) * (1 + asym * Acbox * cos(omega_a*time + (phi + phicbox)));

    }





void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{ 
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(hcomp->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
  
   for(Int_t i=hcomp->FindBin(fit_start); i<=hcomp->FindBin(fit_stop); i++)
    {
      for(Int_t j=hcomp->FindBin(fit_start); j<=hcomp->FindBin(fit_stop); j++)
	{
	  // if(countfcn<=1000)
	  // {
	   if(j<=i+5 && j>=i-5)
	   {
	     ch=ch+((hcomp->GetBinContent(i))-(fuser->Eval(hcomp->GetBinCenter(i))))*cov[i][j]*((hcomp->GetBinContent(j))-(fuser->Eval(hcomp->GetBinCenter(j))));
	    }
	   //  }
	   /*  else
	    {
	      if(j!=i)
		{ ch=ch+((hcomp->GetBinContent(i))-(fuser->Eval(hcomp->GetBinCenter(i))))*cov[i][j]*((hcomp->GetBinContent(j))-(fuser->Eval(hcomp->GetBinCenter(j))));}
		}*/
		}
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}





void run2c_reprocessed_fr()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");

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
   h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   h_sum->Draw();

   hcalo = new TH1D("calosum histogram for fr", "hcalo", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));

   for(Int_t k=1; k<=h_sum->GetNbinsX(); k++)
     {
       hcalo->SetBinContent(k,h_sum->GetBinContent(k));
       hcalo->SetBinError(k,h_sum->GetBinError(k));
     }
   hcalo->Rebin(4);
   
   int iend;
   int flg = 0;
   int mdim = (hcalo->GetNbinsX()/2)-1;
   if(mdim!=2099){cout<<"dimension mismatch"<<endl; return;}
   TArrayD data((mdim)*(mdim));

   c2=new TCanvas("c2","5 parameter wiggle_fit");
   c2->Divide(1,3);
   c2->cd(1);


   
   //##### building covariance matrix and its inverse

     for(Int_t k=hcalo->FindBin(10000); k<=hcalo->GetNbinsX(); k++)
    {
      if(hcalo->GetBinError(k)!=0)
	{ i0=k;
	  break;}
    }
  cout<<"i0 is"<<i0<<" "<<endl;
  
   for (Int_t i = 0; i < (mdim)*(mdim); i++)
      {
	const Int_t ir = (i)/(mdim);
	const Int_t ic = (i)%(mdim);
    data[i] = 0.0;
 
    
    if ( ir == ic && ic>=i0)
      {
	if((2*ir-1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = ((0.25)*(hcalo->GetBinError(2*ir-1))*(hcalo->GetBinError(2*ir-1))+(1)*(hcalo->GetBinError(2*ir))*(hcalo->GetBinError(2*ir))+(0.25)*(hcalo->GetBinError(2*ir+1))*(hcalo->GetBinError(2*ir+1)));
	 cout<<"diagonal "<<ir<<" "<<ic<<endl;
	 // cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<h[1]->GetBinError(2*ir+1)<<" "<<endl;
	  if(data[ir*mdim+ic]==0)
	    {
	      cout<<ir<<" "<<ic<<"breaking"<<endl;
      	      iend=ir;
	      break;
	    }
	  else
	    {
	      iend=mdim;
	    }
       }
      }   
    if ( ic == ir+1 && ic>=i0 && ir>=i0)
      {
        if((2*ir+1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = (0.25)*(hcalo->GetBinError(2*ir+1))*(hcalo->GetBinError(2*ir+1));
	 //	 	cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<((0.25)*(h[1]->GetBinError(2*ir+1))*(h[1]->GetBinError(2*ir+1)))<<" "<<endl;
	 cout<<"diagonal+1  "<<ir<<" "<<ic<<endl;
       }
      }
     if ( ic == ir-1 && ic>=i0 && ir>=i0)
      {
        if((2*ir-1)>=i0)
       {
	 data[(ir)*mdim+(ic)] = (0.25)*(hcalo->GetBinError(2*ir-1))*(hcalo->GetBinError(2*ir-1));
	 //	 	cout<<i<<" "<<(ir)*mdim+(ic)<<" "<<((0.25)*(h[1]->GetBinError(2*ir+1))*(h[1]->GetBinError(2*ir+1)))<<" "<<endl;
	 cout<<"diagonal+1  "<<ir<<" "<<ic<<endl;
       }
      }

    flg=i;
     }
  cout<<"flg "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());

  cov.ResizeTo(i0, iend-1, i0, iend-1, -1);
  
  //   cov.Print();

   Double_t det1;
   cov.Invert(&det1);
   //   cov.Print();
   cout<<"Matrix inverted "<<i0<<endl;

   //######### end building covarince matrix and its inverse 


   h1[1]= new TH1D("180 degree shifted histogram", "h1[1]", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   h1[0]= new TH1D("-180 degree shifted histogram", "h1[0]", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   hcomp= new TH1D("superposed histogram", "h_calo1", hcalo->GetNbinsX(),  hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));
   for(Int_t k=1; k<=h1[1]->GetNbinsX(); k++)
     {
       h1[1]->SetBinContent(k,hcalo->GetBinContent(k+1));
       h1[0]->SetBinContent(k,hcalo->GetBinContent(k-1));
     }
   for(Int_t l=1; l<=hcomp->GetNbinsX(); l++)
     {
       hcomp->SetBinContent(l,(hcalo->GetBinContent(l)+h1[1]->GetBinContent(l))/2);
     }
   h1[1]->SetLineColor(kRed);
   hcomp->SetLineColor(kBlack);

   //  hcalo->GetYaxis()->SetRangeUser(-100., 60000000.);//rebin 0
   hcomp->Rebin(2);
   // hcomp->Fit("fprec","R","",fit_start,fit_stop);
       for(Int_t p=2; p<=hcomp->GetNbinsX(); p++)
     {
       hcomp->SetBinError(p,sqrt((1/4)*(hcalo->GetBinError(2*p-1))*(hcalo->GetBinError(2*p-1))+(1)*(hcalo->GetBinError(2*p))*(hcalo->GetBinError(2*p))+(1/4)*(hcalo->GetBinError(2*p+1))*(hcalo->GetBinError(2*p+1))));

       }

        hcomp->GetYaxis()->SetRangeUser(-100.,  14000000000);
	hcomp->GetYaxis()->SetTitle("ADC counts");
	hcomp->GetXaxis()->SetTitle("time [ns]");

   //cout<<(1/4)*(h[1]->GetBinError(2*p-1))*(h[1]->GetBinError(2*p-1))+(1)*(h[1]->GetBinError(2*p))*(h[1]->GetBinError(2*p))+(1/4)*(h[1]->GetBinError(2*p+1))*(h[1]->GetBinError(2*p+1))<<" "<<endl;



	
	

    fit_func= new TF1("fprec", fprec,  30000,  309000, 22);
    fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo", "A2_cbo","phi2_cbo");
  //, "A_vbo", "Tau_vbo", "omega_vbo", "phi_vbo");
  // );
     fit_func->SetParName(11, "A3_cbo");
     fit_func->SetParName(12, "phi3_cbo");
     fit_func->SetParName(13, "lm_amp");
     fit_func->SetParName(14, "A_vbo");
     fit_func->SetParName(15, "Tau_vbo");
     fit_func->SetParName(16, "omega_vbo");
     fit_func->SetParName(17, "phi_vbo");
     fit_func->SetParName(18, "A_vbo2");
     fit_func->SetParName(19, "Tau_vbo2");
     fit_func->SetParName(20, "omega_vbo2");
     fit_func->SetParName(21, "phi_vbo2");
     fit_func->SetParName(22, "A_pr");
     fit_func->SetParName(23, "Tau_pr");
     fit_func->SetParName(24, "omega_pr");
     fit_func->SetParName(25, "phi_pr");
     fit_func->SetParName(26, "A_vd");
     fit_func->SetParName(27, "Tau_vd");
     fit_func->SetParName(28, "constant_vd");


   
      fit_func->SetNpx(1000000);



 // gStyle->SetOptFit(1111);




      fit_start=30000;
      fit_stop=300000;

   h_sum->Rebin(8);
   fit_func->SetParameters(11601300000, 64426.2, 0.2378, 0, 2.22, -0.002412, 291700, 0.002339, 3.1);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func->SetParameter(9,-0.00062);
     fit_func->SetParameter(10,-2.26);
     fit_func->SetParameter(11, -0.000525);
     fit_func->SetParameter(12,1.52);

     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, -0.003705);
     // fit_func->SetParameter(15, 19012);
     fit_func->SetParameter(15,25283);
     fit_func->SetParameter(16, 0.153509);
     fit_func->SetParameter(17, -2.787);

     fit_func->SetParameter(18, 0.0002308);
     // fit_func->SetParameter(19, 6400);
     fit_func->SetParameter(19,252747);
     fit_func->SetParameter(20, 0.139593);
     fit_func->SetParameter(21, -0.1612);

     fit_func->SetParameter(22, -0.00123);
     fit_func->SetParameter(23,14760.9);
     fit_func->SetParameter(24, 0.01084);
     fit_func->SetParameter(25, -10.4);

     fit_func->SetParameter(26,-0.00842);
     fit_func->SetParameter(27, 248023);
     fit_func->SetParameter(28, 0.00569);

     // fit_func->FixParameter(7,0.00234);    
     //  fit_func->FixParameter(24, 0.010647);
     // fit_func->FixParameter(16,0.1535);
     //   fit_func->FixParameter(20,0.1392);
     //  fit_func->FixParameter(20,0);
     //fit_func->FixParameter(18,0);
     //fit_func->FixParameter(19,0);
     //fit_func->FixParameter(21,0);
     // fit_func->FixParameter(5,0);
     //fit_func->FixParameter(9,0);
     //fit_func->FixParameter(13,0);
     //fit_func->FixParameter(14,0);
     //  fit_func->FixParameter(16,0.1535);
     //fit_func->FixParameter(22,0);
     //fit_func->FixParameter(26,0);
     //fit_func->FixParameter(28,0);

        h_sum->GetYaxis()->SetRangeUser(-100., 30000000000.);
	h_sum->GetYaxis()->SetTitle("ADC counts");
	h_sum->GetXaxis()->SetTitle("time [ns]");
	//fit_func->SetParLimits(8, -3.14, 3.14);
	 // fit_func->SetParLimits(10, -3.14, 3.14);
	//h_sum->Draw();
	  // fit_func->Draw("same");
	//	h_sum->Fit("fprec","R","",fit_start,fit_stop);
         
	   printf("start Q-method analyzer\n");
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
     gFitter = TVirtualFitter::Fitter(hcomp);
     gFitter->SetFCN(chi2);  
     hcomp->Fit("fprec","RE","",fit_start,fit_stop);
     // hcomp->Draw();
     //fit_func->Draw("same");
      gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
  
		h_res= new TH1D("residual histogram", "h_res", hcomp->GetNbinsX(), hcomp->GetBinLowEdge(1), (hcomp->GetBinLowEdge(hcomp->GetNbinsX())+hcomp->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- hcomp->GetBinLowEdge(1))/hcomp->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcomp->GetBinLowEdge(1))/hcomp->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (hcomp->GetBinContent(ibin)- fit_func->Eval( hcomp->GetBinCenter(ibin) ) );
      if(hcomp->GetBinError(ibin)!=0){res=(res/hcomp->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");	
        h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm->SetBins(h_res->GetNbinsX(), 0, 1/h_res->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm->GetXaxis()->SetRangeUser(0,1/(2*h_res->GetBinWidth(1)));
   hm->GetXaxis()->SetTitle("Freq [GHz]");
   hm->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm->Draw();     


     
    	
       
 }