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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo;
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

Double_t wiggle(Double_t *x, Double_t *par)
{
  Double_t f_x=par[0]*exp(-x[0]/par[1])*(1+par[2]*cos(par[3]*x[0]+par[4]));
  return f_x;
}

/*Double_t cbo(Double_t *x, Double_t *par)
{
  Double_t f_x=(1-par[0]*exp(-x[0]/par[1])*(cos(par[2]*x[0]+par[3])));
  return f_x;
  }*/

 /*double paramToFreq(double blindedValue)
{

  double unblindedR = blindedValue - deltaR;
  return 2 * TMath::Pi() * refFreq * (1 + (unblindedR * precisionR));
}
 */

Double_t fprec(Double_t *x, Double_t *par)
{
    double norm0 = par[0];

    double A1 = par[1];

    double tau1 = par[2];

    double A2 = par[3];

    double tau2 = par[4];

    double omega_r = par[5];

    double phi = par[6];

    double time = x[0];


    return norm0*(1+A1*exp(-time/tau1))*(1+A2*exp(-time/tau2)*cos(omega_r*time + phi));
    }

void ped_ring_raw_fit()

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
  dir[1]->GetObject("qHist1D_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzerDB",dir[2]);
  dir[2]->GetObject("qHist1D_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzerDB",dir[3]);
  dir[3]->GetObject("qHist1D_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzerDB",dir[4]);
  dir[4]->GetObject("qHist1D_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzerDB",dir[5]);
  dir[5]->GetObject("qHist1D_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

    _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzerDB",dir[6]);
  dir[6]->GetObject("qHist1D_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzerDB",dir[7]);
  dir[7]->GetObject("qHist1D_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzerDB",dir[8]);
  dir[8]->GetObject("qHist1D_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzerDB",dir[9]);
  dir[9]->GetObject("qHist1D_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzerDB",dir[10]);
  dir[10]->GetObject("qHist1D_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

  _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzerDB",dir[11]);
  dir[11]->GetObject("qHist1D_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzerDB",dir[12]);
  dir[12]->GetObject("qHist1D_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzerDB",dir[13]);
  dir[13]->GetObject("qHist1D_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzerDB",dir[14]);
  dir[14]->GetObject("qHist1D_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzerDB",dir[15]);
  dir[15]->GetObject("qHist1D_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzerDB",dir[16]);
  dir[16]->GetObject("qHist1D_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzerDB",dir[17]);
  dir[17]->GetObject("qHist1D_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzerDB",dir[18]);
  dir[18]->GetObject("qHist1D_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzerDB",dir[19]);
  dir[19]->GetObject("qHist1D_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzerDB",dir[20]);
  dir[20]->GetObject("qHist1D_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzerDB",dir[21]);
  dir[21]->GetObject("qHist1D_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzerDB",dir[22]);
  dir[22]->GetObject("qHist1D_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzerDB",dir[23]);
  dir[23]->GetObject("qHist1D_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
  _file[24]->GetObject("QFillByFillAnalyzerDB",dir[24]);
  dir[24]->GetObject("qHist1D_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
     sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
  
  //return;
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  // h_sum->Sumw2(kTRUE);
   
   for(Int_t i=1; i<=24; i++){h_sum->Add(h[i],1);}

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
   

   
   

 c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);
 
 



  fit_func= new TF1("fprec", fprec,  3000,  45000, 7);
  fit_func->SetParNames("norm0", "A1", "tau1", "A2", "tau2", "omega_r", "phi");
 

   
 fit_func->SetNpx(1000000);



 // gStyle->SetOptFit(1111);
 //  c2->Divide(1,3);
 //  c2->cd(1);
   h_sum->Rebin(2);
   h_sum->Scale(0.5);




      fit_start=3500;
      fit_stop=40000;

      //  h_sum->Rebin(8);
      fit_func->SetParameters(-32200000000,0,9500,0.035,56000,0.01067,0);

      // fit_func->FixParameter(5,0.01061);
   
        h_sum->GetXaxis()->SetRangeUser(fit_start, fit_stop);
	h_sum->GetYaxis()->SetTitle("ADC counts");
	h_sum->GetXaxis()->SetTitle("time [ns]");
	//fit_func->SetParLimits(8, -3.14, 3.14);
	 // fit_func->SetParLimits(10, -3.14, 3.14);
	//		h_sum->Draw();
	// fit_func->Draw("same");
		h_sum->Fit("fprec","R","",fit_start,fit_stop);
     
		h_res= new TH1D("residual histogram", "h_res", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");
	h_res->GetXaxis()->SetRangeUser(3000, 45000);
        h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm->SetBins(h_res->GetNbinsX(), 0, 1000/h_res->GetBinWidth(1));
   //hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm->GetXaxis()->SetRangeUser(0,1000/(2*h_res->GetBinWidth(1)));
   hm->GetXaxis()->SetTitle("Freq [MHz]");
   hm->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm->Draw();     


     
    	
       
 }
