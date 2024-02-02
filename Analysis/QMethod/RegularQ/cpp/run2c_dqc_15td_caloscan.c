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
TH1D *qHist_1D[25], *qHist_2[25], *h[25], *hg[25], *h_sum, *h_sum2, *h_err, *h_err2, *hcalo;
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
    
    double time = x[0];
    
    double omega_a = 1.e-3 * getBlinded.paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)))  * (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo)));

    //    return norm * exp(-time/life) * (1 + asym* (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo)))*cos(omega_a*time + phi)) * (1-A1_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi1_cbo)));
    //* (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) *  (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));

    }

void run2c_dqc_15td_caloscan()

 {
   Double_t chi15[40], chi60[40], life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");

  _file[1]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
  dir[1]->GetObject("qHist1D_sig_1_2",qHist_2[1]);
  hg[1]=(TH1D*)qHist_2[1]->Clone();


  _file[2]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();
   dir[2]->GetObject("qHist1D_sig_2_2",qHist_2[2]);
  hg[2]=(TH1D*)qHist_2[2]->Clone();


   _file[3]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();
   dir[3]->GetObject("qHist1D_sig_3_2",qHist_2[3]);
  hg[3]=(TH1D*)qHist_2[3]->Clone();

   _file[4]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();
 dir[4]->GetObject("qHist1D_sig_4_2",qHist_2[4]);
  hg[4]=(TH1D*)qHist_2[4]->Clone();

   _file[5]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();
 dir[5]->GetObject("qHist1D_sig_5_2",qHist_2[5]);
  hg[5]=(TH1D*)qHist_2[5]->Clone();

   _file[6]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();
 dir[6]->GetObject("qHist1D_sig_6_2",qHist_2[6]);
  hg[6]=(TH1D*)qHist_2[6]->Clone();

   _file[7]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();
 dir[7]->GetObject("qHist1D_sig_7_2",qHist_2[7]);
  hg[7]=(TH1D*)qHist_2[7]->Clone();

   _file[8]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();
 dir[8]->GetObject("qHist1D_sig_8_2",qHist_2[8]);
  hg[8]=(TH1D*)qHist_2[8]->Clone();

   _file[9]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();
 dir[9]->GetObject("qHist1D_sig_9_2",qHist_2[9]);
  hg[9]=(TH1D*)qHist_2[9]->Clone();

   _file[10]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();
 dir[10]->GetObject("qHist1D_sig_10_2",qHist_2[10]);
  hg[10]=(TH1D*)qHist_2[10]->Clone();

   _file[11]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();
 dir[11]->GetObject("qHist1D_sig_11_2",qHist_2[11]);
  hg[11]=(TH1D*)qHist_2[11]->Clone();

   _file[12]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();
 dir[12]->GetObject("qHist1D_sig_12_2",qHist_2[12]);
  hg[12]=(TH1D*)qHist_2[12]->Clone();

   _file[13]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();
 dir[13]->GetObject("qHist1D_sig_13_2",qHist_2[13]);
  hg[13]=(TH1D*)qHist_2[13]->Clone();

   _file[14]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();
 dir[14]->GetObject("qHist1D_sig_14_2",qHist_2[14]);
  hg[14]=(TH1D*)qHist_2[14]->Clone();

   _file[15]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();
 dir[15]->GetObject("qHist1D_sig_15_2",qHist_2[15]);
  hg[15]=(TH1D*)qHist_2[15]->Clone();

   _file[16]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();
 dir[16]->GetObject("qHist1D_sig_16_2",qHist_2[16]);
  hg[16]=(TH1D*)qHist_2[16]->Clone();

   _file[17]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();
 dir[17]->GetObject("qHist1D_sig_17_2",qHist_2[17]);
  hg[17]=(TH1D*)qHist_2[17]->Clone();

   _file[18]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();
 dir[18]->GetObject("qHist1D_sig_18_2",qHist_2[18]);
  hg[18]=(TH1D*)qHist_2[18]->Clone();

   _file[19]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();
 dir[19]->GetObject("qHist1D_sig_19_2",qHist_2[19]);
  hg[19]=(TH1D*)qHist_2[19]->Clone();

   _file[20]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();
 dir[20]->GetObject("qHist1D_sig_20_2",qHist_2[20]);
  hg[20]=(TH1D*)qHist_2[20]->Clone();

   _file[21]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();
 dir[21]->GetObject("qHist1D_sig_21_2",qHist_2[21]);
  hg[21]=(TH1D*)qHist_2[21]->Clone();

   _file[22]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();
 dir[22]->GetObject("qHist1D_sig_22_2",qHist_2[22]);
  hg[22]=(TH1D*)qHist_2[22]->Clone();

   _file[23]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();
 dir[23]->GetObject("qHist1D_sig_23_2",qHist_2[23]);
  hg[23]=(TH1D*)qHist_2[23]->Clone();

   _file[24]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[24]->GetObject("QFillByFillAnalyzer",dir[24]);
  dir[24]->GetObject("qHist1D_sig_24_0",qHist_1D[24]);
  h[24]=(TH1D*)qHist_1D[24]->Clone();
  dir[24]->GetObject("qHist1D_sig_24_2",qHist_2[24]);
  hg[24]=(TH1D*)qHist_2[24]->Clone();

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);

  h_sum2= new TH1D("calo histogram sum2", "h_sum2", hg[1]->GetNbinsX(), 100001, 352001);
  h_sum2->Sumw2(kTRUE);
 
   for(Int_t i=1; i<=24; i++)
    {
      h_sum->Add(h[i],1);
      h_sum2->Add(hg[i],1);
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   //   h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

     double binwidth2=hg[1]->GetBinWidth(1);
   int nbins2=hg[1]->GetNbinsX();
   double binlowedge2=hg[1]->GetBinLowEdge(1);
   double binhighedge2=hg[1]->GetBinLowEdge(nbins2)+binwidth2;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   h_sum2->SetBins(nbins2, rawBinToNs*(binlowedge2-inject_time), rawBinToNs*(binhighedge2-inject_time));
   h_sum->Rebin(8);
   h_sum2->Rebin(2);

   
    h_err= new TH1D("calo sum histogram error", "h_err", h_sum->GetNbinsX(), 100001, 352001);   
    h_err2= new TH1D("calo sum histogram error2", "h_err2", h_sum2->GetNbinsX(), 100001, 352001);

    for(Int_t k=1;k<=h_sum->GetNbinsX();k++)
      {
	h_err->SetBinContent(k,h_sum->GetBinError(k));
      }

     for(Int_t k=1;k<=h_sum2->GetNbinsX();k++)
      {
	h_err2->SetBinContent(k,h_sum2->GetBinError(k));
      }

   
     /*   c1->Divide(1,2);
   c1->cd(1);
   h_sum2->Divide(h_sum);
   h_sum2->Draw();

   c1->cd(2);
     h_err2->Divide(h_err);
     h_err2->Draw();
     
     return;*/
   for(Int_t i=1; i<=24; i++)
     {
       h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       h[i]->Rebin(8);
       h[i]->GetYaxis()->SetRangeUser(-100., 30000000000./24);
       h[i]->GetYaxis()->SetTitle("ADC counts");
       h[i]->GetXaxis()->SetTitle("time [ns]");

       hg[i]->SetBins(nbins2, rawBinToNs*(binlowedge2-inject_time), rawBinToNs*(binhighedge2-inject_time));
       hg[i]->Rebin(2);
       hg[i]->GetYaxis()->SetRangeUser(-100., 30000000000./24);
       hg[i]->GetYaxis()->SetTitle("ADC counts");
       hg[i]->GetXaxis()->SetTitle("time [ns]");
       
     }
   //  h_sum->Rebin(8);
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   //   h[1]->Draw();

   hcalo = new TH1D("calosum histogram for fr", "hcalo", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));

   for(Int_t k=1; k<=h_sum->GetNbinsX(); k++)
     {
       hcalo->SetBinContent(k,h_sum->GetBinContent(k));
       hcalo->SetBinError(k,h_sum->GetBinError(k));
     }
   //   hcalo->Rebin(4);
   


 c2=new TCanvas("c2","5 parameter wiggle_fit");
  c2->Divide(1,3);
 c2->cd(1);
 
 



  fit_func= new TF1("fprec", fprec,  30000,  309000, 18);
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

   
 fit_func->SetNpx(1000000);



 gStyle->SetOptFit(1111);




      fit_start=30000;
      fit_stop=300000;

      cout<<"Calo "<<m<<endl;
   fit_func->SetParameters(1351557083, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 193436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     
     h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     //  cout<<" chi-squared 60 "<<chi60[m]<<endl;
        m=m+1;
     
	/*  		h_res= new TH1D("residual histogram", "h_res", h[1]->GetNbinsX(), h[1]->GetBinLowEdge(1), (h[1]->GetBinLowEdge(h[1]->GetNbinsX())+h[1]->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- h[1]->GetBinLowEdge(1))/h[1]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h[1]->GetBinLowEdge(1))/h[1]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[1]->GetBinContent(ibin)- fit_func->Eval( h[1]->GetBinCenter(ibin) ) );
      if(h[1]->GetBinError(ibin)!=0){res=(res/h[1]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");	
        h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   hm->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hm->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hm->GetXaxis()->SetTitle("Freq [GHz]");
   hm->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hm->Draw();     
     
   return;*/
    cout<<"Calo "<<m<<endl;
    fit_func->SetParameters(1351557083/0.95, 64439.4, 0.234041, 0, 2.22, 0.00261778, 208510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[2]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     

   
     cout<<"Calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.04, 64439.4, 0.234041, 0, 2.22, 0.00261778, 208510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 20000);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     // h[3]->Draw();
     //fit_func->Draw("same");
     h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     //  return;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
      fit_func->SetParameters(1351557083/1.12, 64439.4, 0.234041, 0, 2.22, 0.00261778, 208510, 0.00234020, 3.57360);
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 19000);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
     //hg[m]->Draw();
     //fit_func->Draw("same");
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     // return;


     cout<<"Calo "<<m<<endl;
          fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[4]->Draw();
     //fit_func->Draw("same");
       h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     fit_func->SetParameter(0,1351557083/1.1);
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //  return;

     cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 193436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[5]->Draw();
     //  fit_func->Draw("same");
         h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     fit_func->SetParameter(0,1351557083/1.05);
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     // return;

     cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.08, 64439.4, 0.234041, 0, 2.22, 0.00261778, 228510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 10000);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[6]->Draw();
     //fit_func->Draw("same");
     h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //  return;

     cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.01, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[7]->Draw();
     //  fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //return;
    
     cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[8]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //return;

     
     cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 208510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[9]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     // return;
     

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[10]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[11]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //return;

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[12]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[13]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

     cout<<"calo "<<m<<endl;
     fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[14]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[15]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[16]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     fit_func->SetParameter(0, 1351557083/1.03);
     fit_func->SetParameter(14,0.05);
     //hg[m]->Draw();
     //fit_func->Draw("same");
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //  return;

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[17]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[18]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[19]->Draw();
     //fit_func->Draw("same");
        h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 193436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[20]->Draw();
     //fit_func->Draw("same");
       h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     fit_func->SetParameter(0,1351557083/1.03);
     // fit_func->SetParameter(14,0.05);
     //  fit_func->SetParameter(15,1827);
     //  fit_func->SetParameter(16,0.1535);
     //  hg[m]->Draw();
     //fit_func->Draw("same");
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //  return;

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 218510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.00638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[21]->Draw();
     //fit_func->Draw("same");
     h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     fit_func->SetParameter(0, 1351557083/1.15);
     // hg[m]->Draw();
     //fit_func->Draw("same");
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     //   return;    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 228510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[22]->Draw();
     //  fit_func->Draw("same");
     h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
     // return;

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[23]->Draw();
     //fit_func->Draw("same");
       h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     m=m+1;
    

      cout<<"calo "<<m<<endl;
       fit_func->SetParameters(1351557083/1.02, 64439.4, 0.234041, 0, 2.22, 0.00261778, 238510, 0.00234020, 3.57360);

   
     fit_func->SetParameter(9,0.000485);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11, -0.000311470);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);

     fit_func->SetParameter(14, 0.0638136);
     fit_func->SetParameter(15, 183436);
     fit_func->SetParameter(16, 0.153506);
     fit_func->SetParameter(17, 0);
  
     //h[24]->Draw();
     //fit_func->Draw("same");
       h[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi15[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //n[m]=m;
     cout<<" chi-squared 15 "<<chi15[m]<<endl;
     hg[m]->Fit("fprec","R","",fit_start,fit_stop);
     chi60[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<" chi-squared 60 "<<chi60[m]<<endl;
     // m=m+1;
    
       c4=new TCanvas("c4","chi-squared vs calo#");
     
     gr1=new TGraphErrors(m+1,n,chi15,0,0);
    gr1->SetTitle("chi2/ndf vs calo#");
    gr1->GetXaxis()->SetTitle("calo #");
    gr1->GetYaxis()->SetTitle("chi2/ndf");
     gr1->SetMarkerStyle(20);
     gr1->SetLineColor(kBlue);
    gr1->Draw();

    //       c4=new TCanvas("c4","chi-squared vs calo#");
     
     gr2=new TGraphErrors(m+1,n,chi60,0,0);
     //gr2->SetTitle("chi2/ndf vs calo#");
    //  gr2->GetXaxis()->SetTitle("calo #");
    //gr2->GetYaxis()->SetTitle("chi2/ndf");
     gr2->SetMarkerStyle(24);
     gr2->SetLineColor(kRed);
    gr2->Draw("same");


   auto legend = new TLegend(0.1,0.7);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry("gr1","15 time decimation","lep");
   legend->AddEntry("gr2","60 time decimation","lep");
   //legend->SetTextAlign(11);
   legend->Draw();	
     
 }
