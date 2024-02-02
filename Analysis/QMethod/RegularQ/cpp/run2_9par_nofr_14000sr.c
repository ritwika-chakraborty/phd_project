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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo;
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

    //  double A2_cbo = par[9];

    //double phi2_cbo = par[10];

    /*  double A_vbo = par[11];

    double Tau_vbo = par[12];

    double omega_vbo = par[13];

    double phi_vbo = par[14];

    double A_vbo2 = par[15];

    double Tau_vbo2 = par[16];

    double omega_vbo2 = par[17];

    double phi_vbo2 = par[18];*/
    
    double time = x[0];
    
    double omega_a = rawBinToNs * 1.e-3 * getBlinded.paramToFreq(R);

        return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi)) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo)));

    //    return norm * exp(-time/life) * (1 + asym* (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo)))*cos(omega_a*time + phi)) * (1-A1_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi1_cbo)));
    //* (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) *  (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));

    }

void run2_9par_nofr_14000sr()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[1]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open("run2_14000sr_thresh_1134_EC_p.root");
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


   
   /*TString blindingString("Ritwika's blinding!"); // your blinding string
  TRandom3 *r3 = new TRandom3(Hash(blindingString));
  deltaR = 24.*2*(r3->Rndm() - 0.5); // random number +- 24 ppm
   */
  // cout<<r3->Rndm()<<" "<<endl;
  
 c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);
 
 



  fit_func= new TF1("fprec", fprec,  114000,  330000, 9);
  fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  //, "A2_cbo","phi2_cbo");
  //, "A_vbo", "Tau_vbo", "omega_vbo", "phi_vbo");
  // );
  /*  fit_func->SetParName(11, "A_vbo");
     fit_func->SetParName(12, "Tau_vbo");
     fit_func->SetParName(13, "omega_vbo");
     fit_func->SetParName(14, "phi_vbo");
     fit_func->SetParName(15, "A_vbo2");
     fit_func->SetParName(16, "Tau_vbo2");
     fit_func->SetParName(17, "omega_vbo2");
     fit_func->SetParName(18, "phi_vbo2");
  */

 fit_func->SetNpx(1000000);



 gStyle->SetOptFit(1111);
 //  c2->Divide(1,3);
 //  c2->cd(1);
 // h_sum->Rebin(8);

 // for(Int_t i=1; i<=24; i++)
 // {
     cout<<"calo "<<1<<" "<<endl;
     h[m]->Rebin(8);
     //                     1580000000/2.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4
     //15800000, 51568, 0.2, 136, 3.14/2, 0.06, 148325, 0.0029, 3.14/2
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 158000, 0.0029, (3/4)*3.14);
			     //, 0, 0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
     fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
      // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
       //  fit_func->Draw("same");
     h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
     m=m+1;
     // }

	       /*       h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);
      Double_t fit_start=130000;
      Double_t fit_stop=290000;

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[m]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[m]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[m]->GetBinContent(ibin)- fit_func->Eval( h[m]->GetBinCenter(ibin) ) );
      if(h[m]->GetBinError(ibin)!=0){res=(res/h[m]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   hm->Draw();
	       */
   
     	     cout<<"calo "<<2<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.3, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14/2);
     //, 0,0);
     /*  fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     //  h[m]->Draw();
     //  fit_func->Draw("same");
     h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     
       
      cout<<"calo "<<3<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.5, 51568, 0.2, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14/2);
			     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
      //  fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
     m=m+1;
     
       

     cout<<"calo "<<4<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14/2);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
     // fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     
       

      cout<<"calo "<<5<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14/2);
     //, 0.01,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
    // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
      // fit_func->Draw("same");
        h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     


     cout<<"calo "<<6<<" "<<endl;
    h[m]->Rebin(8);
    fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14);
    //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
      //  fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
      m=m+1;
     


     cout<<"calo "<<7<<" "<<endl;
    h[m]->Rebin(8);
    fit_func->SetParameters(21800000000/3.8, 51568, 0.2, 0, 3.14/2, 0.05, 158325, 0.0029, 3.14/2);
    //, 0.0, 1.57);
    /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
    */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       //  fit_func->Draw("same");
     h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
      m=m+1;
     
     

     cout<<"calo "<<8<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/3.9, 51568, 0.2, 0, 3.14/2, 0.05, 148325, 0.0029, 3.14/2);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     // h[m]->Draw();
      //  fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
      m=m+1;
     

     cout<<"calo "<<9<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4, 51568, 0.2, 0, 3.14/2, 0.05, 158000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /*   fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
     //   h[m]->Draw();
     //  fit_func->Draw("same");
        h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
      m=m+1;

        
          cout<<"calo "<<10<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/3.9, 51568, 0.2, 0, 3.14/2, 0.01, 148000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      //  h[m]->Draw();
      // fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;

       
       cout<<"calo "<<11<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4, 51568, 0.2, 0, 3.14/2, 0.01, 168000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      // fit_func->Draw("same");
       h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       
       
        cout<<"calo "<<12<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/3.4, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //,0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     //  fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       // fit_func->Draw("same");
        h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
      m=m+1;

      
       cout<<"calo "<<13<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 150000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /*  fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      // fit_func->Draw("same");
       h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       

        cout<<"calo "<<14<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     //  fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      // fit_func->Draw("same");
            h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;

       


         cout<<"calo "<<15<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     //  fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      //  fit_func->Draw("same");
       h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       

         cout<<"calo "<<16<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       // fit_func->Draw("same");
         h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
      n[m]=m;
        m=m+1;
	

           cout<<"calo "<<17<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     //  fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      //  fit_func->Draw("same");
         h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       

           cout<<"calo "<<18<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       //  fit_func->Draw("same");
	   h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       

            cout<<"calo "<<19<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /*  fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
   */
     // c2->cd(i);
     //  fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      //  h[m]->Draw();
      //  fit_func->Draw("same");
	 	   h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;
       

      cout<<"calo "<<20<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 140000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       //  fit_func->Draw("same");
        h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;


        cout<<"calo "<<21<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.4, 51568, 0.2, 0, 3.14/2, 0.05, 150000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
      //   fit_func->Draw("same");
      h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
       m=m+1;

       
     cout<<"calo "<<22<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/3.5, 51568, 0.2, 0, 3.14/2, 0.05, 148000, 0.0029, (3*3.14)/4);
     //, 0,0);
     /*  fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      // h[m]->Draw();
       //  fit_func->Draw("same");
           h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
     m=m+1;

     
                     cout<<"calo "<<23<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/3.4, 51568, 0.2, 0, 3.14/2, 0.05, 150000, 0.0029, 0);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      //   h[m]->Draw();
      //  fit_func->Draw("same");
       h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
            m=m+1;


             cout<<"calo "<<24<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(21800000000/4.1, 51568, 0.2, 0, 3.14/2, 0.05, 150000, 0.0029, 0);
     //, 0,0);
     /* fit_func->SetParameter(11,0);
     fit_func->SetParameter(12,0);
     fit_func->SetParameter(13,0);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,0);
     fit_func->SetParameter(16,0);
     fit_func->SetParameter(17,0);
     fit_func->SetParameter(18,0);
     */
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 1500000000.);
      //  h[m]->Draw();
      // fit_func->Draw("same");
     h[m]->Fit("fprec","R","",130000,290000);
     norm[m]=fit_func->GetParameter(5);
     dnorm[m]=fit_func->GetParError(5);
     life[m]=fit_func->GetParameter(6);
     dlife[m]=fit_func->GetParError(6);
     freq[m]=fit_func->GetParameter(7);
     dfreq[m]=fit_func->GetParError(7);
     phi[m]=fit_func->GetParameter(8);
     dphi[m]=fit_func->GetParError(8);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      cout<<"reduced chi "<<chi[m]<<" "<<endl;
     n[m]=m;
      m=m+1;

        hcalo= new TH1D("hcalo","hcalo", h[m]->GetNbinsX(), 100001, 352001);
   for(Int_t j=1;j<=h_sum->GetNbinsX();j++)
     {
       hcalo->SetBinContent(j, h_sum->GetBinContent(j));
       hcalo->SetBinError(j,h_sum->GetBinError(j));
     }

   
      hcalo->Rebin(4);

   h_sum->Rebin(8);
   fit_func->SetParameters(218000000000/1.9, 51568, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
   /*   fit_func->SetParameter(11,10000);
     fit_func->SetParameter(12,10000);
     fit_func->SetParameter(13,0.01);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,10000);
     fit_func->SetParameter(16,2000);
     fit_func->SetParameter(17,0.01);
     fit_func->SetParameter(18,0);
   */

        h_sum->GetYaxis()->SetRangeUser(-100., 30000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h_sum->Draw();
	  // fit_func->Draw("same");
	  h_sum->Fit("fprec","R","",130000,330000);
     
            h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);
      fit_start=130000;
      fit_stop=330000;

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h_sum->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum->GetBinContent(ibin)- fit_func->Eval( h_sum->GetBinCenter(ibin) ) );
      if(h_sum->GetBinError(ibin)!=0){res=(res/h_sum->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   hm->Draw();     


     
     c3=new TCanvas("c3","cbo parameters vs calo#");
     c3->Divide(2,2);
     c3->cd(1);
     gr1=new TGraphErrors(m,n,norm,0,dnorm);
    gr1->SetTitle("A_cbo vs calo#");
    gr1->GetXaxis()->SetTitle("calo #");
    gr1->GetYaxis()->SetTitle("A_cbo");
     gr1->SetMarkerStyle(20);
    gr1->Draw();

    c3->cd(2);
      gr2=new TGraphErrors(m,n,life,0,dlife);
    gr2->SetTitle("Tau_cbo vs calo#");
    gr2->GetXaxis()->SetTitle("calo #");
    gr2->GetYaxis()->SetTitle("Tau_cbo");
     gr2->SetMarkerStyle(20);
    gr2->Draw();

    c3->cd(3);
     gr3=new TGraphErrors(m,n,freq,0,dfreq);
    gr3->SetTitle("omega_cbo vs calo#");
    gr3->GetXaxis()->SetTitle("calo #");
    gr3->GetYaxis()->SetTitle("omega_cbo");
     gr3->SetMarkerStyle(20);
    gr3->Draw();

    c3->cd(4);
     gr4=new TGraphErrors(m,n,phi,0,dphi);
    gr4->SetTitle("phi_cbo vs calo#");
    gr4->GetXaxis()->SetTitle("calo #");
    gr4->GetYaxis()->SetTitle("phi_cbo");
     gr4->SetMarkerStyle(20);
    gr4->Draw();
   
      
     c4=new TCanvas("c4","chi-squared vs calo#");
     
     gr5=new TGraphErrors(m,n,chi,0,0);
    gr5->SetTitle("chi2/ndf vs calo#");
    gr5->GetXaxis()->SetTitle("calo #");
    gr5->GetYaxis()->SetTitle("chi2/ndf");
     gr5->SetMarkerStyle(20);
    gr5->Draw();
	
    /*     c5=new TCanvas("c5","sigma"); 
   hsigma= new TH1D("sigma histogram", "hsigma", 100000, -5, 5);
   for(Int_t u=1; u<=24; u++)
     {
       for(Int_t v=1; v<=hsigma->GetNbinsX(); v++)
	 {
           if((int (chi[u]*100000))==(int (hsigma->GetBinCenter(v)*100000)))
	     {
	       //	       cout<<"voila "<<v<<" "<<hsigma->GetBinCenter(v)<<endl;
	       hsigma->AddBinContent(v,1);
	     }
	 }
     }
   // hsigma->Rebin(10);
   hsigma->Draw();
   // hsigma->Fit("gaus");     	            
   */
    double sigma;
     for(Int_t u=1; u<=24; u++)
     {
       sigma=(chi[u]-1)*(chi[u]-1);
     }

     sigma=sqrt(sigma/24);

     cout<<"standard deviation is"<<sigma;
      
 }
