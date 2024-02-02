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


char root_file_name[128] = "run2e_corr.root";
TH1D *qHist_1D[25],*h[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno,*h_res, *hfftsq;
TH1 *hfft,*hfft2;
TH1F *hlm;
TF1 *fit_func, *fit_func2;
TCanvas *c2;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t T_a = 4200;
Double_t fit_start, fit_stop;

Double_t fprec(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    return asym*cos(omega_a*time + phi);
}

Double_t fprec2(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_cbo_A = par[7];

    double phi_cbo_A= par[8];

    double A_cbo_phi= par[9];

    double phi_cbo_phi=par[10];

    double asym_vw= par[11];

    double tau_vw = par[12];

    double omega_vw = par[13];

    double phi_vw = par[14];

    double lost_muon_amp = par[15];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double Acbob=(1+ asym_cbo_A*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double phicbob=(A_cbo_phi*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double Nlm= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3)));

    double Nlmf= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time + T_a/2) * 1.0e-3)));

    double Nlmb= ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( (time - T_a/2) * 1.0e-3)));

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    double fb=(1+ asym*Acbob*cos(omega_a*(time - T_a/2) + phi + phicbob));
    
    return (2*f*Ncbo*Nvw*Nlm - ff*Ncbof*Nvwf*Nlmf - fb*Ncbob*Nvwb*Nlmb)/(2*f*Ncbo*Nvw*Nlm + ff*Ncbof*Nvwf*Nlmf + fb*Ncbob*Nvwb*Nlmb);
    //return 2*f + ff + fb;
}


void qratio_binshift()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c2=new TCanvas("c2","ratio wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");

  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open(root_file_name);
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open(root_file_name);
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open(root_file_name);
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open(root_file_name);
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

  _file[11]=TFile::Open(root_file_name);
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open(root_file_name);
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open(root_file_name);
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open(root_file_name);
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open(root_file_name);
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open(root_file_name);
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open(root_file_name);
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open(root_file_name);
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open(root_file_name);
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open(root_file_name);
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open(root_file_name);
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open(root_file_name);
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open(root_file_name);
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open(root_file_name);
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

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   //   h_sum->Draw();

   h_sum->Rebin(4);
   h_sum->Scale(0.25);
   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
     {
       hp->SetBinContent(ibin,h_sum->GetBinContent(ibin+28));
     }

   for(int ibin=28;ibin<=h_sum->GetNbinsX();ibin++)
     {
       hm->SetBinContent(ibin,h_sum->GetBinContent(ibin-28));
     }
   h1->Scale(0.25);
   h2->Scale(0.25);
   hp->Scale(0.25);
   hm->Scale(0.25);

     h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

     for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
     {
       h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
     }
   

    h_num=new TH1D("calo histogram sum numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
     h_deno=new TH1D("calo histogram sum denominator", "h_deno", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

     h_num->Add(h1,1);
     h_num->Add(h2,1);
     h_num->Add(hp,-1);
     h_num->Add(hm,-1);

     h_deno->Add(h1,1);
     h_deno->Add(h2,1);
     h_deno->Add(hp,1);
     h_deno->Add(hm,1);

     h_num->Divide(h_deno);

     h_num->Rebin(2);
     h_num->Scale(0.5);


     h_ratio->Rebin(2);
     h_ratio->Scale(0.5);
     
     /*   double t1,t2,t3,t4,t5,t6;
   for(int ibin=1;ibin<=h_num->GetNbinsX();ibin++)
     {
       t1=4*(h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(h_sum->GetBinError(2*ibin-1))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t2=(-4)*(h_sum->GetBinContent(2*ibin-1))*(h_sum->GetBinError(2*ibin-30))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t3=(-4)*(h_sum->GetBinContent(2*ibin-1))*(h_sum->GetBinError(2*ibin+28))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t4=4*(h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(h_sum->GetBinError(2*ibin))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t5=(-4)*(h_sum->GetBinContent(2*ibin))*(h_sum->GetBinError(2*ibin-29))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t6=(-4)*(h_sum->GetBinContent(2*ibin))*(h_sum->GetBinError(2*ibin+29))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t1=t1/2;
       t2=t2/2;
       t3=t3/2;
       t4=t4/2;
       t5=t5/2;
       t6=t6/2;
       h_ratio->SetBinError(ibin,sqrt((t1*t1)+(t2*t2)+(t3*t3)+(t4*t4)+(t5*t5)+(t6*t6)));
     }
     */
       fit_func= new TF1("fprec", fprec,  30000,  309000, 3);
     fit_func->SetParNames("A", "Blind R", "Phi");
     fit_func->SetNpx(1000000);

        fit_func2= new TF1("fprec2", fprec2,  30000,  309000, 16);
	fit_func2->SetParNames("A", "Blind R", "Phi", "A_cbo", "tau_cbo", "omega_cbo", "phi_cbo","asym_cbo_A","asym_cbo_phi","phi_cbo_A","phi_cbo_phi");
	fit_func2->SetParName(11,"A_vw");
	fit_func2->SetParName(12,"tau_vw");
	fit_func2->SetParName(13,"omega_vw");
	fit_func2->SetParName(14,"phi_vw");
	fit_func2->SetParName(15,"lm_amp");

     fit_func2->SetNpx(1000000);
     fit_func2->SetParameters(0.243, 0.0, 2.19, 0.003,225000,0.00234,0.25);
     fit_func2->SetParameter(7,-0.0005);
     fit_func2->SetParameter(8,5.82);
     fit_func2->SetParameter(9,0.000234);
     fit_func2->SetParameter(10,-4.5);
     fit_func2->SetParameter(11,0.00075);
     fit_func2->SetParameter(12,9000);
     fit_func2->SetParameter(13,0.1538);
     fit_func2->SetParameter(14,-3.38);
     fit_func2->SetParameter(15,0);

     // fit_func2->FixParameter(13,0.1535);
     // fit_func2->FixParameter(15,0);
  
     


     //  gStyle->SetOptFit(1111);

     fit_start=30000;
     fit_stop=300000;

     fit_func->SetParameters(0.2, 0.0, 3.14/2);
     c2->Divide(1,3);
     c2->cd(1);

     //h_num->GetYaxis()->SetRangeUser(-100., 30000000000.);
	h_num->GetYaxis()->SetTitle("ADC counts");
	h_num->GetXaxis()->SetTitle("time [ns]");
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//	h_num->Draw();
	//fit_func2->Draw("same");
		h_num->Fit("fprec2","R","",fit_start,fit_stop);

		h_num->Draw();
		h_res= new TH1D("residual histogram", "h_res", h_num->GetNbinsX(), h_num->GetBinLowEdge(1), (h_num->GetBinLowEdge(h_num->GetNbinsX())+h_num->GetBinWidth(1)));


        c2->cd(2);
	for (int ibin = ((fit_start + 1- h_num->GetBinLowEdge(1))/h_num->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_num->GetBinLowEdge(1))/h_num->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_num->GetBinContent(ibin)- fit_func2->Eval( h_num->GetBinCenter(ibin) ) );
      if(h_num->GetBinError(ibin)!=0){res=(res/h_num->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }

	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");	
        h_res->Draw();

  c2->cd(3);
     hfft = h_res->FFT(hfft, "MAG");
   hfft->SetLineColor(kBlack);
   hfft->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hfft->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1)))/2);
   hfft->GetXaxis()->SetTitle("Freq [GHz]");
   hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft->Draw();
     
 }
