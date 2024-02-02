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


TCanvas *c1, *ch1, *ch2, *ch3, *ch4, *ch5, *ch6;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo;
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
Double_t fit_start, fit_stop;

Double_t fVD(Double_t *x, Double_t *par)
{
  double A = par[0];

  double tau_A = par[1];

  double B = par[2];

  double tau_B = par[2];

  double time = x[0];

  return (A * exp(-time/tau_A)) + (B * exp(-time/tau_B));
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

void run2_hslices_early_st()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40],norm1[40],dnorm1[40],norm2[40],dnorm2[40],norm3[40],dnorm3[40],norm4[40],dnorm4[40],norm5[40],dnorm5[40], norm6[40], norm7[40], norm8[40], norm9[40], norm10[40], dnorm6[40], dnorm7[40], dnorm8[40], dnorm9[40], dnorm10[40], diff[40],ratio[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
    hlm->SetName("hlm");
 
  _file[1]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_hslice_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_hslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_hslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_hslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_hslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2c_dqc_thresh_300.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_hslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

  h_res[1]= new TH1D("residual histogram 0", "h_res 0", 2100, 100001, 352001);
  h_res[2]= new TH1D("residual histogram 1", "h_res 1", 2100, 100001, 352001);
  h_res[3]= new TH1D("residual histogram 2", "h_res 2", 2100, 100001, 352001);
  h_res[4]= new TH1D("residual histogram 3", "h_res 3", 2100, 100001, 352001);
  h_res[5]= new TH1D("residual histogram 4", "h_res 4", 2100, 100001, 352001);
  h_res[6]= new TH1D("residual histogram 5", "h_res 5", 2100, 100001, 352001);
  hdiff= new TH1D("difference:top-bottom h-slice","hdiff", 2100, 100001, 352001);


  
 
  
 



  fit_func= new TF1("fprec", fprec,  114000,  330000, 22);
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
     
    
 fit_func->SetNpx(1000000);


 gStyle->SetOptFit(1111);

  h[1]->Rebin(8);

      fit_start=108800;
      fit_stop=275600;

      ch1= new TCanvas("ch1", "ch1");
      ch1->Divide(1,3);
      ch1->cd(1);
fit_func->SetParameters(11192270000, 51004.6, 0.207044, 0, 2.3174, -0.0114425, 351840, 0.00292833, 1.52848);

     fit_func->SetParameter(9,-0.00987335);//A2_cbo
     fit_func->SetParameter(10, 2.46424); //phi2_cbo
     fit_func->SetParameter(11,-0.00914354); //A3_cbo
     fit_func->SetParameter(12, 2.63175); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

  

     //########## My VBO #########     

     fit_func->SetParameter(14, -0.988849); //A_vbo
     fit_func->SetParameter(15,  351099.1); //Tau_vbo
     fit_func->SetParameter(16, 0.0174065); //omega_vbo
     fit_func->SetParameter(17, 1.31441); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.934538); //A2_vbo
     fit_func->SetParameter(19, 17551.6); //Tau2_vbo
     fit_func->SetParameter(20, 0.191891); //omega2_vbo
     fit_func->SetParameter(21, -4.95067); //phi2_vbo
  
     fit_func->FixParameter(13,0); //lost muon
     fit_func->FixParameter(1,51520);//muon life
     fit_func->FixParameter(7,2.92672e-03);//omega_cbo
     fit_func->FixParameter(16,0.0174068);//omega_vbo
     fit_func->FixParameter(20,0.191892);//omega2_vbo
     
     
       h[1]->GetYaxis()->SetRangeUser(-100., 2000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
       //	h[1]->Draw();
       //	  fit_func->Draw("same");
		  	h[1]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        
      for (int ibin = ((fit_start + 1- 100001)/h[1]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[1]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[1]->GetBinContent(ibin)- fit_func->Eval( h[1]->GetBinCenter(ibin) ) );
      if(h[1]->GetBinError(ibin)!=0){res=(res/h[1]->GetBinContent(ibin));}
      h_res[1]->SetBinContent(ibin, (res)  );
     
      }
      ch1->cd(2);
      h_res[1]->Draw();

     hm = h_res[1]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 10);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   ch1->cd(3);
   hm->Draw();
      return;

     hm->Reset();






 gStyle->SetOptFit(1111);

 h[2]->Rebin(8);

      fit_start=128800;
      fit_stop=275600;


fit_func->SetParameters(30675200000, 51619.2, 0.237364, 0, 2.21615, -0.00506999, 203100, 0.00292511, 1.88194);

     fit_func->SetParameter(9, 0.0000868526);//A2_cbo
     fit_func->SetParameter(10, 10.4626); //phi2_cbo
     fit_func->SetParameter(11,-0.0015387); //A3_cbo
     fit_func->SetParameter(12,-3.16144); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

      //  fit_func->SetParameters(170610000000, 51526, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
      /*   fit_func->SetParameter(9,-0.001);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11,-0.0003);
     fit_func->SetParameter(12,0);

     fit_func->SetParameter(13,0);
      */



     //########## My VBO #########     

     /*  fit_func->SetParameter(14, -0.0000603170); //A_vbo
     fit_func->SetParameter(15,-4725250); //Tau_vbo
     fit_func->SetParameter(16,0.0263706); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     */

       fit_func->SetParameter(14, -0.312484); //A_vbo
     fit_func->SetParameter(15,  65991.6); //Tau_vbo
     fit_func->SetParameter(16, 0.0174077); //omega_vbo
     fit_func->SetParameter(17, 1.0902); //phi_vbo
     
   
       fit_func->SetParameter(18, 4602.99); //A2_vbo
     fit_func->SetParameter(19, 9767); //Tau2_vbo
     fit_func->SetParameter(20, 0.191945); //omega2_vbo
     fit_func->SetParameter(21, -4.21671); //phi2_vbo
     
   
     /*  fit_func->SetParameter(18, 0.488191); //A2_vbo
     fit_func->SetParameter(19, 71144.1); //Tau2_vbo
     fit_func->SetParameter(20,0.192032); //omega2_vbo
     fit_func->SetParameter(21,-1.57859); //phi2_vbo
     */
     /*   fit_func->FixParameter(16, 0.0174077); //omega_vbo
      fit_func->FixParameter(15, 65991.6); //Tau_vbo
      fit_func->FixParameter(20, 0.191945); //omega2_vbo
      fit_func->FixParameter(19, 9767); //Tau2_vbo
     */
     //  fit_func->FixParameter(5,0); //A_cbo
     //  fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
        fit_func->FixParameter(1,51520);//muon life
     
     
     //########### Tim's VBO ###########

     /*   fit_func->SetParameter(14, -0.0857001);
     fit_func->SetParameter(15,26441.3/1.25);
     fit_func->SetParameter(16,1.21776125);
     fit_func->SetParameter(17,-91.4);
     
   
     fit_func->SetParameter(18, -0.000416472);
     fit_func->SetParameter(19,92957.4/1.25);
     fit_func->SetParameter(20,1.00613*1.25);
     fit_func->SetParameter(21,-1.6304);
     

      fit_func->FixParameter(16, 1.21776125);
      fit_func->FixParameter(15, 26441.3/1.25);
      fit_func->FixParameter(20, 1.00613*1.25);
      fit_func->FixParameter(19, 92957.4/1.25);
     */
        h[2]->GetYaxis()->SetRangeUser(-100., 4000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h[1]->Draw();
	  // fit_func->Draw("same");
	  	h[2]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        
      for (int ibin = ((fit_start + 1- 100001)/h[2]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[2]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[2]->GetBinContent(ibin)- fit_func->Eval( h[2]->GetBinCenter(ibin) ) );
      if(h[2]->GetBinError(ibin)!=0){res=(res/h[2]->GetBinContent(ibin));}
      h_res[2]->SetBinContent(ibin, (res)  );
     
      }
  

  
     hm = h_res[2]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 300);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   

   //  return;
      
   hm->Reset();
   


 gStyle->SetOptFit(1111);

 h[3]->Rebin(8);

      fit_start=128800;
      fit_stop=275600;


fit_func->SetParameters(46630700000, 51520, 0.236204, 0, 2.10579, -0.00485021, 175756, 0.00292464, 1.87751);

     fit_func->SetParameter(9,-0.000618599);//A2_cbo
     fit_func->SetParameter(10,0); //phi2_cbo
     fit_func->SetParameter(11,-0.00118749); //A3_cbo
     fit_func->SetParameter(12,0); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

      //  fit_func->SetParameters(170610000000, 51526, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
      /*   fit_func->SetParameter(9,-0.001);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11,-0.0003);
     fit_func->SetParameter(12,0);

     fit_func->SetParameter(13,0);
      */



     //########## My VBO #########     

     /*   fit_func->SetParameter(14, -0.0000603170); //A_vbo
     fit_func->SetParameter(15,-4725250); //Tau_vbo
     fit_func->SetParameter(16,0.0266772); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.488191); //A2_vbo
     fit_func->SetParameter(19, 71144.1); //Tau2_vbo
     fit_func->SetParameter(20,0.192030); //omega2_vbo
     fit_func->SetParameter(21,-1.57859); //phi2_vbo
     */
        fit_func->SetParameter(14, 0.0676934); //A_vbo
     fit_func->SetParameter(15,  63960.3); //Tau_vbo
     fit_func->SetParameter(16, 0.0174078); //omega_vbo
     fit_func->SetParameter(17, 0); //phi_vbo
     
   
       fit_func->SetParameter(18, 3.71857); //A2_vbo
     fit_func->SetParameter(19, 18016.4); //Tau2_vbo
     fit_func->SetParameter(20, 0.191889); //omega2_vbo
     fit_func->SetParameter(21, 0); //phi2_vbo
   
  
     //   fit_func->FixParameter(16, 0.0266772); //omega_vbo
     // fit_func->FixParameter(15, 76526); //Tau_vbo
     // fit_func->FixParameter(20, 0.192030); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //   fit_func->FixParameter(5,0); //A_cbo
     //   fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       // fit_func->FixParameter(13,0); //lost muon
       fit_func->FixParameter(1,51520);//muon life
     
     
     //########### Tim's VBO ###########

     /*   fit_func->SetParameter(14, -0.0857001);
     fit_func->SetParameter(15,26441.3/1.25);
     fit_func->SetParameter(16,1.21776125);
     fit_func->SetParameter(17,-91.4);
     
   
     fit_func->SetParameter(18, -0.000416472);
     fit_func->SetParameter(19,92957.4/1.25);
     fit_func->SetParameter(20,1.00613*1.25);
     fit_func->SetParameter(21,-1.6304);
     

      fit_func->FixParameter(16, 1.21776125);
      fit_func->FixParameter(15, 26441.3/1.25);
      fit_func->FixParameter(20, 1.00613*1.25);
      fit_func->FixParameter(19, 92957.4/1.25);
     */
        h[3]->GetYaxis()->SetRangeUser(-100., 7000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h[1]->Draw();
	  // fit_func->Draw("same");
	  	h[3]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        
      for (int ibin = ((fit_start + 1- 100001)/h[3]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[3]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[3]->GetBinContent(ibin)- fit_func->Eval( h[3]->GetBinCenter(ibin) ) );
      if(h[3]->GetBinError(ibin)!=0){res=(res/h[3]->GetBinContent(ibin));}
      h_res[3]->SetBinContent(ibin, (res)  );
     
      }
  
     hm = h_res[3]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 300);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   //   return;
   hm->Reset();
  
   



 gStyle->SetOptFit(1111);

 h[4]->Rebin(8);

      fit_start=128800;
      fit_stop=275600;


fit_func->SetParameters(46492300000, 51491.1, 0.236766, 0, 2.10115, -0.003972515, 202241, 0.00292518, 3.14/4);

     fit_func->SetParameter(9, -0.000588348);//A2_cbo
     fit_func->SetParameter(10,0); //phi2_cbo
     fit_func->SetParameter(11,0.000768829); //A3_cbo
     fit_func->SetParameter(12,0); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

      //  fit_func->SetParameters(170610000000, 51526, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
      /*   fit_func->SetParameter(9,-0.001);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11,-0.0003);
     fit_func->SetParameter(12,0);

     fit_func->SetParameter(13,0);
      */
     //  fit_func->ReleaseParameter(16);
      fit_func->ReleaseParameter(15);
      // fit_func->ReleaseParameter(20);
     fit_func->ReleaseParameter(19);


     //########## My VBO #########     

     fit_func->SetParameter(14, -0.0625487); //A_vbo
     fit_func->SetParameter(15, 65995.3); //Tau_vbo
     fit_func->SetParameter(16,0.0174080); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.0808253); //A2_vbo
     fit_func->SetParameter(19, 32738.5); //Tau2_vbo
     fit_func->SetParameter(20,0.191901); //omega2_vbo
     fit_func->SetParameter(21,0); //phi2_vbo
  
     //   fit_func->FixParameter(16, 0.0266772); //omega_vbo
     // fit_func->FixParameter(15, 76526); //Tau_vbo
     // fit_func->FixParameter(20, 0.192030); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //  fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //   fit_func->FixParameter(13,0); //lost muon
       fit_func->FixParameter(1,51520);//muon life
     
     
     //########### Tim's VBO ###########

     /*   fit_func->SetParameter(14, -0.0857001);
     fit_func->SetParameter(15,26441.3/1.25);
     fit_func->SetParameter(16,1.21776125);
     fit_func->SetParameter(17,-91.4);
     
   
     fit_func->SetParameter(18, -0.000416472);
     fit_func->SetParameter(19,92957.4/1.25);
     fit_func->SetParameter(20,1.00613*1.25);
     fit_func->SetParameter(21,-1.6304);
     

      fit_func->FixParameter(16, 1.21776125);
      fit_func->FixParameter(15, 26441.3/1.25);
      fit_func->FixParameter(20, 1.00613*1.25);
      fit_func->FixParameter(19, 92957.4/1.25);
     */
        h[4]->GetYaxis()->SetRangeUser(-100., 7000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h[1]->Draw();
	  // fit_func->Draw("same");
	  	h[4]->Fit("fprec","R","",fit_start,fit_stop);
     
          

       
      for (int ibin = ((fit_start + 1- 100001)/h[4]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[4]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[4]->GetBinContent(ibin)- fit_func->Eval( h[4]->GetBinCenter(ibin) ) );
      if(h[4]->GetBinError(ibin)!=0){res=(res/h[4]->GetBinContent(ibin));}
      h_res[4]->SetBinContent(ibin, (res)  );
     
      }
  

  
     hm = h_res[4]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 300);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   
   //  return;
  
   hm->Reset();
   





 gStyle->SetOptFit(1111);

 h[5]->Rebin(8);

      fit_start=128800;
      fit_stop=275600;


fit_func->SetParameters(31834600000, 51552.6, 0.239048, 0, 2.21288, 0.00204738, 167854, 0.00292210, -1.01587);

     fit_func->SetParameter(9, 0.00107021);//A2_cbo
     fit_func->SetParameter(10, -24.9388); //phi2_cbo
     fit_func->SetParameter(11, 0.000627496); //A3_cbo
     fit_func->SetParameter(12, 2.36654); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

      //  fit_func->SetParameters(170610000000, 51526, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
      /*   fit_func->SetParameter(9,-0.001);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11,-0.0003);
     fit_func->SetParameter(12,0);

     fit_func->SetParameter(13,0);
      */



    fit_func->ReleaseParameter(16);
      fit_func->ReleaseParameter(15);
       fit_func->ReleaseParameter(20);
     fit_func->ReleaseParameter(19);


     //########## My VBO #########     

     fit_func->SetParameter(14, -0.274595); //A_vbo
     fit_func->SetParameter(15, 68508.7); //Tau_vbo
     fit_func->SetParameter(16,0.0174070); //omega_vbo
     fit_func->SetParameter(17,-1.92499); //phi_vbo
     
   
     fit_func->SetParameter(18, -0.00673255); //A2_vbo
     fit_func->SetParameter(19, 79149.9); //Tau2_vbo
     fit_func->SetParameter(20,0.191965); //omega2_vbo
     fit_func->SetParameter(21,-0.605600); //phi2_vbo
 
     //  fit_func->FixParameter(16, 0.017407); //omega_vbo
     //fit_func->FixParameter(15, 72203.9); //Tau_vbo
     //  fit_func->FixParameter(20, 0.191965); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //   fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
     //fit_func->FixParameter(7,2.91616e-03); //omega_cbo
       fit_func->FixParameter(1,51520);//muon life
     
     
     //########### Tim's VBO ###########

     /*   fit_func->SetParameter(14, -0.0857001);
     fit_func->SetParameter(15,26441.3/1.25);
     fit_func->SetParameter(16,1.21776125);
     fit_func->SetParameter(17,-91.4);
     
   
     fit_func->SetParameter(18, -0.000416472);
     fit_func->SetParameter(19,92957.4/1.25);
     fit_func->SetParameter(20,1.00613*1.25);
     fit_func->SetParameter(21,-1.6304);
     

      fit_func->FixParameter(16, 1.21776125);
      fit_func->FixParameter(15, 26441.3/1.25);
      fit_func->FixParameter(20, 1.00613*1.25);
      fit_func->FixParameter(19, 92957.4/1.25);
     */
        h[5]->GetYaxis()->SetRangeUser(-100., 7000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h[1]->Draw();
	  // fit_func->Draw("same");
	  	h[5]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        
      for (int ibin = ((fit_start + 1- 100001)/h[5]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[5]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[5]->GetBinContent(ibin)- fit_func->Eval( h[5]->GetBinCenter(ibin) ) );
      if(h[5]->GetBinError(ibin)!=0){res=(res/h[5]->GetBinContent(ibin));}
      h_res[5]->SetBinContent(ibin, (res)  );
     
      }


  
     hm = h_res[5]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 300);
   hm->GetXaxis()->SetRangeUser(0, 1050);
 
   //    return;
 
   hm->Reset();





 gStyle->SetOptFit(1111);

 h[6]->Rebin(8);

      fit_start=128800;
      fit_stop=275600;


fit_func->SetParameters(7841460000, 51792.1, 0.210504, 0, 2.32337, 0.00793022, 424647, 0.00292776, -1.63352);

     fit_func->SetParameter(9, -0.00700043);//A2_cbo
     fit_func->SetParameter(10, 1.66952); //phi2_cbo
     fit_func->SetParameter(11, 0.0023177); //A3_cbo
     fit_func->SetParameter(12, 4.68799); //phi3_cbo
     // fit_func->FixParameter(12,0);

     fit_func->SetParameter(13, 0.0); //lost muon

      //  fit_func->SetParameters(170610000000, 51526, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
      /*   fit_func->SetParameter(9,-0.001);
     fit_func->SetParameter(10,0);
     fit_func->SetParameter(11,-0.0003);
     fit_func->SetParameter(12,0);

     fit_func->SetParameter(13,0);
      */



     //########## My VBO #########     

     fit_func->SetParameter(14, -0.469014); //A_vbo
     fit_func->SetParameter(15, 72590.7); //Tau_vbo
     fit_func->SetParameter(16,0.0174071); //omega_vbo
     fit_func->SetParameter(17,-1.93462); //phi_vbo
     
   
     fit_func->SetParameter(18, -5.71193); //A2_vbo
     fit_func->SetParameter(19, 21492.3); //Tau2_vbo
     fit_func->SetParameter(20, 0.1919); //omega2_vbo
     fit_func->SetParameter(21, -4.50907); //phi2_vbo
  
     // fit_func->FixParameter(16, 0.0174071); //omega_vbo
     //fit_func->FixParameter(15, 72203.9); //Tau_vbo
     // fit_func->FixParameter(20, 0.1919); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //   fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
     //fit_func->FixParameter(7,2.91616e-03); //omega_cbo
       fit_func->FixParameter(1,51520);//muon life
     
     
     //########### Tim's VBO ###########

     /*   fit_func->SetParameter(14, -0.0857001);
     fit_func->SetParameter(15,26441.3/1.25);
     fit_func->SetParameter(16,1.21776125);
     fit_func->SetParameter(17,-91.4);
     
   
     fit_func->SetParameter(18, -0.000416472);
     fit_func->SetParameter(19,92957.4/1.25);
     fit_func->SetParameter(20,1.00613*1.25);
     fit_func->SetParameter(21,-1.6304);
     

      fit_func->FixParameter(16, 1.21776125);
      fit_func->FixParameter(15, 26441.3/1.25);
      fit_func->FixParameter(20, 1.00613*1.25);
      fit_func->FixParameter(19, 92957.4/1.25);
     */
        h[6]->GetYaxis()->SetRangeUser(-100., 1500000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h[1]->Draw();
	  // fit_func->Draw("same");
	  	h[6]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        
      for (int ibin = ((fit_start + 1- 100001)/h[6]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[6]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[6]->GetBinContent(ibin)- fit_func->Eval( h[6]->GetBinCenter(ibin) ) );
      if(h[6]->GetBinError(ibin)!=0){res=(res/h[6]->GetBinContent(ibin));}
      h_res[6]->SetBinContent(ibin, (res)  );
     
      }


     hm = h_res[6]->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 300);
   hm->GetXaxis()->SetRangeUser(0, 1050);

   //  c1->Clear();
   c1= new TCanvas("c","fit residuals");
   c1->Divide(1,6);
   for(Int_t k=1;k<=6;k++)
     {
    c1->cd(k);
    //  h_res[k]->Rebin(10);
    h_res[k]->Draw();
     }

         for (int ibin = ((fit_start + 1- 100001)/hdiff->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/hdiff->GetBinWidth(1))+1; ++ibin)
     {
       
       hdiff->SetBinContent(ibin,h_res[1]->GetBinContent(ibin)-h_res[6]->GetBinContent(ibin));
     
      }
   	 c2= new TCanvas("c2","top-bottom");

   	  h_res[1]->Rebin(10);
   

        VD_func= new TF1("fVD", fVD,  110000,  350000, 4);
        VD_func->SetParNames("A", "tau_A", "B", "tau_B");
	VD_func->SetParameters(10000000000, 3600, 0.01, 15300);
	//VD_func->FixParameter(1,(5.07e3)/1.25);
	VD_func->FixParameter(2, 0);
       	VD_func->FixParameter(3,(1.915e5)/1.25);
	h_res[1]->Fit("fVD","R","",fit_start,fit_stop);
	
	//h_res[1]->Draw();
	//VD_func->Draw("same");
 }
