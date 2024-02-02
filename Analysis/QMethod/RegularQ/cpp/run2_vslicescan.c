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
Double_t fit_start, fit_stop;



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

    //    return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3*1.25)))  * (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) * (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));

    return norm * exp(-time/life) * (1 + asym * cos( omega_a*time + phi)) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo)));
    //* ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3*1.25)));

    //     return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo)));

  
   
    }

void run2_vslicescan()

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
 
  _file[1]=TFile::Open("run2c_thresh_300_EC.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_vslice_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2c_thresh_300_EC.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_vslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2c_thresh_300_EC.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_vslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2c_thresh_300_EC.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_vslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2c_thresh_300_EC.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_vslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2c_thresh_300_EC.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_vslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

    _file[7]=TFile::Open("run2c_thresh_300_EC.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_vslice_6",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();


    _file[8]=TFile::Open("run2c_thresh_300_EC.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_vslice_7",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();


    _file[9]=TFile::Open("run2c_thresh_300_EC.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_vslice_8",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

  
  h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);


  
 
  
 



  fit_func= new TF1("fprec", fprec,  114000,  330000, 9);
  fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo");
  //, "A2_cbo","phi2_cbo");
  //, "A_vbo", "Tau_vbo", "omega_vbo", "phi_vbo");
  // );
  /*   fit_func->SetParName(11, "A3_cbo");
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
  */ 
    
 fit_func->SetNpx(1000000);

   c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

   h[1]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


      fit_func->SetParameters(10081300000, 51537.3, 0.137692, 0, 2.13957, -0.000781678, 570423, 0.00292886, 1.71319);

      /*        fit_func->SetParameter(9,-0.00987335);//A2_cbo
     fit_func->SetParameter(10, 2.46424); //phi2_cbo
     fit_func->SetParameter(11,-0.00914354); //A3_cbo
     fit_func->SetParameter(12, 2.63175); //phi3_cbo
      */
      // fit_func->FixParameter(12,0);

      //  fit_func->SetParameter(13, 0.0001); //lost muon
	// fit_func->FixParameter(9, 0);
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

      /*   fit_func->SetParameter(14, -0.488849); //A_vbo
     fit_func->SetParameter(15,  71099.1); //Tau_vbo
     fit_func->SetParameter(16, 0.0174065); //omega_vbo
     fit_func->SetParameter(17, 1.31441); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.534538); //A2_vbo
     fit_func->SetParameter(19, 75551.6); //Tau2_vbo
     fit_func->SetParameter(20, 0.192033); //omega2_vbo
     fit_func->SetParameter(21, -4.95067); //phi2_vbo
      */ 
     //   fit_func->FixParameter(16, 0.0175628); //omega_vbo
     //   fit_func->FixParameter(15, 24003.8); //Tau_vbo
     //    fit_func->FixParameter(20, 0.192033); //omega2_vbo
     //  fit_func->FixParameter(19, 75551.6); //Tau2_vbo
     //fit_func->FixParameter(5,0); //A_cbo
	/*    fit_func->FixParameter(9,0); //A2_cbo
       fit_func->FixParameter(11,0); //A3_cbo
       fit_func->FixParameter(14,0); //A_vbo
       fit_func->FixParameter(18,0); //A2_vbo
       fit_func->FixParameter(13,0); //lost muon
	*/
	//  fit_func->FixParameter(7,0.002926);
     
       
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
        h[1]->GetYaxis()->SetRangeUser(-100., 1500000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	// h[1]->Draw();
	//	fit_func->Draw("same");
      	h[1]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[1]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[1]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[1]->GetBinContent(ibin)- fit_func->Eval( h[1]->GetBinCenter(ibin) ) );
      if(h[1]->GetBinError(ibin)!=0){res=(res/h[1]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 220);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();
   
   h_res->Reset();
   hm->Reset();
   c2->Delete();


    c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

 h[9]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


fit_func->SetParameters(33600000000, 51405.5, 0.237432, 0, 2.21583, -0.00494141, 206796, 0.00292631, 3.14/4);

/*   fit_func->SetParameter(9,-0.00177251);//A2_cbo
     fit_func->SetParameter(10,0); //phi2_cbo
     fit_func->SetParameter(11,-0.00233818); //A3_cbo
     fit_func->SetParameter(12,0); //phi3_cbo
*/   // fit_func->FixParameter(12,0);

     //  fit_func->SetParameter(13, 0.0); //lost muon

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
     fit_func->SetParameter(16,0.0263706); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.488191); //A2_vbo
     fit_func->SetParameter(19, 71144.1); //Tau2_vbo
     fit_func->SetParameter(20,0.192032); //omega2_vbo
     fit_func->SetParameter(21,-1.57859); //phi2_vbo
     */
     //   fit_func->FixParameter(16, 0.0263706); //omega_vbo
     // fit_func->FixParameter(15, 76526); //Tau_vbo
     //  fit_func->FixParameter(20, 0.192032); //omega2_vbo
     // fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //  fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //   fit_func->FixParameter(13,0); //lost muon
     
     
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
        h[9]->GetYaxis()->SetRangeUser(-100., 5000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//	h[9]->Draw();
	//	fit_func->Draw("same");
	  	h[9]->Fit("fprec","R","",fit_start,fit_stop);
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[9]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[9]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[9]->GetBinContent(ibin)- fit_func->Eval( h[9]->GetBinCenter(ibin) ) );
      if(h[9]->GetBinError(ibin)!=0){res=(res/h[9]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 220);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();
   
   return;
      h_res->Reset();
   hm->Reset();
   c2->Delete();


    c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

 h[3]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


fit_func->SetParameters(46450600000, 51517, 0.236198, 0, 2.10495, 0.00405394, 213327, 0.00292536, 3.14/4);

     fit_func->SetParameter(9,-0.000775268);//A2_cbo
     fit_func->SetParameter(10,0); //phi2_cbo
     fit_func->SetParameter(11,0.000473534); //A3_cbo
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

     fit_func->SetParameter(14, -0.0000603170); //A_vbo
     fit_func->SetParameter(15,-4725250); //Tau_vbo
     fit_func->SetParameter(16,0.0266772); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.488191); //A2_vbo
     fit_func->SetParameter(19, 71144.1); //Tau2_vbo
     fit_func->SetParameter(20,0.192030); //omega2_vbo
     fit_func->SetParameter(21,-1.57859); //phi2_vbo
  
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
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[3]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[3]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[3]->GetBinContent(ibin)- fit_func->Eval( h[3]->GetBinCenter(ibin) ) );
      if(h[3]->GetBinError(ibin)!=0){res=(res/h[3]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 520);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();

        h_res->Reset();
   hm->Reset();
   c2->Delete();


    c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

 h[4]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


fit_func->SetParameters(46184400000, 51406.3, 0.236746, 0, 2.1015, -0.00411908, 193719, 0.00292477, 3.14/4);

     fit_func->SetParameter(9,-0.000714869);//A2_cbo
     fit_func->SetParameter(10,0); //phi2_cbo
     fit_func->SetParameter(11,0.000940921); //A3_cbo
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

     fit_func->SetParameter(14, -0.0000603170); //A_vbo
     fit_func->SetParameter(15,-4725250); //Tau_vbo
     fit_func->SetParameter(16,0.0266772); //omega_vbo
     fit_func->SetParameter(17,0); //phi_vbo
     
   
     fit_func->SetParameter(18, 0.488191); //A2_vbo
     fit_func->SetParameter(19, 71144.1); //Tau2_vbo
     fit_func->SetParameter(20,0.192030); //omega2_vbo
     fit_func->SetParameter(21,-1.57859); //phi2_vbo
  
     //   fit_func->FixParameter(16, 0.0266772); //omega_vbo
     // fit_func->FixParameter(15, 76526); //Tau_vbo
     // fit_func->FixParameter(20, 0.192030); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //  fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
     
     
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
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[4]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[4]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[4]->GetBinContent(ibin)- fit_func->Eval( h[4]->GetBinCenter(ibin) ) );
      if(h[4]->GetBinError(ibin)!=0){res=(res/h[4]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 520);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();



   h_res->Reset();
   hm->Reset();
   c2->Delete();


    c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

 h[5]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


fit_func->SetParameters(31812300000, 51544.1, 0.239049, 0, 2.21287, 0.00205411, 167517, 0.00292215, -1.02843);

     fit_func->SetParameter(9, -0.00102843);//A2_cbo
     fit_func->SetParameter(10, -2.95292); //phi2_cbo
     fit_func->SetParameter(11, -0.000629902); //A3_cbo
     fit_func->SetParameter(12, -0.828768); //phi3_cbo
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

     fit_func->SetParameter(14, -0.24166); //A_vbo
     fit_func->SetParameter(15, 72203.9); //Tau_vbo
     fit_func->SetParameter(16,0.017407); //omega_vbo
     fit_func->SetParameter(17,-1.934); //phi_vbo
     
   
     fit_func->SetParameter(18, -0.000289683); //A2_vbo
     fit_func->SetParameter(19, 369567); //Tau2_vbo
     fit_func->SetParameter(20, 0.18991); //omega2_vbo
     fit_func->SetParameter(21, -1.80206); //phi2_vbo
  
     // fit_func->FixParameter(16, 0.017407); //omega_vbo
     //fit_func->FixParameter(15, 72203.9); //Tau_vbo
     // fit_func->FixParameter(20, 0.192); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //   fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
     //fit_func->FixParameter(7,2.91616e-03); //omega_cbo
     
     
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
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[5]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[5]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[5]->GetBinContent(ibin)- fit_func->Eval( h[5]->GetBinCenter(ibin) ) );
      if(h[5]->GetBinError(ibin)!=0){res=(res/h[5]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 520);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();



     h_res->Reset();
   hm->Reset();
   c2->Delete();


    c2=new TCanvas("c2","5 parameter wiggle_fit");
 c2->Divide(1,3);
 c2->cd(1);


 gStyle->SetOptFit(1111);

 h[6]->Rebin(8);

      fit_start=130000;
      fit_stop=330000;


fit_func->SetParameters(7766410000, 51669, 0.210523, 0, 2.32289, 0.00784334, 433956, 0.00292787, -1.65506);

     fit_func->SetParameter(9, -0.00719657);//A2_cbo
     fit_func->SetParameter(10, -1.45541); //phi2_cbo
     fit_func->SetParameter(11, 0.00242554); //A3_cbo
     fit_func->SetParameter(12, 4.75028); //phi3_cbo
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

     fit_func->SetParameter(14, -0.327309); //A_vbo
     fit_func->SetParameter(15, 85198.3); //Tau_vbo
     fit_func->SetParameter(16,0.0174075); //omega_vbo
     fit_func->SetParameter(17,-2.01175); //phi_vbo
     
   
     fit_func->SetParameter(18, -0.000176130); //A2_vbo
     fit_func->SetParameter(19, 136267000); //Tau2_vbo
     fit_func->SetParameter(20, 0.189941); //omega2_vbo
     fit_func->SetParameter(21, -4.97264); //phi2_vbo
  
     // fit_func->FixParameter(16, 0.017407); //omega_vbo
     //fit_func->FixParameter(15, 72203.9); //Tau_vbo
     // fit_func->FixParameter(20, 0.192); //omega2_vbo
     //  fit_func->FixParameter(19, 35578.1); //Tau2_vbo
     //  fit_func->FixParameter(5,0); //A_cbo
     //   fit_func->FixParameter(9,0); //A2_cbo
     //  fit_func->FixParameter(11,0); //A3_cbo
     //  fit_func->FixParameter(14,0); //A_vbo
     //  fit_func->FixParameter(18,0); //A2_vbo
       //  fit_func->FixParameter(13,0); //lost muon
     //fit_func->FixParameter(7,2.91616e-03); //omega_cbo
     
     
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
     
          

        c2->cd(2);
      for (int ibin = ((fit_start + 1- 100001)/h[6]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h[6]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h[6]->GetBinContent(ibin)- fit_func->Eval( h[6]->GetBinCenter(ibin) ) );
      if(h[6]->GetBinError(ibin)!=0){res=(res/h[6]->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
  h_res->Draw();

  c2->cd(3);
     hm = h_res->FFT(hm, "MAG");
   hm->SetLineColor(kBlack);
   //  hm->SetNPX(10000);
   hm->GetYaxis()->SetRangeUser(0, 520);
   hm->GetXaxis()->SetRangeUser(0, 1050);
   hm->Draw();



   
      
 }
