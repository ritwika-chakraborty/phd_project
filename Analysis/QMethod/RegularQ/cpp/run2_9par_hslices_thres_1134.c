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
TH1D *qHist_1D[25], *h[25], *h_sum, *hcalo;
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res;
TH1 *hm;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10;
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

double paramToFreq(double blindedValue)
{

  double unblindedR = blindedValue - deltaR;
  return 2 * TMath::Pi() * refFreq * (1 + (unblindedR * precisionR));
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
    
    double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

        return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi)) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo)));

    //    return norm * exp(-time/life) * (1 + asym* (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo)))*cos(omega_a*time + phi)) * (1-A1_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi1_cbo)));
    //* (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) *  (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));

    }

void run2_9par_hslices_thres_1134()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40],norm1[40],dnorm1[40],norm2[40],dnorm2[40],norm3[40],dnorm3[40],norm4[40],dnorm4[40],norm5[40],dnorm5[40], norm6[40], norm7[40], norm8[40], norm9[40], norm10[40], dnorm6[40], dnorm7[40], dnorm8[40], dnorm9[40], dnorm10[40], diff[40],ratio[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c1=new TCanvas("c1","9 parameter wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

  _file[1]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_hslice_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_hslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_hslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_hslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_hslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2_982sr_thresh_1134.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_hslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();



  
  gStyle->SetOptFit(1111);
 
   

    
 

   
  TString blindingString("Ritwika's blinding!"); // your blinding string
  TRandom3 *r3 = new TRandom3(Hash(blindingString));
  deltaR = 24.*2*(r3->Rndm() - 0.5); // random number +- 24 ppm

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
     fit_func->SetParameters(1580000000/3.3, 51568, 0.2, 136, 3.14/2, 0.05, 158000, 0.0029, (3/4)*3.14);
			     //, 0, 0);

     // c2->cd(i);
      // fit_func->SetParLimits(8, -6.28, 6.28);
      h[m]->GetYaxis()->SetRangeUser(-100., 120000000.);
     // h[m]->Draw();
       //  fit_func->Draw("same");
      //   h[m]->Fit("fprec","R","",130000,330000);
     norm[m]=fit_func->GetParameter(0);
     dnorm[m]=fit_func->GetParError(0);
h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);

     //  diff[m]=norm1[m]-norm1[1];
     //ratio[m]=norm1[m]/norm1[1];
   
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //  h[m]->Fit("fprec","R","",230000,330000);
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
                    m=m+1;
     // }

		    /*       h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);
      Double_t fit_start=130000;
      Double_t fit_stop=330000;

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
     fit_func->SetParameters(15800000000/8.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 158000, 0.0029, 3.14/2);
     //, 0,0);
     
     //  fit_func->SetParameter(11,0);
     // fit_func->SetParameter(12,0);
     //fit_func->SetParameter(13,0);
     //fit_func->SetParameter(14,0);
     // fit_func->SetParameter(15,0);
     //fit_func->SetParameter(16,0);
     //fit_func->SetParameter(17,0);
     //fit_func->SetParameter(18,0);
		      
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 120000000.);
     //    h[m]->Draw();
     //  fit_func->Draw("same");
     //   h[m]->Fit("fprec","R","",130000,330000);
    h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);
     //     diff[m]=norm1[m]-norm1[1];
     // ratio[m]=norm1[m]/norm1[1];
   
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
      m=m+1;

      
 
 
           cout<<"calo "<<3<<" "<<endl;
     h[m]->Rebin(8);
     //   fit_func->SetParameters(15800000000/5.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 150000, 0.0029, 3.14/2);
     fit_func->SetParameters(1580000000/3.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 160000, 0.0029, 3.14/2);
			     //, 0,0);
     //  fit_func->SetParameter(11,0);
     // fit_func->SetParameter(12,0);
     // fit_func->SetParameter(13,0);
     //fit_func->SetParameter(14,0);
     // fit_func->SetParameter(15,0);
     //fit_func->SetParameter(16,0);
     //fit_func->SetParameter(17,0);
     //fit_func->SetParameter(18,0);
     
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     //   h[m]->GetYaxis()->SetRangeUser(-100., 480000000.);
      h[m]->GetYaxis()->SetRangeUser(-100., 80000000.);
      //   h[m]->Draw();
      //  fit_func->Draw("same");
      //	     h[m]->Fit("fprec","R","",130000,330000);
     norm[m]=fit_func->GetParameter(0);
     dnorm[m]=fit_func->GetParError(0);
h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);
     //   diff[m]=norm1[m]-norm1[1];
     // ratio[m]=norm1[m]/norm1[1];
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
        m=m+1;
      

         cout<<"calo "<<4<<" "<<endl;
     h[m]->Rebin(8);
      fit_func->SetParameters(15800000000/5.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 150000, 0.0029, 3.14/2);
     //, 0,0);
      //  fit_func->SetParameter(11,0);
      //fit_func->SetParameter(12,0);
      //fit_func->SetParameter(13,0);
      //fit_func->SetParameter(14,0);
      // fit_func->SetParameter(15,0);
      //fit_func->SetParameter(16,0);
      //fit_func->SetParameter(17,0);
      //fit_func->SetParameter(18,0);
     
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 480000000.);
     // h[m]->Draw();
     // fit_func->Draw("same");
     //   h[m]->Fit("fprec","R","",130000,330000);
     norm[m]=fit_func->GetParameter(0);
     dnorm[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);

     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
      //      diff[m]=norm1[m]-norm1[1];
      // ratio[m]=norm1[m]/norm1[1];
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     


      cout<<"calo "<<5<<" "<<endl;
     h[m]->Rebin(8);
     fit_func->SetParameters(15800000000/8.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 160000, 0.0029, 3.14/2);
     //   fit_func->SetParameter(11,0);
     // fit_func->SetParameter(12,0);
     // fit_func->SetParameter(13,0);
     //fit_func->SetParameter(14,0);
     // fit_func->SetParameter(15,0);
     //fit_func->SetParameter(16,0);
     // fit_func->SetParameter(17,0);
     // fit_func->SetParameter(18,0);
     
    // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 320000000.);
     // h[m]->Draw();
      // fit_func->Draw("same");
     //   h[m]->Fit("fprec","R","",130000,330000);
     norm[m]=fit_func->GetParameter(0);
     dnorm[m]=fit_func->GetParError(0);
    h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     //  diff[m]=norm1[m]-norm1[1];
     // ratio[m]=norm1[m]/norm1[1];
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     


     cout<<"calo "<<6<<" "<<endl;
    h[m]->Rebin(8);
    fit_func->SetParameters(1580000000/3.5, 51568, 0.2, 135.54, 3.14/2, 0.05, 138000, 0.0029, 3.14/2);
    //  fit_func->SetParameter(11,0);
    // fit_func->SetParameter(12,0);
    // fit_func->SetParameter(13,0);
    // fit_func->SetParameter(14,0);
    //  fit_func->SetParameter(15,0);
    // fit_func->SetParameter(16,0);
    // fit_func->SetParameter(17,0);
    // fit_func->SetParameter(18,0);
     
     // c2->cd(i);
     // fit_func->SetParLimits(8, -6.28, 6.28);
     h[m]->GetYaxis()->SetRangeUser(-100., 80000000.);
     // h[m]->Draw();
      //  fit_func->Draw("same");
     //  h[m]->Fit("fprec","R","",130000,330000);
     norm[m]=fit_func->GetParameter(0);
     dnorm[m]=fit_func->GetParError(0);
    h[m]->Fit("fprec","R","",130000,140500);
     norm1[m]=fit_func->GetParameter(0);
     dnorm1[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",140500,151000);
     norm2[m]=fit_func->GetParameter(0);
     dnorm2[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",151000,161500);
     norm3[m]=fit_func->GetParameter(0);
     dnorm3[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",161500,172000);
     norm4[m]=fit_func->GetParameter(0);
     dnorm4[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",172000,182500);
     norm5[m]=fit_func->GetParameter(0);
     dnorm5[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",182500,193000);
     norm6[m]=fit_func->GetParameter(0);
     dnorm6[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",193000,214000);
     norm7[m]=fit_func->GetParameter(0);
     dnorm7[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",214000,235000);
     norm8[m]=fit_func->GetParameter(0);
     dnorm8[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",225000,260000);
     norm9[m]=fit_func->GetParameter(0);
     dnorm9[m]=fit_func->GetParError(0);
     h[m]->Fit("fprec","R","",260000,330000);
     norm10[m]=fit_func->GetParameter(0);
     dnorm10[m]=fit_func->GetParError(0);
     //   diff[m]=norm1[m]-norm1[1];
     //ratio[m]=norm1[m]/norm1[1];
     chi[m]=fit_func->GetChisquare()/fit_func->GetNDF();
     n[m]=m;
     cout<<"reduced chi "<<chi[m]<<" "<<endl;
       m=m+1;
     
     

   
      /*        hcalo= new TH1D("hcalo","hcalo", h[m]->GetNbinsX(), 100001, 352001);
   for(Int_t j=1;j<=h_sum->GetNbinsX();j++)
     {
       hcalo->SetBinContent(j, h_sum->GetBinContent(j));
       hcalo->SetBinError(j,h_sum->GetBinError(j));
     }

   
      hcalo->Rebin(4);
      */
      // h_sum->Rebin(8);
      //fit_func->SetParameters(15800000000/1.5, 51568, 0.2, 135.54, 3.14/2, 0.03, 100000, 0.0029, 0.1);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
   /*   fit_func->SetParameter(11,10000);
     fit_func->SetParameter(12,10000);
     fit_func->SetParameter(13,0.01);
     fit_func->SetParameter(14,0);
      fit_func->SetParameter(15,10000);
     fit_func->SetParameter(16,2000);
     fit_func->SetParameter(17,0.01);
     fit_func->SetParameter(18,0);*/
   

      /*       h_sum->GetYaxis()->SetRangeUser(-100., 2000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//  h_sum->Draw();
	//  fit_func->Draw("same");
	  	h_sum->Fit("fprec","R","",130000,340000);
     
            h_res= new TH1D("residual histogram", "h_res", 2100, 100001, 352001);
      fit_start=130000;
      fit_stop=340000;

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

      */
     
      /*   c3=new TCanvas("c3","cbo parameters vs calo#");
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
      */
      
          c4=new TCanvas("c4","chi-squared vs calo#");
     
    gr1=new TGraphErrors(m,n,norm1,0,dnorm1);
    gr1->SetTitle("norm vs slice#");
    gr1->GetXaxis()->SetTitle("hslice #");
    gr1->GetYaxis()->SetTitle("norm");
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerSize(2);
    gr1->Draw();
	  
    gr2=new TGraphErrors(m,n,norm2,0,dnorm2);
    gr2->SetLineColor(kRed);
     gr2->SetMarkerStyle(8);
    gr2->SetMarkerSize(2);
    gr2->Draw("same");

    gr3=new TGraphErrors(m,n,norm3,0,dnorm3);
    gr3->SetLineColor(kGreen);	
    gr3->SetMarkerStyle(8);
    gr3->SetMarkerSize(2);
    gr3->Draw("same");

    gr4=new TGraphErrors(m,n,norm4,0,dnorm4);
    gr4->SetLineColor(kBlue);
     gr4->SetMarkerStyle(8);
    gr4->SetMarkerSize(2);
   gr4->Draw("same");

    gr5=new TGraphErrors(m,n,norm5,0,dnorm5);
    gr5->SetLineColor(kYellow);
    gr5->SetMarkerStyle(8);
    gr5->SetMarkerSize(2);
    gr5->Draw("same");

    gr6=new TGraphErrors(m,n,norm6,0,dnorm6);
    gr6->SetLineColor(kMagenta);
    gr6->SetMarkerStyle(8);
    gr6->SetMarkerSize(2);
    gr6->Draw("same");

    gr7=new TGraphErrors(m,n,norm7,0,dnorm7);
    gr7->SetLineColor(kCyan);
    gr7->SetMarkerStyle(8);
    gr7->SetMarkerSize(2);
    gr7->Draw("same");

    gr8=new TGraphErrors(m,n,norm8,0,dnorm8);
    gr8->SetLineColor(kTeal);
    gr8->SetMarkerStyle(8);
    gr8->SetMarkerSize(2);
    gr8->Draw("same");

    gr9=new TGraphErrors(m,n,norm9,0,dnorm9);
    gr9->SetLineColor(kViolet);
    gr9->SetMarkerStyle(8);
    gr9->SetMarkerSize(2);
    gr9->Draw("same");

    gr10=new TGraphErrors(m,n,norm10,0,dnorm10);
    gr10->SetLineColor(kGray);
    gr10->SetMarkerStyle(8);
    gr10->SetMarkerSize(2);
    gr10->Draw("same");
      
      c3=new TCanvas("c3","cbo parameters vs calo#");
     c3->Divide(1,2);
     c3->cd(1);
      
      
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

      
 }
