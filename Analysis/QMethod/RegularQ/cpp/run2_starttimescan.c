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



TCanvas *c1, *c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10, *c11, *c12, *c13, *c14;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *h[25], *h_sum[10], *hcalo;
//TH1 *hm;

TH1F *hlm, h0;

TCanvas *c2[10];
TF1 *fit_func[10];
TH1D *h_res[10];
TH1 *hm[10];
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t fit_start[10], fit_stop;



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

    return norm * exp(-time/life) * (1 + asym * (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo))) * cos(omega_a*time + (phi + A3_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi3_cbo))))) * (1-A_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi_cbo))) * ( 1- lost_muon_amp * hlm->GetBinContent(hlm->GetXaxis()->FindBin( time * 1.0e-3*1.25)));

    //    return norm * exp(-time/life) * (1 + asym* (1-A2_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi2_cbo)))*cos(omega_a*time + phi)) * (1-A1_cbo * exp(-time/Tau_cbo) * (cos(omega_cbo*time+phi1_cbo)));
    //* (1-A_vbo * exp(-time/Tau_vbo) * (cos(omega_vbo*time+phi_vbo))) *  (1-A_vbo2 * exp(-time/Tau_vbo2) * (cos(omega_vbo2*time+phi_vbo2)));

    }

void run2_starttimescan()

 {
   Double_t chi[10], dchi[10], R[10], dR[10], lifetime[10],dlifetime[10],freq[10],dfreq[10],norm[10],dnorm[10],phi[10],dphi[10], Acbo[10], dAcbo[10], Tcbo[10],dTcbo[10], phicbo[10], dphicbo[10], lm[10], dlm[10], asym[10], dasym[10], Wcbo[10], dWcbo[10];
   Double_t n[10];
   Int_t ifit=0;
   Double_t N;

 
 TFile *_file[25];
 TDirectoryFile *dir[25];

   _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");
 

  _file[1]=TFile::Open("run2c_thresh_300_EC.root");
  _file[1]->GetObject("QFillByFillAnalyzer",dir[1]);
  dir[1]->GetObject("qHist1D_sig_1_0",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

  _file[2]=TFile::Open("run2c_thresh_300_EC.root");
  _file[2]->GetObject("QFillByFillAnalyzer",dir[2]);
  dir[2]->GetObject("qHist1D_sig_2_0",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open("run2c_thresh_300_EC.root");
  _file[3]->GetObject("QFillByFillAnalyzer",dir[3]);
  dir[3]->GetObject("qHist1D_sig_3_0",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open("run2c_thresh_300_EC.root");
  _file[4]->GetObject("QFillByFillAnalyzer",dir[4]);
  dir[4]->GetObject("qHist1D_sig_4_0",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open("run2c_thresh_300_EC.root");
  _file[5]->GetObject("QFillByFillAnalyzer",dir[5]);
  dir[5]->GetObject("qHist1D_sig_5_0",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

   _file[6]=TFile::Open("run2c_thresh_300_EC.root");
  _file[6]->GetObject("QFillByFillAnalyzer",dir[6]);
  dir[6]->GetObject("qHist1D_sig_6_0",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();

   _file[7]=TFile::Open("run2c_thresh_300_EC.root");
  _file[7]->GetObject("QFillByFillAnalyzer",dir[7]);
  dir[7]->GetObject("qHist1D_sig_7_0",qHist_1D[7]);
  h[7]=(TH1D*)qHist_1D[7]->Clone();

   _file[8]=TFile::Open("run2c_thresh_300_EC.root");
  _file[8]->GetObject("QFillByFillAnalyzer",dir[8]);
  dir[8]->GetObject("qHist1D_sig_8_0",qHist_1D[8]);
  h[8]=(TH1D*)qHist_1D[8]->Clone();

   _file[9]=TFile::Open("run2c_thresh_300_EC.root");
  _file[9]->GetObject("QFillByFillAnalyzer",dir[9]);
  dir[9]->GetObject("qHist1D_sig_9_0",qHist_1D[9]);
  h[9]=(TH1D*)qHist_1D[9]->Clone();

   _file[10]=TFile::Open("run2c_thresh_300_EC.root");
  _file[10]->GetObject("QFillByFillAnalyzer",dir[10]);
  dir[10]->GetObject("qHist1D_sig_10_0",qHist_1D[10]);
  h[10]=(TH1D*)qHist_1D[10]->Clone();

   _file[11]=TFile::Open("run2c_thresh_300_EC.root");
  _file[11]->GetObject("QFillByFillAnalyzer",dir[11]);
  dir[11]->GetObject("qHist1D_sig_11_0",qHist_1D[11]);
  h[11]=(TH1D*)qHist_1D[11]->Clone();

   _file[12]=TFile::Open("run2c_thresh_300_EC.root");
  _file[12]->GetObject("QFillByFillAnalyzer",dir[12]);
  dir[12]->GetObject("qHist1D_sig_12_0",qHist_1D[12]);
  h[12]=(TH1D*)qHist_1D[12]->Clone();

   _file[13]=TFile::Open("run2c_thresh_300_EC.root");
  _file[13]->GetObject("QFillByFillAnalyzer",dir[13]);
  dir[13]->GetObject("qHist1D_sig_13_0",qHist_1D[13]);
  h[13]=(TH1D*)qHist_1D[13]->Clone();

   _file[14]=TFile::Open("run2c_thresh_300_EC.root");
  _file[14]->GetObject("QFillByFillAnalyzer",dir[14]);
  dir[14]->GetObject("qHist1D_sig_14_0",qHist_1D[14]);
  h[14]=(TH1D*)qHist_1D[14]->Clone();

   _file[15]=TFile::Open("run2c_thresh_300_EC.root");
  _file[15]->GetObject("QFillByFillAnalyzer",dir[15]);
  dir[15]->GetObject("qHist1D_sig_15_0",qHist_1D[15]);
  h[15]=(TH1D*)qHist_1D[15]->Clone();

   _file[16]=TFile::Open("run2c_thresh_300_EC.root");
  _file[16]->GetObject("QFillByFillAnalyzer",dir[16]);
  dir[16]->GetObject("qHist1D_sig_16_0",qHist_1D[16]);
  h[16]=(TH1D*)qHist_1D[16]->Clone();

   _file[17]=TFile::Open("run2c_thresh_300_EC.root");
  _file[17]->GetObject("QFillByFillAnalyzer",dir[17]);
  dir[17]->GetObject("qHist1D_sig_17_0",qHist_1D[17]);
  h[17]=(TH1D*)qHist_1D[17]->Clone();

   _file[18]=TFile::Open("run2c_thresh_300_EC.root");
  _file[18]->GetObject("QFillByFillAnalyzer",dir[18]);
  dir[18]->GetObject("qHist1D_sig_18_0",qHist_1D[18]);
  h[18]=(TH1D*)qHist_1D[18]->Clone();

   _file[19]=TFile::Open("run2c_thresh_300_EC.root");
  _file[19]->GetObject("QFillByFillAnalyzer",dir[19]);
  dir[19]->GetObject("qHist1D_sig_19_0",qHist_1D[19]);
  h[19]=(TH1D*)qHist_1D[19]->Clone();

   _file[20]=TFile::Open("run2c_thresh_300_EC.root");
  _file[20]->GetObject("QFillByFillAnalyzer",dir[20]);
  dir[20]->GetObject("qHist1D_sig_20_0",qHist_1D[20]);
  h[20]=(TH1D*)qHist_1D[20]->Clone();

   _file[21]=TFile::Open("run2c_thresh_300_EC.root");
  _file[21]->GetObject("QFillByFillAnalyzer",dir[21]);
  dir[21]->GetObject("qHist1D_sig_21_0",qHist_1D[21]);
  h[21]=(TH1D*)qHist_1D[21]->Clone();

   _file[22]=TFile::Open("run2c_thresh_300_EC.root");
  _file[22]->GetObject("QFillByFillAnalyzer",dir[22]);
  dir[22]->GetObject("qHist1D_sig_22_0",qHist_1D[22]);
  h[22]=(TH1D*)qHist_1D[22]->Clone();

   _file[23]=TFile::Open("run2c_thresh_300_EC.root");
  _file[23]->GetObject("QFillByFillAnalyzer",dir[23]);
  dir[23]->GetObject("qHist1D_sig_23_0",qHist_1D[23]);
  h[23]=(TH1D*)qHist_1D[23]->Clone();

   _file[24]=TFile::Open("run2c_thresh_300_EC.root");
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
 
   
    h_sum[0]= new TH1D("calo histogram sum 0", "h_sum 0", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[0]->Sumw2(kTRUE);
  h_sum[1]= new TH1D("calo histogram sum1", "h_sum1", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[1]->Sumw2(kTRUE);
  h_sum[2]= new TH1D("calo histogram sum2", "h_sum2", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[2]->Sumw2(kTRUE);
 h_sum[3]= new TH1D("calo histogram sum3", "h_sum3", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[3]->Sumw2(kTRUE);
 h_sum[4]= new TH1D("calo histogram sum4", "h_sum4", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[4]->Sumw2(kTRUE);
 h_sum[5]= new TH1D("calo histogram sum5", "h_sum5", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[5]->Sumw2(kTRUE);
 h_sum[6]= new TH1D("calo histogram sum6", "h_sum6", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[6]->Sumw2(kTRUE);
 h_sum[7]= new TH1D("calo histogram sum7", "h_sum7", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[7]->Sumw2(kTRUE);
 h_sum[8]= new TH1D("calo histogram sum8", "h_sum8", h[1]->GetNbinsX(), 100001, 352001);
  h_sum[8]->Sumw2(kTRUE);
   h_sum[9]= new TH1D("calo histogram sum9", "h_sum 9", h[1]->GetNbinsX(), 100001, 352001);
   h_sum[9]->Sumw2(kTRUE);
 
    
  

   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   //  h_sum->Draw();

  // cout<<h_sum->GetBinContent(500)<<" "<<endl;
 c2[0]=new TCanvas("c2_0","13 parameter wiggle_fit");
 c2[1]=new TCanvas("c2_1","13 parameter wiggle_fit");
 c2[2]=new TCanvas("c2_2","13 parameter wiggle_fit");
 c2[3]=new TCanvas("c2_3","13 parameter wiggle_fit");
 c2[4]=new TCanvas("c2_4","13 parameter wiggle_fit");
 c2[5]=new TCanvas("c2_5","13 parameter wiggle_fit");
 c2[6]=new TCanvas("c2_6","13 parameter wiggle_fit");
 c2[7]=new TCanvas("c2_7","13 parameter wiggle_fit");
 c2[8]=new TCanvas("c2_8","13 parameter wiggle_fit");
 c2[9]=new TCanvas("c2_9","13 parameter wiggle_fit");


 h_res[0]= new TH1D("residual histogram 0", "h_res 0", 2100, 100001, 352001);
 h_res[1]= new TH1D("residual histogram 1", "h_res 1", 2100, 100001, 352001);
 h_res[2]= new TH1D("residual histogram 2", "h_res 2", 2100, 100001, 352001);
 h_res[3]= new TH1D("residual histogram 3", "h_res 3", 2100, 100001, 352001);
 h_res[4]= new TH1D("residual histogram 4", "h_res 4", 2100, 100001, 352001);
 h_res[5]= new TH1D("residual histogram 5", "h_res 5", 2100, 100001, 352001);
 h_res[6]= new TH1D("residual histogram 6", "h_res 6", 2100, 100001, 352001);
 h_res[7]= new TH1D("residual histogram 7", "h_res 7", 2100, 100001, 352001);
 h_res[8]= new TH1D("residual histogram 8", "h_res 8", 2100, 100001, 352001);
 h_res[9]= new TH1D("residual histogram 9", "h_res 9", 2100, 100001, 352001);
 

 fit_func[0]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[1]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[2]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[3]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[4]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[5]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[6]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[7]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[8]= new TF1("fprec", fprec,  100000,  350000, 14);
  fit_func[9]= new TF1("fprec", fprec,  100000,  350000, 14);
  
	   for(Int_t i=0; i<10;i++)
	     {
	     fit_func[i]->SetParNames("N_0", "Tau", "A", "Blind R", "Phi", "A_cbo", "Tau_cbo", "omega_cbo", "phi_cbo", "A2_cbo","phi2_cbo");

     fit_func[i]->SetParName(11, "A3_cbo");
     fit_func[i]->SetParName(12, "phi3_cbo");
     fit_func[i]->SetParName(13, "lm_amp");
    
 fit_func[i]->SetNpx(1000000);
	     }
   
	     fit_func[0]->SetParameters(150909000000, 50274, 0.2, 0, 3.14/2, -0.96, 18500, 0.0028, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
       fit_func[0]->SetParameter(9,0.4);
     fit_func[0]->SetParameter(10,0);
     fit_func[0]->SetParameter(11, 0.0001);
     fit_func[0]->SetParameter(12,0);
     fit_func[0]->SetParameter(13,0.0);

           fit_func[1]->SetParameters(170650000000, 51562, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[1]->SetParameter(9,0.4);
     fit_func[1]->SetParameter(10,0);
     fit_func[1]->SetParameter(11, 0.0001);
     fit_func[1]->SetParameter(12,0);
     fit_func[1]->SetParameter(13,0.0);
	   
       fit_func[2]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[2]->SetParameter(9,0.4);
     fit_func[2]->SetParameter(10,0);
     fit_func[2]->SetParameter(11, 0.0001);
     fit_func[2]->SetParameter(12,0);
     fit_func[2]->SetParameter(13,0.0);

       fit_func[3]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[3]->SetParameter(9,0.4);
     fit_func[3]->SetParameter(10,0);
     fit_func[3]->SetParameter(11, 0.0001);
     fit_func[3]->SetParameter(12,0);
     fit_func[3]->SetParameter(13,0.0);

      fit_func[4]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[4]->SetParameter(9,0.4);
     fit_func[4]->SetParameter(10,0);
     fit_func[4]->SetParameter(11, 0.0001);
     fit_func[4]->SetParameter(12,0);
     fit_func[4]->SetParameter(13,0.0);

     fit_func[5]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[5]->SetParameter(9,0.4);
     fit_func[5]->SetParameter(10,0);
     fit_func[5]->SetParameter(11, 0.0001);
     fit_func[5]->SetParameter(12,0);
     fit_func[5]->SetParameter(13,0.0);

     fit_func[6]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[6]->SetParameter(9,0.4);
     fit_func[6]->SetParameter(10,0);
     fit_func[6]->SetParameter(11, 0.0001);
     fit_func[6]->SetParameter(12,0);
     fit_func[6]->SetParameter(13,0.0);

     fit_func[7]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[7]->SetParameter(9,0.4);
     fit_func[7]->SetParameter(10,0);
     fit_func[7]->SetParameter(11, 0.0001);
     fit_func[7]->SetParameter(12,0);
     fit_func[7]->SetParameter(13,0.0);

     fit_func[8]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[8]->SetParameter(9,0.4);
     fit_func[8]->SetParameter(10,0);
     fit_func[8]->SetParameter(11, 0.0001);
     fit_func[8]->SetParameter(12,0);
     fit_func[8]->SetParameter(13,0.0);

         fit_func[9]->SetParameters(170870000000, 51546, 0.2, 0, 3.14/2, 0.03, 200000, 0.0029, 3.14/4);
   //, 0.002, 0.2);
   //, 1000000, 10000, 0.01, 0);
   
     fit_func[9]->SetParameter(9,0.4);
     fit_func[9]->SetParameter(10,0);
     fit_func[9]->SetParameter(11, 0.0001);
     fit_func[9]->SetParameter(12,0);
     fit_func[9]->SetParameter(13,0.0);
     
 for( ifit=0; ifit<10; ifit++) 
     {
            
        for(Int_t i=1; i<=24; i++)
        {
         h_sum[ifit]->Add(h[i],1); 
        }

	h_sum[ifit]->Rebin(8);
	
        c2[ifit]->Divide(1,3);
        c2[ifit]->cd(1);
 
         gStyle->SetOptFit(1111);

 
           if(ifit>0)
          {
           fit_start[ifit]=fit_start[ifit-1]+10000;
          }
          else
          {
           fit_start[ifit]=130000;
          }
          fit_stop=330000;
      

        h_sum[ifit]->GetYaxis()->SetRangeUser(-100., 30000000000.);
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//h_sum->Draw();
        //fit_func->Draw("same");
	h_sum[ifit]->Fit("fprec","R","",fit_start[ifit],fit_stop);
	h_sum[ifit]->Draw();
     
        n[ifit]=ifit;
	chi[ifit]=fit_func[ifit]->GetChisquare()/fit_func[ifit]->GetNDF();
	norm[ifit]=fit_func[ifit]->GetParameter(0);
	dnorm[ifit]=fit_func[ifit]->GetParError(0);
	lifetime[ifit]=fit_func[ifit]->GetParameter(1);
	dlifetime[ifit]=fit_func[ifit]->GetParError(1);
        R[ifit]=fit_func[ifit]->GetParameter(3);
	dR[ifit]=fit_func[ifit]->GetParError(3);
	asym[ifit]=fit_func[ifit]->GetParameter(2);
	dasym[ifit]=fit_func[ifit]->GetParError(2);
	phi[ifit]=fit_func[ifit]->GetParameter(4);
	dphi[ifit]=fit_func[ifit]->GetParError(4);
	Acbo[ifit]=fit_func[ifit]->GetParameter(5);
	dAcbo[ifit]=fit_func[ifit]->GetParError(5);
        Tcbo[ifit]=fit_func[ifit]->GetParameter(6);
	dTcbo[ifit]=fit_func[ifit]->GetParError(6);
        Wcbo[ifit]=fit_func[ifit]->GetParameter(7);
	dWcbo[ifit]=fit_func[ifit]->GetParError(7);
        phicbo[ifit]=fit_func[ifit]->GetParameter(8);
	dphicbo[ifit]=fit_func[ifit]->GetParError(8);
	lm[ifit]=fit_func[ifit]->GetParameter(9);
	dlm[ifit]=fit_func[ifit]->GetParError(9);







        c2[ifit]->cd(2);
      for (int ibin = ((fit_start[ifit] + 1- 100001)/h_sum[ifit]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - 100001)/h_sum[ifit]->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_sum[ifit]->GetBinContent(ibin)- fit_func[ifit]->Eval( h_sum[ifit]->GetBinCenter(ibin) ) );
      if(h_sum[ifit]->GetBinError(ibin)!=0){res=(res/h_sum[ifit]->GetBinError(ibin));}
      h_res[ifit]->SetBinContent(ibin, (res)  );
     
      }
  h_res[ifit]->Draw();

  c2[ifit]->cd(3);
  // if(ifit>0){hm[ifit-1]->Delete();}
  
   hm[ifit] = h_res[ifit]->FFT(hm[ifit], "MAG");
   hm[ifit]->SetLineColor(kBlack);
   hm[ifit]->Draw();
   
     }

  c3=new TCanvas("c3","chi-squared vs fit start-time");
  gr1=new TGraphErrors(ifit,n,chi,0,0);
    gr1->SetTitle("chi-squared vs start-time");
    gr1->GetXaxis()->SetTitle("start-time");
    gr1->GetYaxis()->SetTitle("chi-squared");
     gr1->SetMarkerStyle(20);
    gr1->Draw();

     c4=new TCanvas("c4"," vs fit start-time");
  gr2=new TGraphErrors(ifit,n,norm,0,dnorm);
    gr2->SetTitle("norm vs start-time");
    gr2->GetXaxis()->SetTitle("start-time");
    gr2->GetYaxis()->SetTitle("N_a");
    gr2->SetMarkerStyle(20);
    gr2->Draw();

     c5=new TCanvas("c5"," vs fit start-time");
  gr3=new TGraphErrors(ifit,n,lifetime,0,dlifetime);
    gr3->SetTitle("lifetime vs start-time");
    gr3->GetXaxis()->SetTitle("start-time");
    gr3->GetYaxis()->SetTitle("T_a");
    gr3->SetMarkerStyle(20);
    gr3->Draw();

     c6=new TCanvas("c6"," vs fit start-time");
  gr4=new TGraphErrors(ifit,n,asym,0,dasym);
    gr4->SetTitle("asymmetry vs start-time");
    gr4->GetXaxis()->SetTitle("start-time");
    gr4->GetYaxis()->SetTitle("A_a");
    gr4->SetMarkerStyle(20);
    gr4->Draw();

     c7=new TCanvas("c7"," vs fit start-time");
  gr5=new TGraphErrors(ifit,n,phi,0,dphi);
    gr5->SetTitle("phi vs start-time");
    gr5->GetXaxis()->SetTitle("start-time");
    gr5->GetYaxis()->SetTitle("Phi_a");
     gr5->SetMarkerStyle(20);
    gr5->Draw();

     c8=new TCanvas("c8"," vs fit start-time");
  gr6=new TGraphErrors(ifit,n,R,0,dR);
    gr6->SetTitle("Blinded R vs start-time");
    gr6->GetXaxis()->SetTitle("start-time");
    gr6->GetYaxis()->SetTitle("Blind R");
    gr6->SetMarkerStyle(20);
    gr6->Draw();


   
    

      
 }
