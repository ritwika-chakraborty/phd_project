#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"

char root_file_name1[128] ="starttime_output_avg_constmat_15_part1.root";
char root_file_name2[128] ="starttime_output_avg_constmat_15_part2.root";
char root_file_name3[128] ="starttime_output_avg_constmat_15_part3.root";
char root_file_name4[128] ="starttime_output_avg_constmat_15_part4.root";
char root_file_name5[128] ="starttime_output_avg_constmat_15_part5.root";
char root_file_name6[128] ="starttime_output_avg_constmat_15_part6.root";

TFile *_file[25];

TDirectory *dir[25];

int m=0;
int hpp=3;
TF1 *fit_func23;

TH1D *h[25],*qHist_1D[25],*h_sum, *h1[10], *h2[10], *h3[10],*h4[10], *h5[10],*h6[10],*h7[10],*h8[10],*h9[10],*h10[10],*h11[10],*h12[10],*h13[10],*h14[10],*h15[10],*h16[10],*h17[10],*h18[10],*h19[10],*h20[10],*h21[10],*h22[10],*h23[10],*h24[10], *hA, *hblindR, *hphi_0, *hA_cbo_N, *htau_cbo, *homega_cbo, *hphi_cbo_N, *hA_cbo_A, *hphi_cbo_A, *hA_cbo_phi, *hphi_cbo_phi, *hA_vw, *htau_vw, *homega_vw, *hphi_vw, *hA_y, *htau_y, *homega_y, *hphi_y, *hA_2cbo, *htau_2cbo, *homega_2cbo, *hphi_2cbo, *hn, *hchisq;

TCanvas *c,*c1, *c1a, *c1b, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10;

TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28;

TGraph *gkband1p,*gkband1m,*gkband2p,*gkband2m,*gkband3p,*gkband3m,*gkband4p,*gkband4m,*gkband5p,*gkband5m,*gkband6p,*gkband6m,*gkband7p,*gkband7m,*gkband8p,*gkband8m,*gkband9p,*gkband9m,*gkband10p,*gkband10m,*gkband11p,*gkband11m,*gkband12p,*gkband12m,*gkband13p,*gkband13m,*gkband14p,*gkband14m,*gkband15p,*gkband15m,*gkband16p,*gkband16m,*gkband17p,*gkband17m,*gkband18p,*gkband18m,*gkband19p,*gkband19m,*gkband20p,*gkband20m,*gkband21p,*gkband21m,*gkband22p,*gkband22m,*gkband23p,*gkband23m;

Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi_0[50], dphi_0[50], A_cbo_N[50], tau_cbo[50], omega_cbo[50], phi_cbo_N[50], A_cbo_A[50], phi_cbo_A[50], A_cbo_phi[50], phi_cbo_phi[50], A_vw[50], tau_vw[50], omega_vw[50], phi_vw[50], A_y[50], tau_y[50], omega_y[50], phi_y[50], A_2cbo[50], tau_2cbo[50], omega_2cbo[50], phi_2cbo[50],  dA_cbo_N[50], dtau_cbo[50], domega_cbo[50], dphi_cbo_N[50], dA_cbo_A[50], dphi_cbo_A[50], dA_cbo_phi[50], dphi_cbo_phi[50], dA_vw[50], dtau_vw[50], domega_vw[50], dphi_vw[50], dA_y[50], dtau_y[50], domega_y[50], dphi_y[50], dA_2cbo[50], dtau_2cbo[50], domega_2cbo[50], dphi_2cbo[50], n[50],kband1p[50],kband1m[50],kband2p[50],kband2m[50],kband3p[50],kband3m[50],kband4p[50],kband4m[50],kband5p[50],kband5m[50],kband6p[50],kband6m[50],kband7p[50],kband7m[50],kband8p[50],kband8m[50],kband9p[50],kband9m[50],kband10p[50],kband10m[50],kband11p[50],kband11m[50],kband12p[50],kband12m[50],kband13p[50],kband13m[50],kband14p[50],kband14m[50],kband15p[50],kband15m[50],kband16p[50],kband16m[50],kband17p[50],kband17m[50],kband18p[50],kband18m[50],kband19p[50],kband19m[50],kband20p[50],kband20m[50],kband21p[50],kband21m[50],kband22p[50],kband22m[50],kband23p[50],kband23m[50];


void starttime_scan_kbands()
{
  _file[1]=TFile::Open(root_file_name1);
  _file[2]=TFile::Open(root_file_name2);
  _file[3]=TFile::Open(root_file_name3);
  _file[4]=TFile::Open(root_file_name4);
  _file[5]=TFile::Open(root_file_name5);
  _file[6]=TFile::Open(root_file_name6);

  hA=new TH1D("hA","hA",18,0,18);
  hblindR=(TH1D*)hA->Clone();
  hblindR->Reset();
       hphi_0=(TH1D*)hA->Clone();
       hA_cbo_N=(TH1D*)hA->Clone();
       htau_cbo=(TH1D*)hA->Clone();
       homega_cbo=(TH1D*)hA->Clone();
       hphi_cbo_N=(TH1D*)hA->Clone();
       hA_cbo_A=(TH1D*)hA->Clone();
       hphi_cbo_A=(TH1D*)hA->Clone();
       hA_cbo_phi=(TH1D*)hA->Clone();
       hphi_cbo_phi=(TH1D*)hA->Clone();
       hA_vw=(TH1D*)hA->Clone();
       htau_vw=(TH1D*)hA->Clone();
       homega_vw=(TH1D*)hA->Clone();
       hphi_vw=(TH1D*)hA->Clone();
       hA_y=(TH1D*)hA->Clone();
       htau_y=(TH1D*)hA->Clone();
       homega_y=(TH1D*)hA->Clone();
       hphi_y=(TH1D*)hA->Clone();
       hA_2cbo=(TH1D*)hA->Clone();
       htau_2cbo=(TH1D*)hA->Clone();
       homega_2cbo=(TH1D*)hA->Clone();
       hphi_2cbo=(TH1D*)hA->Clone();
       hn=(TH1D*)hA->Clone();
       hchisq=(TH1D*)hA->Clone();
 
       hphi_0->Reset();
       hA_cbo_N->Reset();
       htau_cbo->Reset();
       homega_cbo->Reset();
       hphi_cbo_N->Reset();
       hA_cbo_A->Reset();
       hphi_cbo_A->Reset();
       hA_cbo_phi->Reset();
       hphi_cbo_phi->Reset();
       hA_vw->Reset();
       htau_vw->Reset();
       homega_vw->Reset();
       hphi_vw->Reset();
       hA_y->Reset();
       htau_y->Reset();
       homega_y->Reset();
       hphi_y->Reset();
       hA_2cbo->Reset();
       htau_2cbo->Reset();
       homega_2cbo->Reset();
       hphi_2cbo->Reset();
       hn->Reset();
       hchisq->Reset();


  for(int i=1;i<=6;i++)
    {	
     h1[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_%d", i)));
     h2[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("blindR_%d", i)));
     h3[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_0_%d", i)));
     h4[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_cbo_N_%d", i)));
     h5[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("tau_cbo_%d", i)));
     h6[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("omega_cbo_%d", i)));
     h7[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_cbo_N_%d", i)));
     h8[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_cbo_A_%d", i)));
     h9[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_cbo_A_%d", i)));
     h10[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_cbo_phi_%d", i)));
     h11[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_cbo_phi_%d", i)));
     h12[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_vw_%d", i)));
     h13[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("tau_vw_%d", i)));
     h14[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("omega_vw_%d", i)));
     h15[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_vw_%d", i)));
     h16[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_y_%d", i)));
     h17[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("tau_y_%d", i)));
     h18[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("omega_y_%d", i)));
     h19[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_y_%d", i)));
     h20[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("A_2cbo_%d", i)));
     h21[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("tau_2cbo_%d", i)));
     h22[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("omega_2cbo_%d", i)));
     h23[i] = (TH1D*)(_file[i]->FindObjectAny(TString::Format("phi_2cbo_%d", i)));
     h24[i]=(TH1D*)(_file[i]->FindObjectAny(TString::Format("n_%d", i)));
    }


  for(int j=0;j<6;j++)
    {
      for(int i=1;i<=hpp;i++)
	{
	 hA->SetBinContent(((j*hpp)+i),h1[j+1]->GetBinContent(i));
	 hblindR->SetBinContent(((j*hpp)+i),h2[j+1]->GetBinContent(i));
	  // cout<<"j "<<j<<" i "<<i<<" j*i "<<((j*3)+i)<<" bincontent "<<hblindR->GetBinContent((j*3)+i)<<endl;
         hphi_0->SetBinContent(((j*hpp)+i),h3[j+1]->GetBinContent(i));
         hA_cbo_N->SetBinContent(((j*hpp)+i),h4[j+1]->GetBinContent(i));
         htau_cbo->SetBinContent(((j*hpp)+i),h5[j+1]->GetBinContent(i));
         homega_cbo->SetBinContent(((j*hpp)+i),h6[j+1]->GetBinContent(i));
         hphi_cbo_N->SetBinContent(((j*hpp)+i),h7[j+1]->GetBinContent(i));
         hA_cbo_A->SetBinContent(((j*hpp)+i),h8[j+1]->GetBinContent(i));
         hphi_cbo_A->SetBinContent(((j*hpp)+i),h9[j+1]->GetBinContent(i));
         hA_cbo_phi->SetBinContent(((j*hpp)+i),h10[j+1]->GetBinContent(i));
         hphi_cbo_phi->SetBinContent(((j*hpp)+i),h11[j+1]->GetBinContent(i));
         hA_vw->SetBinContent(((j*hpp)+i),h12[j+1]->GetBinContent(i));
         htau_vw->SetBinContent(((j*hpp)+i),h13[j+1]->GetBinContent(i));
         homega_vw->SetBinContent(((j*hpp)+i),h14[j+1]->GetBinContent(i));
         hphi_vw->SetBinContent(((j*hpp)+i),h15[j+1]->GetBinContent(i));
         hA_y->SetBinContent(((j*hpp)+i),h16[j+1]->GetBinContent(i));
         htau_y->SetBinContent(((j*hpp)+i),h17[j+1]->GetBinContent(i));
         homega_y->SetBinContent(((j*hpp)+i),h18[j+1]->GetBinContent(i));
         hphi_y->SetBinContent(((j*hpp)+i),h19[j+1]->GetBinContent(i));
         hA_2cbo->SetBinContent(((j*hpp)+i),h20[j+1]->GetBinContent(i));
         htau_2cbo->SetBinContent(((j*hpp)+i),h21[j+1]->GetBinContent(i));
         homega_2cbo->SetBinContent(((j*hpp)+i),h22[j+1]->GetBinContent(i));
         hphi_2cbo->SetBinContent(((j*hpp)+i),h23[j+1]->GetBinContent(i));
         hn->SetBinContent(((j*hpp)+i),h24[j+1]->GetBinContent(i));

         hA->SetBinError(((j*hpp)+i),h1[j+1]->GetBinError(i));
         hblindR->SetBinError(((j*hpp)+i),h2[j+1]->GetBinError(i)) ;
         hphi_0->SetBinError(((j*hpp)+i),h3[j+1]->GetBinError(i));
         hA_cbo_N->SetBinError(((j*hpp)+i),h4[j+1]->GetBinError(i));
         htau_cbo->SetBinError(((j*hpp)+i),h5[j+1]->GetBinError(i));
         homega_cbo->SetBinError(((j*hpp)+i),h6[j+1]->GetBinError(i));
         hphi_cbo_N->SetBinError(((j*hpp)+i),h7[j+1]->GetBinError(i));
         hA_cbo_A->SetBinError(((j*hpp)+i),h8[j+1]->GetBinError(i));
         hphi_cbo_A->SetBinError(((j*hpp)+i),h9[j+1]->GetBinError(i));
         hA_cbo_phi->SetBinError(((j*hpp)+i),h10[j+1]->GetBinError(i));
         hphi_cbo_phi->SetBinError(((j*hpp)+i),h11[j+1]->GetBinError(i));
         hA_vw->SetBinError(((j*hpp)+i),h12[j+1]->GetBinError(i));
         htau_vw->SetBinError(((j*hpp)+i),h13[j+1]->GetBinError(i));
         homega_vw->SetBinError(((j*hpp)+i),h14[j+1]->GetBinError(i));
         hphi_vw->SetBinError(((j*hpp)+i),h15[j+1]->GetBinError(i));
         hA_y->SetBinError(((j*hpp)+i),h16[j+1]->GetBinError(i));
         htau_y->SetBinError(((j*hpp)+i),h17[j+1]->GetBinError(i));
         homega_y->SetBinError(((j*hpp)+i),h18[j+1]->GetBinError(i));
         hphi_y->SetBinError(((j*hpp)+i),h19[j+1]->GetBinError(i));
         hA_2cbo->SetBinError(((j*hpp)+i),h20[j+1]->GetBinError(i));
         htau_2cbo->SetBinError(((j*hpp)+i),h21[j+1]->GetBinError(i));
         homega_2cbo->SetBinError(((j*hpp)+i),h22[j+1]->GetBinError(i));
         hphi_2cbo->SetBinError(((j*hpp)+i),h23[j+1]->GetBinError(i));
         hn->SetBinError(((j*hpp)+i),h24[j+1]->GetBinError(i));
	 
        }
     }
  /*  hchisq->SetBinContent(1,);
  hchisq->SetBinContent(2,);
  hchisq->SetBinContent(3,);
  hchisq->SetBinContent(4,);
  hchisq->SetBinContent(5,);
  hchisq->SetBinContent(6,);
  hchisq->SetBinContent(7,);
  hchisq->SetBinContent(8,);
  hchisq->SetBinContent(9,);
  hchisq->SetBinContent(10,);
  hchisq->SetBinContent(11,);
  hchisq->SetBinContent(12,);
  hchisq->SetBinContent(13,);
  hchisq->SetBinContent(14,);
  hchisq->SetBinContent(15,);
  hchisq->SetBinContent(16,);
  hchisq->SetBinContent(17,);
  hchisq->SetBinContent(18,);
  */
   for(int i=1; i<=18;i++)
       {

       A[m]=hA->GetBinContent(i);
       blindR[m]=hblindR->GetBinContent(i);
       phi_0[m]=hphi_0->GetBinContent(i);
       A_cbo_N[m]=hA_cbo_N->GetBinContent(i);
       tau_cbo[m]=htau_cbo->GetBinContent(i);
       omega_cbo[m]=homega_cbo->GetBinContent(i);
       phi_cbo_N[m]=hphi_cbo_N->GetBinContent(i);
       A_cbo_A[m]=hA_cbo_A->GetBinContent(i);
       phi_cbo_A[m]=hphi_cbo_A->GetBinContent(i);
       A_cbo_phi[m]=hA_cbo_phi->GetBinContent(i);
       phi_cbo_phi[m]=hphi_cbo_phi->GetBinContent(i);
       A_vw[m]=hA_vw->GetBinContent(i);
       tau_vw[m]=htau_vw->GetBinContent(i);
       omega_vw[m]=homega_vw->GetBinContent(i);
       phi_vw[m]=hphi_vw->GetBinContent(i);
       A_y[m]=hA_y->GetBinContent(i);
       tau_y[m]=htau_y->GetBinContent(i);
       omega_y[m]=homega_y->GetBinContent(i);
       phi_y[m]=hphi_y->GetBinContent(i);
       A_2cbo[m]=hA_2cbo->GetBinContent(i);
       tau_2cbo[m]=htau_2cbo->GetBinContent(i);
       omega_2cbo[m]=homega_2cbo->GetBinContent(i);
       phi_2cbo[m]=hphi_2cbo->GetBinContent(i);

       n[m]=hn->GetBinContent(i);

       dA[m]=hA->GetBinError(i);
       dblindR[m]=hblindR->GetBinError(i);
       dphi_0[m]=hphi_0->GetBinError(i);
       dA_cbo_N[m]=hA_cbo_N->GetBinError(i);
       dtau_cbo[m]=htau_cbo->GetBinError(i);
       domega_cbo[m]=homega_cbo->GetBinError(i);
       dphi_cbo_N[m]=hphi_cbo_N->GetBinError(i);
       dA_cbo_A[m]=hA_cbo_A->GetBinError(i);
       dphi_cbo_A[m]=hphi_cbo_A->GetBinError(i);
       dA_cbo_phi[m]=hA_cbo_phi->GetBinError(i);
       dphi_cbo_phi[m]=hphi_cbo_phi->GetBinError(i);
       dA_vw[m]=hA_vw->GetBinError(i);
       dtau_vw[m]=htau_vw->GetBinError(i);
       domega_vw[m]=homega_vw->GetBinError(i);
       dphi_vw[m]=hphi_vw->GetBinError(i);
       dA_y[m]=hA_y->GetBinError(i);
       dtau_y[m]=htau_y->GetBinError(i);
       domega_y[m]=homega_y->GetBinError(i);
       dphi_y[m]=hphi_y->GetBinError(i);
       dA_2cbo[m]=hA_2cbo->GetBinError(i);
       dtau_2cbo[m]=htau_2cbo->GetBinError(i);
       domega_2cbo[m]=homega_2cbo->GetBinError(i);
       dphi_2cbo[m]=hphi_2cbo->GetBinError(i);

       chisq[m]=hchisq->GetBinContent(i);

          if(m==0)
	 {
	   cout<<"m = "<<m<<endl;
	   kband1p[m]=A[m];
	   kband1m[m]=A[m];
	   kband2p[m]=blindR[m];
	   kband2m[m]=blindR[m];
	   kband3p[m]=phi_0[m];
	   kband3m[m]=phi_0[m];
           kband4p[m]=A_cbo_N[m];
	   kband4m[m]=A_cbo_N[m];
           kband5p[m]=tau_cbo[m];
	   kband5m[m]=tau_cbo[m];
           kband6p[m]=omega_cbo[m];
	   kband6m[m]=omega_cbo[m];
           kband7p[m]=phi_cbo_N[m];
	   kband7m[m]=phi_cbo_N[m];
           kband8p[m]=A_cbo_A[m];
	   kband8m[m]=A_cbo_A[m];
           kband9p[m]=phi_cbo_A[m];
	   kband9m[m]=phi_cbo_A[m];
           kband10p[m]=A_cbo_phi[m];
	   kband10m[m]=A_cbo_phi[m];
           kband11p[m]=phi_cbo_phi[m];
	   kband11m[m]=phi_cbo_phi[m];
           kband12p[m]=A_vw[m];
	   kband12m[m]=A_vw[m];
           kband13p[m]=tau_vw[m];
	   kband13m[m]=tau_vw[m];
           kband14p[m]=omega_vw[m];
	   kband14m[m]=omega_vw[m];
           kband15p[m]=phi_vw[m];
	   kband15m[m]=phi_vw[m];
           kband16p[m]=A_y[m];
	   kband16m[m]=A_y[m];
           kband17p[m]=tau_y[m];
	   kband17m[m]=tau_y[m];
           kband18p[m]=omega_y[m];
	   kband18m[m]=omega_y[m];
           kband19p[m]=phi_y[m];
	   kband19m[m]=phi_y[m];
           kband20p[m]=A_2cbo[m];
	   kband20m[m]=A_2cbo[m];
           kband21p[m]=tau_2cbo[m];
	   kband21m[m]=tau_2cbo[m];
           kband22p[m]=omega_2cbo[m];
	   kband22m[m]=omega_2cbo[m];
           kband23p[m]=phi_2cbo[m];
	   kband23m[m]=phi_2cbo[m];

	 }

         else
	 {
	   cout<<"m = "<<m<<endl;
	   kband1p[m]=kband1p[0] + sqrt( TMath::Abs(dA[m]*dA[m] - dA[0]*dA[0]) );
	   kband1m[m]=kband1m[0] - sqrt( TMath::Abs(dA[m]*dA[m] - dA[0]*dA[0]) );
	   kband2p[m]=kband2p[0] + sqrt( TMath::Abs(dblindR[m]*dblindR[m] - dblindR[0]*dblindR[0]) );
					 kband2m[m]=kband2m[0] - sqrt( TMath::Abs(dblindR[m]*dblindR[m] - dblindR[0]*dblindR[0]) );
					 kband3p[m]=kband3p[0] + sqrt( TMath::Abs(dphi_0[m]*dphi_0[m] - dphi_0[0]*dphi_0[0]) );
					 kband3m[m]=kband3m[0] - sqrt( TMath::Abs(dphi_0[m]*dphi_0[m] - dphi_0[0]*dphi_0[0]) );
					 kband4p[m]=kband4p[0] + sqrt( TMath::Abs(dA_cbo_N[m]*dA_cbo_N[m] - dA_cbo_N[0]*dA_cbo_N[0]) );
					 kband4m[m]=kband4m[0] - sqrt( TMath::Abs(dA_cbo_N[m]*dA_cbo_N[m] - dA_cbo_N[0]*dA_cbo_N[0]) );
					 kband5p[m]=kband5p[0] + sqrt( TMath::Abs(dtau_cbo[m]*dtau_cbo[m] - dtau_cbo[0]*dtau_cbo[0]) );
					 kband5m[m]=kband5m[0] - sqrt( TMath::Abs(dtau_cbo[m]*dtau_cbo[m] - dtau_cbo[0]*dtau_cbo[0]) );
					 kband6p[m]=kband6p[0] + sqrt( TMath::Abs(domega_cbo[m]*domega_cbo[m] - domega_cbo[0]*domega_cbo[0]) );
					 kband6m[m]=kband6m[0] - sqrt( TMath::Abs(domega_cbo[m]*domega_cbo[m] - domega_cbo[0]*domega_cbo[0]) );
					 kband7p[m]=kband7p[0] + sqrt( TMath::Abs(dphi_cbo_N[m]*dphi_cbo_N[m] - dphi_cbo_N[0]*dphi_cbo_N[0]) );
					 kband7m[m]=kband7m[0] - sqrt( TMath::Abs(dphi_cbo_N[m]*dphi_cbo_N[m] - dphi_cbo_N[0]*dphi_cbo_N[0]) );
					 kband8p[m]=kband8p[0] + sqrt( TMath::Abs(dA_cbo_A[m]*dA_cbo_A[m] - dA_cbo_A[0]*dA_cbo_A[0]) );
					 kband8m[m]=kband8m[0] - sqrt( TMath::Abs(dA_cbo_A[m]*dA_cbo_A[m] - dA_cbo_A[0]*dA_cbo_A[0]) );
					 kband9p[m]=kband9p[0] + sqrt( TMath::Abs(dphi_cbo_A[m]*dphi_cbo_A[m] - dphi_cbo_A[0]*dphi_cbo_A[0]) );
					 kband9m[m]=kband9m[0] - sqrt( TMath::Abs(dphi_cbo_A[m]*dphi_cbo_A[m] - dphi_cbo_A[0]*dphi_cbo_A[0]) );
					 kband10p[m]=kband10p[0] + sqrt( TMath::Abs(dA_cbo_phi[m]*dA_cbo_phi[m] - dA_cbo_phi[0]*dA_cbo_phi[0]) );
					   kband10m[m]=kband10m[0] - sqrt( TMath::Abs(dA_cbo_phi[m]*dA_cbo_phi[m] - dA_cbo_phi[0]*dA_cbo_phi[0]) );
					   kband11p[m]=kband11p[0] + sqrt( TMath::Abs(dphi_cbo_phi[m]*dphi_cbo_phi[m] - dphi_cbo_phi[0]*dphi_cbo_phi[0]) );
					   kband11m[m]=kband11m[0] - sqrt( TMath::Abs(dphi_cbo_phi[m]*dphi_cbo_phi[m] - dphi_cbo_phi[0]*dphi_cbo_phi[0]) );
					   kband12p[m]=kband12p[0] + sqrt( TMath::Abs(dA_vw[m]*dA_vw[m] - dA_vw[0]*dA_vw[0]) );
					   kband12m[m]=kband12m[0] - sqrt( TMath::Abs(dA_vw[m]*dA_vw[m] - dA_vw[0]*dA_vw[0]) );
					   kband13p[m]=kband13p[0] + sqrt( TMath::Abs(dtau_vw[m]*dtau_vw[m] - dtau_vw[0]*dtau_vw[0]) );
					   kband13m[m]=kband13m[0] - sqrt( TMath::Abs(dtau_vw[m]*dtau_vw[m] - dtau_vw[0]*dtau_vw[0]) );
					   kband14p[m]=kband14p[0] + sqrt( TMath::Abs(domega_vw[m]*domega_vw[m] - domega_vw[0]*domega_vw[0]) );
					   kband14m[m]=kband14m[0] - sqrt( TMath::Abs(domega_vw[m]*domega_vw[m] - domega_vw[0]*domega_vw[0]) );
					   kband15p[m]=kband15p[0] + sqrt( TMath::Abs(dphi_vw[m]*dphi_vw[m] - dphi_vw[0]*dphi_vw[0]) );
					   kband15m[m]=kband15m[0] - sqrt( TMath::Abs(dphi_vw[m]*dphi_vw[m] - dphi_vw[0]*dphi_vw[0]) );
					   kband16p[m]=kband16p[0] + sqrt( TMath::Abs(dA_y[m]*dA_y[m] - dA_y[0]*dA_y[0]) );
					   kband16m[m]=kband16m[0] - sqrt( TMath::Abs(dA_y[m]*dA_y[m] - dA_y[0]*dA_y[0]) );
					   kband17p[m]=kband17p[0] + sqrt( TMath::Abs(dtau_y[m]*dtau_y[m] - dtau_y[0]*dtau_y[0]) );
					   kband17m[m]=kband17m[0] - sqrt( TMath::Abs(dtau_y[m]*dtau_y[m] - dtau_y[0]*dtau_y[0]) );
					   kband18p[m]=kband18p[0] + sqrt( TMath::Abs(domega_y[m]*domega_y[m] - domega_y[0]*domega_y[0]) );
					   kband18m[m]=kband18m[0] - sqrt( TMath::Abs(domega_y[m]*domega_y[m] - domega_y[0]*domega_y[0]) );
					   kband19p[m]=kband19p[0] + sqrt( TMath::Abs(dphi_y[m]*dphi_y[m] - dphi_y[0]*dphi_y[0]) );
					   kband19m[m]=kband19m[0] - sqrt( TMath::Abs(dphi_y[m]*dphi_y[m] - dphi_y[0]*dphi_y[0]) );
					   kband20p[m]=kband20p[0] + sqrt( TMath::Abs(dA_2cbo[m]*dA_2cbo[m] - dA_2cbo[0]*dA_2cbo[0]) );
					   kband20m[m]=kband20m[0] - sqrt( TMath::Abs(dA_2cbo[m]*dA_2cbo[m] - dA_2cbo[0]*dA_2cbo[0]) );
					   kband21p[m]=kband21p[0] + sqrt( TMath::Abs(dtau_2cbo[m]*dtau_2cbo[m] - dtau_2cbo[0]*dtau_2cbo[0]) );
					   kband21m[m]=kband21m[0] - sqrt( TMath::Abs(dtau_2cbo[m]*dtau_2cbo[m] - dtau_2cbo[0]*dtau_2cbo[0]) );
					   kband22p[m]=kband22p[0] + sqrt( TMath::Abs(domega_2cbo[m]*domega_2cbo[m] - domega_2cbo[0]*domega_2cbo[0]) );
					   kband22m[m]=kband22m[0] - sqrt( TMath::Abs(domega_2cbo[m]*domega_2cbo[m] - domega_2cbo[0]*domega_2cbo[0]) );
					   kband23p[m]=kband23p[0] + sqrt( TMath::Abs(dphi_2cbo[m]*dphi_2cbo[m] - dphi_2cbo[0]*dphi_2cbo[0]) );
					   kband23m[m]=kband23m[0] - sqrt( TMath::Abs(dphi_2cbo[m]*dphi_2cbo[m] - dphi_2cbo[0]*dphi_2cbo[0]) );

	 }

	  //n[m]=;
       m=m+1;
     }


      c1=new TCanvas("c1","blind R vs fit start time");
      // c1->Divide(1,2);
      //c1->cd(1);
    /*   gr1=new TGraphErrors(m,n,chisq,0,0);
    gr1->SetTitle("chisq vs fit start time");
    gr1->GetXaxis()->SetTitle("Start time [ns]");
    gr1->GetYaxis()->SetTitle("reduced chi-sq");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw();
    */

      //c1->cd(1);
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle("A vs fit start time");
    gr2->GetXaxis()->SetTitle("Start time [ns]");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw();
    gkband1p=new TGraph(m,n,kband1p);
    gkband1m=new TGraph(m,n,kband1m);
    gkband1m->SetLineColor(kBlue);
    gkband1p->SetLineColor(kBlue);
    gkband1m->SetLineWidth(2);
    gkband1p->SetLineWidth(2);
    gkband1p->Draw("same");
    gkband1m->Draw("same");

    c1a=new TCanvas("c1a","blind R vs fit start time");     
    gr3=new TGraphErrors(m,n,blindR,0,dblindR);
    gr3->SetTitle("blindR vs fit start time");
    gr3->GetXaxis()->SetTitle("Start time [ns]");
    gr3->GetYaxis()->SetTitle("Blind R [ppm]");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw();
    gkband2p=new TGraph(m,n,kband2p);
    gkband2m=new TGraph(m,n,kband2m);
    gkband2m->SetLineColor(kBlue);
    gkband2p->SetLineColor(kBlue);
    gkband2m->SetLineWidth(2);
    gkband2p->SetLineWidth(2);
    gkband2p->Draw("same");
    gkband2m->Draw("same");


     c1b=new TCanvas("c1b","blind R vs fit start time");
    gr4=new TGraphErrors(m,n,phi_0,0,dphi_0);
    gr4->SetTitle("phase vs fit start time");
    gr4->GetXaxis()->SetTitle("Start time [ns]");
    gr4->GetYaxis()->SetTitle("phi_0");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->Draw();
    gkband3p=new TGraph(m,n,kband3p);
    gkband3m=new TGraph(m,n,kband3m);
    gkband3m->SetLineColor(kBlue);
    gkband3p->SetLineColor(kBlue);
    gkband3m->SetLineWidth(2);
    gkband3p->SetLineWidth(2);
    gkband3p->Draw("same");
    gkband3m->Draw("same");


    c4=new TCanvas("c4","cbo vs fit start time");
    c4->Divide(2,2);
    c4->cd(1);
    gr8=new TGraphErrors(m,n,A_cbo_N,0,dA_cbo_N);
    gr8->SetTitle("A_cbo_N vs fit start time");
    gr8->GetXaxis()->SetTitle("Start time [ns]");
    gr8->GetYaxis()->SetTitle("A_cbo_N");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kRed);
    gr8->Draw();
    gkband4p=new TGraph(m,n,kband4p);
    gkband4m=new TGraph(m,n,kband4m);
    gkband4m->SetLineColor(kBlue);
    gkband4p->SetLineColor(kBlue);
    gkband4m->SetLineWidth(2);
    gkband4p->SetLineWidth(2);
    gkband4p->Draw("same");
    gkband4m->Draw("same");


     c4->cd(2);
    gr5=new TGraphErrors(m,n,tau_cbo,0,dtau_cbo);
    gr5->SetTitle("tau_cbo_ vs fit start time");
    gr5->GetXaxis()->SetTitle("Start time [ns]");
    gr5->GetYaxis()->SetTitle("tau_cbo");
    gr5->SetMarkerStyle(20);
    gr5->SetLineColor(kRed);
    gr5->Draw();
    gkband5p=new TGraph(m,n,kband5p);
    gkband5m=new TGraph(m,n,kband5m);
    gkband5m->SetLineColor(kBlue);
    gkband5p->SetLineColor(kBlue);
    gkband5m->SetLineWidth(2);
    gkband5p->SetLineWidth(2);
    gkband5p->Draw("same");
    gkband5m->Draw("same");


     c4->cd(3);
    gr6=new TGraphErrors(m,n,omega_cbo,0,domega_cbo);
    gr6->SetTitle("omega_cbo vs fit start time");
    gr6->GetXaxis()->SetTitle("Start time [ns]");
    gr6->GetYaxis()->SetTitle("omega_cbo");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kRed);
    gr6->Draw();
    gkband6p=new TGraph(m,n,kband6p);
    gkband6m=new TGraph(m,n,kband6m);
    gkband6m->SetLineColor(kBlue);
    gkband6p->SetLineColor(kBlue);
    gkband6m->SetLineWidth(2);
    gkband6p->SetLineWidth(2);
    gkband6p->Draw("same");
    gkband6m->Draw("same");


     c4->cd(4);
    gr7=new TGraphErrors(m,n,phi_cbo_N,0,dphi_cbo_N);
    gr7->SetTitle("phi_cbo_N vs fit start time");
    gr7->GetXaxis()->SetTitle("Start time [ns]");
    gr7->GetYaxis()->SetTitle("phi_cbo_N");
    gr7->SetMarkerStyle(20);
    gr7->SetLineColor(kRed);
    gr7->Draw();
    gkband7p=new TGraph(m,n,kband7p);
    gkband7m=new TGraph(m,n,kband7m);
    gkband7m->SetLineColor(kBlue);
    gkband7p->SetLineColor(kBlue);
    gkband7m->SetLineWidth(2);
    gkband7p->SetLineWidth(2);
    gkband7p->Draw("same");
    gkband7m->Draw("same");


      c6=new TCanvas("c6","A,phi,cbo vs fit start time");
    c6->Divide(2,2);
    c6->cd(1);
    gr13=new TGraphErrors(m,n,A_cbo_A,0,dA_cbo_A);
    gr13->SetTitle("A_cbo_A vs fit start time");
    gr13->GetXaxis()->SetTitle("Start time [ns]");
    gr13->GetYaxis()->SetTitle("A_cbo_A");
    gr13->SetMarkerStyle(20);
    gr13->SetLineColor(kRed);
    gr13->Draw();
    gkband8p=new TGraph(m,n,kband8p);
    gkband8m=new TGraph(m,n,kband8m);
    gkband8m->SetLineColor(kBlue);
    gkband8p->SetLineColor(kBlue);
    gkband8m->SetLineWidth(2);
    gkband8p->SetLineWidth(2);
    gkband8p->Draw("same");
    gkband8m->Draw("same");

     c6->cd(2);
    gr14=new TGraphErrors(m,n,phi_cbo_A,0,dphi_cbo_A);
    gr14->SetTitle("phi_cbo_A vs fit start time");
    gr14->GetXaxis()->SetTitle("Start time [ns]");
    gr14->GetYaxis()->SetTitle("phi_cbo_A");
    gr14->SetMarkerStyle(20);
    gr14->SetLineColor(kRed);
    gr14->Draw();
    gkband9p=new TGraph(m,n,kband9p);
    gkband9m=new TGraph(m,n,kband9m);
    gkband9m->SetLineColor(kBlue);
    gkband9p->SetLineColor(kBlue);
    gkband9m->SetLineWidth(2);
    gkband9p->SetLineWidth(2);
    gkband9p->Draw("same");
    gkband9m->Draw("same");

     c6->cd(3);
    gr15=new TGraphErrors(m,n,A_cbo_phi,0,dA_cbo_phi);
    gr15->SetTitle("A_cbo_phi vs fit start time");
    gr15->GetXaxis()->SetTitle("Start time [ns]");
    gr15->GetYaxis()->SetTitle("A_cbo_phi");
    gr15->SetMarkerStyle(20);
    gr15->SetLineColor(kRed);
    gr15->Draw();
    gkband10p=new TGraph(m,n,kband10p);
    gkband10m=new TGraph(m,n,kband10m);
    gkband10m->SetLineColor(kBlue);
    gkband10p->SetLineColor(kBlue);
    gkband10m->SetLineWidth(2);
    gkband10p->SetLineWidth(2);
    gkband10p->Draw("same");
    gkband10m->Draw("same");

     c6->cd(4);
    gr16=new TGraphErrors(m,n,phi_cbo_phi,0,dphi_cbo_phi);
    gr16->SetTitle("phi_cbo_phi vs fit start time");
    gr16->GetXaxis()->SetTitle("Start time [ns]");
    gr16->GetYaxis()->SetTitle("phi_cbo_phi");
    gr16->SetMarkerStyle(20);
    gr16->SetLineColor(kRed);
    gr16->Draw();
    gkband11p=new TGraph(m,n,kband11p);
    gkband11m=new TGraph(m,n,kband11m);
    gkband11m->SetLineColor(kBlue);
    gkband11p->SetLineColor(kBlue);
    gkband11m->SetLineWidth(2);
    gkband11p->SetLineWidth(2);
    gkband11p->Draw("same");
    gkband11m->Draw("same");


     c7=new TCanvas("c7","2cbo vs fit start time");
     c7->Divide(2,2);
    c7->cd(1);
    gr17=new TGraphErrors(m,n,A_2cbo,0,dA_2cbo);
    gr17->SetTitle("A_2cbo vs fit start time");
    gr17->GetXaxis()->SetTitle("Start time [ns]");
    gr17->GetYaxis()->SetTitle("A_2cbo");
    gr17->SetMarkerStyle(20);
    gr17->SetLineColor(kRed);
    gr17->Draw();
    gkband20p=new TGraph(m,n,kband20p);
    gkband20m=new TGraph(m,n,kband20m);
    gkband20m->SetLineColor(kBlue);
    gkband20p->SetLineColor(kBlue);
    gkband20m->SetLineWidth(2);
    gkband20p->SetLineWidth(2);
    gkband20p->Draw("same");
    gkband20m->Draw("same");


     c7->cd(2);
    gr18=new TGraphErrors(m,n,tau_2cbo,0,dtau_2cbo);
    gr18->SetTitle("tau_2cbo_ vs fit start time");
    gr18->GetXaxis()->SetTitle("Start time [ns]");
    gr18->GetYaxis()->SetTitle("tau_2cbo");
    gr18->SetMarkerStyle(20);
    gr18->SetLineColor(kRed);
    gr18->Draw();
    gkband21p=new TGraph(m,n,kband21p);
    gkband21m=new TGraph(m,n,kband21m);
    gkband21m->SetLineColor(kBlue);
    gkband21p->SetLineColor(kBlue);
    gkband21m->SetLineWidth(2);
    gkband21p->SetLineWidth(2);
    gkband21p->Draw("same");
    gkband21m->Draw("same");


     c7->cd(3);
    gr19=new TGraphErrors(m,n,omega_2cbo,0,domega_2cbo);
    gr19->SetTitle("omega_2cbo vs fit start time");
    gr19->GetXaxis()->SetTitle("Start time [ns]");
    gr19->GetYaxis()->SetTitle("omega_2cbo");
    gr19->SetMarkerStyle(20);
    gr19->SetLineColor(kRed);
    gr19->Draw();
    gkband22p=new TGraph(m,n,kband22p);
    gkband22m=new TGraph(m,n,kband22m);
    gkband22m->SetLineColor(kBlue);
    gkband22p->SetLineColor(kBlue);
    gkband22m->SetLineWidth(2);
    gkband22p->SetLineWidth(2);
    gkband22p->Draw("same");
    gkband22m->Draw("same");


     c7->cd(4);
    gr20=new TGraphErrors(m,n,phi_2cbo,0,dphi_2cbo);
    gr20->SetTitle("phi_2cbo vs fit start time");
    gr20->GetXaxis()->SetTitle("Start time [ns]");
    gr20->GetYaxis()->SetTitle("phi_2cbo");
    gr20->SetMarkerStyle(20);
    gr20->SetLineColor(kRed);
    gr20->Draw();
    gkband23p=new TGraph(m,n,kband23p);
    gkband23m=new TGraph(m,n,kband23m);
    gkband23m->SetLineColor(kBlue);
    gkband23p->SetLineColor(kBlue);
    gkband23m->SetLineWidth(2);
    gkband23p->SetLineWidth(2);
    gkband23p->Draw("same");
    gkband23m->Draw("same");


    c8=new TCanvas("c8","vw vs fit start time");
    c8->Divide(2,2);
    c8->cd(1);
    gr21=new TGraphErrors(m,n,A_vw,0,dA_vw);
    gr21->SetTitle("A_vw vs fit start time");
    gr21->GetXaxis()->SetTitle("Start time [ns]");
    gr21->GetYaxis()->SetTitle("A_vw");
    gr21->SetMarkerStyle(20);
    gr21->SetLineColor(kRed);
    gr21->Draw();
    gkband12p=new TGraph(m,n,kband12p);
    gkband12m=new TGraph(m,n,kband12m);
    gkband12m->SetLineColor(kBlue);
    gkband12p->SetLineColor(kBlue);
    gkband12m->SetLineWidth(2);
    gkband12p->SetLineWidth(2);
    gkband12p->Draw("same");
    gkband12m->Draw("same");


    c8->cd(2);
    gr22=new TGraphErrors(m,n,tau_vw,0,dtau_vw);
    gr22->SetTitle("tau_vw vs fit start time");
    gr22->GetXaxis()->SetTitle("Start time [ns]");
    gr22->GetYaxis()->SetTitle("tau_vw");
    gr22->SetMarkerStyle(20);
    gr22->SetLineColor(kRed);
    gr22->Draw();
    gkband13p=new TGraph(m,n,kband13p);
    gkband13m=new TGraph(m,n,kband13m);
    gkband13m->SetLineColor(kBlue);
    gkband13p->SetLineColor(kBlue);
    gkband13m->SetLineWidth(2);
    gkband13p->SetLineWidth(2);
    gkband13p->Draw("same");
    gkband13m->Draw("same");


    c8->cd(3);
    gr23=new TGraphErrors(m,n,omega_vw,0,domega_vw);
    gr23->SetTitle("omega_vw vs fit start time");
    gr23->GetXaxis()->SetTitle("Start time [ns]");
    gr23->GetYaxis()->SetTitle("omega_vw");
    gr23->SetMarkerStyle(20);
    gr23->SetLineColor(kRed);
    gr23->Draw();
    gkband14p=new TGraph(m,n,kband14p);
    gkband14m=new TGraph(m,n,kband14m);
    gkband14m->SetLineColor(kBlue);
    gkband14p->SetLineColor(kBlue);
    gkband14m->SetLineWidth(2);
    gkband14p->SetLineWidth(2);
    gkband14p->Draw("same");
    gkband14m->Draw("same");


    c8->cd(4);
    gr24=new TGraphErrors(m,n,phi_vw,0,dphi_vw);
    gr24->SetTitle("phi_vw vs fit start time");
    gr24->GetXaxis()->SetTitle("Start time [ns]");
    gr24->GetYaxis()->SetTitle("phi_vw");
    gr24->SetMarkerStyle(20);
    gr24->SetLineColor(kRed);
    gr24->Draw();
    gkband15p=new TGraph(m,n,kband15p);
    gkband15m=new TGraph(m,n,kband15m);
    gkband15m->SetLineColor(kBlue);
    gkband15p->SetLineColor(kBlue);
    gkband15m->SetLineWidth(2);
    gkband15p->SetLineWidth(2);
    gkband15p->Draw("same");
    gkband15m->Draw("same");


    c9=new TCanvas("c9","y vs fit start time");
    c9->Divide(2,2);
    c9->cd(1);
    gr25=new TGraphErrors(m,n,A_y,0,dA_y);
    gr25->SetTitle("A_y vs fit start time");
    gr25->GetXaxis()->SetTitle("Start time [ns]");
    gr25->GetYaxis()->SetTitle("A_y");
    gr25->SetMarkerStyle(20);
    gr25->SetLineColor(kRed);
    gr25->Draw();
    gkband16p=new TGraph(m,n,kband16p);
    gkband16m=new TGraph(m,n,kband16m);
    gkband16m->SetLineColor(kBlue);
    gkband16p->SetLineColor(kBlue);
    gkband16m->SetLineWidth(2);
    gkband16p->SetLineWidth(2);
    gkband16p->Draw("same");
    gkband16m->Draw("same");


    c9->cd(2);
    gr26=new TGraphErrors(m,n,tau_y,0,dtau_y);
    gr26->SetTitle("tau_y vs fit start time");
    gr26->GetXaxis()->SetTitle("Start time [ns]");
    gr26->GetYaxis()->SetTitle("tau_y");
    gr26->SetMarkerStyle(20);
    gr26->SetLineColor(kRed);
    gr26->Draw();
    gkband17p=new TGraph(m,n,kband17p);
    gkband17m=new TGraph(m,n,kband17m);
    gkband17m->SetLineColor(kBlue);
    gkband17p->SetLineColor(kBlue);
    gkband17m->SetLineWidth(2);
    gkband17p->SetLineWidth(2);
    gkband17p->Draw("same");
    gkband17m->Draw("same");


    c9->cd(3);
    gr27=new TGraphErrors(m,n,omega_y,0,domega_y);
    gr27->SetTitle("omega_y vs fit start time");
    gr27->GetXaxis()->SetTitle("Start time [ns]");
    gr27->GetYaxis()->SetTitle("omega_y");
    gr27->SetMarkerStyle(20);
    gr27->SetLineColor(kRed);
    gr27->Draw();
    gkband18p=new TGraph(m,n,kband18p);
    gkband18m=new TGraph(m,n,kband18m);
    gkband18m->SetLineColor(kBlue);
    gkband18p->SetLineColor(kBlue);
    gkband18m->SetLineWidth(2);
    gkband18p->SetLineWidth(2);
    gkband18p->Draw("same");
    gkband18m->Draw("same");


    c9->cd(4);
    gr28=new TGraphErrors(m,n,phi_y,0,dphi_y);
    gr28->SetTitle("phi_y vs fit start time");
    gr28->GetXaxis()->SetTitle("Start time [ns]");
    gr28->GetYaxis()->SetTitle("phi_y");
    gr28->SetMarkerStyle(20);
    gr28->SetLineColor(kRed);
    gr28->Draw();
    gkband19p=new TGraph(m,n,kband19p);
    gkband19m=new TGraph(m,n,kband19m);
    gkband19m->SetLineColor(kBlue);
    gkband19p->SetLineColor(kBlue);
    gkband19m->SetLineWidth(2);
    gkband19p->SetLineWidth(2);
    gkband19p->Draw("same");
    gkband19m->Draw("same");


   
}
