#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28,*gr[50];
TCanvas *c1;

Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi[50], dphi[50], A_cbo_N[50][3], tau_cbo[50][3], omega_cbo[50][3], phi_cbo_N[50][3], A_vw[50][3], tau_vw[50][3], omega_vw[50][3], phi_vw[50][3], A_y[50][3], tau_y[50][3], omega_y[50][3], phi_y[50][3],  dA_cbo_N[50][3], dtau_cbo[50][3], domega_cbo[50][3], dphi_cbo_N[50][3], dA_vw[50][3], dtau_vw[50][3], domega_vw[50][3], dphi_vw[50][3], dA_y[50][3], dtau_y[50][3], domega_y[50][3], dphi_y[50][3], n[50][3];
Double_t run[50][3][15],drun[50][3][15];
Double_t chisq2[50],A2[50],dA2[50],  blindR2[50], dblindR2[50], n1[10],n2[10],n3[10], phi2[50], dphi2[50],chisq3[50],A3[50],dA3[50],  blindR3[50], dblindR3[50], phi3[50], dphi3[50];
Int_t i,j,k;
Int_t m1=0;
Int_t m2=0;
Int_t m3=0;

void benchmark_rm()
{
  double r2_random[15]={ 0.22891822, -44.642105, 4.0108341, 0.0024277649, 256475.25, 0.0023403797, 5.7060738, 0.00022390353, 140446.26, 0.013929809, 5.7741761, 0.0011018196, 35039.870, 0.014038975, 4.3928759 };
  
  double dr2_random[15]= { 9.6696120e-06, 0.52891616, 8.1230354e-05, 6.8517095e-05, 22054.713, 3.3346335e-07, 0.028051768, 4.8663040e-05, 53805.241, 3.0307469e-06, 0.24807596, 0.00032825290, 7293.3331, 4.5319703e-06, 0.22513096 };
  
  double r3btom_random[15]={ 0.23677772, -36.553136, 4.1533796, 0.0019333776, 233240.61, 0.0024037604, 5.9967015, 0.00046285354, 91436.332, 0.014340165, 6.2050396, 0.0032506165, 18043.010, 0.014549783, 3.0436976 };
    
  double dr3btom_random[15]={ 8.0993655e-06, 0.44284751, 6.7849196e-05, 5.8960010e-05, 19687.634, 3.7755463e-07, 0.030829530, 5.1776525e-05, 13317.136, 1.8431945e-06, 0.12375066, 0.0011199286, 2653.9067, 8.7842378e-06, 0.34773253 };

  double r3no_random[15]={ 0.22822197, -33.730830, 3.9974838, 0.0010987376, 245097.42, 0.0023312571, 4.2880445, 0.00058900949, 78806.082, 0.013906808, 5.7031054, 0.0016116001, 31899.538, 0.014042792, 2.1372778 };

  double dr3no_random[15]={ 1.2571352e-05, 0.69079041, 0.00010592169, 9.2586698e-05, 63562.826, 9.7670086e-07, 0.080633209, 0.00012441364, 19097.623, 3.0261655e-06, 0.19622486, 0.00057445496, 6994.6425, 6.1651465e-06, 0.29049156 };

  // ##### my rootfile
  /*  double r2_copy[15]={2.28896703e-01, -4.45657131e+01,  4.01447359e+00,  2.42059790e-03, 2.56445817e+05,  2.34052660e-03, -5.69467481e-01,  2.00065372e-04, 1.66057038e+05,  1.39296003e-02, -4.09496151e-01,  1.51663366e-03, 2.72848136e+04,  1.40381095e-02, -2.08118691e+01};

  double dr2_copy[15]={9.07402041e-06, 5.05812149e-01, 7.68872597e-05, 2.56468788e-05, 8.38848524e+03, 1.27728821e-07, 1.06185843e-02,3.32090764e-05, 5.37309810e+04, 1.95250023e-06, 1.66557619e-01, 2.33626374e-04, 2.60628583e+03, 3.50155308e-06, 1.54279369e-01};

  double r3btom_copy[15]={2.29348767e-01, -3.53960656e+01,  4.02301876e+00,  1.87084739e-03, 2.28222025e+05,  2.32891792e-03, -4.57698457e-01,  4.63898575e-04, 8.57981276e+04,  1.38925126e-02, -2.33276242e-01, -3.89195384e-03, 1.59145941e+04,  1.41000360e-02, -1.88021448e+01};

    double dr3btom_copy[15]={7.15870673e-06, 3.99048568e-01, 6.04882574e-05, 2.09821029e-05, 7.25083722e+03, 1.39997570e-07, 1.13306404e-02, 3.66164284e-05, 8.88817222e+03, 1.20371850e-06, 7.85104086e-02, 6.19017796e-04, 1.07877204e+03, 4.24337948e-06, 1.58323515e-01};

    double r3no_copy[15]={2.28454223e-01, -3.38584842e+01,  3.99678108e+00, -1.33525138e-03, 1.89168098e+05,  2.33129609e-03,  2.51527144e+00, -6.78155275e-04, 6.07425492e+04,  1.38990451e-02,  2.59281192e+00,  2.29604425e-03, 2.37142513e+04,  1.40463433e-02,  1.82530450e+00};

    double dr3no_copy[15]={1.14829732e-05, 6.42560246e-01, 9.74514968e-05, 3.64695300e-05, 1.26819748e+04, 3.55818868e-07, 2.75090088e-02, 9.41849059e-05, 8.46214225e+03, 2.28540490e-06, 1.38052838e-01, 4.08868366e-04, 2.34156457e+03, 4.13193151e-06, 1.76173661e-01};
  */


  // ##### Tim's rootfile
  double r2_copy[15]={2.28918374e-01, -4.47917664e+01,  4.01080398e+00,  2.40902556e-03, 2.62328852e+05,  2.34043548e-03, -5.73836517e-01,  1.90988354e-04,  1.56763402e+05,  1.39287127e-02, -5.36731354e-01,  1.43956492e-03, 2.81320790e+04,  1.40384360e-02, -1.44600867e+01};

  double dr2_copy[15]={9.09564627e-06, 5.07062742e-01, 7.70692459e-05, 2.55110957e-05, 8.73240508e+03, 1.27061250e-07, 1.06122383e-02, 3.41236950e-05, 5.25293323e+04, 2.14221959e-06, 1.79337917e-01, 2.20398747e-04, 2.72494046e+03, 3.44024202e-06, 1.53089640e-01};

  double r3btom_copy[15]={2.29379339e-01, -3.54059657e+01,  4.02358794e+00,  1.87365824e-03, 2.28553814e+05,  2.32873805e-03, -4.68137877e-01,  4.91361824e-04, 7.96875997e+04,  1.38918889e-02, -2.81524137e-01, -4.88593996e-03, 1.45618986e+04,  1.40998745e-02, -1.88370944e+01};

    double dr3btom_copy[15]={7.75357558e-06, 4.32011008e-01, 6.55030973e-05, 2.27160561e-05, 7.85850213e+03, 1.51291030e-07, 1.22477680e-02, 4.32361535e-05, 8.70556717e+03, 1.36625190e-06, 8.74770213e-02, 8.85679537e-04, 1.04795603e+03, 4.92475655e-06, 1.80540779e-01};

    double r3no_copy[15]={2.28222023e-01, -3.37529289e+01,  3.99747592e+00, -1.10788202e-03, 2.27163073e+05,  2.33129710e-03,  2.59577168e+00, -6.79749506e-04, 6.35771822e+04,  1.39007513e-02,  2.72785752e+00,  2.79100978e-03, 2.15810702e+04,  1.40488743e-02,  1.97766521e+00};

    double dr3no_copy[15]={1.21167189e-05, 6.78402336e-01, 1.02922793e-04, 3.55649287e-05, 2.05576311e+04, 4.00223910e-07, 3.23659670e-02, 9.53488801e-05, 9.16886285e+03, 2.26079024e-06, 1.39469752e-01, 5.35043816e-04, 2.14797542e+03, 4.57235199e-06, 1.89420900e-01};

    double r2_q[15]= {0.22892070, -44.752954, 4.0108119, 0.0024142786, 261456.72, 0.0023404561, -0.57193519, 0.0011495125, 32329.863, 0.014034116, 4.1152916, 0.00020461523, 144178.31, 0.013931995, -0.19715734};

    double dr2_q[15]={7.8494464e-06, 0.43894356, 6.7315334e-05, 2.2342062e-05, 7453.9381, 1.0901652e-07, 0.0092113335, 0.00018549282, 3506.3478, 2.3736321e-06, 0.11656488, 2.9068693e-05, 34500.897, 1.9184774e-06, 0.16048126};

    double r3btom_q[15]={0.22938337, -35.548603, 4.0235557, 0.0018710910, 229273.85, 0.0023286677, -0.47613009, 0.0050624326, 14468.270, 0.014091322, 2.7784940, 0.00044674013, 87287.137, 0.013891308, -0.33461725};

    double dr3btom_q[15]={6.7539053e-06, 0.37742592, 5.7740238e-05, 2.0278116e-05, 6931.0350, 1.2961392e-07, 0.010633471, 0.00089389755, 967.39680, 3.9181017e-06, 0.15056707, 3.1629517e-05, 7821.1023, 1.0804943e-06, 0.072812543};

    double r3no_q[15]={0.22822855, -33.996651, 3.9974341, 0.0011269970, 218667.07, 0.0023312317, -0.55255669, 0.0022263735, 24121.815, 0.014042373, 1.6064561, 0.00059034351, 70375.388, 0.013898790, -0.58847719};

    double dr3no_q[15]={9.8988522e-06, 0.55592007, 8.5123673e-05, 3.0887076e-05, 16168.346, 3.2015734e-07, 0.026118342, 0.00047167157, 2799.6418, 3.6228934e-06, 0.16352736, 8.4567609e-05, 10930.260, 2.0322440e-06, 0.13190071};
  

  
    m1=1;
    m2=1;
    m3=1;

    // Get the parameters for each dataset
    for(int i=1;i<=3;i++)
      {
	if(m1==1)
	  {
	   chisq[m1]=1.04;
	   A[m1]=r2_random[0];
	   dA[m1]=dr2_random[0];
	   blindR[m1]=r2_random[1];
	   dblindR[m1]=dr2_random[1];
	   phi[m1]=r2_random[2];
	   dphi[m1]=dr2_random[2];
	   n1[m1]=m1;
	  }
	if(m1==2)
	  {
	    chisq[m1]=1.07;
	    A[m1]=r3btom_random[0];
	    dA[m1]=dr3btom_random[0];
	    blindR[m1]=r3btom_random[1];
	    dblindR[m1]=dr3btom_random[1];
	    phi[m1]=r3btom_random[2];
	    dphi[m1]=dr3btom_random[2];
	    n1[m1]=m1;
	  }
	if(m1==3)
	  {
	    chisq[m1]=1.02;
	    A[m1]=r3no_random[0];
	    dA[m1]=dr3no_random[0];
	    blindR[m1]=r3no_random[1];
	    dblindR[m1]=dr3no_random[1];
	    phi[m1]=r3no_random[2];
	    dphi[m1]=dr3no_random[2];
	    n1[m1]=m1;
	  }
	m1=m1+1;
      }
    
    for(int i=1;i<=3;i++)
      {
	if(m2==1)
	  {
	   chisq2[m2]=1.01;
	   A2[m2]=r2_copy[0];
	   dA2[m2]=dr2_copy[0];
	   blindR2[m2]=r2_copy[1];
	   dblindR2[m2]=dr2_copy[1];
	   phi2[m2]=r2_copy[2];
	   dphi2[m2]=dr2_copy[2];
	   n2[m2]=m2+0.03;
	  }
	if(m2==2)
	  {
	    chisq2[m2]=1.06;
	    A2[m2]=r3btom_copy[0];
	    dA2[m2]=dr3btom_copy[0];
	    blindR2[m2]=r3btom_copy[1];
	    dblindR2[m2]=dr3btom_copy[1];
	    phi2[m2]=r3btom_copy[2];
	    dphi2[m2]=dr3btom_copy[2];
	    n2[m2]=m2+0.03;
	  }
	if(m2==3)
	  {
	    chisq2[m2]=0.99;
	    A2[m2]=r3no_copy[0];
	    dA2[m2]=dr3no_copy[0];
	    blindR2[m2]=r3no_copy[1];
	    dblindR2[m2]=dr3no_copy[1];
	    phi2[m2]=r3no_copy[2];
	    dphi2[m2]=dr3no_copy[2];
	    n2[m2]=m2+0.03;
	  }
	m2=m2+1;
      }

     for(int i=1;i<=3;i++)
      {
	if(m3==1)
	  {
	   chisq3[m3]=1.9;
	   A3[m3]=r2_q[0];
	   dA3[m3]=dr2_q[0];
	   blindR3[m3]=r2_q[1];
	   dblindR3[m3]=dr2_q[1];
	   phi3[m3]=r2_q[2];
	   dphi3[m3]=dr2_q[2];
	   n3[m3]=m3+0.06;
	  }
	if(m3==2)
	  {
	    chisq3[m3]=2.33;
	    A3[m3]=r3btom_q[0];
	    dA3[m3]=dr3btom_q[0];
	    blindR3[m3]=r3btom_q[1];
	    dblindR3[m3]=dr3btom_q[1];
	    phi3[m3]=r3btom_q[2];
	    dphi3[m3]=dr3btom_q[2];
	    n3[m3]=m3+0.06;
	  }
	if(m3==3)
	  {
	    chisq3[m3]=1.87;
	    A3[m3]=r3no_q[0];
	    dA3[m3]=dr3no_q[0];
	    blindR3[m3]=r3no_q[1];
	    dblindR3[m3]=dr3no_q[1];
	    phi3[m3]=r3no_q[2];
	    dphi3[m3]=dr3no_q[2];
	    n3[m3]=m3+0.06;
	  }
	m3=m3+1;
      }


    c1=new TCanvas("c1","run 3no blind R vs random seed");
    c1->Divide(2,2);


    c1->cd(1);
    gr7=new TGraphErrors(m2,n2,chisq2,0,0);
    gr7->SetTitle("chi-sq vs dataset");
    gr7->GetXaxis()->SetTitle("dataset");
    gr7->GetYaxis()->SetTitle("chi-sq");
    gr7->SetMarkerStyle(22);
    gr7->SetLineColor(kBlue);
    gr7->RemovePoint(0);
    
    gr8=new TGraphErrors(m1,n1,chisq,0,0);
    gr8->SetTitle("chi-sq vs dataset");
    gr8->GetXaxis()->SetTitle("dataset");
    gr8->GetYaxis()->SetTitle("chi-sq");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kRed);
    gr8->RemovePoint(0);

    gr9=new TGraphErrors(m3,n3,chisq3,0,0);
    gr9->SetTitle("chi-sq vs dataset");
    gr9->GetXaxis()->SetTitle("dataset");
    gr9->GetYaxis()->SetTitle("chi-sq");
    gr9->SetMarkerStyle(21);
    gr9->SetLineColor(kBlack);
    gr9->RemovePoint(0);


    gr7->Draw("ap");
    gr8->Draw("p");
    gr9->Draw("p");

    c1->cd(2);
    gr1=new TGraphErrors(m2,n2,A2,0,dA2);
    gr1->SetTitle("A vs dataset");
    gr1->GetXaxis()->SetTitle("dataset");
    gr1->GetYaxis()->SetTitle("A");
    gr1->SetMarkerStyle(22);
    gr1->SetLineColor(kBlue);
    gr1->RemovePoint(0);
    
    gr2=new TGraphErrors(m1,n1,A,0,dA);
    gr2->SetTitle("A vs dataset");
    gr2->GetXaxis()->SetTitle("dataset");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->RemovePoint(0);

    gr10=new TGraphErrors(m3,n3,A3,0,dA3);
    gr10->SetTitle("A vs dataset");
    gr10->GetXaxis()->SetTitle("dataset");
    gr10->GetYaxis()->SetTitle("A");
    gr10->SetMarkerStyle(21);
    gr10->SetLineColor(kBlack);
    gr10->RemovePoint(0);


    gr2->Draw("ap");
    gr1->Draw("p");
    gr10->Draw("p");

    c1->cd(3);
    gr3=new TGraphErrors(m2,n2,blindR2,0,dblindR2);
    gr3->SetTitle("R vs dataset");
    gr3->GetXaxis()->SetTitle("dataset");
    gr3->GetYaxis()->SetTitle("R [ppm]");
    gr3->SetMarkerStyle(22);
    gr3->SetLineColor(kBlue);
    gr3->RemovePoint(0);
    
    gr4=new TGraphErrors(m1,n1,blindR,0,dblindR);
    gr4->SetTitle("R vs dataset");
    gr4->GetXaxis()->SetTitle("dataset");
    gr4->GetYaxis()->SetTitle("R [ppm]");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    gr4->RemovePoint(0);

    gr11=new TGraphErrors(m3,n3,blindR3,0,dblindR3);
    gr11->SetTitle("R vs dataset");
    gr11->GetXaxis()->SetTitle("dataset");
    gr11->GetYaxis()->SetTitle("R [ppm]");
    gr11->SetMarkerStyle(21);
    gr11->SetLineColor(kBlack);
    gr11->RemovePoint(0);

    
    gr4->Draw("ap");
    gr3->Draw("p");
    gr11->Draw("p");

     c1->cd(4);
    gr5=new TGraphErrors(m2,n2,phi2,0,dphi2);
    gr5->SetTitle("phase vs dataset");
    gr5->GetXaxis()->SetTitle("dataset");
    gr5->GetYaxis()->SetTitle("phase");
    gr5->SetMarkerStyle(22);
    gr5->SetLineColor(kBlue);
    gr5->RemovePoint(0);
    
    gr6=new TGraphErrors(m1,n1,phi,0,dphi);
    gr6->SetTitle("phase vs dataset");
    gr6->GetXaxis()->SetTitle("dataset");
    gr6->GetYaxis()->SetTitle("phase");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kRed);
    gr6->RemovePoint(0);

    gr12=new TGraphErrors(m3,n3,phi3,0,dphi3);
    gr12->SetTitle("phase vs dataset");
    gr12->GetXaxis()->SetTitle("dataset");
    gr12->GetYaxis()->SetTitle("phase");
    gr12->SetMarkerStyle(21);
    gr12->SetLineColor(kBlack);
    gr12->RemovePoint(0);


    gr6->Draw("ap");
    gr5->Draw("p");
    gr12->Draw("p");


}
