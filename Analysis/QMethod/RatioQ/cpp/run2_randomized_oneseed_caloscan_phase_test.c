#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include <sys/time.h>
#include "Blinders.hh"

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, "Ritwika's new  Blinding" );


//char root_file_name[128] = "highstat_wigsim_output_statfluc_binintegral_e16.root";
char root_file_name1[128] = "r3no_w8_xtalreject.root";
char root_file_name[128] = "r3no_w8_xtalreject.root";
//char root_file_name[128] = "run2C_thresh300_calosum_ratio.root";
//char root_file_name[128] = "run2_thresh300_fbfDQC.root";
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*hr1[25],*hr2[25],*hr3[25],*hr4[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[50], *hr_sum, *hsum1,*hsum2,*hsum3,*hsum4,*hcomp_sum1[25],*hcomp_sum2[25],*hcomp_sum3[25],*hcomp_sum4[25],*hr_shift1,*hr_noshift1,*hr_shift2,*hr_noshift2,*hr_shift3,*hr_noshift3,*hr_shift4,*hr_noshift4,*hcalosum,*h_res_sum,*hcalosum1,*hcalosum2,*hcalosum3,*hcalosum4;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cauto,*cfit;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12, *gr13, *gr14, *gr15, *gr16, *gr17, *gr18, *gr19, *gr20, *gr21, *gr22, *gr23, *gr24, *gr25, *gr26, *gr27, *gr28;
TGraph *gkband1p,*gkband1m,*gkband2p,*gkband2m,*gkband3p,*gkband3m,*gkband4p,*gkband4m,*gkband5p,*gkband5m,*gkband6p,*gkband6m,*gkband7p,*gkband7m,*gkband8p,*gkband8m,*gkband9p,*gkband9m,*gkband10p,*gkband10m,*gkband11p,*gkband11m,*gkband12p,*gkband12m,*gkband13p,*gkband13m,*gkband14p,*gkband14m,*gkband15p,*gkband15m,*gkband16p,*gkband16m,*gkband17p,*gkband17m,*gkband18p,*gkband18m,*gkband19p,*gkband19m,*gkband20p,*gkband20m,*gkband21p,*gkband21m,*gkband22p,*gkband22m,*gkband23p,*gkband23m;
TH1 *hfft[50],*hFFTsq[50],*hAuto[50],*hfftsum;
//TH1 *hfft, *hFFTsq, *hAuto;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=30000;
Double_t fit_stop=300000;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
Int_t m=0;
Int_t corr_coeff;
Int_t countfcn=0;

Double_t chisq[50], A[50], dA[50], blindR[50], dblindR[50], phi_0[50], dphi_0[50], A_cbo_N[50], tau_cbo[50], omega_cbo[50], phi_cbo_N[50], A_cbo_A[50], phi_cbo_A[50], A_cbo_phi[50], phi_cbo_phi[50], A_vw[50], tau_vw[50], omega_vw[50], phi_vw[50], A_y[50], tau_y[50], omega_y[50], phi_y[50], A_2cbo[50], tau_2cbo[50], omega_2cbo[50], phi_2cbo[50],  dA_cbo_N[50], dtau_cbo[50], domega_cbo[50], dphi_cbo_N[50], dA_cbo_A[50], dphi_cbo_A[50], dA_cbo_phi[50], dphi_cbo_phi[50], dA_vw[50], dtau_vw[50], domega_vw[50], dphi_vw[50], dA_y[50], dtau_y[50], domega_y[50], dphi_y[50], dA_2cbo[50], dtau_2cbo[50], domega_2cbo[50], dphi_2cbo[50], n[50],kband1p[50],kband1m[50],kband2p[50],kband2m[50],kband3p[50],kband3m[50],kband4p[50],kband4m[50],kband5p[50],kband5m[50],kband6p[50],kband6m[50],kband7p[50],kband7m[50],kband8p[50],kband8m[50],kband9p[50],kband9m[50],kband10p[50],kband10m[50],kband11p[50],kband11m[50],kband12p[50],kband12m[50],kband13p[50],kband13m[50],kband14p[50],kband14m[50],kband15p[50],kband15m[50],kband16p[50],kband16m[50],kband17p[50],kband17m[50],kband18p[50],kband18m[50],kband19p[50],kband19m[50],kband20p[50],kband20m[50],kband21p[50],kband21m[50],kband22p[50],kband22m[50],kband23p[50],kband23m[50];

//Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25];
TFile *foutput;



TH2D *hcov;
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
TMatrixD cov(mdim,mdim),cov2(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim));

TMatrixD ratio_cov_matrix(TH1D *hr_sum){

  double covar;

  for(Int_t k=hr_sum->FindBin(10000); k<=hr_sum->GetNbinsX(); k++)
    {
      if(hr_sum->GetBinError(k)!=0)
	{ i0=k;
	  break;}
    }
  cout<<"i0 is"<<i0<<" "<<endl;
  
   for (Int_t i = 0; i < (mdim)*(mdim); i++)
      {
	const Int_t ir = (i)/(mdim);
	const Int_t ic = (i)%(mdim);

        
	
        data[i] = 0.0;
 
        if(ir==ic)
	  {
	    //  data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir))/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
           data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	  }

	  if(ir==ic-1)
	  {
           data[i]=(corr_coeff)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
	  }


	  if(ir==ic+1)
	  {
	   data[i]=(corr_coeff)*sqrt(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic)*hr_sum->GetBinError(ic));
      	  }


	
      }
  cout<<"flgs "<<flg<<" "<<endl;
  cov.SetMatrixArray(data.GetArray());


  cov.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
   
  //   cov.Print();
    cov.SetTol(1.e-23);
    Double_t det1;
    cov.Invert(&det1);
     cout<<mdim<<endl;
    //cout<<"Matrix inverted "<<i0<<endl;
  return cov;

}


Double_t fprec3(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    //  return 0.22*cos(omega_a*time + phi);
    return (2*f - ff - fb)/(2*f + ff + fb);

}



TH1D* construct_rhist_rand(TH1D *hsum1, TH1D *hsum2, TH1D *hsum3, TH1D *hsum4)
{
 


   //create the component histograms of the ratio histograms   
   h1=(TH1D*)hsum1->Clone();
   h2=(TH1D*)hsum2->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

      nbinshift=(0.5*T_a_true)/h1->GetBinWidth(1);
      /*  if(double(0.5*T_a_true/h1->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   //   nbinshift=20;
   cout<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=hsum3->GetNbinsX(); ibin++)             
     {
       hp->SetBinContent( ibin, hsum3->GetBinContent( ibin + nbinshift));
       hp->SetBinError(ibin, hsum3->GetBinError( ibin + nbinshift));
     }

   for(int ibin=nbinshift; ibin<=hsum4->GetNbinsX(); ibin++)
     {
       hm->SetBinContent( ibin, hsum4->GetBinContent( ibin - nbinshift));
       hm->SetBinError( ibin, hsum4->GetBinError(ibin - nbinshift));
     }
 
    cout<<hsum1->GetBinWidth(100)<<endl;
   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*hsum1->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   double blife=exp(-(T_a)/(2*lifetime));
   double deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo_histogram_sum_ratio", "h_ratio", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));
       
      for(int ibin=hsum1->FindBin(0);ibin<=hsum1->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo histogram sum numerator", "h_num", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));

      h_num->Add(h1,1);
      h_num->Add(h2,1);
      h_num->Add(hp,-1);
      h_num->Add(hm,-1);

      h_deno->Add(h1,1);
      h_deno->Add(h2,1);
      h_deno->Add(hp,1);
      h_deno->Add(hm,1);
     
      h_num->Divide(h_deno);

     }
    //Assign error bars
   
 if(useerrorbars)
     {
       
  
       long double t1,t2,t3, t4;
      for(int ibin=hsum1->FindBin(0); ibin<=hsum1->GetNbinsX()-nbinshift; ibin++)
      {
	t1=2*(h1->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

	t2=2*(h2->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t3=-2*(hp->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	
	t4=-2*(hm->GetBinError(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));

	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( (t1*t2) + (t2*t2) + (t3*t3) + (t4*t4) ) );
       }
      }

     }

 //send a common histogram name for fitting
 if(usesetbins)
      {
       hcalo=(TH1D*)h_ratio->Clone();
      }
   else
      {
       hcalo=(TH1D*)h_num->Clone();
      }
 cout<<hsum4->GetBinWidth(500)<<endl;


 /*  foutput=new TFile("ratio_4hist_randomized.root","new");
 h_ratio->Write();
 foutput->Close();
 */

 cout<<hcalo->GetBinWidth(500)<<endl; 
 return hcalo;
  }


void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t inf)
{

 
  
    TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  //  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;
  double ch=0;

   if(hr_sum->FindBin(fit_start)<i0)
    {
      cout<<"Wrong Start of Fit!! Returning"<<endl;
      return;
    }
   
   for(Int_t i=hr_sum->FindBin(fit_start); i<=hr_sum->FindBin(fit_stop); i++)
    {
      for(Int_t j=hr_sum->FindBin(fit_start); j<=hr_sum->FindBin(fit_stop); j++)
	{
	  //  for(m=0; m<=105; m++)
	  // {
	  //  if(i==j || i==j+15 || i==j-15 || i==j+30 || i==j-30 || i==j+45 || i==j-45 || i==j+60 || i==j-60 || i==j+75 || i==j-75 || i==j+90 || i==j-90 || i==j+105 || i== j-105 || i==j+120 || i==j-120 || i==j+135 || i==j-135 || i==j+150 || i==j+165 || i==j-165 || i==j+180 || i==j-180)
	   if(cov2[i][j]!=0)
	     {
	      ch=ch+((hr_sum->GetBinContent(i))-(fuser->Eval(hr_sum->GetBinCenter(i))))*cov2[i][j]*((hr_sum->GetBinContent(j))-(fuser->Eval(hr_sum->GetBinCenter(j))));
	      //cout<<m<<endl;
	    
	     }
	     // }
	} 
     }
   
  f=TMath::Abs(ch);
  countfcn++;
  cout<<"fcn call number"<<countfcn<<"  "<<f<<endl;
}


void randomized_oneseed_3par_fft_caloscan()
{

  
  _file[1]=TFile::Open(root_file_name);
  _file[2]=TFile::Open(root_file_name1);
  //  _file[1]->GetObject("QRatioFillByFillAnalyzerDB",dir[1]);
  // dir[1]->GetObject("qHist_1D_sig_sum_1",qHist_1D[1]);
  // _file[1]->GetObject("hwiggle",qHist_1D[1]);
  // hsum1=(TH1D*)qHist_1D[1]->Clone();

    for(int i=1;i<=24;i++)
    {
     hr1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_0_0", i)));
    }

  for(int i=1;i<=24;i++)
    {
     hr2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_1_0", i)));
    }


  for(int i=1;i<=24;i++)
    {
     hr3[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_2_0", i)));
    }


  for(int i=1;i<=24;i++)
    {
     hr4[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_rm_%d_3_0", i)));
    }
  
  hcalosum1= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_0_0")));
  hcalosum2= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_1_0")));
  hcalosum3= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_2_0")));
  hcalosum4= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_rm_sum_3_0")));

  int start_bin=400;
  double binwidth=hr1[1]->GetBinWidth(1);
  int nbins=hr1[1]->GetNbinsX();
  double binlowedge=hr1[1]->GetBinLowEdge(1);
  double binhighedge=hr1[1]->GetBinLowEdge(nbins)+binwidth;
   
  cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

  hcalosum2->Scale(hcalosum1->Integral(start_bin,16800)/hcalosum2->Integral(start_bin,16800));
  hcalosum3->Scale(hcalosum1->Integral(start_bin,16800)/hcalosum3->Integral(start_bin,16800));
  hcalosum4->Scale(hcalosum1->Integral(start_bin,16800)/hcalosum4->Integral(start_bin,16800));
 
    cout<<"Normalized from bin "<<start_bin<<endl;
 
    //h_sum= new TH1D("calo histogram sum", "h_sum", hsum1->GetNbinsX(), 100001, 352001);
  //hrsum= new TH1D("calo histogram sum ratio", "calosum ratio", hsum1->GetNbinsX(), 100001, 352001);

     hcalosum1->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hcalosum2->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hcalosum3->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hcalosum4->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     
   cout<<hcalosum1->GetBinWidth(1)<<" "<<hcalosum1->GetNbinsX()<<" "<<hcalosum1->GetBinLowEdge(1)<<" "<<(hcalosum1->GetBinLowEdge(hcalosum1->GetNbinsX())+hcalosum1->GetBinWidth(1))<<" "<<endl;

    hcalosum1->Rebin(8);
    hcalosum1->Scale(0.125);
    hcalosum2->Rebin(8);
    hcalosum2->Scale(0.125);
    hcalosum3->Rebin(8);
    hcalosum3->Scale(0.125);
    hcalosum4->Rebin(8);
    hcalosum4->Scale(0.125);



 for(int i=1;i<=24;i++)
   {
   
    hr2[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr2[i]->Integral(start_bin,16800));
    hr3[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr3[i]->Integral(start_bin,16800));
    hr4[i]->Scale(hr1[i]->Integral(start_bin,16800)/hr4[i]->Integral(start_bin,16800));
 
    cout<<"Normalized from bin "<<start_bin<<endl;
 
    //h_sum= new TH1D("calo histogram sum", "h_sum", hsum1->GetNbinsX(), 100001, 352001);
  //hrsum= new TH1D("calo histogram sum ratio", "calosum ratio", hsum1->GetNbinsX(), 100001, 352001);

     hr1[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr2[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr3[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hr4[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     
   cout<<hr1[i]->GetBinWidth(1)<<" "<<hr1[i]->GetNbinsX()<<" "<<hr1[i]->GetBinLowEdge(1)<<" "<<(hr1[i]->GetBinLowEdge(hr1[i]->GetNbinsX())+hr1[i]->GetBinWidth(1))<<" "<<endl;

    hr1[i]->Rebin(4);
    hr1[i]->Scale(0.25);
    hr2[i]->Rebin(4);
    hr2[i]->Scale(0.25);
    hr3[i]->Rebin(4);
    hr3[i]->Scale(0.25);
    hr4[i]->Rebin(4);
    hr4[i]->Scale(0.25);
   
   

    hr_noshift1=(TH1D*)hr1[i]->Clone();
    hr_noshift2=(TH1D*)hr2[i]->Clone();
    hr_noshift3=(TH1D*)hr3[i]->Clone();
    hr_noshift4=(TH1D*)hr4[i]->Clone();
	
    hr_shift1=(TH1D*)hr_noshift1->Clone();
    hr_shift1->Reset();
    hr_shift2=(TH1D*)hr_noshift2->Clone();
    hr_shift2->Reset();
    hr_shift3=(TH1D*)hr_noshift3->Clone();
    hr_shift3->Reset();
    hr_shift4=(TH1D*)hr_noshift4->Clone();
    hr_shift4->Reset();


    for(int ibin=1;ibin<=hr_shift1->GetNbinsX();ibin++)
	{
	  hr_shift1->SetBinContent(ibin,hr_noshift1->GetBinContent(ibin+1));
	  hr_shift1->SetBinError(ibin,hr_noshift1->GetBinError(ibin+1));
	  hr_shift2->SetBinContent(ibin,hr_noshift2->GetBinContent(ibin+1));
	  hr_shift2->SetBinError(ibin,hr_noshift2->GetBinError(ibin+1));
          hr_shift3->SetBinContent(ibin,hr_noshift3->GetBinContent(ibin+1));
	  hr_shift3->SetBinError(ibin,hr_noshift3->GetBinError(ibin+1));
          hr_shift4->SetBinContent(ibin,hr_noshift4->GetBinContent(ibin+1));
	  hr_shift4->SetBinError(ibin,hr_noshift4->GetBinError(ibin+1));
	}

	hcomp_sum1[i]=(TH1D*)hr1[i]->Clone();
	hcomp_sum1[i]->Reset();
	hcomp_sum2[i]=(TH1D*)hr2[i]->Clone();
	hcomp_sum2[i]->Reset();
	hcomp_sum3[i]=(TH1D*)hr3[i]->Clone();
	hcomp_sum3[i]->Reset();
	hcomp_sum4[i]=(TH1D*)hr4[i]->Clone();
	hcomp_sum4[i]->Reset();


	for(int ibin=1;ibin<=hr2[i]->GetNbinsX();ibin++)
	  {
	    hcomp_sum1[i]->SetBinContent(ibin,0.5*(hr_noshift1->GetBinContent(ibin)+hr_shift1->GetBinContent(ibin)));
	    hcomp_sum1[i]->SetBinError(ibin,0.5*sqrt((hr_noshift1->GetBinError(ibin)*hr_noshift1->GetBinError(ibin))+(hr_shift1->GetBinError(ibin)*hr_shift1->GetBinError(ibin))));
	    hcomp_sum2[i]->SetBinContent(ibin,0.5*(hr_noshift2->GetBinContent(ibin)+hr_shift2->GetBinContent(ibin)));
	    hcomp_sum2[i]->SetBinError(ibin,0.5*sqrt((hr_noshift2->GetBinError(ibin)*hr_noshift2->GetBinError(ibin))+(hr_shift2->GetBinError(ibin)*hr_shift2->GetBinError(ibin))));
	    hcomp_sum3[i]->SetBinContent(ibin,0.5*(hr_noshift3->GetBinContent(ibin)+hr_shift3->GetBinContent(ibin)));
	    hcomp_sum3[i]->SetBinError(ibin,0.5*sqrt((hr_noshift3->GetBinError(ibin)*hr_noshift3->GetBinError(ibin))+(hr_shift3->GetBinError(ibin)*hr_shift3->GetBinError(ibin))));
	    hcomp_sum4[i]->SetBinContent(ibin,0.5*(hr_noshift4->GetBinContent(ibin)+hr_shift4->GetBinContent(ibin)));
	    hcomp_sum4[i]->SetBinError(ibin,0.5*sqrt((hr_noshift4->GetBinError(ibin)*hr_noshift4->GetBinError(ibin))+(hr_shift4->GetBinError(ibin)*hr_shift4->GetBinError(ibin))));
	  }

	
	hcomp_sum1[i]->Rebin(2);
	hcomp_sum1[i]->Scale(0.5);
	hcomp_sum2[i]->Rebin(2);
	hcomp_sum2[i]->Scale(0.5);
	hcomp_sum3[i]->Rebin(2);
	hcomp_sum3[i]->Scale(0.5);
	hcomp_sum4[i]->Rebin(2);
	hcomp_sum4[i]->Scale(0.5);

	hr_shift1->Reset();
	hr_noshift1->Reset();
	hr_shift2->Reset();
	hr_noshift2->Reset();
	hr_shift3->Reset();
	hr_noshift3->Reset();
        hr_shift4->Reset();
	hr_noshift4->Reset();


   }
  //return;

 //   hr->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   //hsum1->Divide(hsum2);
   
   
 //  hr = construct_rhist_copy(hsum1,hsum2,hsum3,hsum4);
   
 
  


     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A", "Blind R", "Phi");
     fit_func3->SetNpx(1000000);

  
 


       cauto=new TCanvas("auto-correlation", "auto-correlation");
       cauto->Divide(6,4);

        cfit=new TCanvas("23 parameter fit","23 parameter fit");

	hcalosum=construct_rhist_rand(hcalosum1,hcalosum2,hcalosum3,hcalosum4);

  
	hcalosum->GetYaxis()->SetTitle("ADC counts");
	hcalosum->GetXaxis()->SetTitle("time [ns]");


        hcalosum->GetXaxis()->SetRangeUser(fit_start,fit_stop);



        hcalosum->Fit("fprec3","RE","",fit_start,fit_stop);


      		
     //hr_sum=(TH1D*)hcalo->Clone();
       
          for(int i=1; i<=24;i++)
       {

	cout<<"fit start time"<<fit_start<<" ns"<<endl; 
        
	 
	//	fit_func3->SetParameters(0.22, 0.0, 2.22);

		fit_func3->SetParameters(0.23, 0.0, 2.22);

		hr1[i]->Rebin(2);
		hr1[i]->Scale(0.5);
		hr2[i]->Rebin(2);
		hr2[i]->Scale(0.5);
                hr3[i]->Rebin(2);
		hr3[i]->Scale(0.5);
                hr4[i]->Rebin(2);
		hr4[i]->Scale(0.5);
		hcalo=construct_rhist_rand(hr1[i],hr2[i],hr3[i],hr4[i]);


		//	hcalo=construct_rhist_rand(hcomp_sum1[i],hcomp_sum2[i],hcomp_sum3[i],hcomp_sum4[i]);
  
	hcalo->GetYaxis()->SetTitle("ADC counts");
	hcalo->GetXaxis()->SetTitle("time [ns]");


        hcalo->GetXaxis()->SetRangeUser(fit_start,fit_stop);



        hcalo->Fit("fprec3","RE","",fit_start,fit_stop);
	//	hcalo->Draw();
	//	fit_func3->Draw();
	//	hcalo->Draw("same");
	
	
   
   
             gStyle->SetOptFit(1111);

       h_res[i]= new TH1D("residual histogram ", "h_res", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalo->GetBinContent(ibin)- fit_func3->Eval( hcalo->GetBinCenter(ibin) ) );
          if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.05);
        h_res[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(h_res[i]->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.05);
       hfft[i]->GetYaxis()->SetLabelSize(0.05);

       

     
       hFFTsq[i] = (TH1D*)hfft[i]->Clone();
       hFFTsq[i]->Reset();
       for (int ib = 1; ib <= hFFTsq[i]->GetNbinsX(); ib++)
        {
         hFFTsq[i]->SetBinContent( ib, hfft[i]->GetBinContent(ib)*hfft[i]->GetBinContent(ib));
        }


       hAuto[i] = hFFTsq[i]->FFT(hAuto[i], "MAG");
       //hAuto[i]->SetTitle("auto-correlation");
       hAuto[i]->GetXaxis()->SetLabelSize(0.04);
       hAuto[i]->GetXaxis()->SetTitle("bins");
       hAuto[i]->GetXaxis()->SetRangeUser(0,h_res[i]->GetNbinsX()/2);
       hAuto[i]->SetLineColor(kRed);
       hAuto[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       cauto->cd(i);
       hAuto[i]->Draw("HIST");

       corr_coeff=hAuto[i]->GetBinContent(2)/hAuto[i]->GetBinContent(1);

       hfft[i]->Reset();
       h_res[i]->Reset();

       cov2.ResizeTo(hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), hcalo->FindBin(fit_start), hcalo->FindBin(fit_stop), -1);
   
  //   cov.Print();
        cov2.SetTol(1.e-23);


       cov2=ratio_cov_matrix(hcalo);
       
       
       hr_sum=(TH1D*)hcalo->Clone();

     
       printf("start Ratio Fit\n");
       struct timeval t_start, t_end;
       gettimeofday(&t_start, NULL);
       gFitter = TVirtualFitter::Fitter(hr_sum);
       gFitter->SetFCN(chi2);

      
       cfit->cd(1);	
       hr_sum->Fit("fprec3","RE","",fit_start,fit_stop);
       countfcn=0;
       gStyle->SetOptFit(1111);

       gettimeofday(&t_end, NULL);
       printf("QRatio fit duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));
     
       //cout<<" setting ch to zero "<<endl;
       //ch=0;
 
     
            n[m]=i;
       m=m+1;
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr_sum->GetBinContent(ibin)- fit_func3->Eval( hr_sum->GetBinCenter(ibin) ) );
          if(hr_sum->GetBinError(ibin)!=0){res=(res/hr_sum->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.05);
        h_res[i]->GetYaxis()->SetLabelSize(0.05);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hr_sum->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.05);
       hfft[i]->GetYaxis()->SetLabelSize(0.05);


       
 
   }

    h_res_sum= new TH1D("calosum residual histogram ", "h_res_sum", hr_sum->GetNbinsX(), hr_sum->GetBinLowEdge(1), (hr_sum->GetBinLowEdge(hr_sum->GetNbinsX())+hr_sum->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr_sum->GetBinLowEdge(1))/hr_sum->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalosum->GetBinContent(ibin)- fit_func3->Eval( hcalosum->GetBinCenter(ibin) ) );
          if(hcalosum->GetBinError(ibin)!=0){res=(res/hcalosum->GetBinError(ibin));}
          h_res_sum->SetBinContent(ibin, (res)  );
     
        }

	h_res_sum->GetYaxis()->SetTitle("Energy [MeV]");
	h_res_sum->GetXaxis()->SetTitle("time [ns]");
	h_res_sum->GetXaxis()->SetLabelSize(0.05);
        h_res_sum->GetYaxis()->SetLabelSize(0.05);

      

 
       hfftsum = h_res_sum->FFT(hfftsum, "MAG");
       hfftsum->SetLineColor(kBlack);
       hfftsum->SetBins(hcalosum->GetNbinsX(),0,1000/h_res_sum->GetBinWidth(1));
       hfftsum->GetXaxis()->SetRangeUser(0,1000/(2*h_res_sum->GetBinWidth(1)));
       hfftsum->GetXaxis()->SetTitle("Freq [MHz]");
       hfftsum->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfftsum->GetXaxis()->SetLabelSize(0.05);
       hfftsum->GetYaxis()->SetLabelSize(0.05);
       hfftsum->SetLineColor(kRed);


	  
  c2 = new TCanvas("c2","start time residuals");
  c2->Divide(4,6);
  c3 = new TCanvas("c3","start time fft");
  c3->Divide(4,6);
   
 for(int i=1;i<=24;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hfft[i]->GetXaxis()->SetRangeUser(2.1,2.4);
     hfftsum->GetXaxis()->SetRangeUser(2.1,2.4);
     hfftsum->Draw();
     hfft[i]->Draw("same");
   }
 
   
}
