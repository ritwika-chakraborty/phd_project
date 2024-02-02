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
//char root_file_name[128] = "run2g_corrOOF_thresh400.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
char root_file_name[128] = "run2C_thresh300_calosum_ratio.root";
TFile *_file[52];
TDirectory *dir[52];
TH1D *h[50],*qHist_1D[52],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[52], *h_res[52], *hr1[25], *hr2[25], *hr3[25], *hr4[25], *hrsum, *hsum1, *hsum2, *hsum3, *hsum4;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6;
TGraphErrors *gr1, *gr2, *gr3, *gr4, *gr5, *gr6, *gr7, *gr8, *gr9, *gr10, *gr11, *gr12;
TH1 *hfft[52];
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start, fit_stop;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
Int_t m=0;
Double_t A[52], dA[52], blindR[52], dblindR[52], phi[52], dphi[52], A_cbo[52], Tau_cbo[52], w_cbo[52], phi_cbo[52], dA_cbo[52], dTau_cbo[52], dw_cbo[52], dphi_cbo[52], A_vw[52], Tau_vw[52], w_vw[52], phi_vw[52], dA_vw[52], dTau_vw[52], dw_vw[52], dphi_vw[52], n[52];

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

    return (2*f - ff - fb)/(2*f + ff + fb);

}

Double_t fprec7(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo - ff*Ncbof - fb*Ncbob)/(2*f*Ncbo + ff*Ncbof + fb*Ncbob);

}

Double_t fprec11(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_cbo= par[3];

    double tau_cbo = par[4];

    double omega_cbo = par[5];

    double phi_cbo = par[6];

    double asym_vw= par[7];

    double tau_vw = par[8];

    double omega_vw = par[9];

    double phi_vw = par[10];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Ncbob=(1+ asym_cbo*exp(-(time - T_a/2)/tau_cbo)*cos(omega_cbo*(time - T_a/2) + phi_cbo));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvwb=(1+ asym_vw*exp(-(time - T_a/2)/tau_vw)*cos(omega_vw*(time - T_a/2) + phi_vw));

    double  f=(1+ asym*cos(omega_a*time + phi));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    // double fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));

    return (2*f*Ncbo*Nvw - ff*Ncbof*Nvwf - fb*Ncbob*Nvwb)/(2*f*Ncbo*Nvw + ff*Ncbof*Nvwf + fb*Ncbob*Nvwb);

}


TH1D* construct_rhist_copy(TH1D *hsum1, TH1D *hsum2, TH1D *hsum3, TH1D *hsum4, int nbins)
{
 
   h1=(TH1D*)hsum1->Clone();
   h2=(TH1D*)hsum2->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

   //   nbinshift=(0.5*T_a_true)/h1->GetBinWidth(1);
      /*   if(double(0.5*T_a_true/h1->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   //   nbinshift=20;
   nbinshift=nbins;
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
 
 // hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //hcalo->Scale(0.125);

 /*  foutput=new TFile("ratio_4hist_randomized.root","new");
 h_ratio->Write();
 foutput->Close();
 */

 cout<<hcalo->GetBinWidth(500)<<endl;
 //return;
 return hcalo;
  
  }

void ratio_nbinshift_randomized()
{
  
  _file[1]=TFile::Open(root_file_name);


    for(int i=1;i<=24;i++)
    {
     hr1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_0", i)));
    }

  for(int i=1;i<=24;i++)
    {
     hr2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_1", i)));
    }


  for(int i=1;i<=24;i++)
    {
     hr3[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_2", i)));
    }


  for(int i=1;i<=24;i++)
    {
     hr4[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d_3", i)));
    }

  hsum1= new TH1D("calo histogram ratio sum 1", "hsum1", hr1[1]->GetNbinsX(), 100001, 352001);
  hsum2= new TH1D("calo histogram ratio sum 2", "hsum2", hr2[1]->GetNbinsX(), 100001, 352001);
  hsum3= new TH1D("calo histogram ratio sum 3", "hsum3", hr3[1]->GetNbinsX(), 100001, 352001);
  hsum4= new TH1D("calo histogram ratio sum 4", "hsum4", hr4[1]->GetNbinsX(), 100001, 352001);
  
  hsum1->Sumw2(kTRUE);
  for(int i=1;i<=24;i++)
    {
      hsum1->Add(hr1[i],1);
    }
  hsum2->Sumw2(kTRUE);
   for(int i=1;i<=24;i++)
    {
      hsum2->Add(hr2[i],1);
    }
  hsum3->Sumw2(kTRUE);
   for(int i=1;i<=24;i++)
    {
      hsum3->Add(hr3[i],1);
    }
   hsum4->Sumw2(kTRUE);
 for(int i=1;i<=24;i++)
    {
      hsum4->Add(hr4[i],1);
    }

 hsum2->SetLineColor(kRed);
 hsum3->SetLineColor(kBlack);
 hsum4->SetLineColor(kGreen);
 
 int start_bin=400;
 hsum2->Scale(hsum1->Integral(start_bin,16800)/hsum2->Integral(start_bin,16800));
 hsum3->Scale(hsum1->Integral(start_bin,16800)/hsum3->Integral(start_bin,16800));
 hsum4->Scale(hsum1->Integral(start_bin,16800)/hsum4->Integral(start_bin,16800));
  
  cout<<"Normalized from bin "<<start_bin<<endl;
  
  h_sum= new TH1D("calo histogram sum", "h_sum", hsum1->GetNbinsX(), 100001, 352001);
  hrsum= new TH1D("calo histogram sum ratio", "calosum ratio", hsum1->GetNbinsX(), 100001, 352001);

   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   
     hsum1->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum2->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum3->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     hsum4->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

   /*  hsum1->Rebin(8);
   hsum1->Scale(0.125);
   hsum2->Rebin(8);
   hsum2->Scale(0.125);
   hsum3->Rebin(8);
   hsum3->Scale(0.125);
   hsum4->Rebin(8);
   hsum4->Scale(0.125);
   */
   




     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A", "Blind R", "Phi");
     fit_func3->SetNpx(1000000);

     fit_func7= new TF1("fprec7", fprec7,  30000,  309000, 7);
     fit_func7->SetParNames("A", "Blind R", "Phi");
     fit_func7->SetParName(3,"A_cbo");
     fit_func7->SetParName(4,"Tau_cbo");
     fit_func7->SetParName(5,"omega_cbo");
     fit_func7->SetParName(6,"phi_cbo");
     fit_func7->SetNpx(1000000);
     
     fit_func11= new TF1("fprec11", fprec11,  30000,  309000, 11);
     fit_func11->SetParNames("A", "Blind R", "Phi");
     fit_func11->SetParName(3,"A_cbo");
     fit_func11->SetParName(4,"Tau_cbo");
     fit_func11->SetParName(5,"omega_cbo");
     fit_func11->SetParName(6,"phi_cbo");
     fit_func11->SetParName(7,"A_vw");
     fit_func11->SetParName(8,"Tau_vw");
     fit_func11->SetParName(9,"omega_vw");
     fit_func11->SetParName(10,"phi_vw");
     fit_func11->SetNpx(1000000);


 

     fit_start=30000;
     fit_stop=300000;

     nbinshift=95;
     //  h_sum->Rebin(8);
     //h_sum->Scale(0.125);
     
     for(int i=1; i<=50;i++)
       {
        
	cout<<"Calo "<<i<<endl;
	//	hr[i]->Rebin(8);
	//	hr[i]->Scale(0.125);
        hr[i] = construct_rhist_copy(hsum1, hsum2, hsum3, hsum4, nbinshift);
       	hr[i]->Rebin(8);
       	hr[i]->Scale(0.125);
	nbinshift=nbinshift+1;
	 
	fit_func3->SetParameters(0.229, 0.0, 2.22);

	hr[i]->Fit("fprec3","RE","",fit_start,fit_stop);
  
	hr[i]->GetYaxis()->SetTitle("ADC counts");
	hr[i]->GetXaxis()->SetTitle("time [ns]");
        hr[i]->GetXaxis()->SetRangeUser(fit_start,fit_stop);

        fit_func7->SetParameter(0, fit_func3->GetParameter(0));
        fit_func7->SetParameter(1, fit_func3->GetParameter(1));
        fit_func7->SetParameter(2, fit_func3->GetParameter(2));
        fit_func7->SetParameter(3, 0.0024);
        fit_func7->SetParameter(4, 250000);
        fit_func7->SetParameter(5, 0.0023405);
        fit_func7->SetParameter(6, 0.48);

   
        hr[i]->Fit("fprec7","RE","",fit_start,fit_stop);
 

       fit_func11->SetParameter(0, fit_func7->GetParameter(0));
       fit_func11->SetParameter(1, fit_func7->GetParameter(1));
       fit_func11->SetParameter(2, fit_func7->GetParameter(2));
       fit_func11->SetParameter(3, fit_func7->GetParameter(3));
       fit_func11->SetParameter(4, fit_func7->GetParameter(4));
       fit_func11->SetParameter(5, fit_func7->GetParameter(5));
       fit_func11->SetParameter(6, fit_func7->GetParameter(6));
       fit_func11->SetParameter(7, 0.005);
       fit_func11->SetParameter(8, 29000);
       fit_func11->SetParameter(9, 0.01401);
       fit_func11->SetParameter(10, 2.6);


    
       hr[i]->Fit("fprec11","RE","",fit_start,fit_stop);
       
       gStyle->SetOptFit(1111);

       A[m]=fit_func11->GetParameter(0);
       dA[m]=fit_func11->GetParError(0);
       
       blindR[m]=fit_func11->GetParameter(1);
       dblindR[m]=fit_func11->GetParError(1);
       
       phi[m]=fit_func11->GetParameter(2);
       dphi[m]=fit_func11->GetParError(2);

       A_cbo[m]=fit_func11->GetParameter(3);
       dA_cbo[m]=fit_func11->GetParError(3);

       Tau_cbo[m]=fit_func11->GetParameter(4);
       dTau_cbo[m]=fit_func11->GetParError(4);

       w_cbo[m]=fit_func11->GetParameter(5);
       dw_cbo[m]=fit_func11->GetParError(5);

       phi_cbo[m]=fit_func11->GetParameter(6);
       dphi_cbo[m]=fit_func11->GetParError(6);

       A_vw[m]=fit_func11->GetParameter(7);
       dA_vw[m]=fit_func11->GetParError(7);

       Tau_vw[m]=fit_func11->GetParameter(8);
       dTau_vw[m]=fit_func11->GetParError(8);

       w_vw[m]=fit_func11->GetParameter(9);
       dw_vw[m]=fit_func11->GetParError(9);

       phi_vw[m]=fit_func11->GetParameter(10);
       dphi_vw[m]=fit_func11->GetParError(10);
	 
       if(fit_func11->GetParError(1)<0.5)
	 {
	  A[m]=fit_func7->GetParameter(0);
          dA[m]=fit_func7->GetParError(0);

          blindR[m]=fit_func7->GetParameter(1);
          dblindR[m]=fit_func7->GetParError(1);

	  phi[m]=fit_func7->GetParameter(2);
          dphi[m]=fit_func7->GetParError(2);
	  
	  if(fit_func7->GetParError(1)<0.5)
	   {
            A[m]=fit_func3->GetParameter(0);
            dA[m]=fit_func3->GetParError(0);
	     
            blindR[m]=fit_func3->GetParameter(1);
            dblindR[m]=fit_func3->GetParError(1);

	    phi[m]=fit_func3->GetParameter(2);
            dphi[m]=fit_func3->GetParError(2);
	   }

	 }
       n[m]=nbinshift;
       m=m+1;
       
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hr[i]->GetNbinsX(), hr[i]->GetBinLowEdge(1), (hr[i]->GetBinLowEdge(hr[i]->GetNbinsX())+hr[i]->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr[i]->GetBinContent(ibin)- fit_func11->Eval( hr[i]->GetBinCenter(ibin) ) );
          if(hr[i]->GetBinError(ibin)!=0){res=(res/hr[i]->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("ADC counts");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetLabelSize(0.07);
        h_res[i]->GetYaxis()->SetLabelSize(0.07);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hr[i]->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [GHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.07);
       hfft[i]->GetYaxis()->SetLabelSize(0.06);

 
   }

  c2 = new TCanvas("c2","per calo residuals");
  c2->Divide(10,5);
  c3 = new TCanvas("c3","per calo fft");
  c3->Divide(10,5);
   
 for(int i=1;i<=50;i++)
   {
     c2->cd(i);
     h_res[i]->Draw();
     c3->cd(i);
     hfft[i]->Draw();
   }

 
    c1=new TCanvas("c1","Ratio blind R vs caloriemeter #");
    gr1=new TGraphErrors(m,n,dblindR,0,0);
    gr1->SetTitle(" dR vs 15 time decimation bin-shift");
    gr1->GetXaxis()->SetTitle("binshift");
    gr1->GetYaxis()->SetTitle("dR [ppm]");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    gr1->Draw();
    
    c4=new TCanvas("c4","nominal fit parameters vs caloriemeter #");
    c4->Divide(1,3);
    c4->cd(1);
    gr2=new TGraphErrors(m,n,A,0,dA);
    gr2->SetTitle(" Asymmetry vs 15 time decimation bin-shift");
    gr2->GetXaxis()->SetTitle("binshift");
    gr2->GetYaxis()->SetTitle("A");
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    gr2->Draw();

    c4->cd(2);
    gr3=new TGraphErrors(m,n,blindR,0,dblindR);
    gr3->SetTitle(" Blind R vs 15 time decimation bin-shift");
    gr3->GetXaxis()->SetTitle("binshift");
    gr3->GetYaxis()->SetTitle("R");
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    gr3->Draw();

    c4->cd(3);
    gr4=new TGraphErrors(m,n,phi,0,dphi);
    gr4->SetTitle(" phase vs 15 time decimation bin-shift");
    gr4->GetXaxis()->SetTitle("binshift");
    gr4->GetYaxis()->SetTitle("phi");
    gr4->SetMarkerStyle(20);
    gr4->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    gr4->Draw();

      gr5=new TGraphErrors(m,n,A_cbo,0,dA_cbo);
    gr5->SetTitle(" A_cbo vs 15 time decimation bin-shift");
    gr5->GetXaxis()->SetTitle("binshift");
    gr5->GetYaxis()->SetTitle("A_cbo");
    gr5->SetMarkerStyle(20);
    gr5->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
   

  
    gr6=new TGraphErrors(m,n,w_cbo,0,dw_cbo);
    gr6->SetTitle(" omega_cbo vs 15 time decimation bin-shift");
    gr6->GetXaxis()->SetTitle("binshift");
    gr6->GetYaxis()->SetTitle("w_cbo");
    gr6->SetMarkerStyle(20);
    gr6->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);

    gr7=new TGraphErrors(m,n,Tau_cbo,0,dTau_cbo);
    gr7->SetTitle(" Tau_cbo vs 15 time decimation bin-shift");
    gr7->GetXaxis()->SetTitle("binshift");
    gr7->GetYaxis()->SetTitle("Tau_cbo");
    gr7->SetMarkerStyle(20);
    gr7->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);

    gr8=new TGraphErrors(m,n,phi_cbo,0,dphi_cbo);
    gr8->SetTitle(" phase_cbo vs 15 time decimation bin-shift");
    gr8->GetXaxis()->SetTitle("binshift");
    gr8->GetYaxis()->SetTitle("phi_cbo");
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
   

    for(int i=1;i<=50;i++)
      {
	cout<<" i is "<<i<<" binshift is "<<n[i]<<" Tau cbo is "<<Tau_cbo[i]<<endl;
	if(Tau_cbo[i] < 100000 || Tau_cbo[i] > 400000 || dTau_cbo[i]>400000)
	  {
	    cout<<"Removing point "<<i<<"cbo parameters"<<endl;
	    gr5->RemovePoint(i);
	    gr6->RemovePoint(i);
	    gr7->RemovePoint(i);
	    gr8->RemovePoint(i);
	  }
      }

    c5=new TCanvas("c5","cbo fit parameters vs caloriemeter #");
    c5->Divide(2,2);
    c5->cd(1);
    gr5->Draw();
    c5->cd(2);
    gr6->Draw();
    c5->cd(3);
    gr7->Draw();
    c5->cd(4);
    gr8->Draw();



    c6=new TCanvas("c6","vw fit parameters vs caloriemeter #");
    c6->Divide(2,2);
    c6->cd(1);
    gr9=new TGraphErrors(m,n,A_vw,0,dA_vw);
    gr9->SetTitle(" A_vw vs 15 time decimation bin-shift");
    gr9->GetXaxis()->SetTitle("binshift");
    gr9->GetYaxis()->SetTitle("A_vw");
    gr9->SetMarkerStyle(20);
    gr9->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    

    
    gr10=new TGraphErrors(m,n,w_vw,0,dw_vw);
    gr10->SetTitle(" omega_vw vs 15 time decimation bin-shift");
    gr10->GetXaxis()->SetTitle("binshift");
    gr10->GetYaxis()->SetTitle("w_vw");
    gr10->SetMarkerStyle(20);
    gr10->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
    

    
    gr11=new TGraphErrors(m,n,Tau_vw,0,dTau_vw);
    gr11->SetTitle(" Tau_vw vs 15 time decimation bin-shift");
    gr11->GetXaxis()->SetTitle("binshift");
    gr11->GetYaxis()->SetTitle("Tau_vw");
    gr11->SetMarkerStyle(20);
    gr11->SetLineColor(kRed);
    // gr1->GetYaxis()->SetRangeUser(-75,-15);
   

    
    gr12=new TGraphErrors(m,n,phi_vw,0,dphi_vw);
    gr12->SetTitle(" phase_vw vs 15 time decimation bin-shift");
    gr12->GetXaxis()->SetTitle("binshift");
    gr12->GetYaxis()->SetTitle("phi_vw");
    gr12->SetMarkerStyle(20);
    gr12->SetLineColor(kRed);
    

    for(int i=1;i<=50;i++)
     {
      cout<<" i is "<<i<<" binshift is "<<n[i]<<" Tau vw is "<<Tau_vw[i]<<endl;
      if(Tau_vw[i]<10000 || Tau_vw[i]>30000 || dTau_vw[i]>30000)
       {
	cout<<"Removing point "<<i<<"vw parameters"<<endl;
	gr9->RemovePoint(i);
	gr10->RemovePoint(i);
	gr11->RemovePoint(i);
	gr12->RemovePoint(i);
	}
      }

    c6=new TCanvas("c6","vw fit parameters vs caloriemeter #");
    c6->Divide(2,2);
    c6->cd(1);
    gr9->Draw();
    c6->cd(2);
    gr10->Draw();
    c6->cd(3);
    gr11->Draw();
    c6->cd(4);
    gr12->Draw();


   
}

