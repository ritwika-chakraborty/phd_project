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
//char root_file_name[128] = "run2c300MeV.root";
char root_file_name[128] = "run2C_thresh300_calosum_ratio.root";
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr1[25], *hr2[25], *hr3[25], *hr4[25], *hrsum, *h_res, *hsum1, *hsum2, *hsum3, *hsum4, *herror1, *herror2, *herror3, *herror4;
TCanvas *c,*c1, *c2,*c3,*c4, *cauto;
TGraphErrors *gr1,*gr2;
TH1 *hfft, *hFFTsq, *hAuto;
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
Double_t blindR[25],dblindR[25],n[25],phase[25],dphase[25];

void ratio_errorbar_randseed()
{

  
  _file[1]=TFile::Open(root_file_name);
  //  _file[1]->GetObject("QRatioFillByFillAnalyzerDB",dir[1]);
  // dir[1]->GetObject("qHist_1D_sig_sum_1",qHist_1D[1]);
  // _file[1]->GetObject("hwiggle",qHist_1D[1]);
  // hsum1=(TH1D*)qHist_1D[1]->Clone();

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

   hsum1->Rebin(8);
   hsum1->Scale(0.125);
   hsum2->Rebin(8);
   hsum2->Scale(0.125);
   hsum3->Rebin(8);
   hsum3->Scale(0.125);
   hsum4->Rebin(8);
   hsum4->Scale(0.125);

   h1=(TH1D*)hsum1->Clone();
   h2=(TH1D*)hsum2->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

      nbinshift=(0.5*T_a_true)/h1->GetBinWidth(1);
     if(double(0.5*T_a_true/h1->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
   
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
      h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));
       
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

 cout<<hcalo->GetBinWidth(500)<<endl;


   hcalo->GetYaxis()->SetTitle("ADC counts");
   hcalo->GetXaxis()->SetTitle("time [ns]");
   hcalo->GetXaxis()->SetRangeUser(30000,306000);

       

   herror1=new TH1D("calo sum 1", "herror 1", hsum1->GetNbinsX(), hsum1->GetBinLowEdge(1), hsum1->GetBinLowEdge(hsum1->GetNbinsX())+hsum1->GetBinWidth(1));

   herror2=new TH1D("calo sum 2", "herror 2", hsum2->GetNbinsX(), hsum2->GetBinLowEdge(1), hsum2->GetBinLowEdge(hsum2->GetNbinsX())+hsum2->GetBinWidth(1));

   herror3=new TH1D("calo sum 3", "herror 3", hsum3->GetNbinsX(), hsum3->GetBinLowEdge(1), hsum3->GetBinLowEdge(hsum3->GetNbinsX())+hsum3->GetBinWidth(1));

   herror4=new TH1D("calo sum 4", "herror 4", hsum4->GetNbinsX(), hsum4->GetBinLowEdge(1), hsum4->GetBinLowEdge(hsum4->GetNbinsX())+hsum4->GetBinWidth(1));

   for(int ibin=1; ibin<=hcalo->GetNbinsX(); ibin++)
	  {
	    if(ibin>39 && ibin<2099){
	      herror1->SetBinContent(ibin, hsum1->GetBinError(ibin)/hsum1->GetBinContent(ibin));
              herror2->SetBinContent(ibin, hsum2->GetBinError(ibin)/hsum2->GetBinContent(ibin));
	      herror3->SetBinContent(ibin, hsum3->GetBinError(ibin)/hsum3->GetBinContent(ibin));
	      herror4->SetBinContent(ibin, hsum4->GetBinError(ibin)/hsum4->GetBinContent(ibin));
	    }
	    //  cout<<hsum1->GetBinError(ibin)/hsum1->GetBinContent(ibin)<<" "<<ibin<<endl;
	  }
    herror1->Rebin(30);
    herror2->Rebin(30);
    herror3->Rebin(30);
    herror4->Rebin(30);
    
    herror1->Scale(0.0333);
    herror2->Scale(0.0333);
    herror3->Scale(0.0333);
    herror4->Scale(0.0333);
    
    herror1->SetLineColor(kRed);
    herror2->SetLineColor(kGreen);
    herror3->SetLineColor(kMagenta);
    herror4->SetLineColor(kBlack);

    herror1->SetLineWidth(2);
    herror2->SetLineWidth(2);
    herror3->SetLineWidth(2);
    herror4->SetLineWidth(2);

    herror1->Draw("hist");
    herror2->Draw("same hist");
    herror3->Draw("same hist");
    herror4->Draw("same hist");
	  
}
