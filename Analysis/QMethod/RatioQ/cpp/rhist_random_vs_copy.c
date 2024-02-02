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


//char root_file_name[128] = "highstat_wigsim_output_statfluc_binintegral_e16.root";
//char root_file_name[128] = "run2g_corrOOF_thresh400.root";
char root_file_name[128] = "lowstat_wigsim_output_copy.root";
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *h_res;
TCanvas *c, *cnew, *cauto;
TH1 *hfft, *hFFTsq, *hAuto;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t T_a_true=4365.000811;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64425;
bool usesetbins=true;
bool useerrorbars=true;
bool constructcopy=true;
Double_t fit_start, fit_stop;
TF1 *fit_func3;

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


void rhist_random_vs_copy()
{
  if(constructcopy)
    {
     _file[1]=TFile::Open(root_file_name);
     _file[1]->GetObject("hwiggle",qHist_1D[1]);
     h[1]=(TH1D*)qHist_1D[1]->Clone();
    }
  else
    {
     _file[1]=TFile::Open(root_file_name);
     _file[1]->GetObject("hwiggle1",qHist_1D[1]);
     h[1]=(TH1D*)qHist_1D[1]->Clone();

      _file[1]->GetObject("hwiggle2",qHist_1D[2]);
     h[2]=(TH1D*)qHist_1D[2]->Clone();

      _file[1]->GetObject("hwiggle3",qHist_1D[3]);
     h[3]=(TH1D*)qHist_1D[3]->Clone();

      _file[1]->GetObject("hwiggle4",qHist_1D[4]);
     h[4]=(TH1D*)qHist_1D[4]->Clone();
   }
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
  
  if (constructcopy)
    {
     for(Int_t i=1; i<=1; i++)
      {
        h_sum->Add(h[i],1); 
      }
    }

   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   
   if(constructcopy)
   {
     h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   }
   else
     {
       for(int i=1;i<=4;i++)
	 {
	  h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
          h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
	 }
     }
   
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

     h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
     h_sum->Scale(0.125);
   

   //create the component histograms of the ratio histograms
   if(constructcopy)
     {
      h1=(TH1D*)h_sum->Clone();
      h2=(TH1D*)h_sum->Clone();
     }
   else
     {
      h1=(TH1D*)h[1]->Clone();
      h2=(TH1D*)h[2]->Clone();
     }
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

       nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
       /*   if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
       */
   //   nbinshift=400;
   cout<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=h_sum->GetNbinsX(); ibin++)             
     {
       if(constructcopy)
	 {
          hp->SetBinContent( ibin, h_sum->GetBinContent( ibin + nbinshift));
          hp->SetBinError(ibin, h_sum->GetBinError( ibin + nbinshift));
	 }
       else
	 {
          hp->SetBinContent( ibin, h[3]->GetBinContent( ibin + nbinshift));
          hp->SetBinError(ibin, h[3]->GetBinError( ibin + nbinshift));
	 }
     }

   for(int ibin=nbinshift; ibin<=h_sum->GetNbinsX(); ibin++)
     {
       if(constructcopy)
	 {
          hm->SetBinContent( ibin, h_sum->GetBinContent( ibin - nbinshift));
          hm->SetBinError( ibin, h_sum->GetBinError(ibin - nbinshift));
	 }
       else
	 {
          hm->SetBinContent( ibin, h[4]->GetBinContent( ibin - nbinshift));
          hm->SetBinError( ibin, h[4]->GetBinError(ibin - nbinshift));
	 }
     }
 

   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*h_sum->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   double blife=exp(-(T_a)/(2*lifetime));
   double deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo histogram sum ratio", "h_ratio fit", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
       h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
     }

   }
   else
     {
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

     }
    //Assign error bars
 if(useerrorbars)
     {
       
      long double t1,t2,t3;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
	t1=4*(h1->GetBinError(ibin))*(hm->GetBinContent(ibin)+hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	t2=-4*(hp->GetBinError(ibin))*(h1->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	t3=-4*(hm->GetBinError(ibin))*(h2->GetBinContent(ibin))/((h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin))*(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
       if(usesetbins)
       {
	h_ratio->SetBinError( ibin, sqrt( (t1*t1) + (t2*t2) + (t3*t3) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	h_num->SetBinError(ibin,  sqrt( (t1*t2) + (t2*t2) + (t3*t3) ) );
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

 // hcalo->Rebin(8);
 // hcalo->Scale(0.125);

  fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
  fit_func3->SetParNames("A", "Blind R", "Phi");
  fit_func3->SetNpx(1000000);

  fit_start=30000;
  fit_stop=300000;
  cnew=new TCanvas("cnew","cnew");

  cnew->Divide(1,3);
  cnew->cd(1);

 
  

  hcalo->Fit("fprec3","RE","",fit_start, fit_stop);
  hcalo->GetXaxis()->SetRangeUser(fit_start, fit_stop);
  hcalo->GetXaxis()->SetTitleSize(0.05);
  hcalo->GetYaxis()->SetTitleSize(0.05);
  hcalo->GetYaxis()->SetTitle("Energy [MeV]");
  hcalo->GetXaxis()->SetTitle("time [ns]");
  hcalo->GetXaxis()->SetLabelSize(0.06);
  hcalo->GetYaxis()->SetLabelSize(0.06);

   h_res= new TH1D("residual histogram ", "residual histogram", hcalo->GetNbinsX(), hcalo->GetBinLowEdge(1), (hcalo->GetBinLowEdge(hcalo->GetNbinsX())+hcalo->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hcalo->GetBinLowEdge(1))/hcalo->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hcalo->GetBinContent(ibin)- fit_func3->Eval( hcalo->GetBinCenter(ibin) ) );
          if(hcalo->GetBinError(ibin)!=0){res=(res/hcalo->GetBinError(ibin));}
          h_res->SetBinContent(ibin, (res)  );
     
        }

	h_res->GetXaxis()->SetTitleSize(0.05);
	h_res->GetYaxis()->SetTitleSize(0.05);
	h_res->GetYaxis()->SetTitle("Energy [MeV]");
	h_res->GetXaxis()->SetTitle("time [ns]");	
      
    cnew->cd(2);
    h_res->GetXaxis()->SetRangeUser(fit_start, fit_stop);
    h_res->GetXaxis()->SetLabelSize(0.06);
    h_res->GetYaxis()->SetLabelSize(0.06);

    h_res->Draw();	
 
       hfft = h_res->FFT(hfft, "MAG");
       hfft->SetLineColor(kBlack);
       hfft->SetBins(hcalo->GetNbinsX(),0,1000/h_res->GetBinWidth(1));
       hfft->GetXaxis()->SetRangeUser(0,1000/(2*h_res->GetBinWidth(1)));
       hfft->GetXaxis()->SetTitleSize(0.05);
       hfft->GetXaxis()->SetTitle("Freq [MHz]");
       hfft->GetYaxis()->SetTitleSize(0.05);
       hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft->GetXaxis()->SetLabelSize(0.06);
       hfft->GetYaxis()->SetLabelSize(0.06);
       hfft->SetTitle("FFT");

       cnew->cd(3);
       hfft->Draw();

       if(constructcopy)
	 {
	   cauto=new TCanvas("auto-correlation", "auto-correlation");
	   hFFTsq = (TH1D*)hfft->Clone();
           hFFTsq->Reset();
           for (int ib = 1; ib <= hFFTsq->GetNbinsX(); ib++)
	     {
	       hFFTsq->SetBinContent( ib, hfft->GetBinContent(ib)*hfft->GetBinContent(ib));
	     }
           hAuto = hFFTsq->FFT(hAuto, "MAG");
           hAuto->SetTitle("auto-correlation");
           hAuto->GetXaxis()->SetLabelSize(0.04);
	   hAuto->GetXaxis()->SetTitle("bins");
	   hAuto->GetXaxis()->SetRangeUser(0,h_res->GetNbinsX()/2);
	   hAuto->SetLineColor(kRed);
	   hAuto->GetYaxis()->SetTitle("FFT Mag [arb. units]");
           hAuto->Draw("HIST");
 
	 }


  }
