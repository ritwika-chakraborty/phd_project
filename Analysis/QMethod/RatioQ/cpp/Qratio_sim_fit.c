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


char root_file_name[128] = "highstat_wigsim_output_copy.root";
TH1D *qHist_1D[25],*h[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno,*h_res, *hfftsq, *hf_modul, *hb_modul;
TH1 *hfft,*hfft2;
TH1F *hlm;
TF1 *fit_func, *fit_func2, *f_modul, *b_modul;
TCanvas *c2;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t T_a = 4500;
Double_t T_a_true=4365.000811;
Double_t fit_start, fit_stop;
Double_t lifetime = 64425;
Double_t deno,flife,blife,lossfrac;

Double_t fprec(Double_t *x, Double_t *par)//approximated cosine function for fitting
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    return asym*cos(omega_a*time + phi);
}

Double_t fprec2(Double_t *x, Double_t *par)//full ratio function for fitting
{

    Double_t asym = par[0];

    Double_t R = par[1];

    Double_t phi = par[2];

     Double_t life = par[3];

    //  Double_t offset= par[3];
    //Double_t phi_p = par[3];

    //  Double_t phi_m = par[4];

    //  Double_t off= par[3];
    
    Double_t time = x[0];

    //  Double_t t1=h1->GetBinLowEdge(h1->FindBin(time));
    // Double_t t2=h1->GetBinLowEdge(h1->FindBin(time)+1);


    Double_t omega_a =1.e-3 * getBlinded.paramToFreq(R);
    //  Double_t omega_a = R;
    Double_t  f=(1+ asym*cos(omega_a*time + phi));
    //   Double_t f_int = (  exp(-t1/life) - exp(-t2/life)  +  (asym/(1+(life*life*omega_a*omega_a)))  *  ( exp(-t1/life) * (cos(omega_a*t1+phi) - life*omega_a*sin(omega_a*t1+phi)) + exp(-t2/life) * (-cos(omega_a*t2+phi) + life*omega_a*sin(omega_a*t2+phi))));
    //  Double_t f_int=(t2-t1)+(asym/omega_a)*(sin(omega_a*t2+phi)-sin(omega_a*t1+phi));
   
      
    Double_t ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));
    //    Double_t ff_int=(exp(-t1/life)-exp(-t2/life)+(asym/(1+life*life*omega_a*omega_a))*(exp(-t1/life)*(cos(omega_a*t1+phi+(omega_a*T_a)/2)-life*omega_a*sin(omega_a*t1+phi+(omega_a*T_a)/2))+exp(-t2/life)*(-cos(omega_a*t2+phi+(omega_a*T_a)/2)+life*omega_a*sin(omega_a*t2+phi+(omega_a*T_a)/2))));
    //  Double_t ff_int=(t2-t1)+(asym/omega_a)*(sin(omega_a*t2+phi+omega_a*T_a/2)-sin(omega_a*t1+phi+omega_a*T_a/2));
 
    Double_t fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));
    //   Double_t fb_int=(exp(-t1/life)-exp(-t2/life)+(asym/(1+life*life*omega_a*omega_a))*(exp(-t1/life)*(cos(omega_a*t1+phi-(omega_a*T_a)/2)-life*omega_a*sin(omega_a*t1+phi-(omega_a*T_a)/2))+exp(-t2/life)*(-cos(omega_a*t2+phi-(omega_a*T_a)/2)+life*omega_a*sin(omega_a*t2+phi-(omega_a*T_a)/2))));
    //  Double_t fb_int=(t2-t1)+(asym/omega_a)*(sin(omega_a*t2+phi-omega_a*T_a/2)-sin(omega_a*t1+phi-omega_a*T_a/2));

    //      return ((2*f_int - ff_int - fb_int)/(2*f_int + ff_int + fb_int));
       return ((2*f - ff - fb)/(2*f + ff + fb)) ;
    
   
}

void Qratio_sim_fit()

 {
   Double_t chi[40],life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

 c2=new TCanvas("c2","ratio wiggle_fit");  
 TFile *_file[25];
 TDirectoryFile *dir[25];

 //obtain the simulated wiggle histogram from root file
   _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("hwiggle",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();
 


  
  gStyle->SetOptFit(1111);
 
  //h_sum is the caloriemeters summed, for simulation, there is just one calo  
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=1; i++)
    {
      h_sum->Add(h[i],1); 
    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   h_sum->Draw();


   //change the time units from clock ticks to nano seconds
   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;
//   h_sum->Draw();
   h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 75ns
   h_sum->Scale(0.125);
   

   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
	h2=(TH1D*)h_sum->Clone();
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
  
   /*  hf_modul=new TH1D("calo histogram forward modulation", "h_fmodul", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   hb_modul=new TH1D("calo histogram backward modulation", "h_bmodul", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

   
    f_modul= new TF1("helperfunc1", helperfunc1,  30000,  309000, 3);
    f_modul->SetNpx(10000000);
    f_modul->SetParameters(0.2, 0.001438716);

    b_modul= new TF1("helperfunc2", helperfunc2,  30000,  309000, 3);
    b_modul->SetNpx(10000000);
    b_modul->SetParameters(0.2, 0.001438716);
    // return;

      for(int ibin=0;ibin<=h_sum->GetNbinsX();ibin++)             
     {
       hf_modul->SetBinContent(ibin,f_modul->Eval(hf_modul->GetBinCenter(ibin)));
       hb_modul->SetBinContent(ibin,b_modul->Eval(hb_modul->GetBinCenter(ibin)));
     }

   */
   
   //shift of half of g-2 period corresponds to 29 75ns bins
   for(int ibin=0;ibin<=h_sum->GetNbinsX();ibin++)             
     {
       hp->SetBinContent(ibin,hp->GetBinContent(ibin+15));
     }

   for(int ibin=15;ibin<=h_sum->GetNbinsX();ibin++)
     {
       hm->SetBinContent(ibin,hm->GetBinContent(ibin-15));
     }

   // assign the correct weights to the 4 histohgrams
   flife=exp((T_a)/(2*lifetime));
   blife=exp(-(T_a)/(2*lifetime));
   deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);

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

     //     h_num->Rebin(2);
     //h_num->Scale(0.5);


     //  hp->Divide(hb_modul);
     //  hm->Divide(hf_modul);
   
     //   return;
       h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
     for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
      
     }

     
     /*           double t1,t2,t3,t4,t5,t6;
   for(int ibin=1;ibin<=h_ratio->GetNbinsX();ibin++)
     {
       t1=4*(h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(h_sum->GetBinError(2*ibin-1))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t2=(-4)*(h_sum->GetBinContent(2*ibin-1))*(h_sum->GetBinError(2*ibin-30))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t3=(-4)*(h_sum->GetBinContent(2*ibin-1))*(h_sum->GetBinError(2*ibin+28))/((2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28))*(2*h_sum->GetBinContent(2*ibin-1)+h_sum->GetBinContent(2*ibin-30)+h_sum->GetBinContent(2*ibin+28)));
       t4=4*(h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(h_sum->GetBinError(2*ibin))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t5=(-4)*(h_sum->GetBinContent(2*ibin))*(h_sum->GetBinError(2*ibin-29))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t6=(-4)*(h_sum->GetBinContent(2*ibin))*(h_sum->GetBinError(2*ibin+29))/((2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29))*(2*h_sum->GetBinContent(2*ibin)+h_sum->GetBinContent(2*ibin-29)+h_sum->GetBinContent(2*ibin+29)));
       t1=t1/2;
       t2=t2/2;
       t3=t3/2;
       t4=t4/2;
       t5=t5/2;
       t6=t6/2;
       h_ratio->SetBinError(ibin,sqrt((t1*t1)+(t2*t2)+(t3*t3)+(t4*t4)+(t5*t5)+(t6*t6)));
     }
     
   cout<<t1<<"  "<<t2<<"  "<<t3<<"  "<<t4<<"  "<<t5<<"  "<<t6<<endl;
   return;
     */ 
     //  h_ratio->Rebin(2);
     //h_ratio->Scale(0.5);
   
     

   //construct the ratio histogram   
 
     for(int ibin=1;ibin<=h_ratio->GetNbinsX();ibin++)
       {
	 h_ratio->SetBinError(ibin,h_num->GetBinError(ibin)/1.7);
       }

     /*         fit_func= new TF1("fprec", fprec,  30000,  309000, 3);
     fit_func->SetParNames("A", "Blind R", "Phi");
     fit_func->SetNpx(10000000);
     fit_func->SetParameters(0.2, 0.0, 0);
     */
     fit_func2= new TF1("fprec2", fprec2,  30000,  309000, 4);
     fit_func2->SetParNames("A","Blind R", "Phi","life");
     fit_func2->SetNpx(10000000);
     fit_func2->SetParameters(0.2, 0, 0, 64425);
     fit_func2->FixParameter(3,64425);
     // fit_func2->Draw();
     // return;
     //  return;
     //  gStyle->SetOptFit(1111);

     fit_start=30000;
     fit_stop=300000;

     // fit_func->SetParameters(0.2, 0.0, 3.14/2);
     c2->Divide(1,4);
     c2->cd(1);
     h[1]->Draw();
     c2->cd(2);


     //fit the ratio histogram
     //   h_ratio->Rebin(2);
     //  h_ratio->Scale(0.5);
	h_ratio->GetYaxis()->SetTitle("ADC counts");
	h_ratio->GetXaxis()->SetTitle("time [ns]");
       	// fit_func->SetParLimits(8, -3.14, 3.14);
	 //  fit_func->SetParLimits(10, -3.14, 3.14);
	//	h_num->Draw("hist");
        //fit_func2->Draw("same");
	
	  	h_ratio->Fit("fprec2","R","",fit_start,fit_stop);

	//	h_ratio->Draw();
       	//	fit_func2->Draw("same");
		h_res= new TH1D("residual histogram", "h_res", h_ratio->GetNbinsX(), h_ratio->GetBinLowEdge(1), (h_ratio->GetBinLowEdge(h_ratio->GetNbinsX())+h_ratio->GetBinWidth(1)));



       //make the fit residual histogram

		c2->cd(3);
	for (int ibin = ((fit_start + 1- h_ratio->GetBinLowEdge(1))/h_ratio->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_ratio->GetBinLowEdge(1))/h_ratio->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_ratio->GetBinContent(ibin)- fit_func2->Eval( h_ratio->GetBinCenter(ibin) ) );
      //if(h_num->GetBinError(ibin)!=0){res=(res/h_num->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");	
        h_res->Draw();

	
	//Take the FFT
	
  c2->cd(4);
     hfft = h_res->FFT(hfft, "MAG");
   hfft->SetLineColor(kBlack);
   hfft->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hfft->GetXaxis()->SetRangeUser(0,(hfft->GetXaxis()->GetXmax())/2);
   hfft->GetXaxis()->SetTitle("Freq [GHz]");
   hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft->Draw();     
        
 }
