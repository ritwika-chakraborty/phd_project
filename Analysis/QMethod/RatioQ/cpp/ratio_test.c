#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TCanvas *c1, *cfunc;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *hwiggle;
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func, *gaussian, *poisson, *ratio_func;
TH1D *qHist_1D[25], *h[25], *h_sum,*h_ratio, *h_num, *h_deno,*h_res, *hfftsq;
TH1D *hm,*hp,*h1,*h2;
TH1 *hfft,*hfft2;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time= 104800; //in clock-ticks
Int_t nbins = 16800;
Double_t minT = 100001.0;
Double_t maxT = 352001.0;
Double_t T_a = 4500;
Double_t T_a_true=4365.000811;
Double_t fit_start, fit_stop;
Double_t lifetime = 64425;
double deno,flife,blife,lossfrac;




Double_t fprec(Double_t *x, Double_t *par)//The wiggle function
{
    double norm = par[0];

    double life = par[1];

    //  double asym = par[2];

    // double omega_a = par[3];

    //  double phi = par[4];
    
    double time = x[0];

    //return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi));
    // return norm;
     return par[0]*exp(-x[0]/par[1]);
     //    return par[0]*x[0];
}


Double_t fratio_apx(Double_t *x, Double_t *par)//approximated cosine function for fitting
{

    double asym = par[0];

    double omega_a = par[1];

    double phi = par[2];
    
    double time = x[0];

    return asym*cos(omega_a*time + phi);
}

Double_t fratio(Double_t *x, Double_t *par)//full ratio function for fitting
{

  //  Double_t asym = par[0];

  // Double_t omega_a = par[1];

  // Double_t phi = par[2];
    Double_t life = par[0];
    Double_t time = x[0];

    Double_t t1=h1->GetBinLowEdge(h1->FindBin(time));
    Double_t t2=h1->GetBinLowEdge(h1->FindBin(time)+1);

    Double_t tf1=hp->GetBinLowEdge(hp->FindBin(time));
    Double_t tf2=hp->GetBinLowEdge(hp->FindBin(time)+1);

    Double_t tb1=hm->GetBinLowEdge(hm->FindBin(time));
    Double_t tb2=hm->GetBinLowEdge(hm->FindBin(time)+1);

    Double_t fshift = exp(-(T_a/2)/(life));
    Double_t bshift = exp((T_a/2)/(life));
   
    //  Double_t  f=(1+ asym*cos(omega_a*time + phi));
    // Double_t f = time;
    //Double_t f_int=t2*t2-t1*t1;
    Double_t f = exp(-time/life);
    Double_t f_int = exp(-t1/life) - exp(-t2/life);
    
    
    // Double_t ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));
    //Double_t ff = time + T_a/2;
    //Double_t ff_int=(tf2*tf2-tf1*tf1) + (T_a)*(tf2-tf1);
    Double_t ff = bshift*exp(-(time+T_a/2)/life);
    Double_t ff_int = exp(-tf1/life) - exp(-tf2/life);
    //Double_t ff = fshift*exp(-time/life)*bshift;
    
    //Double_t fb=(1+ asym*cos(omega_a*(time - T_a/2) + phi));
    //Double_t fb = time - T_a/2;
    //Double_t fb_int=(tb2*tb2-tb1*tb1) - (T_a)*(tb2-tb1);
    Double_t fb = fshift*exp(-(time-T_a/2)/life);
    Double_t fb_int = exp(-tb1/life) - exp(-tb2/life);
    //Double_t fb = bshift*exp(-time/life)*fshift;
    
    //  return ((2*f - exp(T_a/(2*life))*ff - exp(-T_a/(2*life))*fb)/(2*f + exp(time/(2*life))*ff + exp(-time/(2*life))*fb)) ; 
        return ((2*f_int - ff_int - fb_int)/(2*f_int + ff_int + fb_int)) ;
    //   return  ((2*f - ff - fb)/(2*f + ff + fb));
    //  return ((2*f));
    // return par[0]*x[0];
	return 0;
   
}

void ratio_test()

 {
   double x, y, dy;

  cfunc=new TCanvas("cfunc","5 parameter wiggle function");  
  cfunc->Divide(1,2);
  cfunc->cd(1);
  //Make the wiggle function
  fit_func= new TF1("fprec", fprec,  0,  350000, 2);
  fit_func->SetParNames("N_0", "Tau", "A", "omega_a", "Phi");
  fit_func->SetNpx(10000000);
  //  fit_func->SetParameters(10000000000000000, 51540, 0.2, 0.001798396, 0);
  // fit_func->SetParameter(0,0.00001);
  fit_func->SetParameters(1000000000000000000,64425/1.25);
  //  fit_func->SetParameters(1,64425/1.25);
  fit_func->Draw();
  //  return;
  //Make gaussian and poission ditribution for adding statistical fluctuation
  gaussian = new TF1("gaussian","exp( -x*x/(2.*[0]*[0]) )",-10., 10.);
  gaussian->SetParNames("sigma");
  gaussian->SetParameter( 0, 1.0);
  poisson = new TF1("poisson","TMath::Poisson( x, [0] )", 0.0, 1000.);
  poisson->SetParNames("sigma");
  poisson->SetParameter( 0, 100.0);

  //Fill the wiggle histogram
  hwiggle = new TH1D("hwiggle","hwiggle", nbins, minT, maxT);
    Double_t binwidth= (maxT - minT)/nbins;
   for (Int_t j = 1; j <= hwiggle->GetNbinsX(); j++){
     x = hwiggle->GetBinCenter(j); //Taking the the bincenter value to fill wiggle
     //   y = fit_func->Eval( x );
       y = (fit_func->Integral(x-(binwidth/2),x+(binwidth/2)));
      if (x >= 100.){
	dy = sqrt( y )*gaussian->GetRandom(); //generating gaussian statistical fluctuations
      } else {
	poisson->SetParameter( 0, y);
	dy = poisson->GetRandom(); //generating poisson statistical fluctuation
      }
      hwiggle->Fill( x, y );
      //hwiggle->SetBinError(j,sqrt(y+dy));
    }
   cfunc->cd(2);
   hwiggle->Draw("hist");
    	
   //h_sum is the caloriemeters summed, for simulation, there is just one calo  
  h_sum= new TH1D("calo histogram sum", "h_sum", hwiggle->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
    
  h_sum->Add(hwiggle,1);

   //change the time units from clock ticks to nano seconds
  if(binwidth!=h_sum->GetBinWidth(1))
    {
      cout<<"binwidth discrepancy"<<endl;
      binwidth=h_sum->GetBinWidth(1);
    }
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;
//   h_sum->Draw();
   h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 75ns
   h_sum->Scale(0.125);

   //h_wiggle->Rebin(4);
   //h_wiggle->Scale(0.25);
   
   //  return;
   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

   flife=exp((T_a/2)/(lifetime));
   blife=exp(-(T_a/2)/(lifetime));
   deno=2+flife+blife;

   /*   flife=flife;
   blife=blife;
   
   cout<<(long double) flife<<endl;
   cout<<(long double) blife<<endl;

   if(flife!=1.03434 && blife!=0.966803)
     {
       cout<<flife-1.03434<<endl;
       cout<<blife-0.966803<<endl;
     }

   flife=1.03434-3.47113e-06;
   blife=0.966803+3.29558e-07;

      if(flife!=1.03434 && blife!=0.966803)
     {
       cout<<flife-1.03434<<endl;
       cout<<blife-0.966803<<endl;
     }
   */
   //shift of half of g-2 period corresponds to 29 75ns bins
   for(int ibin=0;ibin<=h_sum->GetNbinsX();ibin++)             
     {
       hp->SetBinContent(ibin,flife*h_sum->GetBinContent(ibin+15));
     }

   for(int ibin=15;ibin<=h_sum->GetNbinsX();ibin++)
     {
       hm->SetBinContent(ibin,blife*h_sum->GetBinContent(ibin-15));
     }
   //  return;
   // assign the correct weights to the 4 histohgrams


   /*  h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   //   return;
   */
   /*   for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
     {
       h1->AddBinContent(ibin,1/deno);
       h2->AddBinContent(ibin,1/deno);
       hp->AddBinContent(ibin,flife/deno);
       hm->AddBinContent(ibin,blife/deno);
     }
   */
   //construct the ratio histogram   
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
     // return;
     //  h_num->Rebin(2);
     //h_num->Scale(0.5);

         h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

     for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
     {
       //   h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
          h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
       //  h_ratio->SetBinContent(ibin,((int)(10000*h1->GetBinContent(ibin))+(int)(10000*h2->GetBinContent(ibin))-(int)(10000*hp->GetBinContent(ibin))-(int)(10000*hm->GetBinContent(ibin)))/((int)(10000*h1->GetBinContent(ibin))+(int)(10000*h2->GetBinContent(ibin))+(int)(10000*hp->GetBinContent(ibin))+(int)(10000*hm->GetBinContent(ibin))));
      
     }
     // return;


     //   h_ratio->Rebin(2);
     //h_ratio->Scale(0.5);
     

     ratio_func= new TF1("fratio", fratio,  30000,  309000, 1);
     //  ratio_func->SetParNames("A", "omega_a", "Phi");
     ratio_func->SetParName(0,"lifetime");
     ratio_func->SetNpx(10000000);
     //  ratio_func->SetNpx(16800/4);
     //  ratio_func->SetParameters(0.2, 0.001438716, 0);
     ratio_func->SetParameter(0,64425);
     ratio_func->FixParameter(0,64425);
     
     gStyle->SetOptFit(1111);

     fit_start=30000;
     fit_stop=300000;

     c2= new TCanvas("c2","c2");
     c2->Divide(1,4);
     c2->cd(1);
     hwiggle->GetXaxis()->SetRangeUser(100001,352001);
     hwiggle->Draw("hist");
     c2->cd(2);


     //fit the ratio histogram
     
	h_num->GetYaxis()->SetTitle("ADC counts");
	h_num->GetXaxis()->SetTitle("time [ns]");
	
	//	h_num->Fit("fratio","R","",fit_start,fit_stop);
	        h_ratio->GetXaxis()->SetRangeUser(30000,300000);
	       	h_ratio->Draw("hist");
	      	ratio_func->Draw("same");
		h_res= new TH1D("residual histogram", "h_res", h_num->GetNbinsX(), h_num->GetBinLowEdge(1), (h_num->GetBinLowEdge(h_num->GetNbinsX())+h_num->GetBinWidth(1)));



       //make the fit residual histogram

		c2->cd(3);
	for (int ibin = ((fit_start + 1- h_num->GetBinLowEdge(1))/h_num->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_num->GetBinLowEdge(1))/h_num->GetBinWidth(1))+1; ++ibin)
     {
       
      double res =  (h_ratio->GetBinContent(ibin)- ratio_func->Eval( h_ratio->GetBinCenter(ibin) ) );
      //if(h_num->GetBinError(ibin)!=0){res=(res/h_num->GetBinError(ibin));}
      h_res->SetBinContent(ibin, (res)  );
     
      }
	h_res->GetYaxis()->SetTitle("ADC counts");
	h_res->GetXaxis()->SetTitle("time [ns]");
	h_res->GetXaxis()->SetRangeUser(30000,300000);
        h_res->Draw();

	
	//Take the FFT
	
    c2->cd(4);
     hfft = h_res->FFT(hfft, "MAG");
   hfft->SetLineColor(kBlack);
   hfft->SetBins(h_res->GetNbinsX(),0,h_res->GetNbinsX()/((h_res->GetBinLowEdge(h_res->GetNbinsX())+h_res->GetBinWidth(1))));
   hfft->GetXaxis()->SetTitle("Freq [GHz]");
   hfft->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   hfft->GetXaxis()->SetRangeUser(0,(hfft->GetXaxis()->GetXmax())/2);
   hfft->Draw();     
     
 }
