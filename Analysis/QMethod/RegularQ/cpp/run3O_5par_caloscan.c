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

//char root_file_name[128] = "run3NO_thresh300_nofbfDQC_wndw_8.root";
//char root_file_name[128] = "run3N_thresh300_nofbfDQC_wndw_8_new.root";
//char root_file_name[128]="run3N_thresh300_nofbfDQC_test.root";
char root_file_name[128]="run3N_part3_partial_4.root";

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *qHist_1D[25], *qHist_2[25], *h[25], *hg[25], *h_sum, *h_sum2, *h_err, *h_err2, *hcalo;
//TH1 *hm;
TH1F *hlm, h0;
TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func;
TH1D *h_res[25];
TH1 *hm[25];
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;
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

 /*double paramToFreq(double blindedValue)
{

  double unblindedR = blindedValue - deltaR;
  return 2 * TMath::Pi() * refFreq * (1 + (unblindedR * precisionR));
}
 */

Double_t fprec(Double_t *x, Double_t *par)
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double R = par[3];

    double phi = par[4];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);
    
    //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R);

    return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi));
}

void run3O_5par_caloscan()

 {
   Double_t chi15[40], chi60[40], life[40],dlife[40],freq[40],dfreq[40],norm[40],dnorm[40],phi[40],dphi[40];
   Double_t n[40];
   Int_t m=1;
   Double_t N;

  
 TFile *_file[25];
 TDirectoryFile *dir[25];

    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");

  _file[1]=TFile::Open(root_file_name);
   for(int i=1;i<=24;i++)
    {
     h[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }

  Double_t sum=0;
  
 for(Int_t k=1; k<=24; k++)
   {
    sum=sum+h[k]->GetBinContent(500);
   }

 cout<<sum<<" "<<endl;
  
  gStyle->SetOptFit(1111);
 
   
  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);


 
   for(Int_t i=1; i<=24; i++)
    {
      h_sum->Add(h[i],1);

    }
   //    h_sum->GetYaxis()->SetRangeUser(-100., 3000000000.);
   //   h_sum->Draw();

   cout<<h_sum->GetBinContent(500)<<" "<<endl;

  double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;

 
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   //  h_sum->Rebin(8);


   
 
   

   for(Int_t i=1; i<=24; i++)
     {
       h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       h[i]->Rebin(8);
       h[i]->GetYaxis()->SetRangeUser(-100., 2000000000);
       h[i]->GetYaxis()->SetTitle("ADC counts");
       h[i]->GetXaxis()->SetTitle("time [ns]");

       
     }
   h_sum->Rebin(8);
   // return;
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;




 c2=new TCanvas("c2","5 parameter wiggle_fit");

 
 
 



  fit_func= new TF1("fprec", fprec,  30000,  309000, 5);
  fit_func->SetParNames("N_0", "Tau", "A", "Blind R", "Phi");

 fit_func->SetNpx(1000000);



  gStyle->SetOptFit(1111);

      fit_start=30000;
      fit_stop=300000;

      //     return;
      for(int icalo=1;icalo<=24;icalo++)
	{
          cout<<"calo "<<icalo<<endl;
	  fit_func->SetParameters(1700000000, 51540*1.25, 0.2, 0.0, 3.14/2);
   
        h[icalo]->GetYaxis()->SetRangeUser(-100., 2000000000.);
	h[icalo]->GetYaxis()->SetTitle("ADC counts");
	h[icalo]->GetXaxis()->SetTitle("time [ns]");
      	//h[1]->Draw();
	//	fit_func->Draw("same");
		 	h[icalo]->Fit("fprec","R","",fit_start,fit_stop);
	//	return;			
		h_res[icalo]= new TH1D("residual histogram", "h_res1", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), (h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)));


        
	  for (int ibin = ((fit_start + 1- h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - h_sum->GetBinLowEdge(1))/h_sum->GetBinWidth(1))+1; ++ibin)
           {
       
      double res =  (h[icalo]->GetBinContent(ibin)- fit_func->Eval( h[icalo]->GetBinCenter(ibin) ) );
      if(h[icalo]->GetBinError(ibin)!=0){res=(res/h[icalo]->GetBinError(ibin));}
      h_res[icalo]->SetBinContent(ibin, (res)  );
     
           }
	//	return;
	h_res[icalo]->GetYaxis()->SetTitle("ADC counts");
	h_res[icalo]->GetXaxis()->SetTitle("time [ns]");	
	//    h_res[1]->Draw();

  
     hm[icalo] = h_res[icalo]->FFT(hm[icalo], "MAG");
     hm[icalo]->SetLineColor(kBlack);
     hm[icalo]->SetBins(h_sum->GetNbinsX(),0,h_sum->GetNbinsX()/((h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))));
     hm[icalo]->GetXaxis()->SetRangeUser(0,h_sum->GetNbinsX()/((h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1)))/2);
     hm[icalo]->GetXaxis()->SetTitle("Freq [GHz]");
     hm[icalo]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
   //hm[1]->Draw();

	}
      c3=new TCanvas("c3","calo residuals");
      c3->Divide(6,4);
      for(int icalo=1;icalo<=24;icalo++)
	{
	 c3->cd(icalo);
	 h_res[icalo]->Draw();
	}

 }
