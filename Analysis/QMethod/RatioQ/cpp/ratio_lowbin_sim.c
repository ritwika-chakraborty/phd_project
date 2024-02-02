#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TCanvas *c1;
//TF1 *fit_func, *fit_func1, *cbo_func;
TH1D *hwiggle, *hwiggle1, *hwiggle2, *hwiggle3, *hwiggle4;
//TH1 *hm;

TCanvas *c2, *c3, *c4, *c5;
TF1 *fit_func, *gaussian, *poisson;
TH1D *h_res;
TH1 *hm;
TGraphErrors *gr2, *gr1, *gr3, *gr4, *gr5;
Double_t deltaR;
Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t fit_start, fit_stop;
Double_t inject_time= 104800; //in clock-ticks
Int_t nbins = 50;
Double_t minT = 100001.0;
Double_t maxT = 352001.0;
TFile *foutput;



Double_t fprec(Double_t *x, Double_t *par)//The wiggle function
{
    double norm = par[0];

    double life = par[1];

    double asym = par[2];

    double omega_a = par[3];

    double phi = par[4];
    
    double time = x[0];

    //  return norm * exp(-time/life) * (1 + asym*cos(omega_a*time + phi));

    return norm;
}

void ratio_lowbin_sim()

 {
   double x, y, dy;

  c1=new TCanvas("c1","5 parameter wiggle function");  

  //Make the wiggle function
  fit_func= new TF1("fprec", fprec,  0,  350000, 1);
  fit_func->SetParNames("N_0");
  fit_func->SetNpx(10000000);
  fit_func->SetParameter(0,1000000);
  fit_func->Draw();

  //Make gaussian and poission ditribution for adding statistical fluctuation
  gaussian = new TF1("gaussian","exp( -x*x/(2.*[0]*[0]) )",-10., 10.);
  gaussian->SetParNames("sigma");
  gaussian->SetParameter( 0, 1.0);
  poisson = new TF1("poisson","TMath::Poisson( x, [0] )", 0.0, 1000.);
  poisson->SetParNames("sigma");
  poisson->SetParameter( 0, 100.0);

  //Fill the wiggle histogram
  hwiggle = new TH1D("hwiggle","hwiggle", nbins, minT, maxT);
  //  hwiggle1 = new TH1D("hwiggle1","hwiggle1", nbins, minT, maxT);
  //hwiggle2 = new TH1D("hwiggle2","hwiggle2", nbins, minT, maxT);
  //hwiggle3 = new TH1D("hwiggle3","hwiggle3", nbins, minT, maxT);
  //hwiggle4 = new TH1D("hwiggle4","hwiggle4", nbins, minT, maxT);
    Double_t binwidth= (maxT - minT)/nbins;
   for (Int_t j = 1; j <= hwiggle->GetNbinsX(); j++){
     x = hwiggle->GetBinCenter(j); //Taking the the bincenter value to fill wiggle
     //  y=fit_func->Integral(hwiggle->GetBinLowEdge(j),hwiggle->GetBinLowEdge(j+1));
       y = binwidth*fit_func->Eval( x );
      if (x >= 100.){
	dy = sqrt( y )*gaussian->GetRandom(); //generating gaussian statistical fluctuations
      } else {
	poisson->SetParameter( 0, y);
	dy = poisson->GetRandom(); //generating poisson statistical fluctuation
      }
      //  hwiggle->SetBinContent(j,y+dy);
      hwiggle->Fill( x, y+dy );
      hwiggle->SetBinError(j,sqrt(y+dy));

      /*     if (x >= 100.){
	dy = sqrt( y )*gaussian->GetRandom(); //generating gaussian statistical fluctuations
      } else {
	poisson->SetParameter( 0, y);
	dy = poisson->GetRandom(); //generating poisson statistical fluctuation
      }
      //  hwiggle->SetBinContent(j,y+dy);
      hwiggle2->Fill( x, y+dy );
      hwiggle2->SetBinError(j,sqrt(y+dy));

       if (x >= 100.){
	dy = sqrt( y )*gaussian->GetRandom(); //generating gaussian statistical fluctuations
      } else {
	poisson->SetParameter( 0, y);
	dy = poisson->GetRandom(); //generating poisson statistical fluctuation
      }
      //  hwiggle->SetBinContent(j,y+dy);
      hwiggle3->Fill( x, y+dy );
      hwiggle3->SetBinError(j,sqrt(y+dy));


       if (x >= 100.){
	dy = sqrt( y )*gaussian->GetRandom(); //generating gaussian statistical fluctuations
      } else {
	poisson->SetParameter( 0, y);
	dy = poisson->GetRandom(); //generating poisson statistical fluctuation
      }
      //  hwiggle->SetBinContent(j,y+dy);
      hwiggle4->Fill( x, y+dy );
      hwiggle4->SetBinError(j,sqrt(y+dy));
      */ 

    }

   //  hwiggle->Rebin(15);
  /*  hwiggle->FillRandom("fprec",1000000000);
  hwiggle->Draw();
  */
  // hwiggle->Draw();
   //Write the wiggle histogram to root outout file
   foutput=new TFile("lowstat_ratio_straight_sim_output_copy.root","new");
   hwiggle->Write();
   // hwiggle2->Write();
   //hwiggle3->Write();
   //hwiggle4->Write();
   foutput->Close();
   
    	
       
 }
