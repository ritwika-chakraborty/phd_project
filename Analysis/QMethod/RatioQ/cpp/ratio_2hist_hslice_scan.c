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
//char root_file_name[128] = "run2_thresh300_fbfDQC.root";
char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res[25];
TCanvas *c,*c1, *c2,*c3;
TGraphErrors *gr1;
TH1 *hfft[25];
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24, *fit_func27;
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
Double_t blindR[25],dblindR[25],n[25];

Double_t fprec3(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));


    return (f - ff)/(f + ff);

}


Double_t fprec7(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_vw= par[3];

    double tau_vw = par[4];

    double omega_vw = par[5];

    double phi_vw = par[6];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));


    double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));



    return (f*Nvw - ff*Nvwf )/(f*Nvw + ff*Nvwf );

}


Double_t fprec11(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_vw= par[3];

    double tau_vw = par[4];

    double omega_vw = par[5];

    double phi_vw = par[6];

    double asym_vbo = par[7];

    double tau_vbo= par[8];

    double omega_vbo= par[9];

    double phi_vbo=par[10];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

 

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));

 
   
    double  f=(1+ asym*cos(omega_a*time + phi));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));


    return (f*Nvw*Nvbo - ff*Nvwf*Nvbof)/(f*Nvw*Nvbo + ff*Nvwf*Nvbof);

}


Double_t fprec15(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_vw= par[3];

    double tau_vw = par[4];

    double omega_vw = par[5];

    double phi_vw = par[6];

    double asym_vbo = par[7];

    double tau_vbo= par[8];

    double omega_vbo= par[9];

    double phi_vbo=par[10];

    double asym_cbo= par[11];

    double tau_cbo = par[12];

    double omega_cbo = par[13];

    double phi_cbo = par[14];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));



    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));



    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));


    double  f=(1+ asym*cos(omega_a*time + phi));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));



    return (f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof )/(f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof );

}


Double_t fprec19(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_vw= par[3];

    double tau_vw = par[4];

    double omega_vw = par[5];

    double phi_vw = par[6];

    double asym_vbo= par[7];

    double tau_vbo = par[8];

    double omega_vbo = par[9];

    double phi_vbo = par[10];

     double asym_cbo= par[11];

    double tau_cbo = par[12];

    double omega_cbo = par[13];

    double phi_cbo = par[14];

    double asym_cbo_A = par[15];

    double phi_cbo_A= par[16];

    double A_cbo_phi= par[17];

    double phi_cbo_phi=par[18];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

 
    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

   
    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

  

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));



    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));

 

    double  f=(1+ asym*Acbo*cos(omega_a*time + phi + phicbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*cos(omega_a*(time + T_a/2) + phi + phicbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));



    return (f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof)/(f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof);

}


Double_t fprec23(Double_t *x, Double_t *par)
{

    double asym = par[0];

    double R = par[1];

    double phi = par[2];

    double asym_vw= par[3];

    double tau_vw = par[4];

    double omega_vw = par[5];

    double phi_vw = par[6];

    double asym_vbo= par[7];

    double tau_vbo = par[8];

    double omega_vbo = par[9];

    double phi_vbo = par[10];

    double asym_cbo= par[11];

    double tau_cbo = par[12];

    double omega_cbo = par[13];

    double phi_cbo = par[14];

    double asym_cbo_A = par[15];

    double phi_cbo_A= par[16];

    double A_cbo_phi= par[17];

    double phi_cbo_phi=par[18];

    double A2_vbo = par[19];

    double phi2_vbo = par[20];

    double A3_vbo = par[21];

    double phi3_vbo = par[22];
    
    double time = x[0];

    double omega_a =1.e-3 * getBlinded.paramToFreq(R);

    double Ncbo=(1+ asym_cbo*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo));

    double Ncbof=(1+ asym_cbo*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo));

    double Acbo=(1+ asym_cbo_A*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_A));

    double Acbof=(1+ asym_cbo_A*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_A));

    double phicbo=(A_cbo_phi*exp(-time/tau_cbo)*cos(omega_cbo*time + phi_cbo_phi));

    double phicbof=(A_cbo_phi*exp(-(time + T_a/2)/tau_cbo)*cos(omega_cbo*(time + T_a/2) + phi_cbo_phi));

    double Nvw=(1+ asym_vw*exp(-time/tau_vw)*cos(omega_vw*time + phi_vw));

    double Nvwf=(1+ asym_vw*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi_vw));

    double Nvbo=(1+ asym_vbo*exp(-time/tau_vbo)*cos(omega_vbo*time + phi_vbo));

    double Nvbof=(1+ asym_vbo*exp(-(time + T_a/2)/tau_vbo)*cos(omega_vbo*(time + T_a/2) + phi_vbo));

    double Avbo=(1+ A2_vbo*exp(-time/tau_vw)*cos(omega_vw*time + phi2_vbo));

    double Avbof=(1+ A2_vbo*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi2_vbo));

    double phivbo=(A3_vbo*exp(-time/tau_vw)*cos(omega_vw*time + phi3_vbo));

    double phivbof=(A3_vbo*exp(-(time + T_a/2)/tau_vw)*cos(omega_vw*(time + T_a/2) + phi3_vbo));

    double  f=(1+ asym*Acbo*Avbo*cos(omega_a*time + phi + phicbo + phivbo));

    // double  f=(1+ asym*cos(omega_a*time + phi));

    double ff=(1+ asym*Acbof*Avbof*cos(omega_a*(time + T_a/2) + phi + phicbof + phivbof));

    // double ff=(1+ asym*cos(omega_a*(time + T_a/2) + phi));

    return (f*Ncbo*Nvw*Nvbo - ff*Ncbof*Nvwf*Nvbof )/(f*Ncbo*Nvw*Nvbo + ff*Ncbof*Nvwf*Nvbof);

}





TH1D* construct_rhist_copy(TH1D *h_sum)
{
 
  //  h_sum->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
  //  h_sum->Scale(0.125);
   

   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

   nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
      /*   if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   //   nbinshift=14;
   cout<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=h_sum->GetNbinsX(); ibin++)             
     {
       hp->SetBinContent( ibin, h_sum->GetBinContent( ibin + nbinshift));
       hp->SetBinError(ibin, h_sum->GetBinError( ibin + nbinshift));
     }

   for(int ibin=nbinshift; ibin<=h_sum->GetNbinsX(); ibin++)
     {
       hm->SetBinContent( ibin, h_sum->GetBinContent( ibin - nbinshift));
       hm->SetBinError( ibin, h_sum->GetBinError(ibin - nbinshift));
     }
 

   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*h_sum->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   //double blife=exp(-(T_a)/(2*lifetime));
   double deno=1+flife;

   h1->Scale(1/deno);
   //h2->Scale(1/deno);
   hp->Scale(flife/deno);
   // hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo histogram sum ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)-hp->GetBinContent(ibin))/(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo histogram sum numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

      h_num->Add(h1,1);
      // h_num->Add(h2,1);
      h_num->Add(hp,-1);
      // h_num->Add(hm,-1);

      h_deno->Add(h1,1);
      //h_deno->Add(h2,1);
      h_deno->Add(hp,1);
      //h_deno->Add(hm,1);
     
      h_num->Divide(h_deno);

     }
    //Assign error bars
   
 if(useerrorbars)
     {
       
       long double t1,t2;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
       	t1=2*(h1->GetBinError(ibin))*(hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+hp->GetBinContent(ibin))*(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));

	t2=-2*(hp->GetBinError(ibin))*(hp->GetBinContent(ibin))/((h1->GetBinContent(ibin)+hp->GetBinContent(ibin))*(h1->GetBinContent(ibin)+hp->GetBinContent(ibin)));
	


 	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt((t1*t1) + (t2*t2) ));
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( (t1*t1) + (t2*t2) ) );
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

   hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
   hcalo->Scale(0.125);

 return hcalo;
  }

void ratio_2hist_hslice_scan()
{
  _file[1]=TFile::Open(root_file_name);
  _file[1]->GetObject("QFillByFillAnalyzerDB",dir[1]);
   dir[1]->GetObject("qHist1D_sig_hslice_0",qHist_1D[1]);
  // _file[1]->GetObject("hwiggle",qHist_1D[1]);
  h[1]=(TH1D*)qHist_1D[1]->Clone();

    _file[2]=TFile::Open(root_file_name);
  _file[2]->GetObject("QFillByFillAnalyzerDB",dir[2]);
  dir[2]->GetObject("qHist1D_sig_hslice_1",qHist_1D[2]);
  h[2]=(TH1D*)qHist_1D[2]->Clone();

   _file[3]=TFile::Open(root_file_name);
  _file[3]->GetObject("QFillByFillAnalyzerDB",dir[3]);
  dir[3]->GetObject("qHist1D_sig_hslice_2",qHist_1D[3]);
  h[3]=(TH1D*)qHist_1D[3]->Clone();

   _file[4]=TFile::Open(root_file_name);
  _file[4]->GetObject("QFillByFillAnalyzerDB",dir[4]);
  dir[4]->GetObject("qHist1D_sig_hslice_3",qHist_1D[4]);
  h[4]=(TH1D*)qHist_1D[4]->Clone();

   _file[5]=TFile::Open(root_file_name);
  _file[5]->GetObject("QFillByFillAnalyzerDB",dir[5]);
  dir[5]->GetObject("qHist1D_sig_hslice_4",qHist_1D[5]);
  h[5]=(TH1D*)qHist_1D[5]->Clone();

    _file[6]=TFile::Open(root_file_name);
  _file[6]->GetObject("QFillByFillAnalyzerDB",dir[6]);
  dir[6]->GetObject("qHist1D_sig_hslice_5",qHist_1D[6]);
  h[6]=(TH1D*)qHist_1D[6]->Clone();



  h_sum= new TH1D("calo histogram sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
  for(Int_t i=1; i<=6; i++)
    {
      h_sum->Add(h[i],1); 
    }

   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;
   for(int i=1;i<=6;i++)
     {
      h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

  






    _file[0]=TFile::Open("r2c.root");
  _file[0]->GetObject("hI",hlm);
  hlm->SetName("hlm");



     fit_func3= new TF1("fprec3", fprec3,  30000,  309000, 3);
     fit_func3->SetParNames("A", "Blind R", "Phi");
     fit_func3->SetNpx(1000000);

       fit_func7= new TF1("fprec7", fprec7,  30000,  309000, 7);
     fit_func7->SetParNames("A", "Blind R", "Phi");
     fit_func7->SetParName(3,"A_vw");
     fit_func7->SetParName(4,"Tau_vw");
     fit_func7->SetParName(5,"omega_vw");
     fit_func7->SetParName(6,"phi_vw");
     fit_func7->SetNpx(1000000);

      fit_func11= new TF1("fprec11", fprec11,  30000,  309000, 11);
     fit_func11->SetParNames("A", "Blind R", "Phi");
     fit_func11->SetParName(3,"A_vw");
     fit_func11->SetParName(4,"Tau_vw");
     fit_func11->SetParName(5,"omega_vw");
     fit_func11->SetParName(6,"phi_vw");
     fit_func11->SetParName(7,"A_vbo");
     fit_func11->SetParName(8,"Tau_vbo");
     fit_func11->SetParName(9,"omega_vbo");
     fit_func11->SetParName(10,"phi_vbo");
     fit_func11->SetNpx(1000000);

      fit_func15= new TF1("fprec15", fprec15,  30000,  309000, 15);
     fit_func15->SetParNames("A", "Blind R", "Phi");
     fit_func15->SetParName(3,"A_vw");
     fit_func15->SetParName(4,"Tau_vw");
     fit_func15->SetParName(5,"omega_vw");
     fit_func15->SetParName(6,"phi_vw");
     fit_func15->SetParName(7,"A_vbo");
     fit_func15->SetParName(8,"Tau_vbo");
     fit_func15->SetParName(9,"omega_vbo");
     fit_func15->SetParName(10,"phi_vbo");
     fit_func15->SetParName(11,"A_cbo");
     fit_func15->SetParName(12,"Tau_cbo");
     fit_func15->SetParName(13,"omega_cbo");
     fit_func15->SetParName(14,"phi_cbo");
     fit_func15->SetNpx(1000000);

     fit_func19= new TF1("fprec19", fprec19,  30000,  309000, 19);
     fit_func19->SetParNames("A", "Blind R", "Phi");
     fit_func19->SetParName(3,"A_vw");
     fit_func19->SetParName(4,"Tau_vw");
     fit_func19->SetParName(5,"omega_vw");
     fit_func19->SetParName(6,"phi_vw");
     fit_func19->SetParName(7,"A_vbo");
     fit_func19->SetParName(8,"Tau_vbo");
     fit_func19->SetParName(9,"omega_vbo");
     fit_func19->SetParName(10,"phi_vbo");
     fit_func19->SetParName(11,"A_cbo");
     fit_func19->SetParName(12,"Tau_cbo");
     fit_func19->SetParName(13,"omega_cbo");
     fit_func19->SetParName(14,"phi_cbo");
     fit_func19->SetParName(15,"A2_cbo");
     fit_func19->SetParName(16,"phi2_cbo");
     fit_func19->SetParName(17,"A3_cbo");
     fit_func19->SetParName(18,"phi3_cbo");
     fit_func19->SetNpx(1000000);


     fit_func23= new TF1("fprec23", fprec23,  30000,  309000, 23);
     fit_func23->SetParNames("A", "Blind R", "Phi");
     fit_func23->SetParName(3,"A_vw");
     fit_func23->SetParName(4,"Tau_vw");
     fit_func23->SetParName(5,"omega_vw");
     fit_func23->SetParName(6,"phi_vw");
     fit_func23->SetParName(7,"A_vbo");
     fit_func23->SetParName(8,"Tau_vbo");
     fit_func23->SetParName(9,"omega_vbo");
     fit_func23->SetParName(10,"phi_vbo");
     fit_func23->SetParName(11,"A_cbo");
     fit_func23->SetParName(12,"Tau_cbo");
     fit_func23->SetParName(13,"omega_cbo");
     fit_func23->SetParName(14,"phi_cbo");
     fit_func23->SetParName(15,"A2_cbo");
     fit_func23->SetParName(16,"phi2_cbo");
     fit_func23->SetParName(17,"A3_cbo");
     fit_func23->SetParName(18,"phi3_cbo");
     fit_func23->SetParName(19,"A2_vbo");
     fit_func23->SetParName(20,"phi2_vbo");
     fit_func23->SetParName(21,"A3_vbo");
     fit_func23->SetParName(22,"phi3_vbo");
     fit_func23->SetNpx(1000000);


     fit_start=30000;
     fit_stop=300000;

     
     for(int i=1; i<=6;i++)
       {

	cout<<"Calo "<<i<<endl; 
        hr[i] = construct_rhist_copy(h[i]);
	 
	fit_func3->SetParameters(0.2, 0.0, 2.4);

  
	hr[i]->GetYaxis()->SetTitle("ADC counts");
	hr[i]->GetXaxis()->SetTitle("time [ns]");
        hr[i]->GetXaxis()->SetRangeUser(fit_start,fit_stop);



        hr[i]->Fit("fprec3","RE","",fit_start,fit_stop);
 
      
		
    
   
   
        fit_func7->SetParameter(0, fit_func3->GetParameter(0));
        fit_func7->SetParameter(1, fit_func3->GetParameter(1));
        fit_func7->SetParameter(2, fit_func3->GetParameter(2));
        fit_func7->SetParameter(3, 0.2);
        fit_func7->SetParameter(4, 85000);
        fit_func7->SetParameter(5, 0.01392);
        fit_func7->SetParameter(6, -2.1);

   
        hr[i]->Fit("fprec7","RE","",fit_start,fit_stop);
      		
   
   
    
       fit_func11->SetParameter(0, fit_func7->GetParameter(0));
       fit_func11->SetParameter(1, fit_func7->GetParameter(1));
       fit_func11->SetParameter(2, fit_func7->GetParameter(2));
       fit_func11->SetParameter(3, fit_func7->GetParameter(3));
       fit_func11->SetParameter(4, fit_func7->GetParameter(4));
       fit_func11->SetParameter(5, fit_func7->GetParameter(5));
       fit_func11->SetParameter(6, fit_func7->GetParameter(6));
       /*    fit_func11->SetParameter(7, 0.0008);
       fit_func11->SetParameter(8, 0.97);
       fit_func11->SetParameter(9, 0.0002);
       fit_func11->SetParameter(10, -2.6);
       */
       fit_func11->SetParameter(7, 0.05);
       fit_func11->SetParameter(8, 50000);
       fit_func11->SetParameter(9, 0.01404);
       fit_func11->SetParameter(10, -1.3);


    
       hr[i]->Fit("fprec11","RE","",fit_start,fit_stop);
      		
       
   
    
       fit_func15->SetParameter(0, fit_func11->GetParameter(0));
       fit_func15->SetParameter(1, fit_func11->GetParameter(1));
       fit_func15->SetParameter(2, fit_func11->GetParameter(2));
       fit_func15->SetParameter(3, fit_func11->GetParameter(3));
       fit_func15->SetParameter(4, fit_func11->GetParameter(4));
       fit_func15->SetParameter(5, fit_func11->GetParameter(5));
       fit_func15->SetParameter(6, fit_func11->GetParameter(6));
       fit_func15->SetParameter(7, fit_func11->GetParameter(7));
       fit_func15->SetParameter(8, fit_func11->GetParameter(8));
       fit_func15->SetParameter(9, fit_func11->GetParameter(9));
       fit_func15->SetParameter(10, fit_func11->GetParameter(10));
       fit_func15->SetParameter(11, 0.01);
       fit_func15->SetParameter(12, 150000);
       fit_func15->SetParameter(13, 0.00234);
       fit_func15->SetParameter(14, 1.6);

   
      
       hr[i]->Fit("fprec15","RE","",fit_start,fit_stop);
          

       fit_func19->SetParameter(0, fit_func15->GetParameter(0));
       fit_func19->SetParameter(1, fit_func15->GetParameter(1));
       fit_func19->SetParameter(2, fit_func15->GetParameter(2));
       fit_func19->SetParameter(3, fit_func15->GetParameter(3));
       fit_func19->SetParameter(4, fit_func15->GetParameter(4));
       fit_func19->SetParameter(5, fit_func15->GetParameter(5));
       fit_func19->SetParameter(6, fit_func15->GetParameter(6));
       fit_func19->SetParameter(7, fit_func15->GetParameter(7));
       fit_func19->SetParameter(8, fit_func15->GetParameter(8));
       fit_func19->SetParameter(9, fit_func15->GetParameter(9));
       fit_func19->SetParameter(10, fit_func15->GetParameter(10));
       fit_func19->SetParameter(11, fit_func15->GetParameter(11));
       fit_func19->SetParameter(12, fit_func15->GetParameter(12));
       fit_func19->SetParameter(13, fit_func15->GetParameter(13));
       fit_func19->SetParameter(14, fit_func15->GetParameter(14));
       fit_func19->SetParameter(15, 0.003);
       fit_func19->SetParameter(16, 0.7);
       fit_func19->SetParameter(17, 0.003);
       fit_func19->SetParameter(18, 1.1);

      
       hr[i]->Fit("fprec19","RE","",fit_start,fit_stop);

       fit_func23->SetParameter(0, fit_func19->GetParameter(0));
       fit_func23->SetParameter(1, fit_func19->GetParameter(1));
       fit_func23->SetParameter(2, fit_func19->GetParameter(2));
       fit_func23->SetParameter(3, fit_func19->GetParameter(3));
       fit_func23->SetParameter(4, fit_func19->GetParameter(4));
       fit_func23->SetParameter(5, fit_func19->GetParameter(5));
       fit_func23->SetParameter(6, fit_func19->GetParameter(6));
       fit_func23->SetParameter(7, fit_func19->GetParameter(7));
       fit_func23->SetParameter(8, fit_func19->GetParameter(8));
       fit_func23->SetParameter(9, fit_func19->GetParameter(9));
       fit_func23->SetParameter(10, fit_func19->GetParameter(10));
       fit_func23->SetParameter(11, fit_func19->GetParameter(11));
       fit_func23->SetParameter(12, fit_func19->GetParameter(12));
       fit_func23->SetParameter(13, fit_func19->GetParameter(13));
       fit_func23->SetParameter(14, fit_func19->GetParameter(14));
       fit_func23->SetParameter(15, fit_func19->GetParameter(15));
       fit_func23->SetParameter(16, fit_func19->GetParameter(16));
       fit_func23->SetParameter(17, fit_func19->GetParameter(17));
       fit_func23->SetParameter(18, fit_func19->GetParameter(18));
       fit_func23->SetParameter(19, 0.002);
       fit_func23->SetParameter(20, 4);
       fit_func23->SetParameter(21, 0.01);
       fit_func23->SetParameter(22, 0.1);


      
       hr[i]->Fit("fprec23","RE","",fit_start,fit_stop);



     
       
        blindR[m]=fit_func23->GetParameter(2);
       // blindR[m]=fit_func27->GetChisquare()/fit_func27->GetNDF();
        dblindR[m]=fit_func23->GetParError(2);
       n[m]=i;
       m=m+1;
       
       gStyle->SetOptFit(1111);
      		
       h_res[i]= new TH1D("residual histogram ", "h_res", hr[i]->GetNbinsX(), hr[i]->GetBinLowEdge(1), (hr[i]->GetBinLowEdge(hr[i]->GetNbinsX())+hr[i]->GetBinWidth(1)));


    
	for (int ibin = ((fit_start + 1- hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ibin <= ((fit_stop + 1 - hr[i]->GetBinLowEdge(1))/hr[i]->GetBinWidth(1))+1; ++ibin)
        {
       
          double res =  (hr[i]->GetBinContent(ibin)- fit_func23->Eval( hr[i]->GetBinCenter(ibin) ) );
          if(hr[i]->GetBinError(ibin)!=0){res=(res/hr[i]->GetBinError(ibin));}
          h_res[i]->SetBinContent(ibin, (res)  );
     
        }

	h_res[i]->GetYaxis()->SetTitle("Energy [MeV]");
	h_res[i]->GetXaxis()->SetTitle("time [ns]");
	h_res[i]->GetXaxis()->SetRangeUser(fit_start, fit_stop);
        h_res[i]->GetXaxis()->SetLabelSize(0.09);
        h_res[i]->GetYaxis()->SetLabelSize(0.09);

      

 
       hfft[i] = h_res[i]->FFT(hfft[i], "MAG");
       hfft[i]->SetLineColor(kBlack);
       hfft[i]->SetBins(hr[i]->GetNbinsX(),0,1000/h_res[i]->GetBinWidth(1));
       hfft[i]->GetXaxis()->SetRangeUser(0,1000/(2*h_res[i]->GetBinWidth(1)));
       hfft[i]->GetXaxis()->SetTitle("Freq [MHz]");
       hfft[i]->GetYaxis()->SetTitle("FFT Mag [arb. units]");
       hfft[i]->GetXaxis()->SetLabelSize(0.09);
       hfft[i]->GetYaxis()->SetLabelSize(0.09);

 
   }

     c2 = new TCanvas("c2","hslice residuals");
     c2->Divide(1,6);
     c3 = new TCanvas("c3","hslice fft");
     c3->Divide(1,6);
   
     for(int i=1;i<=6;i++)
     {
      c2->cd(i);
      h_res[i]->Rebin(2);
      h_res[i]->Draw();
      c3->cd(i);
      hfft[i]->Draw();
     }
    c1=new TCanvas("c1","Phase vs hslice #");
    gr1=new TGraphErrors(m,n,blindR,0,dblindR);
    gr1->SetTitle("g-2 phase vs horizontal slice");
    gr1->GetXaxis()->SetTitle("hslice #");
    gr1->GetYaxis()->SetTitle("phi");
    gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    //gr1->GetYaxis()->SetRangeUser(-55,-30);
    gr1->Draw();

   
}
