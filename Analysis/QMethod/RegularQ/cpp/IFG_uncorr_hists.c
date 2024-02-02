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

char root_file_name1[128] = "run2_NE_no2FCalo18.root";
char root_file_name2[128] = "run3a_NE.root";
char root_file_name3[128] = "run3b_NE.root";

Double_t refFreq = 0.00143;
Double_t precisionR = 1.0;
Double_t rawBinToNs = 1.25;
Double_t inject_time = 104800;

TH1D *hcorr2[30], *hcorr3a[30], *hcorr3b[30], *hcorr2sum, *hcorr3asum, *hcorr3bsum, *h2IFG, *h3aIFG, *h3bIFG, *huncorr2sum, *huncorr3asum, *huncorr3bsum, *hsig2_amp_mult_[50], *hsig3a_amp_mult_[50], *hsig3b_amp_mult_[50], *hsig2_life_mult_[50], *hsig3a_life_mult_[50], *hsig3b_life_mult_[50];
TF1 *func_IFG, *IFG_amp_mult, *IFG_time_mult;

TCanvas *c1, *c2;


void IFG_uncorr_hists()
{

  TFile *_file[10], *foutput;
  
  char *histname = new char[10];
 

  


  _file[1]=TFile::Open(root_file_name1);
  _file[2]=TFile::Open(root_file_name2);
  _file[3]=TFile::Open(root_file_name3);



    for(int i=1;i<=24;i++)
    {
     hcorr2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
     hcorr3a[i]= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
     hcorr3b[i]= (TH1D*)(_file[3]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));
    }
  

    hcorr2sum= new TH1D("run 2 calo histogram sum", "h2calosum", hcorr2[1]->GetNbinsX(), 100001, 352001);
    hcorr2sum->Sumw2(kTRUE);
    hcorr3asum= new TH1D("run 3a calo histogram sum", "h3acalosum", hcorr3a[1]->GetNbinsX(), 100001, 352001);
    hcorr3asum->Sumw2(kTRUE);
    hcorr3bsum= new TH1D("run 3b calo histogram sum", "h3bcalosum", hcorr3b[1]->GetNbinsX(), 100001, 352001);
    hcorr3bsum->Sumw2(kTRUE);
    
   for(Int_t i=1; i<=24; i++)
    {
     hcorr2sum->Add(hcorr2[i],1);
     hcorr3asum->Add(hcorr3a[i],1);
     hcorr3bsum->Add(hcorr3b[i],1);
    }


   double binwidth=hcorr3asum->GetBinWidth(1);
   int nbins=hcorr3asum->GetNbinsX();
   double binlowedge=hcorr3asum->GetBinLowEdge(1);
   double binhighedge=hcorr3asum->GetBinLowEdge(nbins)+binwidth;

   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   

   hcorr2sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   hcorr3asum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   hcorr3bsum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   h2IFG= new TH1D("run 2 IFG", "h2IFG", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   h3aIFG= new TH1D("run 3a IFG", "h3aIFG", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   h3bIFG= new TH1D("run 3b IFG", "h3bIFG", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
 


  
  
  c1= new TCanvas("IFG","IFG");
  c1->Divide(1,3);

  // ###### create functions for pedestal drift and wiggle ditribution for one crystal
   
  func_IFG = new TF1("func_IFG","1-[0]*exp(-x/[1])",rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
  
  func_IFG->SetParameters(0.01938,8635);
  for(int i=1;i<=h2IFG->GetNbinsX();i++)
    {
     h2IFG->SetBinContent(i,func_IFG->Eval(h2IFG->GetBinCenter(i)));
    }

   func_IFG->SetParameters(0.03868,6788);
  for(int i=1;i<=h3aIFG->GetNbinsX();i++)
    {
     h3aIFG->SetBinContent(i,func_IFG->Eval(h3aIFG->GetBinCenter(i)));
    }

   func_IFG->SetParameters(0.03967,6757);
  for(int i=1;i<=h3bIFG->GetNbinsX();i++)
    {
     h3bIFG->SetBinContent(i,func_IFG->Eval(h3bIFG->GetBinCenter(i)));
    }

  c1->cd(1);
  h2IFG->Draw();
  c1->cd(2);
  h3aIFG->Draw();
  c1->cd(3);
  h3bIFG->Draw();


   huncorr2sum= new TH1D("run 2 uncorr", "huncorr2sum", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   huncorr3asum= new TH1D("run 3a uncorr", "huncorr3asum", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   huncorr3bsum= new TH1D("run 3b uncorr", "huncorr3bsum", nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

    for(int i=1;i<=h2IFG->GetNbinsX();i++)
    {
     huncorr2sum->SetBinContent(i,h2IFG->GetBinContent(i)*hcorr2sum->GetBinContent(i));
     huncorr2sum->SetBinError(i,hcorr2sum->GetBinError(i));
    }
    for(int i=1;i<=h3aIFG->GetNbinsX();i++)
    {
     huncorr3asum->SetBinContent(i,h3aIFG->GetBinContent(i)*hcorr3asum->GetBinContent(i));
     huncorr3asum->SetBinError(i,hcorr3asum->GetBinError(i));
    }
    for(int i=1;i<=h3bIFG->GetNbinsX();i++)
    {
     huncorr3bsum->SetBinContent(i,h3bIFG->GetBinContent(i)*hcorr3bsum->GetBinContent(i));
     huncorr3bsum->SetBinError(i,hcorr3bsum->GetBinError(i));
    }

  c2= new TCanvas("IFG uncorrected wiggle","IFG uncorrected wiggle");
  c2->Divide(1,3);  

  c2->cd(1);
  huncorr2sum->Draw();
  c2->cd(2);
  huncorr3asum->Draw();
  c2->cd(3);
  huncorr3bsum->Draw();

   IFG_amp_mult = new TF1("IFG_amp_mult","1-[0]*[2]*exp(-x/[1])",rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   IFG_time_mult = new TF1("IFG_time_mult","1-[0]*exp(-x/([1]*[2]))",rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));

   double amp=0.0;
   double life=0.05;

   for(int k=1;k<=40;k++)
     {
       IFG_amp_mult->SetParameters(0.01938,8635,amp);
       for(int i=1;i<=h2IFG->GetNbinsX();i++)
	 {
          h2IFG->SetBinContent(i,IFG_amp_mult->Eval(h2IFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig2_amp_mult_%d",k);    
       hsig2_amp_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h2IFG->GetNbinsX();i++)
	 {
	   hsig2_amp_mult_[k]->SetBinContent(i,huncorr2sum->GetBinContent(i)/h2IFG->GetBinContent(i));
	   hsig2_amp_mult_[k]->SetBinError(i,huncorr2sum->GetBinError(i));
	 }
       amp=amp+0.05;
       h2IFG->Reset();
       IFG_time_mult->SetParameters(0.01938,8635,life);
       for(int i=1;i<=h2IFG->GetNbinsX();i++)
	 {
          h2IFG->SetBinContent(i,IFG_time_mult->Eval(h2IFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig2_life_mult_%d",k);    
       hsig2_life_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h2IFG->GetNbinsX();i++)
	 {
	   hsig2_life_mult_[k]->SetBinContent(i,huncorr2sum->GetBinContent(i)/h2IFG->GetBinContent(i));
	   hsig2_life_mult_[k]->SetBinError(i,huncorr2sum->GetBinError(i));
	 }
       life=life+0.05;
       h2IFG->Reset();
     }

   amp=0.0;
   life=0.05;

   for(int k=1;k<=40;k++)
     {
       IFG_amp_mult->SetParameters(0.03868,6788,amp);
       for(int i=1;i<=h3aIFG->GetNbinsX();i++)
	 {
          h3aIFG->SetBinContent(i,IFG_amp_mult->Eval(h3aIFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig3a_amp_mult_%d",k);    
       hsig3a_amp_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h3aIFG->GetNbinsX();i++)
	 {
	   hsig3a_amp_mult_[k]->SetBinContent(i,huncorr3asum->GetBinContent(i)/h3aIFG->GetBinContent(i));
	   hsig3a_amp_mult_[k]->SetBinError(i,huncorr3asum->GetBinError(i));
	 }
       amp=amp+0.05;
       h3aIFG->Reset();
       IFG_time_mult->SetParameters(0.03868,6788,life);
       for(int i=1;i<=h3aIFG->GetNbinsX();i++)
	 {
          h3aIFG->SetBinContent(i,IFG_time_mult->Eval(h3aIFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig3a_life_mult_%d",k);    
       hsig3a_life_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h3aIFG->GetNbinsX();i++)
	 {
	   hsig3a_life_mult_[k]->SetBinContent(i,huncorr3asum->GetBinContent(i)/h3aIFG->GetBinContent(i));
	   hsig3a_life_mult_[k]->SetBinError(i,huncorr3asum->GetBinError(i));
	 }
       life=life+0.05;
       h3aIFG->Reset();
     }

   amp=0.0;
   life=0.05;

   for(int k=1;k<=40;k++)
     {
       IFG_amp_mult->SetParameters(0.03967,6757,amp);
       for(int i=1;i<=h3bIFG->GetNbinsX();i++)
	 {
          h3bIFG->SetBinContent(i,IFG_amp_mult->Eval(h3bIFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig3b_amp_mult_%d",k);    
       hsig3b_amp_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h2IFG->GetNbinsX();i++)
	 {
	   hsig3b_amp_mult_[k]->SetBinContent(i,huncorr3bsum->GetBinContent(i)/h3bIFG->GetBinContent(i));
	   hsig3b_amp_mult_[k]->SetBinError(i,huncorr3bsum->GetBinError(i));
	 }
       amp=amp+0.05;
       h3bIFG->Reset();
       IFG_time_mult->SetParameters(0.03967,6757,life);
       for(int i=1;i<=h3bIFG->GetNbinsX();i++)
	 {
          h3bIFG->SetBinContent(i,IFG_time_mult->Eval(h3bIFG->GetBinCenter(i)));
	 }
       sprintf(histname,"hsig3b_life_mult_%d",k);    
       hsig3b_life_mult_[k]= new TH1D(histname, histname, nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
       for(int i=1;i<=h3bIFG->GetNbinsX();i++)
	 {
	   hsig3b_life_mult_[k]->SetBinContent(i,huncorr3bsum->GetBinContent(i)/h3bIFG->GetBinContent(i));
	   hsig3b_life_mult_[k]->SetBinError(i,huncorr3bsum->GetBinError(i));
	 }
       life=life+0.05;
       h3bIFG->Reset();
     }

   foutput=new TFile("gain_systematics_scan.root","new");
   for(int k=1;k<=40;k++)
     {
       hsig2_amp_mult_[k]->Write();
       hsig2_life_mult_[k]->Write();
       hsig3a_amp_mult_[k]->Write();
       hsig3a_life_mult_[k]->Write();
       hsig3b_amp_mult_[k]->Write();
       hsig3b_life_mult_[k]->Write();
     }
     foutput->Close();



}
