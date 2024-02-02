#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TFile *f = TFile::Open("RootOutput_35026.root");

TH1D *h[21132], *hsum, *hx[54];
TH1 *hm1, *hm2, *hm3, *hnew;

void daq_test_hist()
{

  hsum=new TH1D("hsum","hsum",54000,0,54000);

 for(int i=1;i<=21132;i++)
   {
     h[i]= (TH1D*)(f->FindObjectAny(TString::Format("qHist_2;%d", i)));
     hsum->Add(h[i],1);
   }
 hsum->Draw();
}

TFile *f = TFile::Open("RootOutput_35033.root");
TH1D *h[21132];
TH1D *h15[20000];
TH1D *hx[54];
TH1D *hx15[54];
TH1D *hsum=new TH1D("hsum","hsum",54000,0,54000);
for(int i=1;i<=18000;i++){h[i]= (TH1D*)(f->FindObjectAny(TString::Format("qHist_2;%d", i)));hsum->Add(h[i],1);}
TH1D *hsum15=new TH1D("hsum15","hsum15",3564,0,3564);
for(int i=1;i<=18000;i++){h15[i]= (TH1D*)(f->FindObjectAny(TString::Format("qHist_0;%d", i)));hsum15->Add(h15[i],1);}
hx[17]=new TH1D("hx","hx",1000,0,1000);
hx15[17]=new TH1D("hx15","hx15",67,0,67);
for(int i=1;i<=1000;i++){hx[17]->SetBinContent(i,hsum->GetBinContent(17*1000+i));}
for(int i=1;i<=67;i++){hx15[17]->SetBinContent(i,hsum15->GetBinContent(17*66+i));}
TH1 *hm15;
hm15=hx15[17]->FFT(hm15,"MAG");
hm15->SetBins(hx15[17]->GetNbinsX(),0,1/(1.25*hx15[17]->GetBinWidth(1)));
hm15->Draw()
hsum->Rebin(15);
hx[17]->Rebin(15);
TH1D *hdiff=new TH1D("hdiff","hdiff",3600,0,3600);
TH1D *hdiff_xtal1=new TH1D("hdiff_xtal1","hdiff_xtal1",67,0,67);
for(int i=1;i<=67;i++){hdiff_xtal1->SetBinContent(i,hx15[17]->GetBinContent(i)-hx[17]->GetBinContent(i));}
for(int i=1;i<=3600;i++){hdiff->SetBinContent(i,hsum15->GetBinContent(i)-hsum->GetBinContent(i));}


for(int i=1;i<=54;i++)
  {
   hx[i]=new TH1D("hx","hx",1000,0,1000);
   for(j=1;j<=1000;j++)
    {
    hx[i]->SetBinContent(j,hsum->GetBinContent(((i-1)*1000)+1+j));
    }
  }




 //hsum->SetBins(54000,0,67500);
 //hsum->Draw();   
 /*   hnew=new TH1D("h26","h26",1080,0,1080);
     
   for(int i=26000;i<=27081;i++)
     {
       hnew->SetBinContent(i-25999,hsum->GetBinContent(i));
     }
 */
 //hsum->Rebin(30);
   
 /*   hm1 = hsum->FFT(hm1, "MAG");
   hm1->SetBins(hsum->GetNbinsX(),0,1/hsum->GetBinWidth(1));
   hm1->GetXaxis()->SetTitle("Freq [GHz]");
   hm1->SetTitle("30 time decimation");
   hm1->Draw();*/
   // hnew->Draw();
   
 //}
