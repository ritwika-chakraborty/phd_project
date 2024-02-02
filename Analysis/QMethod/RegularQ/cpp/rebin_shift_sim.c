#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

TH1D *h1,*h2,*hm,*hp, *h_ratio;
int nshift_ratio;
int nshift_fr;
void rebin_shift_sim()
{

 h1= new TH1D("h1", "h1", 16800, -5998.75, 309001.25);
 h2= new TH1D("h2", "h2", 16800, -5998.75, 309001.25);
 hp= new TH1D("hm", "hm", 16800, -5998.75, 309001.25);
 hm= new TH1D("hp", "hp", 16800, -5998.75, 309001.25);

  for(int i=1;i<=h1->GetNbinsX();i++)
   {    

    h1->SetBinContent(i,1);
    h2->SetBinContent(i,0.99);
    h2->SetLineColor(kBlack);
    hp->SetBinContent(i,0.98);
    hm->SetBinContent(i,0.97);
   }
  h1->Rebin(4);
  h2->Rebin(4);
  hp->Rebin(4);
  hm->Rebin(4);
  
  nshift_ratio=29;
  
   for(int i=1;i<=h1->GetNbinsX();i++)
   {    
     hp->SetBinContent(i,h1->GetBinContent(i+nshift_ratio));
     hp->SetLineColor(kRed);
     hm->SetBinContent(i,h1->GetBinContent(i-nshift_ratio));
     hm->SetLineColor(kGreen);
   }

 h1->Draw();
 h2->Draw("same");
 hp->Draw("same");
 hm->Draw("same");

 h_ratio= new TH1D("h_ratio", "h_ratio", h1->GetNbinsX(), -5998.75, 309001.25);

   for(int i=1;i<=h1->GetNbinsX();i++)
   {    
     h_ratio->SetBinContent(i,(h1->GetBinContent(i)+h2->GetBinContent(i)-hp->GetBinContent(i)-hm->GetBinContent(i))/(h1->GetBinContent(i)+h2->GetBinContent(i)+hp->GetBinContent(i)+hm->GetBinContent(i)));
     h_ratio->SetLineColor(kMagenta);
   }

   h_ratio->Draw("same");
}
