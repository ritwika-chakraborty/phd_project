#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

//char root_file_name1[128] = "run2losses.root";
//char root_file_name2[128] = "run3losses.root";

char root_file_name1[128] = "run2-Wmod.root";
char root_file_name2[128] = "run3a-Wmod.root";
char root_file_name3[128] = "run3b-Wmod.root";





TFile *_file[50], *foutput;
TDirectory *dir[50];
TH1F *h1,*h2, *h3;
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cauto,*cfit;

void muonloss_setunits()
{


  _file[1]=TFile::Open(root_file_name1);
  _file[1]->GetObject("hI",h1);
   h1->SetBins(h1->GetNbinsX(), h1->GetBinLowEdge(1)*1000,  ((h1->GetBinLowEdge(h1->GetNbinsX()))+h1->GetBinWidth(1))*1000); 
  
  _file[2]=TFile::Open(root_file_name2);
  _file[2]->GetObject("hI",h2);
   h2->SetBins(h2->GetNbinsX(), h2->GetBinLowEdge(1)*1000,  ((h2->GetBinLowEdge(h2->GetNbinsX()))+h2->GetBinWidth(1))*1000);
   
  _file[3]=TFile::Open(root_file_name3);
  _file[3]->GetObject("hI",h3);
   h3->SetBins(h3->GetNbinsX(), h3->GetBinLowEdge(1)*1000,  ((h3->GetBinLowEdge(h3->GetNbinsX()))+h3->GetBinWidth(1))*1000); 
 
    
    foutput=new TFile("reconWmuonloss2_mod.root","new");
    h1->Write();
    foutput->Close();

    foutput=new TFile("reconWmuonloss3a_mod.root","new");
    h2->Write();
    foutput->Close();

    foutput=new TFile("reconWmuonloss3b_mod.root","new");
    h3->Write();
    foutput->Close();
    
}
