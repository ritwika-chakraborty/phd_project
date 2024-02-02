#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

char root_file_name1[128] = "run2C_final.root";
char root_file_name2[128] = "run2D_final.root";
char root_file_name3[128] = "run2E_final.root";
char root_file_name4[128] = "run2F_final.root";
char root_file_name5[128] = "run2G_final.root";
char root_file_name6[128] = "run2H_final.root";


TFile *_file[50], *foutput;
TDirectory *dir[50];
TH1D *h[25],*hr1[50],*hr2[50],*hr3[50],*hr4[50], *hr5[50], *hr6[50];
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cauto,*cfit;

void run2_root_file_wo_2F_calo18()
{


  _file[1]=TFile::Open(root_file_name1);
  _file[2]=TFile::Open(root_file_name2);
  _file[3]=TFile::Open(root_file_name3);
  _file[4]=TFile::Open(root_file_name4);
  _file[5]=TFile::Open(root_file_name5);
  _file[6]=TFile::Open(root_file_name6);


    for(int i=1;i<=24;i++)
    {
     cout<<"Calo "<<i<<endl;
      
     hr1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr2[i]= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr3[i]= (TH1D*)(_file[3]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr4[i]= (TH1D*)(_file[4]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr5[i]= (TH1D*)(_file[5]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr6[i]= (TH1D*)(_file[6]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     cout<<hr1[i]->GetBinContent(1000)+hr2[i]->GetBinContent(1000)+hr3[i]->GetBinContent(1000)+hr4[i]->GetBinContent(1000)+hr5[i]->GetBinContent(1000)+hr6[i]->GetBinContent(1000)<<endl;

     hr1[i]->Sumw2(kTRUE);

     hr1[i]->Add(hr2[i]);
     hr1[i]->Add(hr3[i]);
     if(i!=18){ hr1[i]->Add(hr4[i]);}
     hr1[i]->Add(hr5[i]);
     hr1[i]->Add(hr6[i]);

     cout<<hr1[i]->GetBinContent(1000)<<endl;

    }

    foutput=new TFile("run2_final_no2FCalo18.root","new");
    for(int i=1;i<=24;i++)
    {
     hr1[i]->Write();
    }
    foutput->Close();
    
}
