#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

char root_file_name[128] = "run3O_thresh300_nofbfDQC_wndw_8.root";
char root_file_name2[128] = "run3N_thresh300_nofbfDQC_wndw_8.root";
TFile *_file[50], *foutput;
TDirectory *dir[50];
TH1D *h[25],*hr1[50],*hr2[50],*hr3[50],*hr4[50];
TCanvas *c,*c1, *c2,*c3, *c4, *c5, *c6, *c7, *c8, *c9, *c10,*cauto,*cfit;

void negative_spike_test()
{


  _file[1]=TFile::Open(root_file_name);
  _file[2]=TFile::Open(root_file_name2);

    for(int i=1;i<=24;i++)
    {
     cout<<"Calo "<<i<<endl;
      
     hr1[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     hr2[i]= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_sig_%d", i)));

     cout<<hr1[i]->GetBinContent(1000)+hr2[i]->GetBinContent(1000)<<endl;

     hr1[i]->Sumw2(kTRUE);

     if(i!=23){ hr1[i]->Add(hr2[i]);}

     cout<<hr1[i]->GetBinContent(1000)<<endl;

    }

    foutput=new TFile("run3NO_thesh300_nofbfDQC_temp.root","new");
    for(int i=1;i<=24;i++)
    {
     hr1[i]->Write();
    }
    foutput->Close();
    
}
