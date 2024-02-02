#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TString.h"


//char root_file_name_1[128] = "run2_thresh300_fbfDQC_wndw_8.root";
char root_file_name_1[128] = "run3N_thresh300_nofbfDQC_wndw_8_new_2.root";
//char root_file_name_2[128] = "run3BtoM_thresh300_nofbfDQC_wndw_8_2.root";
char root_file_name_2[128] = "run3O_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_3[128] = "run3NO_thresh300_nofbfDQC_wndw_8_2.root";

TH1D *hr2[25], *hr3a[25], *hr3b[25];

TCanvas *cenergy;

TFile *_file[10];


void energydist_comp()
{

  
    _file[1]=TFile::Open(root_file_name_1);
    _file[2]=TFile::Open(root_file_name_2);
    _file[3]=TFile::Open(root_file_name_3);
    
    
    for(int i=1;i<=24;i++)
    {
     hr2[i]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_ydiff_%d", i)));
     hr2[i]->SetLineColor(kRed);
     hr2[i]->SetLineWidth(3);
     
     hr3a[i]= (TH1D*)(_file[2]->FindObjectAny(TString::Format("qHist1D_ydiff_%d", i)));
     hr3a[i]->Scale(hr2[i]->GetBinContent(223)/hr3a[i]->GetBinContent(223));
     hr3a[i]->SetLineColor(kGreen);
     hr3a[i]->SetLineWidth(3);
     
     hr3b[i]= (TH1D*)(_file[3]->FindObjectAny(TString::Format("qHist1D_ydiff_%d", i)));
     hr3b[i]->Scale(hr2[i]->GetBinContent(223)/hr3b[i]->GetBinContent(223));
     hr3b[i]->SetLineWidth(3);
    }

    cenergy=new TCanvas("cEnergy","cEnergy");
    cenergy->Divide(6,4);
     
    for(int i=1;i<=24;i++)
    {
      cenergy->cd(i);
      hr2[i]->GetXaxis()->SetRangeUser(0,7000);
      hr2[i]->Draw("hist");
      hr3a[i]->GetXaxis()->SetRangeUser(0,7000);
      hr3a[i]->Draw("same hist");
      hr3b[i]->GetXaxis()->SetRangeUser(0,7000);
      hr3b[i]->Draw("same hist");
    }
    

}
