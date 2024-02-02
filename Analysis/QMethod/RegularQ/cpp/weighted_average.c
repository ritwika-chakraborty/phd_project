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

//char root_file_name[128] = "r2_w8_calos.root";
//char root_file_name[128] = "r3no_w8_xtalreject.root";
/*char root_file_name_0[128] = "run3B_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_1[128] = "run3C_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_2[128] = "run3D_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_3[128] = "run3E_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_4[128] = "run3F_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_5[128] = "run3G_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_6[128] = "run3I_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_7[128] = "run3J_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_8[128] = "run3K_thresh300_fbfDQC_wndw_8_new.root";
char root_file_name_9[128] = "run3L_thresh300_nofbfDQC_wndw_8.root";
char root_file_name_10[128] = "run3M_thresh300_nofbfDQC_wndw_8_new.root";
char root_file_name_11[128] = "run3N_thresh300_nofbfDQC_wndw_8_new_2.root";
char root_file_name_12[128] = "run3O_thresh300_nofbfDQC_wndw_8.root";
*/

char root_file_name_0[128] = "r3b_w8_xtalreject.root";
char root_file_name_1[128] = "r3c_w8_xtalreject.root";
char root_file_name_2[128] = "r3d_w8_xtalreject.root";
char root_file_name_3[128] = "r3e_w8_xtalreject.root";
char root_file_name_4[128] = "r3f_w8_xtalreject.root";
char root_file_name_5[128] = "r3g_w8_xtalreject.root";
char root_file_name_6[128] = "r3ij_w8_xtalreject.root";
//char root_file_name_7[128] = "r3ij_w8_xtalreject.root";
char root_file_name_7[128] = "r3k_w8_xtalreject.root";
char root_file_name_8[128] = "r3l_w8_xtalreject.root";
char root_file_name_9[128] = "r3m_w8_xtalreject.root";
char root_file_name_10[128] = "r3n_w8_xtalreject.root";
char root_file_name_11[128] = "r3o_w8_xtalreject.root";

TH1D *h[20];
Double_t weight[20],frac_weight[20];
Double_t sum_weight=0.0;
TFile *_file[50];

void weighted_average()
{
   _file[1]=TFile::Open(root_file_name_0);
    _file[2]=TFile::Open(root_file_name_1);
    _file[3]=TFile::Open(root_file_name_2);
    _file[4]=TFile::Open(root_file_name_3);
    _file[5]=TFile::Open(root_file_name_4);
    _file[6]=TFile::Open(root_file_name_5);
    _file[7]=TFile::Open(root_file_name_6);
    _file[8]=TFile::Open(root_file_name_7);
    _file[9]=TFile::Open(root_file_name_8);
    _file[10]=TFile::Open(root_file_name_9);
    _file[11]=TFile::Open(root_file_name_10);
    _file[12]=TFile::Open(root_file_name_11);
    // _file[13]=TFile::Open(root_file_name_12);

    
    for(int idataset=1;idataset<=10;idataset++)
      {
       h[idataset]= (TH1D*)(_file[idataset]->FindObjectAny(TString::Format("Hseqnumstats")));
       weight[idataset]=h[idataset]->Integral();
       sum_weight=sum_weight+weight[idataset];
      }
      for(int idataset=1;idataset<=11;idataset++)
      {
       frac_weight[idataset]=weight[idataset]/sum_weight;
      }

    
}
