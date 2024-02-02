#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

char root_file_name[128] = "run2d_laser_26066_241.root";
TH1D *h[24][54];

void laser_filter()
{

  TFile *_file,*file_,*f, *fxtal;
  TDirectoryFile *dir, *dirxtal;

  int count_xtal=0;

  _file=TFile::Open(root_file_name);
  _file->GetObject("QFillByFillAnalyzer",dir);
  // dir->GetObject("qHist1D_1_0",qHist_data);
  // hdata=(TH1D*)qHist_data->Clone();

  //TH1 *h;
  //SomeRootFile->GetObject(TString::Format("histo%d", i), h);

  for(int i=0; i<24; i++)
    {
      for(int j=0; j<54; j++)
	{
	  dir->GetObject(TString::Format("qHist1D_sig_xtal_%d_%d", i,j),h[i][j]);
	}
    }

  for(int i=0; i<24; i++)
    {
      for(int j=0; j<54; j++)
	{
	  if(h[i][j]->GetBinContent(694)!=0)
	    {
	      count_xtal++;
	      cout<<"In calo "<<i<<endl;
	      cout<<"xtal "<<j<<endl;
	      
	    }
	}
    }
  cout<<count_xtal<<endl;
}
