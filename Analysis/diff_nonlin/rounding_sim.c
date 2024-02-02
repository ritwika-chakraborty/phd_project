#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

char root_file_name[128] = "run2c_newconfig_4fills.root";
char root_file_name1[128] = "RootOutput_35029.root";
char root_file_name2[128] = "RootOutput_35041.root";

TH1D *hdata, *qHist_data, *hpeddrift[54], *hpedconst, *hnoise, *h[1010], *hsum[253], *hx[1010], *hfillsum, *hdfillsum, *hfill[54], *hdfill[54], *hrebin, *hrebin_true, *qHist_simped, *hsimped, *hsimnoise, *hsimnoise_test, *hwiggle[54], *hwigfunc[54], *h1[253], *h2[253], *h3[253], *h4[253];
TF1 *func_peddrift, *func_pedconst, *func_noise, *func_wiggle;
TCanvas *c1, *c2;

double xtal_ped[54];

void rounding_sim()
{

  TFile *_file,*file_,*f;
 TDirectoryFile *dir;

  _file=TFile::Open(root_file_name);
  _file->GetObject("QFillByFillAnalyzer",dir);
  dir->GetObject("qHist1D_1_0",qHist_data);
  hdata=(TH1D*)qHist_data->Clone();
  
  c1= new TCanvas("function","function");



  file_=TFile::Open(root_file_name1);

  f=TFile::Open(root_file_name2);

  f->GetObject("qHist_2_10",qHist_simped);
  hsimped=(TH1D*)qHist_simped->Clone();

  for(int i=1;i<=15;i++)
    {
     h[i]= (TH1D*)(file_->FindObjectAny(TString::Format("qHist_0;%d", i)));
     h[i]->Scale(0.0667);
     //  h[i]->Scale(0.25);
    }
  

   
  hsimnoise=new TH1D("hsimnoise","hsimnoise",13608000,0,13608000);

  int j=1;
  int p=1;
  
  for(int k=1;k<=907200;k++)
    {
      for(int i=p;i<=13608000;i++)
	{
	  hsimnoise->SetBinContent(i, (h[j]->GetBinContent(k)));
	 j=j+1;
	 p=p+1;
	 if(i%15==0)
	   {
	     j=1;
	     break;
	   }
	}
      
    }
  
  
  for(int i=0;i<=53;i++)
    {
      xtal_ped[i]=0.0;
    }
  
  int m=0; 
  
  for(int i=1; i<=54270; i++)
    {
      xtal_ped[m]+=hsimped->GetBinContent(i);
      if(i%1005==0)
	{
	  xtal_ped[m]=xtal_ped[m]/1005;
	  // cout<<xtal_ped[m]<<endl;
	  m=m+1;
	}
    }

    m=0;

   for(int i=1; i<=13608000; i++)
    {
      hsimnoise->SetBinContent(i, hsimnoise->GetBinContent(i)+ (xtal_ped[m]));
      if(i%252000==0)
	{
	  m=m+1;
	}
    }
   // hsimnoise->Draw();
  
   
   
   
  

  

  
  func_peddrift= new TF1("func_peddrift","[0]*exp(-x/[1])+[2]",100001, 352001);
  func_peddrift->SetParameter(0,-14000000);
  func_peddrift->SetParameter(1,14500);
  func_peddrift->SetParameter(2,-400);

  func_wiggle= new TF1("func_wiggle","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",100001, 352001);
  func_wiggle->SetParameter(0,60000);
  func_wiggle->SetParameter(1,52000);
  func_wiggle->SetParameter(2,0.2);
  func_wiggle->SetParameter(3,0.001432);
  func_wiggle->SetParameter(4,0);
  func_wiggle->SetNpx(100000);
  //func_peddrift->Draw("same");

  /*  func_pedconst= new TF1("func_pedconst","[0]",100001, 352001);
  func_pedconst->SetParameter(0,-1711.4567);
  func_pedconst->SetLineColor(kGreen);
  //  func_pedconst->Draw("same");

  func_noise= new TF1("func_noise","[0]",100001, 352001);
  func_noise->SetParameter(0,0.1);
  func_noise->SetLineColor(kGreen);
  //func_noise->Draw("same");
  */

  char *histname = new char[10];
  char *wigname =new char[10];
  char *wighist= new char[10];

  double v;

  for(int i=0;i<=53;i++)
    {
      /*  sprintf(wighist,"hwigfunc_%d",i);
     hwigfunc[i]= new TH1D(wighist, wighist, 252000, 100001, 352001);
     for(int u=1;u<=252000;u++)
       {
	 hwigfunc[i]->SetBinContent(u, func_wiggle->Eval(hwigfunc[i]->GetBinCenter(u)));
	 }*/
     sprintf(histname,"hpeddrift_%d",i);    
     hpeddrift[i]= new TH1D(histname, histname, 252000, 100001, 352001);
     sprintf(wigname,"hwiggle_%d",i);
     hwiggle[i]= new TH1D(wigname, wigname, 252000, 100001, 352001);
     for(int u=1;u<=10;u++)
       {
	 v=func_wiggle->GetRandom();
         hwiggle[i]->SetBinContent(hwiggle[i]->FindBin(v), func_wiggle->Eval(hwiggle[i]->GetBinCenter(v)));
       }
    }
  //  hpedconst= new TH1D("hpedconst", "hpedconst", 252000, 100001, 352001);
  //hnoise= new TH1D("hnoise", "hnoise", 252000, 100001, 352001);
  //  hfill= new TH1D("hfill", "hfill", 252000, 100001, 352001);
  
  
  for(int n=0;n<=53;n++)
    {
     for(int i=1;i<=252000;i++)
      {
        hpeddrift[n]->SetBinContent(i,func_peddrift->Eval(hpeddrift[n]->GetBinCenter(i)));
      }
    }
  
  for(int i=0;i<=53;i++)
    {
     hpeddrift[i]->Scale(0.06667);
     hpeddrift[i]->Scale(0.0185);
     //  hpeddrift[i]->Scale(0.25);
    }
  
  int l=1;
  
  for(int i=0;i<=53;i++)
    {
      for(int k=1;k<=252000;k++)
	{
	    hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l)+hwiggle[i]->GetBinContent(k));
	    //  hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l));
	  l=l+1;
        }
    }
  
  
  char *newname1 = new char[10];
  char *newname2 = new char[10];
  for(int i=0;i<=53;i++)
    {
     sprintf(newname1,"hfill_%d",i);    
     hfill[i]= new TH1D(newname1, newname1, 16800, 100001, 352001);
     sprintf(newname2,"hdfill_%d",i);    
     hdfill[i]= new TH1D(newname2, newname2, 16800, 100001, 352001);
    }

  int adc=0;
  double fadc=0;
  int k=1;
  int t;
  
  for(int i=0;i<=53;i++)
    {
      t=1;
      k=1;
      while(k<=252000-15+1)
	{
	  for(int r=k; r<k+15; r++)
	    {
	      adc += int (hpeddrift[i]->GetBinContent(r));
	      fadc += int (hpeddrift[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+15;
          hfill[i]->SetBinContent(t,adc);
	  hdfill[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      //   hfill[i]->Scale(4);
      //   hdfill[i]->Scale(4);
      hdfill[i]->SetLineColor(kRed);
    }
  
   hfillsum= new TH1D("hfillsum", "hfillsum", 16800, 100001, 352001);
   hdfillsum= new TH1D("hdfillsum", "hdfillsum", 16800, 100001, 352001);

   for(int i=0;i<=53;i++)
     {
       hfillsum->Add(hfill[i],1);
       hdfillsum->Add(hdfill[i],1);
     }
   hdfillsum->SetLineColor(kRed);
   hdata->SetLineColor(kBlack);
   hfillsum->Draw("hist");
   hdfillsum->Draw("same hist");
   hdata->Draw("same hist");
}
