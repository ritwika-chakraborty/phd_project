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

char root_file_name[128] = "run2c_newconfig_4fills.root";
char root_file_name1[128] = "RootOutput_35029.root";
char root_file_name2[128] = "RootOutput_35041.root";
char root_file_name3[128] = "run2d_dqc_26013_1flush_xtals.root";
char root_file_name4[128] = "run2CDEFGH.root";

TH1D *hdata, *hdataxtal, *qHist_data, *qHist_data_xtal, *hpeddrift1[54], *hpeddrift2[54], *hpeddrift3[54], *hpeddrift4[54], *hpedconst, *hnoise, *h[1010], *hsum[253], *hx[1010], *hfillsum[54], *hdfillsum[54], *hfill1[54], *hfill2[54], *hfill3[54], *hfill4[54], *hdfill1[54], *hdfill2[54], *hdfill3[54], *hdfill4[54], *hrebin, *hrebin_true, *qHist_simped, *hsimped, *hsimnoise, *hsimnoise_test, *hwiggle1[54], *hwiggle2[54], *hwiggle3[54], *hwiggle4[54], *hwigfunc[54], *h1[253], *h2[253], *h3[253], *h4[253], *hsimnoise1, *hsimnoise2, *hsimnoise3, *hsimnoise4, *hsig[54], *hdsig[54], *hrawdiff[54], *hdiff[54], *hxtalsum, *hdxtalsum, *hpedsum[54], *hdpedsum[54], *hpeddiff[54], *hxtalsumdiff, *hedist, *qHist_edist, *hnoisedist, *hratio, *hxtalsum0, *hdxtalsum0, *hxtalsumdiff0;
TH1D *h60fill1[54], *h60dfill1[54], *h60fill2[54], *h60dfill2[54],  *h60fill3[54],*h60dfill3[54], *h60fill4[54],*h60dfill4[54], *h60fillsum[54], *h60dfillsum[54], *h60sig[54], *h60dsig[54], *h60xtalsum, *h60dxtalsum, *h60xtalsumdiff;
TF1 *func_peddrift, *func_pedconst, *func_noise, *func_wiggle;
TCanvas *c1, *c2;

double xtal_ped[54];

int nflush=1250;

void rounding_sim_corr()
{

  TFile *_file,*file_,*f, *fxtal, *fedist, *foutput;
  TDirectoryFile *dir, *dirxtal, *diredist;
  
  char *histname = new char[10];
  char *wigname =new char[10];
  char *wighist= new char[10];
  char *newname1 = new char[10];
  char *newname2 = new char[10];
  char *newname3 = new char[10];
  char *newname4 = new char[10];
  char *newname5 = new char[10];
  char *newname6 = new char[10];
  char *newname7 = new char[10];
  char *newname8 = new char[10];
  char *newname9 = new char[10];
  char *newname10 = new char[10];
  char *newname11 = new char[10];


  


  _file=TFile::Open(root_file_name);
  _file->GetObject("QFillByFillAnalyzer",dir);
  dir->GetObject("qHist1D_1_0",qHist_data);
  hdata=(TH1D*)qHist_data->Clone();

   fedist=TFile::Open(root_file_name4);
   fedist->GetObject("QFillByFillAnalyzer",diredist);
   diredist->GetObject("qHist1D_ydiff_1",qHist_edist);
   hedist=(TH1D*)qHist_edist->Clone();


  fxtal=TFile::Open(root_file_name3);
  fxtal->GetObject("QFillByFillAnalyzer",dirxtal);
  dirxtal->GetObject("qHist1D_xtal_9",qHist_data_xtal);
  hdataxtal=(TH1D*)qHist_data_xtal->Clone();
  
  c1= new TCanvas("function","function");

  // ###### create functions for pedestal drift and wiggle ditribution for one crystal
   
  func_peddrift= new TF1("func_peddrift","[0]*exp(-x/[1])+[2]",100001, 352001);
  func_peddrift->SetParameter(0,-7000*2);
  func_peddrift->SetParameter(1,14500);
  func_peddrift->SetParameter(2,-1.5);

  func_wiggle= new TF1("func_wiggle","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",100001, 352001);
  func_wiggle->SetParameter(0,10);
  func_wiggle->SetParameter(1,52000);
  func_wiggle->SetParameter(2,0.2);
  func_wiggle->SetParameter(3,0.001432);
  func_wiggle->SetParameter(4,0);
  func_wiggle->SetNpx(100000);
  //  func_peddrift->Draw();
  //return;



  file_=TFile::Open(root_file_name1);

  f=TFile::Open(root_file_name2);

  f->GetObject("qHist_2_10",qHist_simped);
  hsimped=(TH1D*)qHist_simped->Clone();

 
  // #### get daq test histograms for noise ###
  
  for(int i=1;i<=252;i++)
    {
     h1[i]= (TH1D*)(file_->FindObjectAny(TString::Format("qHist_2;%d", i)));
    }
  for(int i=1;i<=252;i++)
    {
     h2[i]= (TH1D*)(file_->FindObjectAny(TString::Format("qHist_2;%d", i+252)));
    }
  for(int i=1;i<=252;i++)
    {
     h3[i]= (TH1D*)(file_->FindObjectAny(TString::Format("qHist_2;%d", i+504)));
    }
  for(int i=1;i<=252;i++)
    {
     h4[i]= (TH1D*)(file_->FindObjectAny(TString::Format("qHist_2;%d", i+756)));
    }
  
  hsimnoise1=new TH1D("hsimnoise_1","hsimnoise_1",13608000,0,13608000);
   hsimnoise2=new TH1D("hsimnoise_2","hsimnoise_2",13608000,0,13608000);
    hsimnoise3=new TH1D("hsimnoise_3","hsimnoise_3",13608000,0,13608000);
     hsimnoise4=new TH1D("hsimnoise_4","hsimnoise_4",13608000,0,13608000);

    for(int i=0;i<=53;i++)
    {
     sprintf(histname,"hpeddrift1_%d",i);    
     hpeddrift1[i]= new TH1D(histname, histname, 252000, 100001, 352001);
     sprintf(wigname,"hwiggle1_%d",i);
     hwiggle1[i]= new TH1D(wigname, wigname, 252000, 100001, 352001);
     /*sprintf(wighist,"hwigfunc_%d",i);
     hwigfunc[i]= new TH1D(wighist, wighist, 252000, 100001, 352001);*/
     sprintf(histname,"hpeddrift2_%d",i);    
     hpeddrift2[i]= new TH1D(histname, histname, 252000, 100001, 352001);
     sprintf(wigname,"hwiggle2_%d",i);
     hwiggle2[i]= new TH1D(wigname, wigname, 252000, 100001, 352001);
     sprintf(histname,"hpeddrift3_%d",i);    
     hpeddrift3[i]= new TH1D(histname, histname, 252000, 100001, 352001);
     sprintf(wigname,"hwiggle3_%d",i);
     hwiggle3[i]= new TH1D(wigname, wigname, 252000, 100001, 352001);
     sprintf(histname,"hpeddrift4_%d",i);    
     hpeddrift4[i]= new TH1D(histname, histname, 252000, 100001, 352001);
     sprintf(wigname,"hwiggle4_%d",i);
     hwiggle4[i]= new TH1D(wigname, wigname, 252000, 100001, 352001);
     
     sprintf(newname1,"hfill1_%d",i);    
     hfill1[i]= new TH1D(newname1, newname1, 16800, 100001, 352001);
     sprintf(newname2,"hdfill1_%d",i);    
     hdfill1[i]= new TH1D(newname2, newname2, 16800, 100001, 352001);
     sprintf(newname1,"hfill2_%d",i);    
     hfill2[i]= new TH1D(newname1, newname1, 16800, 100001, 352001);
     sprintf(newname2,"hdfill2_%d",i);    
     hdfill2[i]= new TH1D(newname2, newname2, 16800, 100001, 352001);
     sprintf(newname1,"hfill3_%d",i);    
     hfill3[i]= new TH1D(newname1, newname1, 16800, 100001, 352001);
     sprintf(newname2,"hdfill3_%d",i);    
     hdfill3[i]= new TH1D(newname2, newname2, 16800, 100001, 352001);
     sprintf(newname1,"hfill4_%d",i);    
     hfill4[i]= new TH1D(newname1, newname1, 16800, 100001, 352001);
     sprintf(newname2,"hdfill4_%d",i);    
     hdfill4[i]= new TH1D(newname2, newname2, 16800, 100001, 352001);

     sprintf(newname1,"h60_fill1_%d",i);    
     h60fill1[i]= new TH1D(newname1, newname1, 4200, 100001, 352001);
     sprintf(newname2,"h60_dfill1_%d",i);    
     h60dfill1[i]= new TH1D(newname2, newname2, 4200, 100001, 352001);
     sprintf(newname1,"h60_fill2_%d",i);    
     h60fill2[i]= new TH1D(newname1, newname1, 4200, 100001, 352001);
     sprintf(newname2,"h60_dfill2_%d",i);    
     h60dfill2[i]= new TH1D(newname2, newname2, 4200, 100001, 352001);
     sprintf(newname1,"h60_fill3_%d",i);    
     h60fill3[i]= new TH1D(newname1, newname1, 4200, 100001, 352001);
     sprintf(newname2,"h60_dfill3_%d",i);    
     h60dfill3[i]= new TH1D(newname2, newname2, 4200, 100001, 352001);
     sprintf(newname1,"h60_fill4_%d",i);    
     h60fill4[i]= new TH1D(newname1, newname1, 4200, 100001, 352001);
     sprintf(newname2,"h60_dfill4_%d",i);    
     h60dfill4[i]= new TH1D(newname2, newname2, 4200, 100001, 352001);
 
     
     sprintf(newname3,"hfillsum_%d",i);    
     hfillsum[i]= new TH1D(newname3, newname3, 16800, 100001, 352001);
     sprintf(newname4,"hdfillsum_%d",i);    
     hdfillsum[i]= new TH1D(newname4, newname4, 16800, 100001, 352001);

     sprintf(newname3,"h60_fillsum_%d",i);    
     h60fillsum[i]= new TH1D(newname3, newname3, 4200, 100001, 352001);
     sprintf(newname4,"h60_dfillsum_%d",i);    
     h60dfillsum[i]= new TH1D(newname4, newname4, 4200, 100001, 352001);

     
     sprintf(newname5,"hsig_%d",i);    
     hsig[i]= new TH1D(newname5, newname5, 16800, 100001, 352001);
     sprintf(newname6,"hdsig_%d",i);    
     hdsig[i]= new TH1D(newname6, newname6, 16800, 100001, 352001);

     sprintf(newname5,"h60_sig_%d",i);    
     h60sig[i]= new TH1D(newname5, newname5, 4200, 100001, 352001);
     sprintf(newname6,"h60_dsig_%d",i);    
     h60dsig[i]= new TH1D(newname6, newname6, 4200, 100001, 352001);

     
     sprintf(newname9,"hpedsum_%d",i);    
     hpedsum[i]= new TH1D(newname9, newname9, 16800, 100001, 352001);
     sprintf(newname10,"hdpedsum_%d",i);    
     hdpedsum[i]= new TH1D(newname10, newname10, 16800, 100001, 352001);
     
     sprintf(newname7,"hdiff_%d",i);    
     hdiff[i]= new TH1D(newname7, newname7, 16800, 100001, 352001);

     sprintf(newname11,"hpeddiff_%d",i);    
     hpeddiff[i]= new TH1D(newname11, newname11, 16800, 100001, 352001);
     
     sprintf(newname8,"hrawdiff_%d",i);    
     hrawdiff[i]= new TH1D(newname8, newname8, 16800, 100001, 352001);
    }

      hnoisedist= new TH1D("hnoisedist","hnoisedist",30,-15,15);
      hxtalsum= new TH1D("hxtalsum", "hxtalsum", 16800, 100001, 352001);
      hdxtalsum= new TH1D("hdxtalsum", "hdxtalsum", 16800, 100001, 352001);
      hxtalsumdiff= new TH1D("hxtalsumdiff", "hxtalsumdiff", 16800, 100001, 352001);

      h60xtalsum= new TH1D("h60_xtalsum", "h60_xtalsum", 4200, 100001, 352001);
      h60dxtalsum= new TH1D("h60_dxtalsum", "h60_dxtalsum", 4200, 100001, 352001);
      h60xtalsumdiff= new TH1D("h60_xtalsumdiff", "h60_xtalsumdiff", 4200, 100001, 352001);




     // ### fill noise histograms using noise distribution from DAQ test



 for(int i=1;i<=252;i++)
  {  
  for(int itbin=1;itbin<54000;itbin++)
    {
     hnoisedist->Fill(h1[i]->GetBinContent(itbin));
    }
  }


 
//######## BEGIN FLUSH LOOP

cout<<"Begin flush loop"<<endl;
struct timeval t_start, t_end;
gettimeofday(&t_start, NULL);

for(int iflush=1; iflush<=nflush; iflush++) 
{

  for(int i=0;i<=53;i++)
      {
       
      hpeddrift1[i]->Reset();
      hwiggle1[i]->Reset();
      hpeddrift2[i]->Reset();
      hwiggle2[i]->Reset();
      hpeddrift3[i]->Reset();
      hwiggle3[i]->Reset();
      hpeddrift4[i]->Reset();
      hwiggle4[i]->Reset();
         
      hfill1[i]->Reset();
      hdfill1[i]->Reset();
      hfill2[i]->Reset();
      hdfill2[i]->Reset();   
      hfill3[i]->Reset();
      hdfill3[i]->Reset();
      hfill4[i]->Reset();    
      hdfill4[i]->Reset();

      hfillsum[i]->Reset();
      hdfillsum[i]->Reset();
       
      hsig[i]->Reset();   
      hdsig[i]->Reset();

      h60fill1[i]->Reset();
      h60dfill1[i]->Reset();
      h60fill2[i]->Reset();
      h60dfill2[i]->Reset();   
      h60fill3[i]->Reset();
      h60dfill3[i]->Reset();
      h60fill4[i]->Reset();    
      h60dfill4[i]->Reset();

      h60fillsum[i]->Reset();
      h60dfillsum[i]->Reset();
       
      h60sig[i]->Reset();   
      h60dsig[i]->Reset(); 

     
      hpedsum[i]->Reset();
      hdpedsum[i]->Reset(); 
     
      hdiff[i]->Reset();

      hpeddiff[i]->Reset(); 
     
      hrawdiff[i]->Reset(); 
    }
      hsimnoise1->Reset();
      hsimnoise2->Reset();
      hsimnoise3->Reset();
      hsimnoise4->Reset();

  
 for(int i=1; i<=13608000; i++)
   {
  hsimnoise1->SetBinContent(i,int (hnoisedist->GetRandom()));
  hsimnoise2->SetBinContent(i,int (hnoisedist->GetRandom()));
  hsimnoise3->SetBinContent(i,int (hnoisedist->GetRandom()));
  hsimnoise4->SetBinContent(i,int (hnoisedist->GetRandom()));
   }
 
 

 // ### fill noise histograms by knitting together histograms from DAQ test from DAQ test
 
 /*  int j=1;
  int p=1;
  int s=1;
  int c=0;
  while(s<=13608000)
  { 
  for(int i=1;i<=252000;i++)
    {
      hsimnoise1->SetBinContent(s,h1[p]->GetBinContent((c*1000)+j));
      j=j+1;
      s=s+1;
      if(i%1000==0)
	{
	  p=p+1;
	  j=1;
	}
    }
  c=c+1;
  p=1;
}

  j=1;
  p=1;
  s=1;
  c=0;
  while(s<=13608000)
  { 
  for(int i=1;i<=252000;i++)
    {
      hsimnoise2->SetBinContent(s,h2[p]->GetBinContent((c*1000)+j));
      j=j+1;
      s=s+1;
      if(i%1000==0)
	{
	  p=p+1;
	  j=1;
	}
    }
  c=c+1;
  p=1;
}
  
 j=1;
  p=1;
  s=1;
  c=0;
  while(s<=13608000)
  { 
  for(int i=1;i<=252000;i++)
    {
      hsimnoise3->SetBinContent(s,h3[p]->GetBinContent((c*1000)+j));
      j=j+1;
      s=s+1;
      if(i%1000==0)
	{
	  p=p+1;
	  j=1;
	}
    }
  c=c+1;
  p=1;
}
  
 j=1;
  p=1;
  s=1;
  c=0;
  while(s<=13608000)
  { 
  for(int i=1;i<=252000;i++)
    {
      hsimnoise4->SetBinContent(s,h4[p]->GetBinContent((c*1000)+j));
      j=j+1;
      s=s+1;
      if(i%1000==0)
	{
	  p=p+1;
	  j=1;
	}
    }
  c=c+1;
  p=1;
}
 */

 // ###### calculate avaerage pedestal from DAQ test histogram
 
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


  // ##### Add the average pedestal to the noise histograms


  m=0;

   for(int i=1; i<=13608000; i++)
    {
      hsimnoise1->SetBinContent(i, hsimnoise1->GetBinContent(i)+ (xtal_ped[m]));
      if(i%252000==0)
	{
	  m=m+1;
	}
    }
     m=0;

   for(int i=1; i<=13608000; i++)
    {
      hsimnoise2->SetBinContent(i, hsimnoise2->GetBinContent(i)+ (xtal_ped[m]));
      if(i%252000==0)
	{
	  m=m+1;
	}
    }
     m=0;

   for(int i=1; i<=13608000; i++)
    {
      hsimnoise3->SetBinContent(i, hsimnoise3->GetBinContent(i)+ (xtal_ped[m]));
      if(i%252000==0)
	{
	  m=m+1;
	}
    }
     m=0;

   for(int i=1; i<=13608000; i++)
    {
      hsimnoise4->SetBinContent(i, hsimnoise4->GetBinContent(i)+ (xtal_ped[m]));
      if(i%252000==0)
	{
	  m=m+1;
	}
    }



  // ###### making wiggle histograms for one crystal for 4 fills

  

  double v;
  int npulse = 5;

  for(int i=0;i<=53;i++)
    {
      /* 
     for(int u=1;u<=252000;u++)
       {
	 hwigfunc[i]->SetBinContent(u, func_wiggle->Eval(hwigfunc[i]->GetBinCenter(u)));
	 }*/
       for(int u=1;u<=npulse;u++)
       {
	 v=func_wiggle->GetRandom();
         hwiggle1[i]->SetBinContent(hwiggle1[i]->FindBin(v), hedist->GetRandom());
	 //* func_wiggle->Eval(hwiggle1[i]->GetBinCenter(v)));
       }
    }
  
   for(int i=0;i<=53;i++)
    {
      /*  
     for(int u=1;u<=252000;u++)
       {
	 hwigfunc[i]->SetBinContent(u, func_wiggle->Eval(hwigfunc[i]->GetBinCenter(u)));
	 }*/
      for(int u=1;u<=npulse;u++)
       {
	 v=func_wiggle->GetRandom();
         hwiggle2[i]->SetBinContent(hwiggle2[i]->FindBin(v), hedist->GetRandom());
	 //* func_wiggle->Eval(hwiggle2[i]->GetBinCenter(v)));
       }
    }

      for(int i=0;i<=53;i++)
    {
      /* 
     for(int u=1;u<=252000;u++)
       {
	 hwigfunc[i]->SetBinContent(u, func_wiggle->Eval(hwigfunc[i]->GetBinCenter(u)));
	 }*/
     for(int u=1;u<=npulse;u++)
       {
	 v=func_wiggle->GetRandom();
         hwiggle3[i]->SetBinContent(hwiggle3[i]->FindBin(v), hedist->GetRandom());
	 //* func_wiggle->Eval(hwiggle3[i]->GetBinCenter(v)));
       }
    }

         for(int i=0;i<=53;i++)
    {
      /*  
     for(int u=1;u<=252000;u++)
       {
	 hwigfunc[i]->SetBinContent(u, func_wiggle->Eval(hwigfunc[i]->GetBinCenter(u)));
	 }*/
     for(int u=1;u<=npulse;u++)
       {
	 v=func_wiggle->GetRandom();
         hwiggle4[i]->SetBinContent(hwiggle4[i]->FindBin(v), hedist->GetRandom());
	 //* func_wiggle->Eval(hwiggle4[i]->GetBinCenter(v)));
       }
    }
	 

	 // ###### Create pedestal drift histograms using the shape from pedestal drift function
    // ##### Add together pedestal drift histogram, noise hostogram and wiggle histograms for each crystals for 4 fills
	 
	 
    for(int n=0;n<=53;n++)
    {
     for(int i=1;i<=252000;i++)
      {
        hpeddrift1[n]->SetBinContent(i,func_peddrift->Eval(hpeddrift1[n]->GetBinCenter(i)));
      }
    }
    
  int l=1;
  for(int i=0;i<=53;i++)
    {
      for(int k=1;k<=252000;k++)
	{
	  hpeddrift1[i]->SetBinContent(k,int (hpeddrift1[i]->GetBinContent(k)+hsimnoise1->GetBinContent(l)+hwiggle1[i]->GetBinContent(k)));
	    //  hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l));
	  l=l+1;
        }
    }

  
    for(int n=0;n<=53;n++)
    {
     for(int i=1;i<=252000;i++)
      {
        hpeddrift2[n]->SetBinContent(i,func_peddrift->Eval(hpeddrift2[n]->GetBinCenter(i)));
      }
    }
  l=1;
  for(int i=0;i<=53;i++)
    {
      for(int k=1;k<=252000;k++)
	{
	  hpeddrift2[i]->SetBinContent(k,int (hpeddrift2[i]->GetBinContent(k)+hsimnoise2->GetBinContent(l)+hwiggle2[i]->GetBinContent(k)));
	    //  hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l));
	  l=l+1;
        }
    }
  

    for(int n=0;n<=53;n++)
    {
     for(int i=1;i<=252000;i++)
      {
        hpeddrift3[n]->SetBinContent(i,func_peddrift->Eval(hpeddrift3[n]->GetBinCenter(i)));
      }
    }
  l=1;
  for(int i=0;i<=53;i++)
    {
      for(int k=1;k<=252000;k++)
	{
	  hpeddrift3[i]->SetBinContent(k,int (hpeddrift3[i]->GetBinContent(k)+hsimnoise3->GetBinContent(l)+hwiggle3[i]->GetBinContent(k)));
	    //  hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l));
	  l=l+1;
        }
    }


    for(int n=0;n<=53;n++)
    {
     for(int i=1;i<=252000;i++)
      {
        hpeddrift4[n]->SetBinContent(i,func_peddrift->Eval(hpeddrift4[n]->GetBinCenter(i)));
      }
    }
  l=1;
  for(int i=0;i<=53;i++)
    {
      for(int k=1;k<=252000;k++)
	{
	  hpeddrift4[i]->SetBinContent(k,int (hpeddrift4[i]->GetBinContent(k)+hsimnoise4->GetBinContent(l)+hwiggle4[i]->GetBinContent(k)));
	    //  hpeddrift[i]->SetBinContent(k,hpeddrift[i]->GetBinContent(k)+hsimnoise->GetBinContent(l));
	  l=l+1;
        }
    }


  // ###### Subtracting average pedestal like and unlike the way done in DAQ with 15 time decimation

  

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
	     adc += int (hpeddrift1[i]->GetBinContent(r));
	     fadc += int (hpeddrift1[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+15;
          hfill1[i]->SetBinContent(t,adc);
	  hdfill1[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      hdfill1[i]->SetLineColor(kRed);
    }


  adc=0;
  fadc=0;
  k=1;
 
  for(int i=0;i<=53;i++)
    {
      t=1;
      k=1;
      while(k<=252000-15+1)
	{
	  for(int r=k; r<k+15; r++)
	    {
	     adc += int (hpeddrift2[i]->GetBinContent(r));
	     fadc += int (hpeddrift2[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+15;
          hfill2[i]->SetBinContent(t,adc);
	  hdfill2[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      hdfill2[i]->SetLineColor(kRed);
    }


  adc=0;
  fadc=0;
  k=1;
  
  for(int i=0;i<=53;i++)
    {
      t=1;
      k=1;
      while(k<=252000-15+1)
	{
	  for(int r=k; r<k+15; r++)
	    {
	     adc += int (hpeddrift3[i]->GetBinContent(r));
	     fadc += int (hpeddrift3[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+15;
          hfill3[i]->SetBinContent(t,adc);
	  hdfill3[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      hdfill3[i]->SetLineColor(kRed);
    }



  adc=0;
  fadc=0;
  k=1;
  
  for(int i=0;i<=53;i++)
    {
      t=1;
      k=1;
      while(k<=252000-15+1)
	{
	  for(int r=k; r<k+15; r++)
	    {
	     adc += int (hpeddrift4[i]->GetBinContent(r));
	     fadc += int (hpeddrift4[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+15;
          hfill4[i]->SetBinContent(t,adc);
	  hdfill4[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      hdfill4[i]->SetLineColor(kRed);
    }
  

  for(int i=0;i<=53;i++)
    {
 /*   hfillsum[i]->Add(hfill1[i],1);
       hfillsum[i]->Add(hfill2[i],1);
        hfillsum[i]->Add(hfill3[i],1);
	 hfillsum[i]->Add(hfill4[i],1);
       hdfillsum[i]->Add(hdfill1[i],1);
        hdfillsum[i]->Add(hdfill2[i],1);
         hdfillsum[i]->Add(hdfill3[i],1);
	  hdfillsum[i]->Add(hdfill4[i],1);
      */
      for(int in=1; in<=16800; in++)
	{
	  hfillsum[i]->SetBinContent(in,hfill1[i]->GetBinContent(in)+hfill2[i]->GetBinContent(in)+hfill3[i]->GetBinContent(in)+hfill4[i]->GetBinContent(in));
	  hdfillsum[i]->SetBinContent(in,hdfill1[i]->GetBinContent(in)+hdfill2[i]->GetBinContent(in)+hdfill3[i]->GetBinContent(in)+hdfill4[i]->GetBinContent(in));
     	}
     
    }

  /*  hdataxtal->SetLineColor(kRed);
  hdataxtal->Draw();
  hfillsum[10]->Draw("same");*/



 // ###### Subtracting average pedestal like and unlike the way done in DAQ with 60 time decimation

  

  adc=0;
  fadc=0;
  k=1;
  for(int i=0;i<=53;i++)
    {
      t=1;
      k=1;
      while(k<=252000-60+1)
	{
	  for(int r=k; r<k+60; r++)
	    {
	     adc += int (hpeddrift1[i]->GetBinContent(r));
	     fadc += int (hpeddrift1[i]->GetBinContent(r));
	     adc -= xtal_ped[i];
	     fadc -= xtal_ped[i];
	    }
	  k=k+60;
          h60fill1[i]->SetBinContent(t,adc);
	  h60dfill1[i]->SetBinContent(t,fadc);
	  t=t+1;
	  adc=0;
	  fadc=0.0;
	}
      h60dfill1[i]->SetLineColor(kRed);
    }


  

  for(int i=0;i<=53;i++)
    {
 /*   hfillsum[i]->Add(hfill1[i],1);
       hfillsum[i]->Add(hfill2[i],1);
        hfillsum[i]->Add(hfill3[i],1);
	 hfillsum[i]->Add(hfill4[i],1);
       hdfillsum[i]->Add(hdfill1[i],1);
        hdfillsum[i]->Add(hdfill2[i],1);
         hdfillsum[i]->Add(hdfill3[i],1);
	  hdfillsum[i]->Add(hdfill4[i],1);
      */
      for(int in=1; in<=4200; in++)
	{
	  h60fillsum[i]->SetBinContent(in,h60fill1[i]->GetBinContent(in));
	  h60dfillsum[i]->SetBinContent(in,h60dfill1[i]->GetBinContent(in));
     	}
     
    }

  
  


  double pedsum=0;
  double dpedsum=0;
  double diff, fdiff;
  
  for(int i=0; i<=53; i++)
    {
      for(int ib=0; ib<=16800; ib++)
	{
	  for(int ic=ib-(4+1);ic<=ib+(4+1);ic++)
	    {
	      if(ic<ib-1 || ic>ib+1)
		{
		  pedsum+=hfillsum[i]->GetBinContent(ic);
		  dpedsum+=hdfillsum[i]->GetBinContent(ic);
		}
	    }
	  diff=hfillsum[i]->GetBinContent(ib)-(pedsum/8);
	  fdiff=hdfillsum[i]->GetBinContent(ib)-(dpedsum/8);
	  hpedsum[i]->SetBinContent(ib,pedsum/8);
	  hdpedsum[i]->SetBinContent(ib,dpedsum/8);
	  if(diff>=300)
	    {
	     hsig[i]->SetBinContent(ib,diff);
	    }
	  if(fdiff>=300)
	    {
	     hdsig[i]->SetBinContent(ib,fdiff);
	    }
	  pedsum=0;
	  dpedsum=0;
	}
      hdsig[i]->SetLineColor(kRed);
      hdpedsum[i]->SetLineColor(kRed);
    }
   for(int i=0; i<=53; i++)
    {


      for(int ibin=1;ibin<=16800;ibin++)
	{
	  hdiff[i]->SetBinContent(ibin,hdsig[i]->GetBinContent(ibin)-hsig[i]->GetBinContent(ibin));
	  hpeddiff[i]->SetBinContent(ibin, hdpedsum[i]->GetBinContent(ibin)-hpedsum[i]->GetBinContent(ibin));
	}
      hdiff[i]->SetLineWidth(3);
      
    }


   for(int i=0; i<=53; i++)
    {
      for(int ibin=1;ibin<=16800;ibin++)
	{
	  hrawdiff[i]->SetBinContent(ibin,hdfillsum[i]->GetBinContent(ibin)-hfillsum[i]->GetBinContent(ibin));
	}
      //  hrawdiff[i]->SetLineWidth(3);
    }



 pedsum=0;
 dpedsum=0;
  
  for(int i=0; i<=53; i++)
    {
      for(int ib=0; ib<=4200; ib++)
	{
	  for(int ic=ib-(4+1);ic<=ib+(4+1);ic++)
	    {
	      if(ic<ib-1 || ic>ib+1)
		{
		  pedsum+=h60fillsum[i]->GetBinContent(ic);
		  dpedsum+=h60dfillsum[i]->GetBinContent(ic);
		}
	    }
	  diff=h60fillsum[i]->GetBinContent(ib)-(pedsum/8);
	  fdiff=h60dfillsum[i]->GetBinContent(ib)-(dpedsum/8);
	  if(diff>=300)
	    {
	     h60sig[i]->SetBinContent(ib,diff);
	    }
	  if(fdiff>=300)
	    {
	     h60dsig[i]->SetBinContent(ib,fdiff);
	    }
	  pedsum=0;
	  dpedsum=0;
	}
      h60dsig[i]->SetLineColor(kRed);
    }



   
    for(int i=0;i<=53;i++)
     {
       hxtalsum->Add(hsig[i],1);
       hdxtalsum->Add(hdsig[i],1);
     }

     for(int i=0;i<=53;i++)
     {
       h60xtalsum->Add(h60sig[i],1);
       h60dxtalsum->Add(h60dsig[i],1);
     }

    
    iflush++;
 } // END FLUSH LOOP
gettimeofday(&t_end, NULL);
 cout<<"End of flush loop. Number of flushes is "<<nflush<<" and time to finish is "<< (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0)<<" seconds"<<endl; 
   hdxtalsum->SetLineColor(kRed);
   h60dxtalsum->SetLineColor(kRed);


   for(int ibin=1; ibin<=16800;ibin++)
     {
       hxtalsumdiff->SetBinContent(ibin, hdxtalsum->GetBinContent(ibin)-hxtalsum->GetBinContent(ibin));
     }

    for(int ibin=1; ibin<=4200;ibin++)
     {
       h60xtalsumdiff->SetBinContent(ibin, h60dxtalsum->GetBinContent(ibin)-h60xtalsum->GetBinContent(ibin));
     }

   c1->Divide(1,3);
   c1->cd(1);
   for(int i=1;i<=53;i++){hfillsum[0]->Add(hfillsum[i],1);}
    for(int i=1;i<=53;i++){h60fillsum[0]->Add(h60fillsum[i],1);}
   h60fillsum[0]->GetYaxis()->SetRangeUser(-50000,2000);
   h60fillsum[0]->Draw();
   c1->cd(2);
   h60dxtalsum->Draw();
   h60xtalsum->Draw("same");
   c1->cd(3);
   //  hxtalsumdiff->Divide(hxtalsum);
   h60xtalsumdiff->Draw();
   
  foutput=new TFile("RoundSimOutput_ped7.root","new");

  hfillsum[0]->Write();
  hxtalsum->Write();
  hdxtalsum->Write();
  hxtalsumdiff->Write();

  h60fillsum[0]->Write();
  h60xtalsum->Write();
  h60dxtalsum->Write();
  h60xtalsumdiff->Write();

  foutput->Close();

   
}
