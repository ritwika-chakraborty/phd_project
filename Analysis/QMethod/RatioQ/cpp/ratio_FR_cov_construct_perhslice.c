#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TFormula.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include <sys/time.h>



float toddiff(struct timeval *tod2, struct timeval *tod1)
{
  float fdt, fmudt;
  long long t1, t2, mut1, mut2;
  long long dt, mudt;
  t1 = tod1->tv_sec;
  mut1 = tod1->tv_usec;
  t2 = tod2->tv_sec;
  mut2 = tod2->tv_usec;
  dt = (t2 - t1);
  mudt = (mut2 - mut1);
  fdt = (float)dt;
  fmudt = (float)mudt;
  return 1.0e6 * fdt + fmudt;
}



//char root_file_name[128] = "highstat_wigsim_output_statfluc_binintegral_e16.root";
//char root_file_name[128] = "run2g_corrOOF_thresh400.root";
//char root_file_name[128] = "run2C_thresh300_reprocessed_new.root";
char root_file_name[128] = "run2_final.root";
char *histname = new char[10];
char *histname2 = new char[10];
TVirtualFitter *gFitter;
TFile *_file[25];
TDirectory *dir[25];
TH1D *h[25],*qHist_1D[25],*h_sum, *h1, *h2, *hp, *hm, *h_ratio, *h_num, *h_deno, *hcalo, *hfit, *hr[25], *h_res, *hr_sum, *hcov_1_p, *hcov_1_m,*hcov_13_p,*hcov_13_m,*hcov_14_p, *hcov_14_m, *hcov_15_p, *hcov_15_m, *hcov_27_p, *hcov_27_m,*hcov_28_p, *hcov_28_m, *hcov_29_p, *hcov_29_m, *hr_noshift, *hr_shift, *hcomp[25],*hcomp_sum,*htemp;
TCanvas *c,*c1, *c2,*c3, *cauto;
TGraphErrors *gr1;
TH1 *hfft, *hFFTsq, *hAuto;
TF1 *fit_func3, *fit_func7, *fit_func11, *fit_func15, *fit_func19, *fit_func23, *fit_func24;
//TH1 *hfft3, *hfft7, *hfft11, *hfft15, *hfft19, *hfft23, *hfft24;
TH1F *hlm;
Double_t fit_start=30000;
Double_t fit_stop=300000;
//Double_t fit_stop=306076.25;
Double_t rawBinToNs = 1.25;
Double_t inject_time;
Double_t T_a_true=4365.411;
Int_t nbinshift;
Double_t T_a;
Double_t lifetime=64440;
bool usesetbins=true;
bool useerrorbars=true;
Int_t m=0;
Double_t blindR[25],dblindR[25],n[25];
Int_t countfcn=0;

TFile *foutput;
TH2D *hcov[75];
int i0, iend;
int flg = 0;
int mdim = 2100;
// int mdim = 10;
// TMatrixD cov(hcomp->GetNbinsX(),hcomp->GetNbinsX());
//TMatrixD cov2(mdim,mdim);
TMatrixD cov(mdim,mdim);
  // TArrayD  data(hcomp->GetNbinsX()*hcomp->GetNbinsX());
TArrayD data((mdim)*(mdim)),data2((mdim)*(mdim));

TArrayD ratio_cov_matrix(TH1D *hr_sum, TH1D *h_sum, Double_t fit_start)
{
  

  /* for (Int_t i = 0; i < mdim*mdim; i++) {
    const Int_t ir = i/mdim;
    const Int_t ic = i%mdim;
    data[i] = 0.0;
    if ( ir == ic ) data[i] = 1;
    if ( ir == ic-1 || ir == ic+1 ) data[i] = 7;
  }
  h.SetMatrixArray(data.GetArray());

  h.Print();

  Double_t det1;
  h.Invert(&det1);

  h.Print();
  cout<<h[1][2]<<" "<<endl;*/

  double covar, dprod1, dprod2, dprod3, dprod4, dprod5, dprod6, dprod7, dprod8, dprod9, dprod10, dprod11,dprod12, dprod13, dprod14, dprod15, dprod16, dprod17,dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI21, dI22, dI23, dI24, dI25, dI26, dI27, dI28, dI29;

  double a,bp,bm,cp,cm,dp,dm,ep,em,fp,fm,gp,gm,hp,hm,ip,im,jp,jm,kp,km,lp,lm,mp,mm,np,nm,op,om,pp,pm,qp,qm,rp,rm,sp,sm,tp,tm,up,um,vp,vm,wp,wm,xp,xm,yp,ym;

   double da,dbp,dbm,dcp,dcm,ddp,ddm,dep,dem,dfp,dfm,dgp,dgm,dhp,dhm,dip,dim,djp,djm,dkp,dkm,dlp,dlm,dmp,dmm,dnp,dnm,dop,dom,dpp,dpm,dqp,dqm,drp,drm,dsp,dsm,dtp,dtm,dup,dum,dvp,dvm,dwp,dwm,dxp,dxm,dyp,dym;

  for(Int_t k=hr_sum->FindBin(10000); k<=hr_sum->GetNbinsX(); k++)
    {
      if(hr_sum->GetBinError(k)!=0)
	{ i0=k;
	  break;}
    }
  cout<<"i0 is"<<i0<<" "<<endl;
   cout<<"check "<<h_sum->GetBinError(500)<<endl;
   for (Int_t i = 0; i < (mdim)*(mdim); i++)
      {
	const Int_t ir = (i)/(mdim);
	const Int_t ic = (i)%(mdim);
        
          a=h_sum->GetBinContent(2*ir);
	  da=h_sum->GetBinError(2*ir);
	  
          bp=h_sum->GetBinContent((2*ir)+1);
	  dbp=h_sum->GetBinError((2*ir)+1);
	  bm=h_sum->GetBinContent((2*ir)-1);
	  dbm=h_sum->GetBinError((2*ir)-1);
	  
	  cp=h_sum->GetBinContent((2*ir)+2);
	  dcp=h_sum->GetBinError((2*ir)+2);
	  cm=h_sum->GetBinContent((2*ir)-2);
	  dcm=h_sum->GetBinError((2*ir)-2);

	  dp=h_sum->GetBinContent((2*ir)+3);
	  ddp=h_sum->GetBinError((2*ir)+3);
	  dm=h_sum->GetBinContent((2*ir)-3);
	  ddm=h_sum->GetBinError((2*ir)-3);

	  ep=h_sum->GetBinContent(2*(ir+nbinshift)-3);
	  dep=h_sum->GetBinError(2*(ir+nbinshift)-3);
	  em=h_sum->GetBinContent(2*(ir-nbinshift)+3);
	  dem=h_sum->GetBinError(2*(ir-nbinshift)+3);

	  fp=h_sum->GetBinContent(2*(ir+nbinshift)-2);
	  dfp=h_sum->GetBinError(2*(ir+nbinshift)-2);
	  fm=h_sum->GetBinContent(2*(ir-nbinshift)+2);
	  dfm=h_sum->GetBinError(2*(ir-nbinshift)+2);

	  gp=h_sum->GetBinContent(2*(ir+nbinshift)-1);
	  dgp=h_sum->GetBinError(2*(ir+nbinshift)-1);
	  gm=h_sum->GetBinContent(2*(ir-nbinshift)+1);
	  dgm=h_sum->GetBinError(2*(ir-nbinshift)+1);

	  hp=h_sum->GetBinContent(2*(ir+nbinshift));
	  dhp=h_sum->GetBinError(2*(ir+nbinshift));
	  hm=h_sum->GetBinContent(2*(ir-nbinshift));
	  dhm=h_sum->GetBinError(2*(ir-nbinshift));

          ip=h_sum->GetBinContent(2*(ir+nbinshift)+1);
	  dip=h_sum->GetBinError(2*(ir+nbinshift)+1);
	  im=h_sum->GetBinContent(2*(ir-nbinshift)-1);
	  dim=h_sum->GetBinError(2*(ir-nbinshift)-1);

	  jp=h_sum->GetBinContent(2*(ir+nbinshift)+2);
	  djp=h_sum->GetBinError(2*(ir+nbinshift)+2);
	  jm=h_sum->GetBinContent(2*(ir-nbinshift)-2);
	  djm=h_sum->GetBinError(2*(ir-nbinshift)-2);

	  kp=h_sum->GetBinContent(2*(ir+nbinshift)+3);
	  dkp=h_sum->GetBinError(2*(ir+nbinshift)+3);
	  km=h_sum->GetBinContent(2*(ir-nbinshift)-3);
	  dkm=h_sum->GetBinError(2*(ir-nbinshift)-3);

	  lp=h_sum->GetBinContent(2*(ir+2*nbinshift)-3);
	  dlp=h_sum->GetBinError(2*(ir+2*nbinshift)-3);
	  lm=h_sum->GetBinContent(2*(ir-2*nbinshift)+3);
	  dlm=h_sum->GetBinError(2*(ir-2*nbinshift)+3);

	  mp=h_sum->GetBinContent(2*(ir+2*nbinshift)-2);
	  dmp=h_sum->GetBinError(2*(ir+2*nbinshift)-2);
	  mm=h_sum->GetBinContent(2*(ir-2*nbinshift)+2);
	  dmm=h_sum->GetBinError(2*(ir-2*nbinshift)+2);

	  np=h_sum->GetBinContent(2*(ir+2*nbinshift)-1);
	  dnp=h_sum->GetBinError(2*(ir+2*nbinshift)-1);
	  nm=h_sum->GetBinContent(2*(ir-2*nbinshift)+1);
	  dnm=h_sum->GetBinError(2*(ir-2*nbinshift)+1);

	  op=h_sum->GetBinContent(2*(ir+2*nbinshift));
	  dop=h_sum->GetBinError(2*(ir+2*nbinshift));
	  om=h_sum->GetBinContent(2*(ir-2*nbinshift));
	  dom=h_sum->GetBinError(2*(ir-2*nbinshift));

	  pp=h_sum->GetBinContent(2*(ir+2*nbinshift)+1);
	  dpp=h_sum->GetBinError(2*(ir+2*nbinshift)+1);
	  pm=h_sum->GetBinContent(2*(ir-2*nbinshift)-1);
	  dpm=h_sum->GetBinError(2*(ir-2*nbinshift)-1);

	  qp=h_sum->GetBinContent(2*(ir+2*nbinshift)+2);
	  dqp=h_sum->GetBinError(2*(ir+2*nbinshift)+2);
	  qm=h_sum->GetBinContent(2*(ir-2*nbinshift)-2);
	  dqm=h_sum->GetBinError(2*(ir-2*nbinshift)-2);

	  rp=h_sum->GetBinContent(2*(ir+2*nbinshift)+3);
	  drp=h_sum->GetBinError(2*(ir+2*nbinshift)+3);
	  rm=h_sum->GetBinContent(2*(ir-2*nbinshift)-3);
	  drm=h_sum->GetBinError(2*(ir-2*nbinshift)-3);

	  sp=h_sum->GetBinContent(2*(ir+3*nbinshift)-3);
	  dsp=h_sum->GetBinError(2*(ir+3*nbinshift)-3);
	  sm=h_sum->GetBinContent(2*(ir-3*nbinshift)+3);
	  dsm=h_sum->GetBinError(2*(ir-3*nbinshift)+3);

	  tp=h_sum->GetBinContent(2*(ir+3*nbinshift)-2);
	  dtp=h_sum->GetBinError(2*(ir+3*nbinshift)-2);
	  tm=h_sum->GetBinContent(2*(ir-3*nbinshift)+2);
	  dtm=h_sum->GetBinError(2*(ir-3*nbinshift)+2);

	  up=h_sum->GetBinContent(2*(ir+3*nbinshift)-1);
	  dup=h_sum->GetBinError(2*(ir+3*nbinshift)-1);
	  um=h_sum->GetBinContent(2*(ir-3*nbinshift)+1);
	  dum=h_sum->GetBinError(2*(ir-3*nbinshift)+1);

	  vp=h_sum->GetBinContent(2*(ir+3*nbinshift));
	  dvp=h_sum->GetBinError(2*(ir+3*nbinshift));
	  vm=h_sum->GetBinContent(2*(ir-3*nbinshift));
	  dvm=h_sum->GetBinError(2*(ir-3*nbinshift));

	  wp=h_sum->GetBinContent(2*(ir+3*nbinshift)+1);
	  dwp=h_sum->GetBinError(2*(ir+3*nbinshift)+1);
	  wm=h_sum->GetBinContent(2*(ir-3*nbinshift)-1);
	  dwm=h_sum->GetBinError(2*(ir-3*nbinshift)-1);

	  xp=h_sum->GetBinContent(2*(ir+3*nbinshift)+2);
	  dxp=h_sum->GetBinError(2*(ir+3*nbinshift)+2);
	  xm=h_sum->GetBinContent(2*(ir-3*nbinshift)-2);
	  dxm=h_sum->GetBinError(2*(ir-3*nbinshift)-2);

	  yp=h_sum->GetBinContent(2*(ir+3*nbinshift)+3);
	  dyp=h_sum->GetBinError(2*(ir+3*nbinshift)+3);
	  ym=h_sum->GetBinContent(2*(ir-3*nbinshift)-3);
	  dym=h_sum->GetBinError(2*(ir-3*nbinshift)-3);
        
	
        data[i] = 0.0;
 
        if(ir==ic)
	  {

	    data[i]=(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ir));
	    /*  if(data[i]<1.0000e-23)
	      {
		cout<<"zero diag elements"<<data[i]<<" "<<ir<<" "<<ic<<" "<<endl;
		data[i]=1.0000e-23;
	      }
	    */
	    //data[i]=1;
	    //cout<<"i "<<i<<"data "<<data[i]<<endl;
	  }

		if(ir==ic-1) //diag+1
	  {
	    dprod1=((32*(2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bp - 4*cp - 2*dp + em + 2*fm + gm + ip + 2*jp + kp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp)))*da*da;
	    dprod3=((32*(2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp)))*dhp*dhp;
	    dprod4=((8*(2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp)))*dim*dim;
	    dprod5=((8*((bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) + 4*(2*a + bm + bp)*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp) + (2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bp - 4*cp - 2*dp + em +  2*fm + gm + ip + 2*jp + kp))/(pow((4*a + 2*bm + 2*bp + gm + gp +  2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm +  ip + 2*jp + kp)))*dbm*dbm;
	    dprod7=-((8*(2*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(em + 2*fm + gm + ip + 2*jp + kp) - 4*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(em + 2*fm + gm + ip + 2*jp + kp)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp) - 2*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bp - 4*cp - 2*dp + em + 2*fm + gm + ip + 2*jp + kp)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dbp*dbp;
	    dprod8=((8*(2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp)))*dgp*dgp;
	    dprod9=((8*((bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) + 4*(2*a + bm + bp)*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp) + (2*a + bm + bp)*(2*bp + 4*cp + 2*dp - em - 2*fm - gm - ip - 2*jp - kp)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dip*dip;
	    dprod10=((32*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dfm*dfm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(em + 2*fm + gm + ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dcp*dcp;
	    dprod12=((32*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*djp*djp;
	    dprod13=((8*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im -  ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dem*dem;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(em + 2*fm + gm + ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*ddp*ddp;
	    dprod15=((8*(bp + 2*cp + dp)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))))*dkp*dkp;

	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dfm*dfm;
	    dI22=-((64*(em + 2*fm + gm + ip + 2*jp + kp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dcp*dcp;
	    dI23=((32*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip +  2*jp + kp),3))*djp*djp;
	    dI24=((8*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dgm*dgm;
	    dI25=((8*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dem*dem;
	    dI26=-((16*(em + 2*fm + gm + ip + 2*jp + kp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dbp*dbp;
	    dI27=-((16*(em + 2*fm + gm + ip + 2*jp + kp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*ddp*ddp;
	    dI28=((8*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dip*dip;
	    dI29=((8*(bp + 2*cp + dp))/pow((2*bp + 4*cp + 2*dp + em + 2*fm + gm + ip + 2*jp + kp),3))*dkp*dkp;

            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
	  }

       	if(ir==ic+1) //diag-1
	  {
            dprod1=((32*(2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bm - 4*cm - 2*dm + ep + 2*fp + gp + im + 2*jm + km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*da*da;
	    dprod3=((32*(2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*dhp*dhp;
	    dprod4=((8*((bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) + 4*(2*a + bm + bp)*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km) + (2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dim*dim;
	    dprod5=((8*(2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*dgm*dgm;
	    dprod6=-((8*(2*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(ep + 2*fp + gp + im + 2*jm + km) - 4*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(ep + 2*fp + gp + im + 2*jm + km)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km) - 2*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bm - 4*cm - 2*dm + ep + 2*fp + gp + im + 2*jm + km)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(-2*bm - 4*cm - 2*dm + ep + 2*fp + gp + im + 2*jm + km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*dbp*dbp;
	    dprod8=((8*((bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) + 4*(2*a + bm + bp)*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km) + (2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dgp*dgp;
	    dprod9=((8*(2*a + bm + bp)*(2*bm + 4*cm + 2*dm - ep - 2*fp - gp - im - 2*jm - km))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km)))*dip*dip;
	    dprod10=((32*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*djm*djm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ep + 2*fp + gp + im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dcm*dcm;
	    dprod12=((32*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dfp*dfp;
	    dprod13=((8*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dkm*dkm;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ep + 2*fp + gp + im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*ddm*ddm;
	    dprod15=((8*(bm + 2*cm + dm)*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))))*dep*dep;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*djm*djm;
	    dI22=-((64*(ep + 2*fp + gp + im + 2*jm + km))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dcm*dcm;
	    dI23=((32*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dfp*dfp;
	    dI24=((8*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dkm*dkm;
	    dI25=((8*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dim*dim;
	    dI26=-((16*(ep + 2*fp + gp + im + 2*jm + km))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*ddm*ddm;
	    dI27=-((16*(ep + 2*fp + gp + im + 2*jm + km))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dbm*dbm;
	    dI28=((8*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dep*dep;
	    dI29=((8*(bm + 2*cm + dm))/pow((2*bm + 4*cm + 2*dm + ep + 2*fp + gp + im + 2*jm + km),3))*dgp*dgp;
	    
	    // covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //  data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	      data[i]=covar;
	  }

       	 if(ir==ic-(nbinshift-1)) //diag+13
	  {
            dprod1=-((32*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dgm*dgm;
	    dprod6=((2*(-4*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(ep + 2*fp + gp)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dbp*dbp;
	    dprod8=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(bm + 2*cm + dm + lp + 2*mp + np) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm + lp + 2*mp + np)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np) + 4*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*ep - 4*fp - 2*gp + lp + 2*mp + np))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np)))*dip*dip;
	    dprod10=-((32*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dcm*dcm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bm + 2*cm + dm + lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dfp*dfp;
	    dprod12=-((32*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dmp*dmp;
	    dprod13=-((8*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*ddm*ddm;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bm + 2*cm + dm + lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dep*dep;
	    dprod15=-((8*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dlp*dlp;
	    dprod16=-((8*(ep + 2*fp + gp)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))))*dnp*dnp;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dcm*dcm;
	    dI22=-((64*(bm + 2*cm + dm + lp + 2*mp + np))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dfp*dfp;
	    dI23=((32*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dmp*dmp;
	    dI24=((8*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*ddm*ddm;
	    dI25=((8*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dbm*dbm;
	    dI26=-((16*(bm + 2*cm + dm + lp + 2*mp + np))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dep*dep;
	    dI27=-((16*(bm + 2*cm + dm + lp + 2*mp + np))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dgp*dgp;
	    dI28=((8*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dlp*dlp;
	    dI29=((8*(ep + 2*fp + gp))/pow((bm + 2*cm + dm + 2*ep + 4*fp + 2*gp + lp + 2*mp + np),3))*dnp*dnp;

	    //covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
       	  }

       	if(ir==ic+(nbinshift-1)) //diag-13
	  {
            dprod1=-((32*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dim*dim;
	    dprod5=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(bp + 2*cp + dp + lm + 2*mm + nm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp + lm + 2*mm + nm)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm) + 4*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm +  2*mm + nm)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dbm*dbm;
	    dprod7=((2*(-4*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(em + 2*fm + gm)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*em - 4*fm - 2*gm + lm + 2*mm + nm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm)))*dip*dip;
	    dprod10=-((32*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dmm*dmm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bp + 2*cp + dp + lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dfm*dfm;
	    dprod12=-((32*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dcp*dcp;
	    dprod13=-((8*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dnm*dnm;
	    dprod14=-((8*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dlm*dlm;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bp + 2*cp + dp + lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*dem*dem;
	    dprod16=-((8*(em + 2*fm + gm)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))))*ddp*ddp;
	    
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))*dmm*dmm;
	    dI22=-((64*(bp + 2*cp + dp + lm + 2*mm + nm))/pow((bp + 2*cp + dp + 2*em +  4*fm + 2*gm + lm + 2*mm + nm),3))*dfm*dfm;
	    dI23=((32*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))*dcp*dcp;
	    dI24=((8*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm +  2*mm + nm),3))*dnm*dnm;
	    dI25=((8*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm +  2*mm + nm),3))*dlm*dlm;
	    dI26=-((16*(bp + 2*cp + dp + lm + 2*mm + nm))/pow((bp + 2*cp + dp + 2*em +  4*fm + 2*gm + lm + 2*mm + nm),3))*dgm*dgm;
	    dI27=-((16*(bp + 2*cp + dp + lm + 2*mm + nm))/pow((bp + 2*cp + dp + 2*em +  4*fm + 2*gm + lm + 2*mm + nm),3))*dem*dem;
	    dI28=((8*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))*dbp*dbp;
	    dI29=((8*(em + 2*fm + gm))/pow((bp + 2*cp + dp + 2*em + 4*fm + 2*gm + lm + 2*mm + nm),3))*ddp*ddp;

	    // covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
      	  }

       	if(ir==ic-(nbinshift)) //diag+14
	  {
            dprod1=-((32*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp))/((pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp)))*dhm*dhm;
	    dprod2=((8*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gp + 2*hp + ip)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*da*da;
	    dprod3=-((8*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + np + 2*op + pp)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op +pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp)))*dgm*dgm;
	    dprod6=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gp + 2*hp + ip)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dbm*dbm;
	    dprod7=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gp + 2*hp + ip)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dbp*dbp;
	    dprod8=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + np + 2*op + pp)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dgp*dgp;
	    dprod9=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + np + 2*op + pp)*(2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gp - 4*hp - 2*ip + np + 2*op + pp)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dop*dop;
	    dprod11=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dnp*dnp;
	    dprod12=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))))*dpp*dpp;
	  
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*da*da;
	    dI22=-((64*(2*a + bm + bp + np + 2*op + pp))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dhp*dhp;
	    dI23=((32*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dop*dop;
	    dI24=((8*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dbm*dbm;
	    dI25=((8*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dbp*dbp;
	    dI26=-((16*(2*a + bm + bp + np + 2*op + pp))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dgp*dgp;
	    dI27=-((16*(2*a + bm + bp + np + 2*op + pp))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dip*dip;
	    dI28=((8*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dnp*dnp;
	    dI29=((8*(gp + 2*hp + ip))/pow((2*a + bm + bp + 2*gp + 4*hp + 2*ip + np + 2*op + pp),3))*dpp*dpp;
	    
	    //covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));

	    covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
       	  }

      	if(ir==ic+(nbinshift)) //diag-14
	  {
            dprod1=-((8*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + nm + 2*om + pm)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm +  2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/((pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dhm*dhm;
	    dprod2=((8*(-4*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gm + 2*hm + im)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/((pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm)))*dhp*dhp;
	    dprod4=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + nm + 2*om + pm)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dim*dim;
	    dprod5=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(2*a + bm + bp + nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + nm + 2*om + pm)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 4*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dgm*dgm;
	    dprod6=((2*(-4*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gm + 2*hm + im)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dbm*dbm;
	    dprod7=((2*(-4*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2)) - 16*(gm + 2*hm + im)*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(2*a + bm + bp - 2*gm - 4*hm - 2*im + nm + 2*om + pm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm)))*dip*dip;
	    dprod10=-((32*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dom*dom;
	    dprod11=-((8*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dpm*dpm;
	    dprod12=-((8*(gm + 2*hm + im)*(-4*a - 2*bm - 2*bp + gm + gp + 2*hm + 2*hp + im + ip))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))))*dnm*dnm;

	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dom*dom;
	    dI22=-((64*(2*a + bm + bp + nm + 2*om + pm))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dhm*dhm;
	    dI23=((32*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*da*da;
	    dI24=((8*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dpm*dpm;
	    dI25=((8*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dnm*dnm;
	    dI26=-((16*(2*a + bm + bp + nm + 2*om + pm))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dim*dim;
	    dI27=-((16*(2*a + bm + bp + nm + 2*om + pm))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dgm*dgm;
	    dI28=((8*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dbm*dbm;
	    dI29=((8*(gm + 2*hm + im))/pow((2*a + bm + bp + 2*gm + 4*hm + 2*im + nm + 2*om + pm),3))*dbp*dbp;


	    covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //  data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
      	  }

		if(ir==ic-(nbinshift+1)) //diag+15
	  {
            dprod1=-((32*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dbm*dbm;
	    dprod7=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(ip + 2*jp + kp) - 16*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(ip + 2*jp + kp)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp)))*dgp*dgp;
	    dprod9=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(bp + 2*cp + dp + pp + 2*qp + rp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bp + 2*cp + dp + pp + 2*qp + rp)*(bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp) + 4*(2*a + bm + bp)*(bp + 2*cp + dp - 2*ip - 4*jp - 2*kp + pp + 2*qp + rp)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dcp*dcp;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bp + 2*cp + dp + pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*djp*djp;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dqp*dqp;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*ddp*ddp;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bp + 2*cp + dp + pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dkp*dkp;
	    dprod15=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*dpp*dpp;
	    dprod16=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))))*drp*drp;
	   
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dcp*dcp;
	    dI22=-((64*(bp + 2*cp + dp + pp + 2*qp + rp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*djp*djp;
	    dI23=((32*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dqp*dqp;
	    dI24=((8*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dbp*dbp;
	    dI25=((8*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*ddp*ddp;
	    dI26=-((16*(bp + 2*cp + dp + pp + 2*qp + rp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dip*dip;
	    dI27=-((16*(bp + 2*cp + dp + pp + 2*qp + rp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dkp*dkp;
	    dI28=((8*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*dpp*dpp;
	    dI29=((8*(ip + 2*jp + kp))/pow((bp + 2*cp + dp + 2*ip + 4*jp + 2*kp + pp + 2*qp + rp),3))*drp*drp;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //  data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
       	  }

       	if(ir==ic+(nbinshift+1)) //diag-15
	  {
            dprod1=-((32*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dhp*dhp;
	    dprod4=-((2*(8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(bm + 2*cm + dm + pm + 2*qm + rm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm + pm + 2*qm + rm)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm) + 4*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dgm*dgm;
	    dprod6=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(im + 2*jm + km) - 16*(gm + gp + 2*hm + 2*hp + im + ip)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(im + 2*jm + km)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm) + 8*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(bm + 2*cm + dm - 2*im - 4*jm - 2*km + pm + 2*qm + rm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm)))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dqm*dqm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bm + 2*cm + dm + pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*djm*djm;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dcm*dcm;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*drm*drm;
	    dprod14=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dpm*dpm;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(bm + 2*cm + dm + pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*dkm*dkm;
	    dprod16=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))))*ddm*ddm;
	   
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dqm*dqm;
	    dI22=-((64*(bm + 2*cm + dm + pm + 2*qm + rm))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*djm*djm;
	    dI23=((32*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dcm*dcm;
	    dI24=((8*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*drm*drm;
	    dI25=((8*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dpm*dpm;
	    dI26=-((16*(bm + 2*cm + dm + pm + 2*qm + rm))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dkm*dkm;
	    dI27=-((16*(bm + 2*cm + dm + pm + 2*qm + rm))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dim*dim;
	    dI28=((8*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*ddm*ddm;
	    dI29=((8*(im + 2*jm + km))/pow((bm + 2*cm + dm + 2*im + 4*jm + 2*km + pm + 2*qm + rm),3))*dbm*dbm;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
      	  }

		 if(ir==ic-((2*nbinshift)-1)) //diag+27
	  {
            dprod1=-((32*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dbp*dbp;
	    dprod8=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(lp + 2*mp + np) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(lp + 2*mp + np)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up) - 4*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(ep + 2*fp + gp - 2*lp - 4*mp - 2*np + sp + 2*tp + up))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up)))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dfp*dfp;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ep + 2*fp + gp + sp + 2*tp + up))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dmp*dmp;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dtp*dtp;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dep*dep;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ep + 2*fp + gp + sp + 2*tp + up))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dlp*dlp;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ep + 2*fp + gp + sp + 2*tp + up))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dnp*dnp;
	    dprod16=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lp + 2*mp +  np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dsp*dsp;
	    dprod17=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lp + 2*mp + np))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))))*dup*dup;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dfp*dfp;
	    dI22=-((64*(ep + 2*fp + gp + sp + 2*tp + up))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dmp*dmp;
	    dI23=((32*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dtp*dtp;
	    dI24=((8*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dep*dep;
	    dI25=((8*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dgp*dgp;
	    dI26=-((16*(ep + 2*fp + gp + sp + 2*tp + up))/pow((ep + 2*fp + gp + 2*lp +  4*mp + 2*np + sp + 2*tp + up),3))*dlp*dlp;
	    dI27=-((16*(ep + 2*fp + gp + sp + 2*tp + up))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dnp*dnp;
	    dI28=((8*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dsp*dsp;
	    dI29=((8*(lp + 2*mp + np))/pow((ep + 2*fp + gp + 2*lp + 4*mp + 2*np + sp + 2*tp + up),3))*dup*dup;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16+dprod17 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	     data[i]=covar;
      	  }

	if(ir==ic+((2*nbinshift)-1))//diag-27
	  {

	    dprod1=-((32*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dim*dim;
	    dprod5=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(lm + 2*mm + nm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(lm + 2*mm + nm)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um) -  4*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(em + 2*fm + gm - 2*lm - 4*mm - 2*nm + sm + 2*tm + um))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um)))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dtm*dtm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(em + 2*fm + gm + sm + 2*tm + um))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dmm*dmm;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dfm*dfm;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dum*dum;
	    dprod14=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dsm*dsm;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(em + 2*fm + gm + sm + 2*tm + um))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dnm*dnm;
	    dprod16=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(em + 2*fm + gm + sm + 2*tm + um))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dlm*dlm;
	    dprod17=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(lm + 2*mm + nm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))))*dem*dem;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;


	    dI21=((32*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dtm*dtm;
	    dI22=-((64*(em + 2*fm + gm + sm + 2*tm + um))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dmm*dmm;
	    dI23=((32*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dfm*dfm;
	    dI24=((8*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dum*dum;
	    dI25=((8*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dsm*dsm;
	    dI26=-((16*(em + 2*fm + gm + sm + 2*tm + um))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dnm*dnm;
	    dI27=-((16*(em + 2*fm + gm + sm + 2*tm + um))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dlm*dlm;
	    dI28=((8*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dgm*dgm;
	    dI29=((8*(lm + 2*mm + nm))/pow((em + 2*fm + gm + 2*lm + 4*mm + 2*nm + sm + 2*tm + um),3))*dem*dem;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16+dprod17 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	     data[i]=covar;		
	  }

	   	if(ir==ic-(2*nbinshift))//diag+28
	  {
            dprod1=-((32*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*da*da;
	    dprod3=((8*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(np + 2*op + pp)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp) - 4*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp)))*dbp*dbp;
	    dprod8=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(np + 2*op + pp)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp) - 4*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dgp*dgp;
	    dprod9=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(np + 2*op + pp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(np + 2*op + pp)*(gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp) - 4*(2*a + bm + bp)*(gp + 2*hp + ip - 2*np - 4*op - 2*pp + up + 2*vp + wp)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dip*dip;
	    dprod10=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip + up + 2*vp + wp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dop*dop;
	    dprod11=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(np + 2*op + pp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dvp*dvp;
	    dprod12=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip + up + 2*vp + wp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dnp*dnp;
	    dprod13=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gp + 2*hp + ip + up + 2*vp + wp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dpp*dpp;
	    dprod14=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(np + 2*op + pp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dup*dup;
	    dprod15=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(np + 2*op + pp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))))*dwp*dwp;

	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;


	    dI21=((32*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dhp*dhp;
	    dI22=-((64*(gp + 2*hp + ip + up + 2*vp + wp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dop*dop;
	    dI23=((32*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dvp*dvp;
	    dI24=((8*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dgp*dgp;
	    dI25=((8*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dip*dip;
	    dI26=-((16*(gp + 2*hp + ip + up + 2*vp + wp))/pow((gp + 2*hp + ip + 2*np +  4*op + 2*pp + up + 2*vp + wp),3))*dnp*dnp;
	    dI27=-((16*(gp + 2*hp + ip + up + 2*vp + wp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dpp*dpp;
	    dI28=((8*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dup*dup;
	    dI29=((8*(np + 2*op + pp))/pow((gp + 2*hp + ip + 2*np + 4*op + 2*pp + up + 2*vp + wp),3))*dwp*dwp;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	     data[i]=covar;
      	  }

	if(ir==ic+(2*nbinshift))//diag-28
	  {

	    dprod1=((8*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(nm + 2*om + pm)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm) - 4*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*dhp*dhp;
	    dprod4=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(nm + 2*om + pm)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm) - 4*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dim*dim;
	    dprod5=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(nm + 2*om + pm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(nm + 2*om + pm)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm) - 4*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(gm + 2*hm + im - 2*nm - 4*om - 2*pm + um + 2*vm + wm))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm)))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(nm + 2*om + pm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dvm*dvm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gm + 2*hm + im + um + 2*vm + wm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dom*dom;
	    dprod12=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(nm + 2*om + pm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dwm*dwm;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(nm + 2*om + pm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dum*dum;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gm + 2*hm + im + um + 2*vm + wm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dpm*dpm;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(gm + 2*hm + im + um + 2*vm + wm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))))*dnm*dnm;

	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;


	    dI21=((32*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dvm*dvm;
	    dI22=-((64*(gm + 2*hm + im + um + 2*vm + wm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dom*dom;
	    dI23=((32*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dhm*dhm;
	    dI24=((8*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dwm*dwm;
	    dI25=((8*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dum*dum;
	    dI26=-((16*(gm + 2*hm + im + um + 2*vm + wm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dpm*dpm;
	    dI27=-((16*(gm + 2*hm + im + um + 2*vm + wm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dnm*dnm;
	    dI28=((8*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dim*dim;
	    dI29=((8*(nm + 2*om + pm))/pow((gm + 2*hm + im + 2*nm + 4*om + 2*pm + um + 2*vm + wm),3))*dgm*dgm;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //  data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;		
	  }

       	if(ir==ic-((2*nbinshift)+1))//diag+29
	  {
            dprod1=-((32*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dhp*dhp;
	    dprod4=-((8*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp)))*dgp*dgp;
	    dprod9=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(pp + 2*qp + rp) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pp + 2*qp + rp)*(ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp) - 4*(2*a + bm + bp)*(ip + 2*jp + kp - 2*pp - 4*qp - 2*rp + wp + 2*xp + yp)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*djp*djp;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp + wp + 2*xp + yp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dqp*dqp;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pp + 2*qp +rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dxp*dxp;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dkp*dkp;
	    dprod14=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp + wp + 2*xp + yp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dpp*dpp;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(ip + 2*jp + kp + wp + 2*xp + yp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*drp*drp;
	    dprod16=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dwp*dwp;
	    dprod17=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pp + 2*qp + rp))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))))*dyp*dyp;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;


	    dI21=((32*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*djp*djp;
	    dI22=-((64*(ip + 2*jp + kp + wp + 2*xp + yp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dqp*dqp;
	    dI23=((32*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dxp*dxp;
	    dI24=((8*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dip*dip;
	    dI25=((8*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dkp*dkp;
	    dI26=-((16*(ip + 2*jp + kp + wp + 2*xp + yp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dpp*dpp;
	    dI27=-((16*(ip + 2*jp + kp + wp + 2*xp + yp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*drp*drp;
	    dI28=((8*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dwp*dwp;
	    dI29=((8*(pp + 2*qp + rp))/pow((ip + 2*jp + kp + 2*pp + 4*qp + 2*rp + wp + 2*xp + yp),3))*dyp*dyp;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16+dprod17 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    // data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	    data[i]=covar;
      	  }

	if(ir==ic+((2*nbinshift)+1))//diag-29
	  {

	    dprod1=-((32*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dhm*dhm;
	    dprod2=((64*(gm + gp + 2*hm + 2*hp + im + ip)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*da*da;
	    dprod3=-((32*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm +  2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dhp*dhp;
	    dprod4=((2*(4*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),2))*(pm + 2*qm + rm) + 16*(2*a + bm + bp)*(4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pm + 2*qm + rm)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym) - 4*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),2))))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dim*dim;
	    dprod5=-((8*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dgm*dgm;
	    dprod6=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dbm*dbm;
	    dprod7=((16*(gm + gp + 2*hm + 2*hp + im + ip)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dbp*dbp;
	    dprod8=-((8*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dgp*dgp;
	    dprod9=-((8*(2*a + bm + bp)*(im + 2*jm + km - 2*pm - 4*qm - 2*rm + wm + 2*xm + ym))/(pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3)*(im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym)))*dip*dip;
	    dprod10=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dxm*dxm;
	    dprod11=-((64*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km + wm + 2*xm + ym))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dqm*dqm;
	    dprod12=((32*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*djm*djm;
	    dprod13=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dym*dym;
	    dprod14=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dwm*dwm;
	    dprod15=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km + wm + 2*xm + ym))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*drm*drm;
	    dprod16=-((16*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(im + 2*jm + km + wm + 2*xm + ym))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dpm*dpm;
	    dprod17=((8*(4*a + 2*bm + 2*bp - gm - gp - 2*hm - 2*hp - im - ip)*(pm + 2*qm + rm))/((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip)*(pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))))*dkm*dkm;
	    
	    dI11=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhm*dhm;
	    dI12=-((64*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*da*da;
	    dI13=((32*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dhp*dhp;
	    dI14=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im +  ip),3))*dim*dim;
	    dI15=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgm*dgm;
	    dI16=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbm*dbm;
	    dI17=-((16*(gm + gp + 2*hm + 2*hp + im + ip))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dbp*dbp;
	    dI18=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dgp*dgp;
	    dI19=((8*(2*a + bm + bp))/pow((4*a + 2*bm + 2*bp + gm + gp + 2*hm + 2*hp + im + ip),3))*dip*dip;

	    dI21=((32*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dxm*dxm;
	    dI22=-((64*(im + 2*jm + km + wm + 2*xm + ym))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dqm*dqm;
	    dI23=((32*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*djm*djm;
	    dI24=((8*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dym*dym;
	    dI25=((8*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dwm*dwm;
	    dI26=-((16*(im + 2*jm + km + wm + 2*xm + ym))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*drm*drm;
	    dI27=-((16*(im + 2*jm + km + wm + 2*xm + ym))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dpm*dpm;
	    dI28=((8*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dkm*dkm;
	    dI29=((8*(pm + 2*qm + rm))/pow((im + 2*jm + km + 2*pm + 4*qm + 2*rm + wm + 2*xm + ym),3))*dim*dim;


	    
            covar=( (hr_sum->GetBinContent(ir)*hr_sum->GetBinContent(ic)) + 0.5*( dprod1+dprod2+dprod3+dprod4+dprod5+dprod6+dprod7+dprod8+dprod9+dprod10+dprod11+dprod12+dprod13+dprod14+dprod15+dprod16+dprod17 ) ) - (( hr_sum->GetBinContent(ir) + 0.5*( dI11+dI12+dI13+dI14+dI15+dI16+dI17+dI18+dI19 ) )*( hr_sum->GetBinContent(ic) + 0.5*( dI21+dI22+dI23+dI24+dI25+dI26+dI27+dI28+dI29 ) ));
	    
	    //  data[i]=covar/(hr_sum->GetBinError(ir)*hr_sum->GetBinError(ic));
	     data[i]=covar;		
	  }
				 
	
      }
  cout<<"flgs "<<flg<<" "<<endl;
  // cout<<"before cov[1000][1000] "<<cov[1000][1000]<<endl;
  //cov.SetMatrixArray(data.GetArray());
  //cout<<"after cov[1000][1000] "<<cov[1000][1000]<<endl;

  //cov.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
   
  //   cov.Print();
  //  cov.SetTol(1.e-50);
    // Double_t det1;
    //cov.Invert(&det1);
    //   cov.Print();


    cout<<mdim<<endl;
    

    // return cov;
    return data;

}

TH1D* fast_rotation_smoothing(TH1D *h_sum)
{
    hr_noshift=(TH1D*)h_sum->Clone();
	
    hr_shift=(TH1D*)hr_noshift->Clone();
    hr_shift->Reset();

    for(int ibin=1;ibin<=hr_shift->GetNbinsX();ibin++)
       {
	hr_shift->SetBinContent(ibin,hr_noshift->GetBinContent(ibin+1));
	hr_shift->SetBinError(ibin,hr_noshift->GetBinError(ibin+1));
       }

    hcomp_sum=(TH1D*)h_sum->Clone();
    hcomp_sum->Reset();

    for(int ibin=1;ibin<=h_sum->GetNbinsX();ibin++)
       {
	hcomp_sum->SetBinContent(ibin,0.5*(hr_noshift->GetBinContent(ibin)+hr_shift->GetBinContent(ibin)));
	hcomp_sum->SetBinError(ibin,0.5*sqrt((hr_noshift->GetBinError(ibin)*hr_noshift->GetBinError(ibin))+(hr_shift->GetBinError(ibin)*hr_shift->GetBinError(ibin))));
       }

    return hcomp_sum;

}


TH1D* construct_rhist_copy(TH1D *h_sum, TH1D *ho_sum)
{
 
  //  h_sum->Rebin(4);//This brings the bin width of the calo summed histograms to 150 ns
  // h_sum->Scale(0.25);
   
  cout<<"check "<<h_sum->GetBinError(500)<<endl;
   //create the component histograms of the ratio histograms   
   h1=(TH1D*)h_sum->Clone();
   h2=(TH1D*)h_sum->Clone();
   
   hp=new TH1D("calo histogram sum plus", "h_sum plus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
   
   hm=new TH1D("calo histogram sum minus", "h_sum minus", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

     nbinshift=(0.5*T_a_true)/h_sum->GetBinWidth(1);
      /*     if(double(0.5*T_a_true/h_sum->GetBinWidth(1))-nbinshift>0.5)
     {
       nbinshift=nbinshift+1;
     }
      */
   //  nbinshift=14;
   cout<<nbinshift<<endl;
   
   
   for(int ibin=0; ibin<=h_sum->GetNbinsX(); ibin++)             
     {
       hp->SetBinContent( ibin, h_sum->GetBinContent( ibin + nbinshift));
       hp->SetBinError(ibin, h_sum->GetBinError( ibin + nbinshift));
     }

   for(int ibin=nbinshift; ibin<=h_sum->GetNbinsX(); ibin++)
     {
       hm->SetBinContent( ibin, h_sum->GetBinContent( ibin - nbinshift));
       hm->SetBinError( ibin, h_sum->GetBinError(ibin - nbinshift));
     }
 

   // assign the correct weights to the 4 histohgrams
   T_a=2*nbinshift*h_sum->GetBinWidth(1);
   double flife=exp((T_a)/(2*lifetime));
   double blife=exp(-(T_a)/(2*lifetime));
   double deno=2+flife+blife;

   h1->Scale(1/deno);
   h2->Scale(1/deno);
   hp->Scale(flife/deno);
   hm->Scale(blife/deno);
   
   if(usesetbins)
     {
      h_ratio=new TH1D("calo_histogram_sum_ratio", "h_ratio", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
       
      for(int ibin=h_sum->FindBin(0);ibin<=h_sum->GetNbinsX();ibin++)
     {
         h_ratio->SetBinContent(ibin,(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)-hp->GetBinContent(ibin)-hm->GetBinContent(ibin))/(h1->GetBinContent(ibin)+h2->GetBinContent(ibin)+hp->GetBinContent(ibin)+hm->GetBinContent(ibin)));
	 //	 cout<<h_ratio->GetBinContent(ibin)<<endl;
       //  h_ratio->SetBinContent(ibin,((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)-(float)hp->GetBinContent(ibin)-(float)hm->GetBinContent(ibin))/((float)h1->GetBinContent(ibin)+(float)h2->GetBinContent(ibin)+(float)hp->GetBinContent(ibin)+(float)hm->GetBinContent(ibin)));
     }

   }
   else
     {
      h_num=new TH1D("calo_histogram_sum_numerator", "h_num", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));
      h_deno=new TH1D("calo histogram sum denominator", "h_deno", h_sum->GetNbinsX(), h_sum->GetBinLowEdge(1), h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1));

      h_num->Add(h1,1);
      h_num->Add(h2,1);
      h_num->Add(hp,-1);
      h_num->Add(hm,-1);

      h_deno->Add(h1,1);
      h_deno->Add(h2,1);
      h_deno->Add(hp,1);
      h_deno->Add(hm,1);
     
      h_num->Divide(h_deno);

     }
    //Assign error bars
   
 if(useerrorbars)
     {
       long double et1,det1,et2,det2,et3,det3,et4,det4,et5,det5,et6,det6,et7,det7,et8,det8,et9,det9;
       long double t1,t2,t3,t4,t5,t6,t7,t8,t9;
      for(int ibin=h_sum->FindBin(0); ibin<=h_sum->GetNbinsX()-nbinshift; ibin++)
      {
	et1=ho_sum->GetBinContent(2*ibin);
        det1=ho_sum->GetBinError(2*ibin);
	et2=ho_sum->GetBinContent((2*ibin)+1);
	det2=ho_sum->GetBinError((2*ibin)+1);
	et3=ho_sum->GetBinContent((2*ibin)-1);
	det3=ho_sum->GetBinError((2*ibin)-1);
	et4=ho_sum->GetBinContent((2*ibin)+27);
	det4=ho_sum->GetBinError((2*ibin))+27;
	et5=ho_sum->GetBinContent((2*ibin)-27);
	det5=ho_sum->GetBinError((2*ibin)-27);
	et6=ho_sum->GetBinContent((2*ibin)+28);
	det6=ho_sum->GetBinError((2*ibin)+28);
	et7=ho_sum->GetBinContent((2*ibin)-28);
	det7=ho_sum->GetBinError((2*ibin)-28);
	et8=ho_sum->GetBinContent((2*ibin)+29);
	det8=ho_sum->GetBinError((2*ibin)+29);
	et9=ho_sum->GetBinContent((2*ibin)-29);
	det9=ho_sum->GetBinError((2*ibin)-29);

	
	t1=((8*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det1;

	t2=((4*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det2;
	
	t3=((4*(et4 + et5 + 2*et6 + 2*et7 + et8 + et9))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det3;
	
	t4=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det4;

	t5=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det5;

        t6=-((8*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det6;

	t7=-((8*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det7;

	t8=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det8;

	t9=-((4*(2*et1 + et2 + et3))/pow((4*et1 + 2*et2 + 2*et3 + et4 + et5 + 2*et6 + 2*et7 + et8 + et9),2))*det9;
 	
       if(usesetbins)
       {
	 h_ratio->SetBinError( ibin, sqrt( ( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9)) ) );
	//	cout<<( (t1*t1) + (t2*t2) + (t3*t3) )<<endl;
       }
       else
       {
	 h_num->SetBinError(ibin,  sqrt( ( (t1*t1) + (t2*t2) + (t3*t3) + (t4*t4) + (t5*t5) + (t6*t6) + (t7*t7) + (t8*t8) + (t9*t9))  ) );
       }
      }
     }

 //send a common histogram name for fitting
 if(usesetbins)
      {
       hcalo=(TH1D*)h_ratio->Clone();
      }
   else
      {
       hcalo=(TH1D*)h_num->Clone();
      }

 // hcalo->Rebin(8);//This brings the bin width of the calo summed histograms to 150 ns
 //hcalo->Scale(0.125);

 /* foutput=new TFile("ratio_4hist_copy.root","new");
 h_ratio->Write();
 foutput->Close();
 */
 
 return hcalo;
  }


void ratio_FR_cov_construct_perhslice()
{

    _file[1]=TFile::Open(root_file_name);

  for(int i=0;i<=5;i++)
    {
     h[i+1]= (TH1D*)(_file[1]->FindObjectAny(TString::Format("qHist1D_sig_hslice_%d", i)));
    }
  /*h[24]->Fit("gaus","","",104750,104850);
  TF1 *fit = h[24]->GetFunction("gaus");
  inject_time=fit->GetParameter(1);
  */
  inject_time=104800;
  
  h_sum= new TH1D("calo_histogram_sum", "h_sum", h[1]->GetNbinsX(), 100001, 352001);
  h_sum->Sumw2(kTRUE);
  for(Int_t i=1; i<=6; i++)
    {
      h_sum->Add(h[i],1); 
    }

   double binwidth=h_sum->GetBinWidth(1);
   int nbins=h_sum->GetNbinsX();
   double binlowedge=h_sum->GetBinLowEdge(1);
   double binhighedge=h_sum->GetBinLowEdge(nbins)+binwidth;
   
   cout<<binwidth<<" "<<nbins<<" "<<binlowedge<<" "<<binhighedge<<" "<<endl;

   
   for(int i=1;i<=6;i++)
     {
      h[i]->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
     }
   h_sum->SetBins(nbins, rawBinToNs*(binlowedge-inject_time), rawBinToNs*(binhighedge-inject_time));
   cout<<h_sum->GetBinWidth(1)<<" "<<h_sum->GetNbinsX()<<" "<<h_sum->GetBinLowEdge(1)<<" "<<(h_sum->GetBinLowEdge(h_sum->GetNbinsX())+h_sum->GetBinWidth(1))<<" "<<endl;

    
        h_sum->Rebin(4);
        h_sum->Scale(0.25);

	for(int i=1;i<=6;i++)
	  {
	    h[i]->Rebin(4);
	    h[i]->Scale(0.25);

     
	    hcomp[i]=fast_rotation_smoothing(h[i]);
	
            hcomp[i]->Rebin(2);
            hcomp[i]->Scale(0.5);

	    	  
            sprintf(histname2,"hcalo_ratio_%d",i);
	
           hr[i]= new TH1D(histname2, histname2, hcomp[i]->GetNbinsX(),hcomp[i]->GetBinLowEdge(1),hcomp[i]->GetBinLowEdge(hcomp[i]->GetNbinsX())+hcomp[i]->GetBinWidth(1));

	   htemp = construct_rhist_copy(hcomp[i], h[i]);

	   for(int ibin=1;ibin<=htemp->GetNbinsX();ibin++)
	     {
	       hr[i]->SetBinContent(ibin,htemp->GetBinContent(ibin));
	       hr[i]->SetBinError(ibin,htemp->GetBinError(ibin));
	     }
           htemp->Reset();
           cout<<"check "<<h[i]->GetBinError(500)<<endl;

	  
	   hr[i]->GetYaxis()->SetTitle("ADC counts");
	   hr[i]->GetXaxis()->SetTitle("time [ns]");
           hr[i]->GetXaxis()->SetRangeUser(fit_start,fit_stop);

  
 
	   // h_sum->Rebin(2);
	   //h_sum->Scale(0.5);

	   //	   cov2.ResizeTo(hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), -1);
   

	   // cov2.SetTol(1.e-23);

           data2=ratio_cov_matrix(hr[i],h[i],fit_start);

	   TMatrixD cov2(mdim,mdim);
		
           cout<<"before cov[1000][1000] "<<cov2[1000][1000]<<endl;
           cov2.SetMatrixArray(data2.GetArray());
           cout<<"after cov[1000][1000] "<<cov2[1000][1000]<<endl;

           //cov2.ResizeTo(hr_sum->FindBin(fit_start),hr_sum->FindBin(fit_stop), hr_sum->FindBin(fit_start), hr_sum->FindBin(fit_stop), -1);
	   cov2.ResizeTo(hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), hr[i]->FindBin(fit_start), hr[i]->FindBin(fit_stop), -1);
   
	   // cov2.Print();
           cov2.SetTol(1.e-50);


	   cout<<"calo "<<i<<endl; 


	   sprintf(histname,"hcov_%d",i);
	
           hcov[i]= new TH2D(histname, histname, 2100, 0.0, 2100.0, 2100, 0.0, 2100.0);


	    	      
           for(int irow=hr[i]->FindBin(fit_start); irow<=hr[i]->FindBin(fit_stop); irow++)
           {
	    for(int icol=hr[i]->FindBin(fit_start); icol<=hr[i]->FindBin(fit_stop); icol++)
	    {
	     hcov[i]->SetBinContent(irow,icol,cov2(irow,icol));
	    }
           }
  
  
	   cov2.Clear();
	// cov.Clear("C");
        //cov.UnitMatrix();
 
   }
  
   
    foutput=new TFile("run2_4hcopy_ratioFR_perhslice_cov_mat_final.root","new");
   for(int i=1;i<=6;i++)
     {
      hr[i]->Write();
      hcov[i]->Write();
     }
   foutput->Close();
     
   
}
