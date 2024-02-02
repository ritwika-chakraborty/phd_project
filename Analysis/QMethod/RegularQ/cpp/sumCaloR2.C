#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TStyle.h"
#include "TMath.h" 
#include "TSpline.h"
#include "TMinuit.h"

#include "Blinders/Blinders.hh"
char blindname[128] = "Tim's 9day";
//char blindname[128] = "Tim's 60hr";

blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;

blinding::Blinders getBlinded( ftype, blindname );
//blinding::Blinders getBlinded( ftype, "daisydog" );
//blinding::Blinders getBlinded( ftype, "Tim's 60hr" );
//blinding::Blinders getBlinded( ftype, "Tim's 9day" );
//blinding::Blinders getBlinded( ftype, "Tim get blinded!" );
//blinding::Blinders getBlinded( ftype, "Tim is having fun getting his omega_a blinded!" );

TH1F *hmuonloss;
TH1F *hmlcorrection, *higcorrection, *hlostmuontimedist, *hdmcorrection;
TH1D *hFRint;

TH1D *hpedestaldrift;
TF1 *p0f;
TF1 *caloslowtermsf[24];

TFile *fpedestaldrift;
TFile *fmuonloss;
TFile *filefr;

// for contour plot
TGraph *gr_infillgain_muonloss1, *gr_infillgain_muonloss2, *gr_infillgain_muonloss3, *gr_infillgain_muonloss4, *gr_infillgain_muonloss5;

// 800 MSPS to ns conversion *(39.998/40.000) from docdb #9599
//double rawBinToNs = 1./(0.8*(39.998/40.000));
// 800 MSPS to ns conversion
double rawBinToNs = 1./0.8;

bool dbg = 0;

double	minTglobal;
double	maxTglobal;

const double precisionR = 1e-6;
const double refFreq = 0.2291; // MHz (from Kim and analysis blinding in bnl 821)
//const double refFreq = 0.2290739; // MHz (figure 38, page 81 of bnl  final report)
// above is defn used for ref. frequency in E821 final report page 82
// where omega_a = 2 pi  0.2291 MHz ( 1 + (R - deltaR)x10^-6)
double deltaR = 0;

double omega_factor = 1.25;

ROOT::Fit::Fitter fitter;
TVirtualFitter *gFitter = 0;
//TFitResultPtr fitres;

int igap = 11; // left-right gap
//int igap = 00; // left-right gap
int inoise = 10; // noise rejection parameter
int iswrb = 1; // software rebining factor - used to select root histograms  - 1,2,4
int iwind = 4; // left / right pedestal window - used to select root histograms - 4,8,12
int ithreshold = 15; // threshold setting - used to select root file - 8,12,16,20,24,28,32
int icalo = 0; // calorimeter index for "full" to select all calo (=0) versus individual calos (>0)
double infillgainscaling = 0.0; // scaling factor for multiplicative infill gain correction
bool fastRotationCorrection = 0; // switch on / off the fast rotation correction
bool fastRotationAnalysis = 0; // switch on / off the fast rotation (75 ns binning) analysis
double scalefastrotation = 1.0; // for scaling the fast roation correction 
bool FFTon = 1; // switch on / off FFT of residuals
int firstcaloinslice = 1;
int lastcaloinslice = 24;
int swrbglobal = -1; // global use of software rebinning, initialized to invalid rebinning
bool fracResiUnits = 1; // fractional or absolute units on residuals histogram

bool run_2 = 0; // true for run 2
bool fixbug = 0; // fix the bug with binning of earlier sorts of 60hr dataset
bool doPhaseShift = 0; // true = add phase shidt from vertical drift analysis
bool doExpRelax = 1; // type of relation - exponential of clone of N(t) vertical drift
bool doPUErrorBarCorr = 0; // apply the simulation calculated corrections for error bars
bool doBinInts = 0; // fit to intergrals of function over bin width
bool doShiftBins = 0; // to shift bins of time distribution by -1, 75ns bin
bool doBlendBins = 0; // rebinned bins are y_n + 0.5( y_n-1 + y_n+1 ) not y_n+y_n+1
bool doUserCovariance = 0; // user calculation of covariance matrix for chi2
bool doMinos = 0; // user calculation of covariance matrix for chi2
bool doFRfit = 0; // adding fast rotation to precessf
bool doDraw = 1; // =0 avoid drawing fitting histograms for remote speed
bool doDraw2 = 1; // =0 avoid drawing vertical drift histograms for remote speed
bool doFRfitWithInts = 0; // doing fast rotation fit with bin integrals
bool doTimeShift = 0; /// do the time shift of the graph
double timeShiftMultiplier = 1.0; /// do the time shift of the graph
bool doConvoluted = 0; // if true fit the a 0.5 + 1.0 + 0.5 convoluted histogram  that cancels-out of Tcyc signa
bool doDriftCorrection = 0; // replay htmp with htmpdriftcorrection in fitting
TH1D *htmpdriftcorrection;


char starttimescantype[32] = "all"; // also near, far, hot, cold, top bottom
char thresholdscantype[32] = ""; // also core2, core4 for vertical centers of calos,  also puLH, puA for PU studies

//char aname[32] = "1"; // point rejection identifier - used to select root file - "1","1e9" (point rejecton on / off)
char aname[32] = "1e0"; // point rejection identifier - used to select root file - "1","1e9" (point rejecton on / off)
//char aname[32] = "1e9"; // point rejection identifier - used to select root file - "1","1e9" (point rejecton on / off)
char dname[32] = "QFillByFillAnalyzer"; // folder name for calo-by-calo histograms, e.g.  QFillByFillAnalyzer, CQBankAnalyzerTim
//char dname[32] = "CQBankAnalyzerTim"; // folder name for calo-by-calo histograms, e.g.  QFillByFillAnalyzer, CQBankAnalyzerTim
//char fname[128] = "noxtaldata/9day-DQC-noinfillgain-w4t15g11.root"; // filename, 9day testing, no global gain, no infill gain
//char fname[128] = "noxtaldata/v4w4n10t32prd1withgaincorrection-gap00.root"; // filename, 9day testing, no global gain, no infill gain
char fname[512] = "noxtaldata/v9.17-w4n10t15prd1e0infillwithstdpnoOOF-gap11-60hr-pufix-gold-withpederr.root"; // filename, 9day testing, no global gain, no infill gain
char fnamepedestal[512] = "noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-9d-pufix-gold-withpederr-gainbugfix-newenergyscale.root"; // filename, 9day testing, no global gain, no infill gain
//char fnamefastrotation[128] = "hcyc-9day-vbo-pw.root"; // fast roation correction file
//char fnamefastrotation[128] = "hcyc-60hr-vbo-pw.root"; // fast roation correction file
char fnamefastrotation[128] = "hcyc-60hr-vbo-pw-gold.root"; // fast roation correction file
//char fnamefastrotation[128] = "hcyc-9day-13par-pw.root"; // fast rotation correction file
//char fnamefastrotation[128] = "hcyc-60hr-13par-pw.root"; // fast rotation correction file
char fnamemuonloss[128] = "muonloss/60hours-DQC.root"; // muon loss correction file
//char fnamemuonloss[128] = "muonloss/60hrs-DQC-double-time.root"; // muon loss correction file
//char fnamemuonloss[128] = "muonloss/60hrs-DQC-triple-time.root"; // muon loss correction file
//char fitname[32] = "5par"; // type of fit to perform on time distribution, e.g 5par, wiggle only
//char fitname[32] = "9par"; // type of fit to perform on time distribution, e.g 5par, wiggle only, CBO norm amp, frq, phs, tau
//char fitname[32] = "13par"; // type of fit to perform on time distribution, e.g 5par, wiggle only, CBO norm amp, frq, phs, tau
//char fitname[32] = "vbo"; // type of fit to perform on time distribution, e.g 5par, wiggle only, CBO norm amp, frq, phs, tau
char fitname[32] = "vbo2"; // type of fit to perform on time distribution, e.g 5par, wiggle only, CBO norm amp, frq, phs, tau
char dataname[32] = "60hr"; // dataset identifier
//char dataname[32] = "9day"; // dataset identifier
char datasubset[128] = "all"; // dataset identifier

double timeMin = 500., timeMax = 520.; // for horizontal,vertical time ranges

/*
60hr calo sum
  17  A_vbo        2.19762e+06   6.58512e+05   1.01605e+02   1.71871e-10
  18  omega_vbo    1.77073e-02   1.14774e-05   4.11235e-08  -4.20249e+00
  19  phi_vbo      4.87942e+00   3.74993e-01   5.14291e-04   3.80564e-04
  23  tau_vbo      1.83632e+04   3.03015e+03   1.35854e+01  -3.54618e-08
  
  26  A_vbo2       1.01492e+07   4.51644e+07   4.29279e+03  -9.99802e-12
  27  omega_vbo2   1.33392e-02   7.76617e-05   1.53648e-07  -8.11118e-01
  28  tau_vbo2     6.82332e+03   7.48268e+03   7.18415e+00  -2.51850e-09
  29  phi_vbo2     1.27120e+00   2.14377e+00   4.23702e-03   4.88432e-05

9day calo sum
  17  A_vbo        2.13546e+06   5.08278e+05   4.72747e+00  -8.27205e-12
  18  omega_vbo    1.57447e-02   6.90120e-06   2.17757e-11  -9.47922e-02
  19  phi_vbo      1.87486e+00   2.39318e-01  -7.35474e-07   1.26059e-06
  23  tau_vbo      2.23763e+04   3.43926e+03  -3.25575e-02   4.90719e-10

  26  A_vbo2       1.99357e+05   5.20210e+04  -9.63130e-03  -3.16034e-11
  27  omega_vbo2   1.83318e-02   3.45980e-06   2.83700e-12   2.65050e-01
  28  tau_vbo2     5.29970e+05   1.04836e+06   3.19322e-01  -1.56849e-12
  29  phi_vbo2     5.88579e+00   2.45007e-01  -2.12323e-07   1.09946e-06

*/

//  magic momentum,  3.094 GeV/c,
//  muon mass,  0.1056583715  GeV
//  muon lifetime, 2.1969811 us
//  lorentz factor, sqrt(1.+3.094*3.094/0.1056583715/0.1056583715) = 29.3001
//  time-dilated lifetime, 29.300*2.1969811 = 64.3715 us
//  for 3.094 *1.001 = 3.09709 GeV/c, gamma = 29.3294, gamma*tau = 64.4361 us
//  for 3.094 *0.999 = 3.09709 GeV/c, gamma = 29.2709, gamma*tau = 64.3075 us
//  64.4361/64.3715 =1.001, 0.1% increase in momentum is 0.1% increase in gamma*tau

bool use_verticaldrift_correction = 0;
bool inflate_errors = 0; // for results of blind E

bool fix_R = 0; // fix blinded R
double set_R = -26.71; // set blinded R term
// -26.71 eg
// - 26.03 9d
// -24.70 hk
// 60hr -34.05

bool fix_bkd = 1; // fix bkd term
double set_bkd = 0.0; // set bkd term

bool fix_asym = 0; // fix ampl for asymmetry
double set_asym = 0.235; // set ampl for asymmetry

bool fix_asymtau = 1; // fix relaxation term for asymmetry
double set_asymtau = 1.e12; // set relaxation term for asymmetry

bool fix_2omega_a = 1; // fix 2omega_a term for asymmetry
double set_A_2omega_a = 0.0; // set 2omega_a term for asymmetry
double set_phi_2omega_a = 0.0; // set 2omega_a term for asymmetry

bool fix_gammatau = 0;
double set_gammatau = 6.4427e4;

bool fix_delta_omega_a = 1;
double set_delta_omega_a = 0.0;

bool fix_cbo_amp = 0;
double set_cbo_amp = 0.00;

bool fix_2cbo_amp = 1;
double set_2cbo_amp = 0.00;
bool fix_2cbo_phi = 1;
double set_2cbo_phi = 0.00;

bool fix_anomasym_amp = 0;
double set_anomasym_amp = 0.0;

bool fix_anomphi_amp = 1;
double set_anomphi_amp = 0.00;
bool fix_anomphi_phi = 1;
double set_anomphi_phi = 0.00;

bool fix_vbo_amp = 0;
bool fix_vbo_phi = 0;
bool fix_vbo_tau = 1;
bool fix_vbo_omega = 1;

bool fix_vbo_amp2 = 0;
bool fix_vbo_phi2 = 0;
bool fix_vbo_tau2 = 1;
bool fix_vbo_omega2 = 1;

bool fix_vbo_Aasym2 = 1;
bool fix_vbo_phiasym2 = 1;

bool fix_delta_omega_cbo = 1;
bool fix_omega_cbo = 0;

bool fix_pileup = 1;
double set_pileup = 0.0;

double set_p0 = 0.0, set_p1 = 0.0, set_p2 = 0.0, set_p3 = 0.0, set_p4 = 0.0, set_p5 = 0.0;

double set_amp_vbo_60hr = 3.38057e-03;
double set_amp_vbo_9day = -1.79829e-03;
double set_amp_vbo_endgame = 2.01815e-03;
double set_amp_vbo_hk = -1.79829e-03;
double set_amp_vbo_lk = -1.79829e-03;

double set_amp_vbo2_60hr = 3.60241e-04;
double set_amp_vbo2_9day = -2.63167e-04;
double set_amp_vbo2_endgame = -8.40008e-05;
double set_amp_vbo2_hk = -2.63167e-04;
double set_amp_vbo2_lk = -2.63167e-04;

double set_phi_vbo_60hr = 1.40407e+00;
double set_phi_vbo_9day = 1.55275e+00;
double set_phi_vbo_endgame = 1.13780e+00;
double set_phi_vbo_hk = 1.55275e+00;
double set_phi_vbo_lk = 1.55275e+00;

double set_phi_vbo2_60hr = 4.28849e+00;
double set_phi_vbo2_9day = 4.26609e+00;
double set_phi_vbo2_endgame = 2.61285e+00;
double set_phi_vbo2_hk = 4.26609e+00;
double set_phi_vbo2_lk = 4.26609e+00;

double set_omega_vbo_60hr = 1.; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
double set_omega_vbo_9day = 1.; // 0.0126? based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
double set_omega_vbo_endgame = 1.; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
double set_omega_vbo_hk = 1.; // 0.0126? based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
double set_omega_vbo_lk = 1.; // 0.0126? based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724

double set_omega_vbo2_60hr = 1.; // f_vw ~ f_bo for 60hr field index
double set_omega_vbo2_9day = 1.; // based on FFT of absolute residuals and peak at 556 -> T_vbo = 190020./556. = 341.763 ct, omega_cbo = 2pi/T_vbo = 0.018
double set_omega_vbo2_endgame = 1.; // f_vw ~ f_bo for 60hr field index
double set_omega_vbo2_hk = 1.; // based on FFT of absolute residuals and peak at 556 -> T_vbo = 190020./556. = 341.763 ct, omega_cbo = 2pi/T_vbo = 0.018
double set_omega_vbo2_lk = 1.; // based on FFT of absolute residuals and peak at 556 -> T_vbo = 190020./556. = 341.763 ct, omega_cbo = 2pi/T_vbo = 0.018

double set_tau_vbo_60hr = 1.8e4; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
double set_tau_vbo_9day = 1.6e4; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
double set_tau_vbo_endgame = 1.8e4; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
double set_tau_vbo_hk = 1.6e4; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
double set_tau_vbo_lk = 1.6e4; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724

double set_tau_vbo2_60hr = 1.8e4; // f_vw ~ f_bo for 60hr field index
double set_tau_vbo2_9day = 1.6e4; // based on FFT of absolute residuals and peak at 556 -> T_vbo = 190020./556. = 341.763 ct, omega_cbo = 2pi/T_vbo = 0.018
double set_tau_vbo2_endgame = 1.8e4; // f_vw ~ f_bo for 60hr field index
double set_tau_vbo2_hk = 1.6e4; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
double set_tau_vbo2_lk = 1.6e4; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724

// freq change parameters

// pre 1/9/2020 values
/*
double set_AExpTermCBO_60hr = 2.90;
double set_BExpTermCBO_60hr = 5.12;
double set_AExpTermCBO_9day = 2.86;
double set_BExpTermCBO_9day = 5.50;
double set_AExpTermCBO_hk = 2.86;
double set_BExpTermCBO_hk = 5.50;
double set_AExpTermCBO_lk = 2.86;
double set_BExpTermCBO_lk = 5.50;
double set_AExpTermCBO_eg = 6.626; // best fit calo data
double set_BExpTermCBO_eg = 6.821; // best fit calo data
*/
// post 1/9/2020 values
double set_AExpTermCBO_60hr = 2.79; // from tracker
double set_BExpTermCBO_60hr = 5.63; // from tracker
double set_AExpTermCBO_9day = 2.80;
double set_BExpTermCBO_9day = 6.18;
double set_AExpTermCBO_hk = 2.86; // TO DO
double set_BExpTermCBO_hk = 5.50; // TO DO
double set_AExpTermCBO_lk = 2.86; // TO DO
double set_BExpTermCBO_lk = 5.50; // TO DO
double set_AExpTermCBO_eg = 6.82; // from tracker
double set_BExpTermCBO_eg = 5.42; // from tracker

// pre 1/9/2020 values
/*
double set_AExpTauTermCBO_60hr = 2.90; // from tracker
double set_BExpTauTermCBO_60hr = 5.12; // from tracker
double set_AExpTauTermCBO_9day = 2.86;
double set_BExpTauTermCBO_9day = 5.50;
double set_AExpTauTermCBO_hk = 2.86;
double set_BExpTauTermCBO_hk = 5.50;
double set_AExpTauTermCBO_lk = 2.86;  // TO DO
double set_BExpTauTermCBO_lk = 5.50;  // TO DO
double set_AExpTauTermCBO_eg = 6.626; // best fit calo data
double set_BExpTauTermCBO_eg = 6.821; // best fit calo data
*/
 // post 1/9/2020 values
double set_AExpTauTermCBO_60hr = 2.79; // from tracker
double set_BExpTauTermCBO_60hr = 5.63; // from tracker
double set_AExpTauTermCBO_9day = 2.80;
double set_BExpTauTermCBO_9day = 6.18;
double set_AExpTauTermCBO_hk = 2.86; // TO DO
double set_BExpTauTermCBO_hk = 5.50; // TO DO
double set_AExpTauTermCBO_lk = 2.86; // TO DO
double set_BExpTauTermCBO_lk = 5.50; // TO DO
double set_AExpTauTermCBO_eg = 6.82; // from tracker
double set_BExpTauTermCBO_eg = 5.42; // from tracker



double set_delta_omega_cbo = 5.60678e-08;
double set_omega_cbo = 2.27430e-03;

bool use_cbo_gauss = 0;
bool fix_tau_cbo = 0;
double set_tau_cbo = 2.e5;

bool fix_fitring_amp = 1;
double set_fitring_amp = 0.0;
bool fix_fitring_tau = 1;
double set_fitring_tau = 2.654e+04;
bool fix_fitring_omega = 1;
double set_fitring_omega = 1.28e-02
;
bool fix_fitring2_amp = 1;
double set_fitring2_amp = 0.0;
bool fix_fitring2_tau = 1;
double set_fitring2_tau = 2.654e+04;
bool fix_fitring2_omega = 1;
double set_fitring2_omega = 1.49e-02;

bool fix_gain_amp = 1;
bool fix_gain_tau = 1;
double set_gain_amp = 0.03;
double set_gain_tau = 2.12812e+04;

bool fix_muonloss_amp = 1;
double set_muonloss_amp = 3.0e-10;

bool fix_detectmuon_amp = 1;
double set_detectmuon_amp = 0.0;

bool fix_pedestaldrift_amp = 1;
double set_pedestaldrift_amp = 0.0;

bool fix_AExpTermCBO = 1;
bool fix_BExpTermCBO = 1;
bool use_new_cbomodel = 1;

//#define Nfit 3 // new years + 11589etc + 11169 etc datasets
//#define Nfit 9 // vertical slices
//#define Nfit 6 // horizontal slices
//#define Nfit 8 // sequence number
//#define Nfit 48 // pileup scan - three threshold settings
//#define Nfit 36 // pileup scan - three threshold settings
//#define Nfit 12 // pileup scan - only one threshold
#define Nfit 40 // start time scan
//#define Nfit 9 // threshold, rebin, window scan
//#define Nfit 12 // opposite pair scan
//#define Nfit 24 // calorimeter scan
//#define Nfit 5 // threshold-only scan
//define Nfit 4 // window-only scan
//#define Nfit 7 // 60hr threshold scan, width=4
//#define Nfit 5 // 60hr noise scan, moise rms multipliers 4, 6, 8, 10, 12
//#define Nfit 62 // 60hr run scan, width=4

// no. of data points in fit scan, set Nfit, Ifit
int Ifit = Nfit;
int Nfstarttimescan = Nfit;

// global variables
int ICalo;
char hname[64], hsname[64], hpname[64];
char foutname[64];
TH1D *h1[54], *h2[54], *hSum, *htmp, *hRMS, *hCnts;
TH2D *s1[54], *s2[54], *sSum, *stmp, *sRMS, *sCnts;
TH2D *htop, *hbot;
TH1D *htop1D, *hbot1D;
TH1D *hpair[12];
TH1D *htmpshift;
TH1D *htmp1, *htmp2, *htmp3, *htmp4, *htmp5;
TH1D *htmp1b, *htmp2b, *htmp3b, *htmp4b, *htmp5b;
TH2D *h2tmp, *h2tmp2;
TH1D *hres, *hres2;
TH1 *hmcalo, *hmcalosum;
TH1D *hfastrotation;
TH1 *hm2 = 0; // for FFT

TFile *file0, *file1, *file2, *file3, *file4;

TGraphErrors *gWiggle;
TGraphErrors *gChifit, *gOmega_afit, *gTaufit, *gAfit, *gArelaxfit, *gPhifit, *gNfit, *gBfit, *gAGAINfit, *gtauGAINfit, *gmuonlossfit,  *gpeddriftfit, *gRPILEUPfit, *gACBOfit, *gtauCBOfit, *gfreqCBOfit, *gdfreqCBOfit, *gphaseCBOfit, *gAVWfit, *gphaseVWfit, *gAVW2fit, *gphaseVW2fit,  *gtauVWfit,  *gtauVW2fit, *gomegaVWfit,  *gomegaVW2fit; 
TGraphErrors *gACBO2fit, *gphase2CBOfit, *gACBO3fit, *gphase3CBOfit, *g2ACBOfit, *g2phaseCBOfit;
TGraphErrors *gACBOmeanX, *gACBOmeanY; 
TGraph *gCBOamp, *gCBOomga;
TGraph *gPed, *g1;
TF1 *precessf, *cbof, *frf, *pol0f, *ringf, *ring2f, *slowtermsf;
TF1 *muonlossf,  *ffastrotation;

//double paramToFreq(double blindedValue);

//double paramToFreq(double blindedValue){
//
//  double unblindedR = blindedValue - deltaR;
//
//  return 2 * TMath::Pi() * refFreq * (1 + (unblindedR * precisionR));
//}

Double_t delta_omega_vbo = 0.0;


Double_t fprec(Double_t *x, Double_t *par);

Double_t fprec(Double_t *x, Double_t *par)
  {
     // precession function with sinosoidal amplitude, freqeuncy and phase
     //
     // par[0] - N0
     // par[1] - tau_mu
     // par[2] - A_a
     // par[3] - R -> omega_a
     // par[4] - phi_a
     // par[5] - time-independent background
     // par[6] - A_cbo
     // par[7] - tau_cbo
     // par[8] - omega_cbo
     // par[9] - phi_cbo
     // par[10] - A_gain
     // par[11] - tau_gain
     
     Double_t xx =x[0];

     // injection start time
     Double_t t0 = 0.0;

     double R = par[3];
     double omega_a = 1.00e-3 * getBlinded.paramToFreq(R); // conversions from bins to ns and ns to us as Qmethod fit in units of ticks

     double delta_omega_a = par[29]; // omega_a linear time dependence
     omega_a = omega_a*( 1.0 + delta_omega_a * ( xx-t0 ) );

     // CBO linear time dependence (tracker team parameterization)
     double Acbo = par[23];
     double Bcbo = par[24];
     double omega_cbo, tauAcbo, tauBcbo;
     double delta_omega_cbo = par[20]; 

     if (use_new_cbomodel) { 
       //new model
      if ( strcmp( "60hr", dataname) == 0 ) tauAcbo = set_AExpTauTermCBO_60hr; // convert tracker parameter to ns
      if ( strcmp( "60hr", dataname) == 0 ) tauBcbo = set_BExpTauTermCBO_60hr; // convert tracker parameter to ns
      if ( strcmp( "9day", dataname) == 0 ) tauAcbo = set_AExpTauTermCBO_9day; // convert tracker parameter to ns
      if ( strcmp( "9day", dataname) == 0 ) tauBcbo = set_BExpTauTermCBO_9day; // convert tracker parameter to ns
      if ( strcmp( "hk", dataname) == 0 ) tauAcbo = set_AExpTauTermCBO_hk; // convert tracker parameter to ns
      if ( strcmp( "hk", dataname) == 0 ) tauBcbo = set_BExpTauTermCBO_hk; // convert tracker parameter to ns
      if ( strcmp( "lk", dataname) == 0 ) tauAcbo = set_AExpTauTermCBO_lk; // convert tracker parameter to ns
      if ( strcmp( "lk", dataname) == 0 ) tauBcbo = set_BExpTauTermCBO_lk; // convert tracker parameter to ns
      if ( strcmp( "endgame", dataname) == 0 ) tauAcbo = set_AExpTauTermCBO_eg; // convert tracker parameter to ns
      if ( strcmp( "endgame", dataname) == 0 ) tauBcbo = set_BExpTauTermCBO_eg; // convert tracker parameter to ns
       omega_cbo =  par[8] * ( 1.0 + Acbo*exp( -( xx-t0 ) / tauAcbo ) / ( par[8] * (xx - t0) ) + Bcbo*exp( -( xx-t0 ) / tauBcbo ) / ( par[8] * (xx - t0) ) );
     } else {
       // old model
       tauAcbo = 73.5e3; // convert tracker parameter to ns
       tauBcbo = 16.6e3; // convert tracker parameter to ns
       omega_cbo =  par[8] * ( 1.0 + delta_omega_cbo * ( xx-t0 ) + Acbo*exp( -( xx-t0 ) / tauAcbo ) + Bcbo*exp( -( xx-t0 ) / tauBcbo ) );
     }
     if (run_2) omega_cbo = par[8];


     double fig, figa;
     if (use_verticaldrift_correction) {

       // vertical drift from fit to difference of residuals from top / bottom of calo's w/o muloss, gain terms and fixed gamma*tau
       double p0, p1, p2, p3, p4, p5;
       if ( strcmp( "60hr", dataname) == 0 )  {
	 p0 = -0.0118174;
	 p1 =  7.37384e-07;
	 p2 = -1.44961e-11;
	 p3 =  1.25639e-16;
	 p4 = -5.13631e-22;
	 p5 =  8.03945e-28;
       }
       if ( strcmp( "9day", dataname) == 0 )  {
	 p0 =   -0.0115946;
	 p1 =  7.21943e-07;
	 p2 = -1.41151e-11;
	 p3 =  1.21412e-16;
	 p4 =  -4.9225e-22;
	 p5 =  7.63965e-28;
       }
       if ( strcmp( "endgame", dataname) == 0 )  {
	 p0 =  -0.00659176;
	 p1 =  5.03643e-07;
	 p2 = -9.70229e-12;
	 p3 =  7.75849e-17;
	 p4 = -2.97651e-22;
	 p5 =  4.45124e-28;
       }
       if ( strcmp( "hk", dataname) == 0 )  {
	 p0 =   -0.0115946;
	 p1 =  7.21943e-07;
	 p2 = -1.41151e-11;
	 p3 =  1.21412e-16;
	 p4 =  -4.9225e-22;
	 p5 =  7.63965e-28;
       }
       if ( strcmp( "lk", dataname) == 0 )  {
	 p0 =   -0.0115946;
	 p1 =  7.21943e-07;
	 p2 = -1.41151e-11;
	 p3 =  1.21412e-16;
	 p4 =  -4.9225e-22;
	 p5 =  7.63965e-28;
       }
       fig = 1.0 + par[10] * ( p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx + p4*xx*xx*xx*xx + p5*xx*xx*xx*xx*xx);
       figa = 1.0 + par[44] * ( p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx + p4*xx*xx*xx*xx + p5*xx*xx*xx*xx*xx);
     }

     double phaseshift = 0.0;
     if (doPhaseShift) {
       double p0 = set_p0, p1 = set_p1, p2 = set_p2, p3 = set_p3, p4 = set_p4, p5 = set_p5;

       /* 
       p0                        =  0.0; // p0 = 2.13544;
       p1                        =  5.55299e-07;
       p2                        = -1.38405e-11;
       p3                        =  1.54076e-16;
       p4                        = -7.82636e-22;
       p5                        =  1.47964e-27;
       */
      
       phaseshift =  p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx + p4*xx*xx*xx*xx + p5*xx*xx*xx*xx*xx;
     }

     double omega_cyc = 2.*TMath::Pi()/149.3;
     double omega_VW = par[17] * ( omega_cyc - 2 * omega_cbo * sqrt( ( 2.*omega_cyc/omega_cbo ) - 1. )); // from joe price
     double omega_y = par[26] * omega_cbo * sqrt( ( 2.*omega_cyc/omega_cbo ) - 1. ); // from joe price

     // time dependent horizontal cbo asymmetry term + vertical cbo asymmetry with exponential envelope
     Double_t A_a = par[2];
     if (!use_cbo_gauss) A_a *= (1.0 + par[12] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 ) + par[13] ) ); // horizontal CBO with freq dependence
     if (use_cbo_gauss) A_a *= (1.0 + par[12] * exp( -( xx - t0 )*( xx - t0 ) / par[7] / par[7] /2. ) * cos( omega_cbo * ( xx - t0 ) + par[13] ) ); // horizontal CBO with freq dependence
     A_a *= (1.0 + par[38] * exp( -( xx - t0 ) / par[27] ) * cos( omega_y * ( xx - t0 ) + par[39] ) ); // vertical CBO f_y (FIX ME!!!)
     //A_a = par[2]; // avoid CBO on assymmetry

     // phase of omega_a wiggle with addition of cbo term (2 parameters amplitude, phase)
     // phase shift of omega_a from vertical drift analysis
     Double_t phi_a = par[4] + phaseshift + par[14] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 ) + par[15] ); // with freq dependence

     //phi_a = par[4]; // avoid CBO on phase

     // 6-parameter wiggle with time-independent background (5-par w/o bkd)
     Double_t f = 0.0;
     f = par[0] * exp( -( xx - t0 ) / par[1] ) * ( 1.0 + A_a * exp( -( xx - t0 ) / par[44] ) * ( cos ( omega_a * ( xx - t0 ) + phi_a ) + par[45] * cos ( 2.*omega_a * ( xx - t0 ) + par[46] ) ) ) + par[5];
     if (!doExpRelax && use_verticaldrift_correction) f = par[0] * exp( -( xx - t0 ) / par[1] ) * ( 1.0 + A_a * figa * cos ( omega_a * ( xx - t0 ) + phi_a ) ) + par[5];

     //if (dbg) printf("f = %f\n", f);

     // time dependent cbo leading-order nomalization term - exponential envelope
     f *= (1.0 + par[6] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 )   +  par[9] ) );
     //if (dbg) printf("f = %f\n", f);

     // time dependent 2*CB0 waist term, freq  2*omega_cbo, tau same as CBO, + amp, phase parameters
     f *= (1.0 + par[36] * exp( -( xx - t0 ) / par[7] ) * cos( 2*omega_cbo * ( xx - t0 ) + par[37] ) );
     //if (dbg) printf("f = %f\n", f);

     // vertical BO modifier, VW = fc-2fy (has much smaller effect from change in field index)
     f *= ( 1.0 + par[16] * exp( -( xx - t0 ) / par[22] ) * cos( omega_VW * ( xx - t0 ) + par[18] ) );
     //if (dbg) printf("f = %f\n", f);

     // vertical BO modifier2, BO = fy, same tau as VW (has much larger effect from change in field index)
     f *= ( 1.0 + par[25] * exp( -( xx - t0 ) / par[27] ) * cos( omega_y * ( xx - t0 ) + par[28] ) ); // tau_VB = tau_VW
     //if (dbg) printf("par[25], par[26], par[27], par[28] %f, %f, %f, %f, f = %f\n", par[25], par[26], par[27], par[28], f);

     // injection induced, pedestal ringing 1, or proxy for oscillations
     f *= (1.0 + par[31] * exp( -( xx - t0 ) / par[33] ) * cos( par[32] * ( xx - t0 ) + par[34] ) );
     //if (dbg) printf("par[31], par[32], par[33], par[34] %f, %f, %f, %f, f = %f\n", par[31], par[32], par[33], par[34], f);

     // injection induced, pedestal ringing 2, or proxy for oscillations
     f *= (1.0 + par[40] * exp( -( xx - t0 ) / par[42] ) * cos( par[41] * ( xx - t0 ) + par[43] ) );
     //if (dbg) printf("par[39], par[40], par[41], par[42] %f, %f, %f, %f, f = %f\n", par[40], par[41], par[42], par[43], f);

     //  use empirical gain multiplier, Tim's  is just an exponential, Aaron's is also a wiggle
     // if (!use_verticaldrift_correction) fig = ( 1.0 + par[10] * exp( -( xx - t0 ) / par[11] ) ); // Tim's empirical infill gain 
     if (!use_verticaldrift_correction) fig = ( 1.0 + par[10] * exp( -( xx - t0 ) / par[1] ) * ( 1.0 + par[11] * cos ( omega_a * ( xx - t0 ) + phi_a ) ) ); // Aaron's empirical infill gain 

     // empirical vertical drift effect
     // June 18, I think that there's a bug here,  the fit to dift was done in clock ticks not ns.
     // this explains diff between tracker and calo in vertical position

     higcorrection->SetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fig ); // diagnostic     
     f *= fig;
     //if (dbg) printf("f = %f\n", f);

     // get pedestal drift histogram
     double fpd = ( 1.0 - par[35] * hpedestaldrift->GetBinContent( hpedestaldrift->GetXaxis()->FindBin( xx) ) ); 
     f *= fpd; // correct for pedestal drift / ringing

     // muon loss multiplier (corrects for effect on positrons of lost muons)- new
     // should i worry about discontinuities in histogram, I dont think so, tested with function as alternative
     //double fml = 1.0 - par[21] * muonlossf->Eval( 1.00e-3*(xx - t0) ); // replace histogram with smooth function, 1.00e-3 ns  to microsec's in Sudeshna muonloss histogram
     double fml = 1.0 - par[21] * hmuonloss->GetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ) ); // 1.00e-3 converts ns to microsec's in Sudeshna muonloss histogram
     hmlcorrection->SetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fml ); // diagnostic
     f *= fml;

     // muon addition multiplier (corrects for effect of detecting the lost muons)
     double fmlfordetectedmuons = par[30] * hlostmuontimedist->GetBinContent( hlostmuontimedist->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ) );
     hdmcorrection->SetBinContent( hlostmuontimedist->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fmlfordetectedmuons );
     f += fmlfordetectedmuons;

     // pileup term
     f += par[19]*f*f;
     //if (dbg) printf("f = %f\n", f);

     return f;
   }

Double_t fcn_R, fcn_dR; // for passing the R-value
Double_t fcn_chi2; // for passing the chisq;
int fcn_ndf; // for passing the ndf;

void chsq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t ) {

  //  minimization using a chisquare method with correlated errors  

  TF1 *fuser   = (TF1*)gFitter->GetUserFunc();
  TH1D *hfit = (TH1D*)gFitter->GetObjectFit();

  Int_t np = fuser->GetNpar();
  fuser->SetParameters( par);
  f = 0;

  double chi2 = 0.;
  Double_t chisq0 = 0., chisq1 = 0., chisq2 = 0.;
  Double_t xi, yi, yerri, fvali, xj, yj, yerrj, fvalj;
  Double_t Mij0 = 1.06, Mij1 = -0.181, Mij2 = 0.0309; // obtained from simulation (smoothing.C) and matrix inversion (matrix.C)

  for (int ib = htmp->FindBin(minTglobal); ib < htmp->FindBin(maxTglobal); ib++) {
    
    chisq0 = 0.;
    chisq1 = 0.; 
    chisq2 = 0.;
    
    xi = htmp->GetBinCenter(ib);
    yi = htmp->GetBinContent(ib);
    yerri = htmp->GetBinError(ib);
    fvali = fuser->EvalPar( &xi, par);
    
    for (int jb = ib-2; jb <= ib+2; jb++) {
      
      xj = htmp->GetBinCenter(jb);
      yj = htmp->GetBinContent(jb);
      yerrj = htmp->GetBinError(jb);
      fvalj = fuser->EvalPar( &xj, par);
      
      if ( ib == jb ) chisq0 += Mij0 * (yi - fvali) * (yi - fvali) / ( yerri * yerri ); // coeff from inversion of covariance matrix
      if ( ib == jb-1 || ib == jb+1 ) chisq1 += Mij1 * (yi - fvali) * (yj - fvalj) / (yerri * yerrj); // coeff from inversion of covariance matrix
      if ( ib == jb-2 || ib == jb+2 ) chisq2 += Mij2 * (yi - fvali) * (yj - fvalj) / (yerri * yerrj); // coeff from inversion of covariance matrix
      
    }

    chi2 += chisq0 + chisq1 + chisq2;
  }

  if (dbg) printf("x %f, y %f +/- %f, f %f, par[0] %f, par[1] %f, chi2 %f\n", xi, yi, yerri, fvali, par[0], par[1], chi2);
  fcn_chi2 = chi2;
  fcn_ndf = fuser->GetNDF();
  fcn_R = fuser->GetParameter(3);
  fcn_dR = fuser->GetParError(3);

  f = chi2;
}

Double_t fprecconv(Double_t *x, Double_t *par);

Double_t fprecconv(Double_t *x, Double_t *par)
  {
     // precession function with sinosoidal amplitude, freqeuncy and phase
     //
     // par[0] - N0
     // par[1] - tau_mu
     // par[2] - A_a
     // par[3] - R -> omega_a
     // par[4] - phi_a
     // par[5] - time-independent background
     // par[6] - A_cbo
     // par[7] - tau_cbo
     // par[8] - omega_cbo
     // par[9] - phi_cbo
     // par[10] - A_gain
     // par[11] - tau_gain
     
     Double_t xx=x[0];

     // injection start time
     //Double_t t0 = 2.49e4;
     Double_t t0 = 0.0;

     double R = par[3];
     //double omega_a = rawBinToNs * 1.e-3 * paramToFreq(R); // conversions from bins to ns and ns to us as Qmethod fit in units of ticks
     double omega_a = 1.00e-3 * getBlinded.paramToFreq(R); // conversions from bins to ns and ns to us as Qmethod fit in units of ticks

     double delta_omega_a = par[29]; // omega_a linear time dependence
     omega_a = omega_a*( 1.0 + delta_omega_a * ( xx-t0 ) );

     // CBO linear time dependence (tracker team parameterization)
     double Acbo = par[23];
     double Bcbo = par[24];
     double omega_cbo, tauAcbo, tauBcbo;
     delta_omega_cbo = par[20]; 

     Double_t f, fig, fpd, fml, fmlfordetectedmuons;
     Double_t omega_cyc, omega_VW, omega_y;
     Double_t A_a, phi_a;
     Double_t fc = 0.0;

     for (int ic = 0; ic <= 1; ic++) {
     
       xx = x[0] + ic*150.; // 150ns bin width for handling convoluted fit

       if (use_new_cbomodel) { 
	 //new model
	 if ( strcmp( "60hr", dataname) == 0 ) tauAcbo = 81.8e3; // convert tracker parameter to ns
	 if ( strcmp( "60hr", dataname) == 0 ) tauBcbo = 7.7e3; // convert tracker parameter to ns
	 if ( strcmp( "9day", dataname) == 0 ) tauAcbo = 72.8e3; // convert tracker parameter to ns
	 if ( strcmp( "9day", dataname) == 0 ) tauBcbo = 8.5e3; // convert tracker parameter to ns
	 if ( strcmp( "endgame", dataname) == 0 ) tauAcbo = 81.8e3; // convert tracker parameter to ns
	 if ( strcmp( "endgame", dataname) == 0 ) tauBcbo = 7.7e3; // convert tracker parameter to ns
	 if ( strcmp( "hk", dataname) == 0 ) tauAcbo = 72.8e3; // convert tracker parameter to ns
	 if ( strcmp( "hk", dataname) == 0 ) tauBcbo = 8.5e3; // convert tracker parameter to ns
	 if ( strcmp( "lk", dataname) == 0 ) tauAcbo = 72.8e3; // convert tracker parameter to ns
	 if ( strcmp( "lk", dataname) == 0 ) tauBcbo = 8.5e3; // convert tracker parameter to ns
	 omega_cbo =  par[8] * ( 1.0 + Acbo*exp( -( xx-t0 ) / tauAcbo ) / ( par[8] * (xx - t0) ) + Bcbo*exp( -( xx-t0 ) / tauBcbo ) / ( par[8] * (xx - t0) ) );
       } else {
	 // old model
	 tauAcbo = 73.5e3; // convert tracker parameter to ns
	 tauBcbo = 16.6e3; // convert tracker parameter to ns
	 omega_cbo =  par[8] * ( 1.0 + delta_omega_cbo * ( xx-t0 ) + Acbo*exp( -( xx-t0 ) / tauAcbo ) + Bcbo*exp( -( xx-t0 ) / tauBcbo ) );
       }
       
       // time dependent horizontal cbo asymmetry term + vertical cbo asymmetry with exponential envelope
       A_a = par[2];
       A_a *= (1.0 + par[12] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 ) + par[13] ) ); // horizontal CBO with freq dependence
       A_a *= (1.0 + par[38] * exp( -( xx - t0 ) / par[27] ) * cos( par[26] * ( xx - t0 ) + par[39] ) ); // vertical CBO f_y 
       //A_a = par[2]; // avoid CBO on assymmetry
       
       // phase of omega_a wiggle with addition of cbo term (2 parameters amplitude, phase)
       phi_a = par[4] + par[14] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 ) + par[15] ); // with freq dependence
       //phi_a = par[4]; // avoid CBO on phase
       
       // 6-parameter wiggle with time-independent background (5-par w/o bkd)
       f = par[0] * exp( -( xx - t0 ) / par[1] ) * ( 1.0 + A_a*cos ( omega_a * ( xx - t0 ) + phi_a ) ) + par[5];
       //if (dbg) printf("f = %f\n", f);
       
       // time dependent cbo leading-order nomalization term - exponential envelope
       f *= (1.0 + par[6] * exp( -( xx - t0 ) / par[7] ) * cos( omega_cbo * ( xx - t0 )   +  par[9] ) );
       //if (dbg) printf("f = %f\n", f);
       
       // time dependent 2*CB0 waist term, freq  2*omega_cbo, tau same as CBO, + amp, phase parameters
       f *= (1.0 + par[36] * exp( -( xx - t0 ) / par[7] ) * cos( 2*omega_cbo * ( xx - t0 ) + par[37] ) );
       //if (dbg) printf("f = %f\n", f);
       
       omega_cyc = 2.*TMath::Pi()/149.3;
       omega_VW = par[17] * ( omega_cyc - 2 * omega_cbo * sqrt( ( 2.*omega_cyc/omega_cbo ) - 1. )); // from joe price
       omega_y = par[26] * omega_cbo * sqrt( ( 2.*omega_cyc/omega_cbo ) - 1. ); // from joe price
       
       // vertical BO modifier, VW = fc-2fy (has much smaller effect from change in field index)
       //f *= ( 1.0 + par[16] * exp( -( xx - t0 ) / par[22] ) * cos( par[17] * ( xx - t0 ) + par[18] ) );
       f *= ( 1.0 + par[16] * exp( -( xx - t0 ) / par[22] ) * cos( omega_VW * ( xx - t0 ) + par[18] ) );
       //if (dbg) printf("f = %f\n", f);
       
       // vertical BO modifier2, BO = fy, same tau as VW (has much larger effect from change in field index)
       //f *= ( 1.0 + par[25] * exp( -( xx - t0 ) / par[22] ) * cos( par[26] * ( xx - t0 ) + par[28] ) ); // tau_VB = tau_VW
       f *= ( 1.0 + par[25] * exp( -( xx - t0 ) / par[27] ) * cos( omega_y * ( xx - t0 ) + par[28] ) ); // tau_VB = tau_VW
       //if (dbg) printf("par[25], par[26], par[27], par[28] %f, %f, %f, %f, f = %f\n", par[25], par[26], par[27], par[28], f);
       
       // injection induced, pedestal ringing 1, or proxy for oscillations
       f *= (1.0 + par[31] * exp( -( xx - t0 ) / par[33] ) * cos( par[32] * ( xx - t0 ) + par[34] ) );
       //if (dbg) printf("par[31], par[32], par[33], par[34] %f, %f, %f, %f, f = %f\n", par[31], par[32], par[33], par[34], f);
       
       // injection induced, pedestal ringing 2, or proxy for oscillations
       f *= (1.0 + par[40] * exp( -( xx - t0 ) / par[42] ) * cos( par[41] * ( xx - t0 ) + par[43] ) );
       //if (dbg) printf("par[39], par[40], par[41], par[42] %f, %f, %f, %f, f = %f\n", par[40], par[41], par[42], par[43], f);
       
       //  use empirical gain multiplier
       if (!use_verticaldrift_correction) fig = ( 1.0 + par[10] * exp( -( xx - t0 ) / par[11] ) ); // empirical infill gain 

       // empirical vertical drift effect
       if (use_verticaldrift_correction) {
	 if (xx <= 43300) 
	   {
	     fig = 2.646-(0.89718-0.897092)/(-3.00230e-02)+3.589e-5*xx-4.135e-10*xx*xx; // from fit to vertical mean energy position
	   } else 
	   {
	     fig = 3.429-2.91e-8*xx; // from fit to vertical mean energy position
	   }
	 fig -= par[11]; // set zero of effect
	 fig = 1.0 + par[10]*fig;
       }

       higcorrection->SetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fig ); // diagnostic     
       f *= fig;
       //if (dbg) printf("f = %f\n", f);
       
       // pileup effect
       f += par[19]*f*f;
       //if (dbg) printf("f = %f\n", f);
       
       // get pedestal drift histogram
       fpd = ( 1.0 - par[35] * hpedestaldrift->GetBinContent( hpedestaldrift->GetXaxis()->FindBin( xx) ) ); 
       f *= fpd; // correct for pedestal drift / ringing
       
       // muon loss multiplier (corrects for effect on positrons of lost muons)- new
       // should i worry about discontinuities in histogram, I dont think so, tested with function as alternative
       //double fml = 1.0 - par[21] * muonlossf->Eval( 1.00e-3*(xx - t0) ); // replace histogram with smooth function, 1.00e-3 ns  to microsec's in Sudeshna muonloss histogram
       fml = 1.0 - par[21] * hmuonloss->GetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ) ); // 1.00e-3 converts ns to microsec's in Sudeshna muonloss histogram
       hmlcorrection->SetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fml ); // diagnostic
       f *= fml;
       
       // muon addition multiplier (corrects for effect of detecting the lost muons)
       fmlfordetectedmuons = par[30] * hlostmuontimedist->GetBinContent( hlostmuontimedist->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ) );
       hdmcorrection->SetBinContent( hlostmuontimedist->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ), fmlfordetectedmuons );
       f += fmlfordetectedmuons;
       
       fc += f;
     }
     
     return fc;
   }

Double_t fvbterms(Double_t *x, Double_t *par);

Double_t fvbterms(Double_t *x, Double_t *par)
  {

     Double_t xx =x[0];

     // injection start time
     Double_t t0 = 0.0;

     // vertical BO modifier, VW = fc-2fy
     Double_t f = 1.0;
     f += ( par[16] * exp( -( xx - t0 ) / par[22] ) * cos( par[17] * ( xx - t0 ) + par[18] ) );
     //if (dbg) printf("f = %f\n", f);

     // vertical BO modifier2, BO = fy, same tau as VW
     f += ( par[25] * exp( -( xx - t0 ) / par[22] ) * cos( par[26] * ( xx - t0 ) + par[28] ) ); // tau_VB = tau_VW

     return f;

  }

Double_t fslowterms(Double_t *x, Double_t *par);


Double_t fslowterms(Double_t *x, Double_t *par)
  {
     // precession function with sinosoidal amplitude, freqeuncy and phase
     //
     // par[0] - N0
     // par[1] - tau_mu
     // par[2] - A_gain
     // par[3] - tau_gain
     // par[4] - muon loss
     
     Double_t xx =x[0];

     // injection start time
     Double_t t0 = 0.0;

     // muon lifetime exponential
     f = par[0] * exp( -( xx - t0 ) / par[1] );

     // Gain multiplier
     f *= ( 1.0 + par[2] * exp( -( xx - t0 ) / par[3] ) );

     // muon loss multiplier (corrects for effect on positrons of lost muons)- 
     //double f *= ( 1.0 - par[21] * muonlossf->Eval( 1.00e-3*(xx - t0) ) ); // replace histogram with smooth function, 1.00e-3 ns  to microsec's in Sudeshna muonloss histogram
     f *= ( 1.0 - par[4] * hmuonloss->GetBinContent( hmuonloss->GetXaxis()->FindBin( 1.00e-3*(xx - t0) ) ) ); // 1.00e-3 converts ns to microsec's in Sudeshna muonloss histogram

     return f;
   }

//  cbof = new TF1("cbof", "[0] * sin( [1]*x + [2] ) + [3] * sin( [4]*x + [5] ) + [6] * sin( [7]*x + [8] )", minT, maxT);

Double_t fcbo(Double_t *x, Double_t *par)
  {
     // sinosoidal amplitude, freqeuncy and phase
     //
     // par[0] - A_cbo1
     // par[1] - omega_cbo1
     // par[2] - phi_cbo1
     // etc

     Double_t xx =x[0];
     
     Double_t f = par[0] * sin ( par[1]*xx + par[2] );
     f += par[3] * sin ( par[4]*xx + par[5] );
     f += par[6] * sin ( par[7]*xx + par[8] );
     f += par[9] * sin ( par[10]*xx + par[11] );
     f += par[12] * sin ( par[13]*xx + par[14] );
     
     return f;
   }

Double_t fvo(Double_t *x, Double_t *par);

Double_t fvo(Double_t *x, Double_t *par)
  {
     //  vvw,vo fit function with exponetial envelope
     //
     // par[0] - A_vo1
     // par[1] - tau_vo1
     // par[1] - omega_vo1
     // par[2] - phi_vo1
     // and same vor vo2

     Double_t xx =x[0];
     
     Double_t f = par[0] *exp(-(xx-27.4e3)/par[1]) * sin ( par[2]*xx + par[3] );
     f += par[4] *exp(-(xx-27.4e3)/par[5]) * sin ( par[6]*xx + par[7] );
     f *= 1.0+par[10]*sin(par[11]*(xx-27.4e3)+par[12]); // omega multiplier
     f += par[8]*exp(-(xx-27.4e3)/par[9]); // bkd 

     return f;
   }

Double_t ffr(Double_t *x, Double_t *par)
  {
     // fast rotation correction function
     //
     // par[0] - A_fr
     // par[1] - tau_fr
     // par[2] - omega_fr
     // par[3] - phi_fr

     Double_t xx =x[0];

     Double_t f = par[0] * exp( -(xx-5.e4)/par[1] ) * sin( par[2]*xx + par[3] );
     return f;
   }

Double_t fring(Double_t *x, Double_t *par);
Double_t fring(Double_t *x, Double_t *par)
{
  // fast rotation correction function
  //
  // par[0] - A_ring
  // par[1] - tau_ring
  // par[2] - omega_ring
  // par[3] - phi_ring
  
  // injection start time
  Double_t t0 = 2.49e4;
  
  Double_t xx =x[0];
  
  Double_t f = par[0] * exp( -( xx - t0 ) / par[1] ) * sin( par[2] * ( xx - t0 ) + par[3] );
  
  return f;
}


// user functions
int fitCBO(int tmn, int tmx);
int goSumCBO();
int sumCBO(int ic);
int caloSumSlice(char* cRun);
int sumXtal(int iRun, int iCal, char* cHist);
int sumXtal2D(int iRun, int iCal, char* cHist);
int sumCalo(int iRun, char* cHist);
int sumCalo2D(int iRun, char* cHist);
int sumCaloXtal(int iRun, char* cHist);
int goXtal(int iRun, char* );
int goXtal2D(int iRun, char* );
int goSumHist(int iRun);
int goSumXtalHist(int iRun); // adds calos for each xtal, i.e. xtal 1, calo 1 +  xtal 1, calo 2, ...
int goSumScat(int iRun);
int goHist(int iRun);
int goScat(int iRun);
int goAll();

int BuildFrInt();

int goDrawSync(int iRun);
int drawSync(int iRun, int iCal);

int goDrawFlash(int iRun);
int drawFlash(int iRun, int iCal);

int goDrawThrSplash(int iRun);
int drawThrSplash(int iRun, int iCal);

int goDrawThrStore(int iRun);
int drawThrStore(int iRun, int iCal);

int goDrawPed1(int iRun);
int drawPed1(int iRun, int iCal);

int goDrawPed2(int iRun);
int drawPed2(int iRun, int iCal);

int goDrawPedvSample(int iRun);
int drawPedvSample(int iRun, int iCal);

int goFitCalo( int iS, int iF);
int fitCalo( char* cFit, int iC, int iS, int iF);
//int viewCalo( char* cRun, double startTime, double endTime);
int viewCalo( char* cFit, char* cRun, char* cHis, int swrb, double startTime, double endTime);  // in 800 MSPS units

int energyPlots();

double GPar_QM[Nfit], dGPar_QM[Nfit], NF_QM[Nfit], dNF_QM[Nfit], tauF_QM[Nfit], err_tauF_QM[Nfit], dtauF_QM[Nfit], AF_QM[Nfit], dAF_QM[Nfit], omega_aF_QM[Nfit], err_omega_aF_QM[Nfit], domega_aF_QM[Nfit], phaseF_QM[Nfit], err_phaseF_QM[Nfit], dphaseF_QM[Nfit], BF_QM[Nfit], dBF_QM[Nfit], chisq_QM[Nfit], dchisq_QM[Nfit], ndf_QM[Nfit];

double ACBOF_QM[Nfit], dACBOF_QM[Nfit], tauCBOF_QM[Nfit], dtauCBOF_QM[Nfit], freqCBOF_QM[Nfit], dfreqCBOF_QM[Nfit], deltafreqCBOF_QM[Nfit], ddeltafreqCBOF_QM[Nfit], phaseCBOF_QM[Nfit], dphaseCBOF_QM[Nfit];
double ACBO2F_QM[Nfit], dACBO2F_QM[Nfit], phaseCBO2F_QM[Nfit], dphaseCBO2F_QM[Nfit];
double ACBO3F_QM[Nfit], dACBO3F_QM[Nfit], phaseCBO3F_QM[Nfit], dphaseCBO3F_QM[Nfit];
double A2CBOF_QM[Nfit], dA2CBOF_QM[Nfit], phase2CBOF_QM[Nfit], dphase2CBOF_QM[Nfit];
double AGAINF_QM[Nfit], dAGAINF_QM[Nfit];
double tauGAINF_QM[Nfit], dtauGAINF_QM[Nfit];
double muonlossF_QM[Nfit], dmuonlossF_QM[Nfit];
double peddriftF_QM[Nfit], dpeddriftF_QM[Nfit];
double RPILEUPF_QM[Nfit], dRPILEUPF_QM[Nfit];
double AVWF_QM[Nfit], dAVWF_QM[Nfit], omegaVWF_QM[Nfit], domegaVWF_QM[Nfit], tauVWF_QM[Nfit], dtauVWF_QM[Nfit], phaseVWF_QM[Nfit], dphaseVWF_QM[Nfit];
double AVW2F_QM[Nfit], dAVW2F_QM[Nfit], omegaVW2F_QM[Nfit], domegaVW2F_QM[Nfit], tauVW2F_QM[Nfit], dtauVW2F_QM[Nfit], phaseVW2F_QM[Nfit], dphaseVW2F_QM[Nfit];
double AFrelax_QM[Nfit], dAFrelax_QM[Nfit];

double A_cbo_save, tau_cbo_save, omega_cbo_save, phi_cbo_save, delta_omega_cbo_save, A2_cbo_save, phi2_cbo_save, A_vbo_save, phi_vbo_save;

double noise1, noise2, noise4, noise8;
double minTCBO = 60000, maxTCBO = 190000; // used for fitting the CBO term

TCanvas *c1 = new TCanvas("c1","c1", 1000, 700);
TGraphErrors *gtmp, *gtmp1, *gtmp2, *gtmp3;
TGraph  *gybandm, *gybandp;

int goCaloSum(char* iRun);
int caloSum(char* cRun, char* cHist);
int caloSumHot(char* cRun, char* cHist);
int caloSumCold(char* cRun, char* cHist);
int caloSum2D(char* cRun, char* cHist);
int caloSumRb(char* cRun, char* cHist, int rbmax);
int caloSumRbHot(char* cRun, char* cHist, int rbmax);
int caloSumRbCold(char* cRun, char* cHist, int rbmax);
int caloSumRbNear(char* cRun, char* cHist, int rbmax);
int caloSumRbFar(char* cRun, char* cHist, int rbmax);
int caloSumRbPairs(char* cRun, char* cHist, int rbmax);

int starttimewindowscan(int jt, int jg, int jrb);
int thresholdscan(int jt, int jg, int jrb);
int kawallBands(char* scanname, char* histname, int jt, int jg, int jrb);
int bintest();

int BuildFrInt(){
  //
  // function does "by-hand" integral of fast rotation function over bin bin width.
  // Had problem with aut bin-integral thru "I" option in for function and this method
  // allow changing of tolerence on bin integral
  //
  hFRint = (TH1D*)htmp->Clone();
  hFRint->Reset();
  hFRint->SetName("hFRint");
  hFRint->SetTitle("hFRint");
  
  for (int ib = 1; ib <= htmp->GetNbinsX(); ib++){
    if ( htmp->GetBinLowEdge( ib) + htmp->GetBinWidth( ib) < 30000. || htmp->GetBinLowEdge( ib) > 70000. ) continue;
    hFRint->SetBinContent( ib, ffastrotation->Integral( htmp->GetBinLowEdge( ib), htmp->GetBinLowEdge( ib) + htmp->GetBinWidth( ib), 1.e-3 )/htmp->GetBinWidth( ib) );
  }

  return 0;
}


TH1D  *h1test;
TH1D  *h2test;
TH1D  *h3test;
int bintest(){

  h1test = new TH1D("h1test","h1test",190020,20001.0,210021.);
  h2test = new TH1D("h2test","h2test",3168,20001.0,210021.);
  h3test = new TH1D("h3test","h3test",3167,20001.0,210021.);
   
  for (int i = 20001; i <= 210021; i++){
    h1test->Fill(i);
    h2test->Fill(i);
    h3test->Fill(i);
  }
  return 0; 
}

/* how to run a few scans and avoid too many files error

root -l
.L sumCalo.C
thresholdscan(10,0,1)
.qqqq

.L sumCalo.C
starttimewindowscan(12,11,1)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(15,11,1)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(18,11,1)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(24,11,1)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(12,11,2)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(15,11,2)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(18,11,2)
.qqqq
root -l
.L sumCalo.C
starttimewindowscan(24,11,2)
.qqqq

*/

/*
iwind = 12
iswrb = 1
igap = 0
ithreshold = 32
*/

int starttimewindowscan(int jt, int jg, int jrb){

  ithreshold = jt;
  iswrb = jrb;
  igap=jg;

  iwind=2;
  fitCalo("starttimescan",0,50000,200000);
  iwind=4;
  fitCalo("starttimescan",0,50000,200000);
  iwind=8;
  fitCalo("starttimescan",0,50000,200000);
  iwind=12;    
  fitCalo("starttimescan",0,50000,200000);
  return 0;
  }

int thresholdscan(int jt, int jg, int jrb){

  ithreshold = jt;
  iswrb = jrb;
  igap=jg;

  iwind=2;
  fitCalo("thresholdscan",0,50000,200000);
  iwind=4;
  fitCalo("thresholdscan",0,50000,200000);
  iwind=8;
  fitCalo("thresholdscan",0,50000,200000);
  iwind=12;    
  fitCalo("thresholdscan",0,50000,200000);
  return 0;
  }

int windowscan(int jt, int jg, int jrb){

  ithreshold = jt;
  iswrb = jrb;
  igap=jg;

  ithreshold=12;
  fitCalo("windowscan",0,50000,200000);
  ithreshold=18;
  fitCalo("windowscan",0,50000,200000);
  ithreshold=24;
  fitCalo("windowscan",0,50000,200000);
  ithreshold=32;    
  fitCalo("windowscan",0,50000,200000);
  return 0;
  }

/*

kawallBands("thresholdscan","gOmega_afit",32,0,1);
c1->SaveAs("w12-t32-g00-swrb1-prd1-gOmega_afit-thresholdscan.png")


kawallBands("gOmega_afit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gOmega_afit-starttimescan.png")
kawallBands("gPhifit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gPhifit-starttimescan.png")
kawallBands("gAfit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gAfit-starttimescan.png")
kawallBands("gTaufit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gTaufit-starttimescan.png")
kawallBands("gmuonlossfit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gmuonlossfit-starttimescan.png")
kawallBands("gChifit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gChifit-starttimescan.png")

kawallBands("gOmega_afit",12,11,1);
c1->SaveAs("w12-t12-g11-swrb1-prd1-gOmega_afit-starttimescan.png")

kawallBands("gOmega_afit",18,11,1);
c1->SaveAs("w12-t18-g11-swrb1-prd1-gOmega_afit-starttimescan.png")

kawallBands("gOmega_afit",24,11,1);
c1->SaveAs("w12-t24-g11-swrb1-prd1-gOmega_afit-starttimescan.png")

kawallBands("gOmega_afit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gOmega_afit-starttimescan.png")

kawallBands("gAfit",12,11,1);
c1->SaveAs("w12-t12-g11-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",18,11,1);
c1->SaveAs("w12-t18-g11-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",24,11,1);
c1->SaveAs("w12-t24-g11-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gPhifit",12,11,1);
c1->SaveAs("w12-t12-g11-swrb1-prd1-gPhifit-starttimescan.png")

kawallBands("gPhifit",18,11,1);
c1->SaveAs("w12-t18-g11-swrb1-prd1-gPhifit-starttimescan.png")

kawallBands("gPhifit",24,11,1);
c1->SaveAs("w12-t24-g11-swrb1-prd1-gPhifit-starttimescan.png")

kawallBands("gPhifit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gPhifit-starttimescan.png")


kawallBands("gTaufit",12,11,1);
c1->SaveAs("w12-t12-g11-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",18,11,1);
c1->SaveAs("w12-t18-g11-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",24,11,1);
c1->SaveAs("w12-t24-g11-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gTaufit-starttimescan.png")


kawallBands("gmuonlossfit",12,11,1);
c1->SaveAs("w12-t12-g11-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",18,11,1);
c1->SaveAs("w12-t18-g11-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",24,11,1);
c1->SaveAs("w12-t24-g11-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gChifit",15,11,1);
c1->SaveAs("w12-t15-g11-swrb1-prd1-gChifit-starttimescan.png")


kawallBands("gA_afit",12,0,1);
c1->SaveAs("w12-t12-g00-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",24,0,1);
c1->SaveAs("w12-t24-g00-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gAfit",32,0,1);
c1->SaveAs("w12-t32-g00-swrb1-prd1-gAfit-starttimescan.png")

kawallBands("gTaufit",12,0,1);
c1->SaveAs("w12-t12-g00-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",24,0,1);
c1->SaveAs("w12-t24-g00-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gTaufit",32,0,1);
c1->SaveAs("w12-t32-g00-swrb1-prd1-gTaufit-starttimescan.png")

kawallBands("gmuonlossfit",12,0,1);
c1->SaveAs("w12-t12-g00-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",18,0,1);
c1->SaveAs("w12-t18-g00-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",24,0,1);
c1->SaveAs("w12-t24-g00-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gmuonlossfit",32,0,1);
c1->SaveAs("w12-t32-g00-swrb1-prd1-gmuonlossfit-starttimescan.png")

kawallBands("gChifit",12,0,1);
c1->SaveAs("w12-t12-g00-swrb1-prd1-gChifit-starttimescan.png")

 */




int nCaloFrStudy = 2;
char fnamefrstudy[128];
bool doFRStudy = 1;
bool frStoreHistograms =1;

int frshiftwrite();
int frshiftwrite(double tb = 30000, double te = 215000){

  TH1D *hs[24], *hns[24], *ds[24], *ss[24];

  for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    doShiftBins = 0;
    fitCalo( "singlecalo", iCalo, tb, te);
    hns[iCalo-1] = (TH1D*)hres->Clone();

    doShiftBins = 1;
    fitCalo( "singlecalo", iCalo, tb, te);
    hs[iCalo-1]=(TH1D*)hres->Clone();
    ds[iCalo-1]=(TH1D*)hres->Clone();
    ss[iCalo-1]=(TH1D*)hres->Clone();
    ds[iCalo-1]->Add(hns[iCalo-1],-1.0);
    ss[iCalo-1]->Add(hns[iCalo-1],+1.0);

    for (int ib = 1; ib <= hs[iCalo-1]->GetNbinsX(); ib++){
      hs[iCalo-1]->SetBinError(ib, 0.0);
      hns[iCalo-1]->SetBinError(ib, 0.0);
      ds[iCalo-1]->SetBinError(ib, 0.0);
      ss[iCalo-1]->SetBinError(ib, 0.0);
    }

    sprintf(hname,"hs%i",iCalo);
    hs[iCalo-1]->SetName(hname);
    hs[iCalo-1]->SetTitle(hname);
    sprintf(hname,"hns%i",iCalo);
    hns[iCalo-1]->SetName(hname);
    hns[iCalo-1]->SetTitle(hname);
    sprintf(hname,"ds%i",iCalo);
    ds[iCalo-1]->SetName(hname);
    ds[iCalo-1]->SetTitle(hname);
    sprintf(hname,"ss%i",iCalo);
    ss[iCalo-1]->SetName(hname);
    ss[iCalo-1]->SetTitle(hname);

    hs[iCalo-1]->SetLineColor(kBlue);
    hns[iCalo-1]->SetLineColor(kRed);
    ds[iCalo-1]->SetLineColor(kBlack);
    ss[iCalo-1]->SetLineColor(kMagenta);
    hs[iCalo-1]->SetLineWidth(2);
    hns[iCalo-1]->SetLineWidth(2);
    ds[iCalo-1]->SetLineWidth(2);
    ss[iCalo-1]->SetLineWidth(2);

    //hs[iCalo-1]->Rebin(8);
    //hns[iCalo-1]->Rebin(8);
    //ds[iCalo-1]->Rebin(8);
  }

  sprintf( fnamefrstudy , "frstudy.root");
  if (frStoreHistograms){

    TFile *ffr;
    ffr = TFile::Open( fnamefrstudy, "UPDATE");
    for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {
      ds[iCalo-1]->Write();
      ss[iCalo-1]->Write();
      hs[iCalo-1]->Write();
      hns[iCalo-1]->Write();
    }
    ffr->Close();
  }

  return 0;
}

int frshiftplot();
int frshiftplot(){

  TH1D hs[24], hns[24], ds[24], ss[24];
  TFile *ffr;
  sprintf( fnamefrstudy , "frstudy.root");
  ffr = TFile::Open( fnamefrstudy, "UPDATE");

  c1->Clear();
  c1->Divide(4,6);

  for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    c1->cd(iCalo);
    
    sprintf( hname, "ds%i", iCalo);
    ffr->GetObject( hname, htmp2);
    htmp2->Rebin(8);
    htmp2->Draw("hist");
  }
  c1->SaveAs("ds.png")
;
 for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    c1->cd(iCalo);
    sprintf( hname, "ss%i", iCalo);
    ffr->GetObject( hname, htmp2);
    htmp2->Rebin(8);
    htmp2->Draw("hist");
  }
  c1->SaveAs("ss.png");

  ffr->Close();

  return 0;
}

int frcorrplot(){

  TH1D hs[24], hns[24], ds[24], ss[24];
  TFile *ffr1, *ffr2;

  sprintf( fnamefrstudy , "hcyc-9day-13par-pw-gold-gainbugfix-thr30-GainON-LossOFF-RingOFF-convoluted-noBinShift.root");
  ffr1 = TFile::Open( fnamefrstudy, "UPDATE");
  sprintf( fnamefrstudy , "hcyc-9day-13par-pw-gold-gainbugfix-thr30-GainON-LossOFF-RingOFF-convoluted-BinShift.root");
  ffr2 = TFile::Open( fnamefrstudy, "UPDATE");
  sprintf( fnamefrstudy , "frstudy.root");
  ffr3 = TFile::Open( fnamefrstudy, "UPDATE");

  c1->Clear();
  c1->Divide(4,6);

  for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    c1->cd(iCalo);
    
    sprintf( hname, "hfull%i", iCalo);
    ffr1->GetObject( hname, htmp2);
    htmp2->Rebin(8);
    htmp2->SetLineColor(kBlue);
    htmp2->SetMaximum(0.02);
    htmp2->SetMinimum(-0.02);
    htmp2->Draw("hist");
  }

 for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    c1->cd(iCalo);
    sprintf( hname, "hfull%i", iCalo);
    ffr2->GetObject( hname, htmp2);
    htmp2->Rebin(8);
    htmp2->SetLineColor(kRed);
    htmp2->Draw("histsame");
  }

  c1->SaveAs("fits-with-without-binshift.png");

  c1->Clear();
  c1->Divide(4,6);

 for (iCalo = 1; iCalo <= nCaloFrStudy; iCalo++) {

    c1->cd(iCalo);
    sprintf( hname, "ds%i", iCalo);
    ffr3->GetObject( hname, htmp3);
    sprintf( hname, "hfull%i", iCalo);
    ffr2->GetObject( hname, htmp2);
    sprintf( hname, "hfull%i", iCalo);
    ffr1->GetObject( hname, htmp1);
    htmp2->Add(htmp1,-1.0);
    htmp2->SetLineColor(kMagenta);
    htmp2->Draw("hist");
    htmp3->Rebin(8);
    htmp3->Draw("histsame");
  }

  c1->SaveAs("fits-diff-with-without-binshift.png");


  ffr1->Close();
  ffr2->Close();
  ffr3->Close();

  return 0;
}

int frcorrectionplot(){
  
  char fname1[128], fname2[128];
  char histname[64];
  
  double xval1[50]={0}, yval1[50]={0};
  double xval2[50]={0}, yval2[50]={0}, dy[50]={0}, err[50]={0};
   
  sprintf( fname1, "60hr-FR.root"); 
  sprintf( fname2, "60hr-noFR.root"); 
  sprintf( fname1, "9day-FR.root"); 
  sprintf( fname2, "9day-noFR.root"); 
  sprintf( fname1, "60hr-FR-55000.root"); 
  sprintf( fname2, "60hr-noFR-55000.root"); 
  sprintf( fname1, "9day-FR-55000.root"); 
  sprintf( fname2, "9day-noFR-55000.root"); 
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-shiftbins-noAmplShift-noTimeShift-ts30000.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-noTimeShift-noAmplShift-ts30000.root");
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-amplFRC-BinShift.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-amplFRC-noBinShift.root");
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-oppamplFRC-BinShift.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-oppamplFRC-noBinShift.root");
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-timeamplFRC-BinShift.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-timeamplFRC-noBinShift.root");
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-timeoppamplFRC-BinShift.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-timeoppamplFRC-noBinShift.root");
  sprintf( fname1, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-noFRC-Blend1.root");
  sprintf( fname2, "graphs-w4-t18-g11-swrb1-prd1e0-caloscan-vbo2-9day-all-noFRC-Blend2.root");

  TFile *file1 = TFile::Open(fname1);
  printf("got file0 %s\n", fname1); 
  TFile *file2 = TFile::Open(fname2);
  printf("got file1 %s\n", fname2); 
  
  sprintf( histname, "gOmega_afit"); 
  file1->GetObject(histname, gtmp1);
  printf("got graph1\n");
  file2->GetObject(histname, gtmp2);
  printf("got graph2\n");

  for (int ip = 1; ip <= gtmp1->GetN(); ip++ ) {
    gtmp1->GetPoint(ip-1, xval1[ip-1], yval1[ip-1]);
    gtmp2->GetPoint(ip-1, xval2[ip-1], yval2[ip-1]);
    dy[ip-1] = yval1[ip-1] - yval2[ip-1];
    printf("x, y1, y2 %f %f %f diff %f\n", xval1[ip-1], yval1[ip-1], yval2[ip-1], dy[ip-1]); 
  }

  gtmp3 = new TGraphErrors(24, xval1, dy, err, err);
  gtmp3->SetMarkerColor(kRed);
  gtmp3->SetMarkerStyle(20);
  gtmp3->SetLineWidth(2);
  gtmp3->SetTitle("frcorrection");
  gtmp3->SetName("frcorrection");
  gtmp3->GetXaxis()->SetTitle("calorimeter index");
  gtmp3->GetYaxis()->SetTitle("FR correction to R[ppm]");

  c1->Clear();
  c1->Divide(1,3);

  c1->cd(1);
  gtmp3->Draw("ap");
  c1->cd(2);
  gtmp1->Draw("ap");
  c1->cd(3);
  gtmp2->Draw("ap");

  file1->Close();
  file2->Close();

  return 0;
}

int kawallBands(char* scanname, char* histname, int jt, int jg, int jrb){

  /*
  KEY: TGraphErrors	gNfit;1	Q-Method normalization
  KEY: TGraphErrors	gTaufit;1	Q-Method tau
  KEY: TGraphErrors	gAfit;1	Q-Method amplitude
  KEY: TGraphErrors	gOmega_afit;1	Q-Method blind R[ppm]
  KEY: TGraphErrors	gPhifit;1	Q-Method phase
  KEY: TGraphErrors	gACBOfit;1	Q-Method CBO ampl
  KEY: TGraphErrors	gfreqCBOfit;1	Q-Method CBO freq
  KEY: TGraphErrors	gphaseCBOfit;1	Q-Method CBO phase
  KEY: TGraphErrors	gACBO2fit;1	Q-Method CBO2 ampl (for A_a)
  KEY: TGraphErrors	gAGAINfit;1	Q-Method GAIN ampl
  KEY: TGraphErrors	gmuonlossfit;1	Q-Method muon loss ampl.
  KEY: TGraphErrors	gChifit;1	Q-Method chi-squared/pdf
  KEY: TGraphErrors	gAVWfit;1	Q-Method VW ampl
  KEY: TGraphErrors	gphaseVWfit;1	Q-Method VW ampl
  KEY: TGraphErrors	gAVWfit2;1	Q-Method VW2 ampl
  KEY: TGraphErrors	gphaseVWfit2;1	Q-Method VW2 ampl
  */

  double xval[50]={0}, yval[50]={0};
  double xerr[50]={0}, yerr[50]={0};
  double ybandm[50], ybandp[50]; 
  char fname[64];
  double xp, yp;
  ithreshold=jt;
  iswrb = jrb;
  igap=jg;

  iwind=4;
  ithreshold=15;
  igap=0;
  sprintf( fname, "graphs-w%i-t%i-g%i-swrb%i-prd%s-%s.root", iwind, ithreshold, igap, iswrb, aname, scanname); 
  printf("file %s\n", fname);
  TFile *file0 = TFile::Open(fname);
  printf("got file %s\n", fname); 
  file0->GetObject(histname, gtmp);
  printf("got graph\n");
   
  for (int ip = 1; ip <= gtmp->GetN(); ip++ ) {
    gtmp->GetPoint(ip-1, xval[ip-1], yval[ip-1]);
    yerr[ip-1] = gtmp->GetErrorY(ip-1);
    printf("x, y, ey %f %f %f \n", xval[ip-1], yval[ip-1], yerr[ip-1]); 
    if ( ( ip > 1 ) && ( yerr[ip-1] <  yerr[ip-2] ) )  yerr[ip-1] =  yerr[ip-2];
  }

#if 1 //for start time scan 

  ybandm[0] =  yval[0];
  ybandp[0] =  yval[0];
  for (int ip = 2; ip <= gtmp->GetN(); ip++ ) {
    ybandm[ip-1] =  yval[0] - sqrt( yerr[ip-1]*yerr[ip-1] - yerr[0]*yerr[0] );
    ybandp[ip-1] =  yval[0] + sqrt( yerr[ip-1]*yerr[ip-1] - yerr[0]*yerr[0] );
  }
  
  gybandm = new TGraph( 50, xval, ybandm); 
  gybandp = new TGraph( 50, xval, ybandp); 
  gybandm->SetLineWidth(2);
  gybandp->SetLineWidth(2);

#endif

#if 0 //for threshold scan - kawall bands from simulation scaled for 2.2 ppm stat errors
   
  ybandm[0] =  yval[0];
  ybandp[0] =  yval[0];
  ybandm[1] =  yval[0]-2.2*5.8/65.7;
  ybandp[1] =  yval[0]+2.2*5.8/65.7;
  ybandm[2] =  yval[0]-2.2*12.8/65.7;
  ybandp[2] =  yval[0]+2.2*12.8/65.7;
  ybandm[3] =  yval[0]-2.2*17.9/65.7;
  ybandp[3] =  yval[0]+2.2*17.9/65.7;
  ybandm[4] =  yval[0]-2.2*27.2/65.7;
  ybandp[4] =  yval[0]+2.2*27.2/65.7;

  gybandm = new TGraph( 5, xval, ybandm); 
  gybandp = new TGraph( 5, xval, ybandp); 
  gybandm->SetLineWidth(2);
  gybandp->SetLineWidth(2);

#endif

  gtmp->SetMarkerColor(kRed);
  gtmp->SetMarkerStyle(20);
  gtmp->SetLineColor(kRed);
  gtmp->SetLineWidth(2);
  //gtmp->SetMaximum(0.25);
  //gtmp->SetMinimum(0.21);
  gtmp->Draw("AP");
  gybandm->Draw("L");
  gybandp->Draw("L");
  
  return 0;

  //TFile *file1 = TFile::Open("graphs-w4-t18-g11-swrb1-prd1-starttimescan.root");
  //iwind=4;
  ithreshold=18;
  sprintf( fname, "graphs-w%i-t%i-g%i-swrb%i-prd%s-%s.root", iwind, ithreshold, igap, iswrb, aname, scanname); 
  TFile *file1 = TFile::Open(fname);
  file1->GetObject(histname, gtmp);
  gtmp->SetMarkerColor(kBlack);
  gtmp->SetMarkerStyle(21);
  gtmp->SetLineColor(kBlack);
  gtmp->SetLineWidth(2);
  gtmp->Draw("P");
  //TFile *file2 = TFile::Open("graphs-w8-t18-g11-swrb1-prd1-starttimescan.root");
  //iwind=8;
  ithreshold=24;
  sprintf( fname, "graphs-w%i-t%i-g%i-swrb%i-prd%s-%s.root", iwind, ithreshold, igap, iswrb, aname, scanname); 
  TFile *file2 = TFile::Open(fname);
  file2->GetObject(histname, gtmp);
  gtmp->SetMarkerColor(kMagenta);
  gtmp->SetMarkerStyle(22);
  gtmp->SetLineColor(kMagenta);
  gtmp->Draw("P");
  //TFile *file3 = TFile::Open("graphs-w12-t18-g11-swrb1-prd1-starttimescan.root");
  //iwind=12;
  ithreshold=32;
  sprintf( fname, "graphs-w%i-t%i-g%i-swrb%i-prd%s-%s.root", iwind, ithreshold, igap, iswrb, aname, scanname); 
  TFile *file3 = TFile::Open(fname);
  file3->GetObject(histname, gtmp);
  gtmp->SetMarkerColor(kOrange+4);
  gtmp->SetMarkerStyle(23);
  gtmp->SetLineColor(kOrange+2);
  gtmp->Draw("P");

  sprintf( fname, "w%i-t%i-g%02i-swrb%i-prd%s-%s-threshold.png", iwind, ithreshold, igap, iswrb, aname, histname); 
  c1->SaveAs(fname);

  return 0;
}


int singlefills( int calo, int offset){

  char fname[64], hname[64];

  int nfills = 8;
  c1->Clear();
  c1->Divide(4,2);

  for (int ifill = 1; ifill <= nfills; ifill++){
    
    sprintf( fname,"noxtaldata/7762_f%i.root", ifill+offset);
    sprintf( fname,"noxtaldata/7844_f%i.root", ifill+offset);
    sprintf( fname,"7847_f%i.root", ifill+offset);
    TFile *file = TFile::Open( fname);
    printf("got file\n");

    TDirectoryFile *dir;
    file->GetObject("CQBankAnalyzerTim",dir);
    printf("got directory\n");
    
    sprintf( hname,"qHist1D_%i_%i", calo, calo); // old name
    sprintf( hname,"qHist1D_%i_%i", calo); // new name
    dir->GetObject( hname, htmp1);
    printf("got histogram\n");

    c1->cd(ifill);
    htmp1->GetXaxis()->SetRangeUser(24000,30000);
    htmp1->SetMaximum(+1000.);
    htmp1->SetMinimum(-1000.);
    htmp1->Draw();

    //file->Delete();
    //dir->Delete();
    
  }
  sprintf( fname,"7762_group%i_calo%i.png", offset, calo);
  c1->SaveAs(fname);

  return 1;
}

int fitCBO(int ic, int tmn, int tmx){

  // fit cbo histograms from residuals
  //
  char fname[64], hname[64];

  //TFile *file1 = TFile::Open("cbo-5par-residuals.root");
  //TFile *file1 = TFile::Open("cbo-5par-residuals-highrate.root");
  TFile *file1 = TFile::Open("cbo-6par-pileupterm-residuals-highrate.root");
  //TFile *file1 = TFile::Open("cbo-6par-pileupterm-residuals-lowrate.root");
  //TFile *file1 = TFile::Open("cbo5par.root");
  //sprintf( hname, "cbo_sumi");
  sprintf( hname, "cbo_calo_%02i", ic-1);
  file1->GetObject(hname, htmp1);

  //TFile *file1 = TFile::Open("hotcold.root");
  //file1->GetObject("hot", htmp1);

  // omega_cbo = 3.22904e-03 from fit in raw bin units
  // omega_cbo-omega_a = 3.22904e-03-1.79982e-03 = 0.00142922 from fit in raw bin units
  // omega_cbo+omega_a = 3.22904e-03+1.79982e-03 = 0.00502886 from fit in raw bin units

  double omega_cbo1 = 3.22904e-03; // omega_cbo
  double omega_cbo2 = 0.00142922; // omega_cbo-omega_a
  double omega_cbo3 = 0.00502886; // omega_cbo+omega_a
  double omega_cbo4 = 121./30325.; // ??? 30325 is FFT bins to omega in raw 1.25ns bins conversion factor
  double omega_cbo5 = 73./30325; // ??? 30325 is FFT bins to omega in raw 1.25ns bins conversion factor

  double A_cbo1 = 0.002;
  double A_cbo2 = 0.05;
  double A_cbo3 = 0.05;
  double A_cbo4 = 0.05;
  double A_cbo5 = 0.05;

  double phi_cbo1 = 0.0;
  double phi_cbo2 = 0.0;
  double phi_cbo3 = 0.0;
  double phi_cbo4 = 0.0;
  double phi_cbo5 = 0.0;

  double minT = 58000.;
  double maxT = 105000.;
  minT= tmn;
  maxT= tmx;

  A_cbo1 = 0.4;
  //cbof = new TF1("cbof", "[0] * sin( [1]*x + [2] ) + [3] * sin( [4]*x + [5] ) + [6] * sin( [7]*x + [8] )", minT, maxT);

  cbof = new TF1("cbof", fcbo, 0.0, 560000.0, 15);
  cbof->SetNpx(10000);
  cbof->SetParameters( A_cbo1, omega_cbo1, phi_cbo1, A_cbo2, omega_cbo2, phi_cbo2, A_cbo3, omega_cbo3, phi_cbo3);
  cbof->SetParNames( "A_cbo1", "omega_cbo1", "phi_cbo1", "A_cbo2", "omega_cbo2", "phi_cbo2", "A_cbo3", "omega_cbo3", "phi_cbo3");
  cbof->SetParameter(9,A_cbo4);
  cbof->SetParameter(10,omega_cbo4);
  cbof->SetParameter(11,phi_cbo4);
  cbof->SetParameter(12,A_cbo5);
  cbof->SetParameter(13,omega_cbo5);
  cbof->SetParameter(14,phi_cbo5);

  cbof->SetParLimits(2,-3.14,3.14);
  cbof->SetParLimits(5,-3.14,3.14);
  cbof->SetParLimits(8,-3.14,3.14);
  cbof->SetParLimits(11,-3.14,3.14);
  cbof->SetParLimits(14,-3.14,3.14);

  //cbof->FixParameter(0,A_cbo1);
  //cbof->FixParameter(3,0.0);
  //cbof->FixParameter(6,0.0);
  cbof->FixParameter(9,0.0);
  cbof->FixParameter(12,0.0);

  cbof->FixParameter(1,omega_cbo1);
  cbof->FixParameter(4,omega_cbo2);
  cbof->FixParameter(7,omega_cbo3);
  cbof->FixParameter(10,omega_cbo3);
  cbof->FixParameter(13,omega_cbo3);

  //cbof->FixParameter(1,phi_cbo1);
  //cbof->FixParameter(5,phi_cbo2);
  //cbof->FixParameter(8,phi_cbo3);
  cbof->FixParameter(11,phi_cbo3);
  cbof->FixParameter(14,phi_cbo3);

  cbof->SetLineWidth(2);
  cbof->SetLineColor(kRed);

  htmp1->SetLineWidth(2);
  htmp1->GetXaxis()->SetRangeUser(50000, 200000);

  cbof->SetParameter(1,2.6e-3);
  htmp1->Fit( cbof, "", "R", minT, maxT);

  // TH1 *hm =0;TVirtualFFT::SetTransform(0);hm = hres->FFT(hm, "MAG"); hm->SetTitle("discrete fourier transform");hm->Draw();
  //Compute the transform and look at the magnitude of the output
  TVirtualFFT::SetTransform(0);
  TH1 *hm = 0;
  hm = htmp1->FFT(hm, "MAG");
  hm->SetTitle("discrete fourier transform");
  hm->SetTitle("FFT of calo");
  hm->SetName("caloFFT");

  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  htmp1->GetXaxis()->SetRangeUser(50000,200000);
  htmp1->Draw();
  htmp1->Draw("histsame");
  //cbof->Draw("same");
  c1->cd(2);
  hm->GetXaxis()->SetRangeUser(0,800);
  hm->Draw();

  sprintf( fname, "cbo_residuals-7par-calo_%02i.png", ic);
  c1->SaveAs(fname);

  return 1;
}

int goSumCBO(){

  int ic;
  for (ic = 1; ic <= 24; ic++){
    sumCBO( ic);
    c1->Update();
    sleep(2);
  }

  return 1;
}

int sumCBO(int iCalo){

  // sum per-calo cbo histograms  with offset in time corresponding to cbo period/24
  //
  char fname[64], hname[64];

  TFile *file1 = TFile::Open("cbo-5par-residuals-highrate.root");
  //TFile *file1 = TFile::Open("cbo5par.root");
  //TFile *file1 = TFile::Open("cbo-6par-pileupterm-residuals-highrate.root");
  //TFile *file1 = TFile::Open("cbo-6par-pileupterm-residuals-lowrate.root");
  file1->GetObject("cbo_sumi", htmp1);
  //file1->GetObject("cbo_calo_01", htmp1);
  htmp2 = (TH1D*)htmp1->Clone();
  //htmp2->Reset();
  htmp2->SetName("CBO Sum");
  htmp2->SetTitle("CBO Sum");

  // omega_cbo = 3.22904e-03 from fit in raw bin units
  // omega_cbo-omega_a = 3.22904e-03-1.79982e-03 = 0.00142922 from fit in raw bin units
  // omega_cbo+omega_a = 3.22904e-03+1.79982e-03 = 0.00502886 from fit in raw bin units

  //double CBOPeriod = 0.0; // turn off shidting
  double CBOPeriod = 2.*3.14159265/3.22904e-03/htmp1->GetBinWidth(1); // omega_cbo
  //double CBOPeriod = 2.*3.14159265/0.00142922/htmp1->GetBinWidth(1); // omega_cbo-omega_a
  //double CBOPeriod = 2.*3.14159265/0.00502886/htmp1->GetBinWidth(1); // omega_cbo+omega_a

  double CBOPeriodOver24 = CBOPeriod/24.;
  printf(" CBOPeriod, CBOPeriodOver24 = %f, %f\n", CBOPeriod, CBOPeriodOver24);

  /*
  for (int ic = 1; ic <= 24; ic++){
    
    if (iCalo == 0 || iCalo == ic){

      sprintf( hname, "cbo_calo_%02i", ic-1);
      file1->GetObject(hname, htmp1);
      printf(" calo %i, offset %i, got histogram %s\n", ic, (int)(ic*CBOPeriodOver24), hname);
      
      for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++){
	
	int ibprime;
	if (iCalo == 0) {
	  ibprime = ib - (int)(ic*CBOPeriodOver24);
	} else {
	  ibprime = ib;
	}

	htmp2->SetBinContent( ibprime, htmp2->GetBinContent(ibprime)+htmp1->GetBinContent(ib));      
	htmp2->SetBinError( ibprime, htmp2->GetBinError(ibprime)*htmp2->GetBinError(ibprime)+htmp1->GetBinError(ib)*htmp1->GetBinError(ib));      
      }
      
      htmp2->SetLineColor(kBlue);
      htmp2->SetLineWidth(2.0);
    }
  }
  */

  htmp2->Draw("hist");

  // /* parameter for high rate data
  double A_cbo = 0.2951; // fractional residuals
  //double A_cbo = 4.0e6; // absolute residuals
  
  double tau_cbo = 1.433e5; // reactional residuals
  //double tau_cbo = 1.264e10; // make infinitely long
  //double tau_cbo = 5.0e4; // reactional residuals

  double omega_cbo = 2.320e-03;
  double phi_cbo = 0.0;
  double p0_cbo = 0.0;
  double delta_omega_cbo = -0.00457;
  // */

  delta_omega_cbo = 0.0;

  /* parameter for low rate data
  double A_cbo = 0.2951;
  double tau_cbo = 1.264e5;
  //double tau_cbo = 1.264e10; // make infinitely long
  double omega_cbo = 2.9023e-03;
  double phi_cbo = 0.0;
  double p0_cbo = 0.0;
  double delta_omega_cbo = -0.005;
  */

  double minT = 60000.;
  double maxT = 192000.; // reduced range from 195000 due to problems fitting

  // reference omega_cbo(t), tau_cbo to sample 58000
  cbof = new TF1("cbof", "[0] *  exp( -( x - 5.8e4 ) / [1] )  * sin( [2]*(1-[5]*(x-5.8e4)/(192000.-58000.))*(x-5.8e4) + [3] ) + [4]", minT, maxT);
  cbof->SetParameters( A_cbo, tau_cbo, omega_cbo, phi_cbo, p0_cbo, delta_omega_cbo);
  cbof->SetParNames( "A_cbo", "tau_cbo", "omega_cbo", "phi_cbo", "p0_cbo", "delta_omega_cbo" );
  cbof->SetParLimits(3,-3.14,3.14);
  //cbof->FixParameter(0,A_cbo);
  //cbof->FixParameter(1,tau_cbo);
  //cbof->FixParameter(2,omega_cbo);
  cbof->FixParameter(4,p0_cbo);
  cbof->FixParameter(5,delta_omega_cbo);
  //cbof->FixParameter(5,0.0);

  cbof->SetLineWidth(2);
  cbof->SetLineColor(kRed);
  cbof->SetNpx(10000);
  htmp2->SetLineWidth(2);
  htmp2->SetMaximum(0.04);
  htmp2->SetMinimum(-0.04);
  htmp2->GetXaxis()->SetRangeUser(30000, 205000);
  htmp2->Fit( cbof, "", "RIE", minT, maxT);
  htmp2->Fit( cbof, "", "RIE", minT, maxT);
  htmp2->Draw("HIST");
  cbof->Draw("same");

  //return 1;

  cbof->FixParameter(1,cbof->GetParameter(1));

  int nfit = 12;
  double DT = (maxT - minT) / nfit;
  double offsetDT = minT;
  double Xcbo[12], Ycbo[12], Ycbo2[12];
  for (int ifit = 1; ifit <= nfit; ifit++){
    htmp2->Fit( cbof, "", "WWRIE", offsetDT, offsetDT+DT);
    Xcbo[ifit-1] = offsetDT + DT/2.;
    Ycbo[ifit-1] = cbof->GetParameter(0);
    Ycbo2[ifit-1] = cbof->GetParameter(2);
    htmp2->Draw("hist");
    cbof->Draw("same");
    c1->Update();
    sleep(1);
    offsetDT += DT; 
  }

  c1->Clear();
  c1->Divide(1,2);

  c1->cd(1);
  gCBOamp = new TGraph(nfit, Xcbo, Ycbo);
  gCBOamp->SetName("gCBOamp");
  gCBOamp->SetTitle("cbo amplitude modifier versus time window");
  gCBOamp->GetXaxis()->SetTitle("time (samples)");
  gCBOamp->GetYaxis()->SetTitle("amplitude multiplier");
  gCBOamp->GetYaxis()->SetTitleOffset(0.5);
  gCBOamp->GetYaxis()->SetTitleSize(0.07);
  gCBOamp->GetXaxis()->SetTitleOffset(0.5);
  gCBOamp->GetXaxis()->SetTitleSize(0.07);
  gCBOamp->SetMarkerStyle(8);
  gCBOamp->SetMarkerSize(1.0);
  gCBOamp->SetMarkerColor(kRed+1);
  gCBOamp->SetLineWidth(2);
  gCBOamp->SetLineColor(kRed+1);
  gCBOamp->SetMinimum(0.0);
  gCBOamp->SetMaximum(1.0);
  //c1->Clear();
  gCBOamp->Draw("AP");

  c1->cd(2);
  gCBOomga = new TGraph(nfit, Xcbo, Ycbo2);
  gCBOomga->SetName("gCBOomga");
  gCBOomga->SetTitle("cbo omega versus time window");
  gCBOomga->GetXaxis()->SetTitle("time (samples)");
  gCBOomga->GetYaxis()->SetTitle("cbo omega");
  gCBOomga->GetYaxis()->SetTitleOffset(0.5);
  gCBOomga->GetYaxis()->SetTitleSize(0.07);
  gCBOomga->GetXaxis()->SetTitleOffset(0.5);
  gCBOomga->GetXaxis()->SetTitleSize(0.07);
  gCBOomga->SetMarkerStyle(8);
  gCBOomga->SetMarkerSize(1.0);
  gCBOomga->SetMarkerColor(kRed+1);
  gCBOomga->SetLineWidth(2);
  gCBOomga->SetLineColor(kRed+1);
  gCBOomga->SetMinimum(0.00310);
  gCBOomga->SetMaximum(0.00330);
  //c1->Clear();
  gCBOomga->Draw("AP");

  return 1;

  // TH1 *hm =0;TVirtualFFT::SetTransform(0);hm = hres->FFT(hm, "MAG"); hm->SetTitle("discrete fourier transform");hm->Draw();
  //Compute the transform and look at the magnitude of the output
  TVirtualFFT::SetTransform(0);

  TH1 *hm2 = 0;
  hm2 = htmp2->FFT(hm2, "MAG");
  hm2->SetTitle("discrete fourier transform");

  hmcalo = (TH1*)hm2->Clone();
  hmcalo->Clear();
  hmcalo->SetTitle("FFT of calo");
  hmcalo->SetName("caloFFT");

  hmcalosum = (TH1*)hm2->Clone();
  hmcalosum->Clear();
  hmcalosum->SetTitle("sum FFTs of calos");
  hmcalosum->SetName("sumcaloFFT");

  // sum of FFT's for each calo
  for (int ic = 1; ic <= 24; ic++){
    
    sprintf( hname, "cbo_calo_%02i", ic-1);
    file1->GetObject(hname, htmp1);
    hmcalo = htmp1->FFT(hmcalo, "MAG");
    hmcalo->Draw();
    hmcalosum->Add(hmcalo,1.0);

  }
  
  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  htmp2->Draw("hist");
  //cbof->Draw("same");
  c1->cd(2);
  //hmcalosum->Draw();
  hmcalosum->Draw();
 
  return 1;
}


int binning(){

  TFile *file1 = TFile::Open("noxtaldata/1912_00_1.root");
  TFile *file2 = TFile::Open("noxtaldata/1912_00_2.root");
  file1->GetObject("qHist1D_sum", htmp1);
  file2->GetObject("qHist1D_sum", htmp2);
  htmp1->Rebin(2);
  htmp1->Draw();
  htmp2->SetLineColor(kRed);
  htmp2->Draw("same");
  printf("rebin 1, integral %f\n", htmp1->Integral());
  printf("rebin 2, integral %f\n", htmp2->Integral());

  return 1;
}

int plotRatioNew(){

  char fname[64], hname[64];
  sprintf( fname,"noxtaldata/786x_w4_t16.root");
  TFile *file1 = TFile::Open( fname);

  sprintf( hname,"qHist1D_sig_sum_%i", 0);
  sprintf( hname,"qHist1D_ped_ent_sum_%i", 0);
  file1->GetObject( hname, htmp1);
  htmp1->Rebin(4);
  htmp1->SetName("sig_reb1");
  htmp1->SetLineColor(kRed);
  sprintf( hname,"qHist1D_sig_sum_%i", 1);
  sprintf( hname,"qHist1D_ped_ent_sum_%i", 1);
  file1->GetObject( hname, htmp2);
  htmp2->Rebin(2);
  htmp2->SetName("sig_reb2");
  htmp2->SetLineColor(kBlue);
  sprintf( hname,"qHist1D_sig_sum_%i", 2);
  sprintf( hname,"qHist1D_ped_ent_sum_%i", 2);
  file1->GetObject( hname, htmp3);
  htmp3->SetName("sig_rb4");
  htmp3->SetLineColor(kViolet);

  //htmp2->Divide(htmp1);
  //htmp3->Divide(htmp1);

  htmp1->SetLineWidth(2);
  htmp1->Draw("HIST");
  htmp2->SetLineWidth(2);
  htmp2->Draw("sameHIST");
  htmp3->SetLineWidth(2);
  htmp3->Draw("sameHIST");

  return 1;
}

int pedcorrection(int swrbindex);
int pedcorrection(int swrbindex){

  char fname[64], hname[64];

  sprintf( fname, "noxtaldata/v4w%it%iprd%s-60hr.root", iwind, ithreshold, aname);
  TFile *file1 = TFile::Open( fname, "UPDATE");
  printf("got file, %s\n", fname);

  sprintf( hname,"qHist1D_ped_usd_sum_%i", swrbindex);
  file1->GetObject( hname, htmp1);
  sprintf( hname,"pedestal_%i", swrbindex);
  htmp1->SetName(hname);
  htmp1->SetTitle(hname);
  printf("got histogram, %s\n", htmp1->GetName());

  htmp2 = (TH1D*) htmp1->Clone();
  sprintf( hname,"pedestalcorrection_%i", swrbindex);
  htmp2->SetName(hname);
  htmp2->SetTitle(hname);
  htmp2->Reset();
  printf("cloned histogram, %s\n", htmp2->GetName());
  
  int bookend = 4;
  for (int ibin = 1 + bookend; ibin <= htmp1->GetNbinsX() - bookend; ibin++) {
 
    printf("ibin %i\n", ibin);

    double ysum = 0.0;
    for (int jbin = ibin - bookend; jbin <= ibin + bookend; jbin++) {
      if (jbin != ibin) ysum += htmp1->GetBinContent(jbin);
    }
    htmp2->SetBinContent( ibin, htmp1->GetBinContent(ibin) - ysum / (2.*bookend) );
  }

  sprintf( hname,"qHist1D_sig_sum_%i", swrbindex);
  file1->GetObject( hname, htmp3);
  sprintf( hname,"pedestalcorrectedsignal_%i", swrbindex);
  htmp3->SetName(hname);
  htmp3->SetTitle(hname);
  printf("got histogram, %s\n", htmp3->GetName());
 
  htmp3->Add(htmp2, 1.0); // make pedestal correction

  htmp3->Write(); 
  printf("wrote histogram %s into file %s\n", htmp3->GetName(), fname);
  file1->Close();
  
  return 0;
}

int plotvertical();
int plotvertical(){

  TFile *file1 = TFile::Open( fname);
  printf("got file, %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got directory, %s\n", dname);

  sprintf( hname,"qHist1D_sig_vslice_%i", 0);
  dir1->GetObject( hname, htmp1);
  htmp1->Clear();
  htmp1->SetName("energysum");
  htmp1->SetTitle("energysum");
  printf("got histogram, %s\n", htmp1->GetName());

  sprintf( hname,"qHist1D_sig_vslice_%i", 0);
  dir1->GetObject( hname, htmp2);
  htmp2->Clear();
  htmp2->SetName("xcoordenergysum");
  htmp2->SetTitle("xcoordenergysum");
  printf("got histogram, %s\n", htmp2->GetName());

  for (int islice = 0; islice <= 8; islice++) {

    sprintf( hname,"qHist1D_sig_vslice_%i", islice);
    dir1->GetObject( hname, htmp3);
    printf("got slice vertical histogram, %s\n", htmp3->GetName());

    for (int ibin = 1; ibin <= htmp3->GetNbinsX(); ibin++) {
      double energy = htmp3->GetBinContent(ibin);
      //printf(" bin %i, energy %f, energy slice %f\n", ibin, energy, energy*(islice+1));
      htmp1->AddBinContent( ibin, energy );
      htmp2->AddBinContent( ibin, energy*(islice+1) );
    }
  }

  htmp2->Divide(htmp1);
  htmp2->Rebin(4);
  htmp2->Scale(0.25);
  htmp2->GetXaxis()->SetRangeUser(25.e3,200.e3);
  htmp2->SetLineColor(kRed);
  htmp2->SetLineWidth(2.0);
  htmp2->Draw("HIST");

  return 1;

}

TH1D* hhorizontalprofile;
int plotvertical2D(int iCaloMin, int iCaloMax);
int plotvertical2D(int iCaloMin, int iCaloMax){

  TFile *file1 = TFile::Open( fname);
  printf("got file, %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got directory, %s\n", dname);

  sprintf( hname,"qHist2D_sig_vslice_%i", 0);
  dir1->GetObject( hname, h2tmp);
  htmp1 = h2tmp->ProjectionX("htmp1",iCaloMin,iCaloMax,"");
  htmp1->Clear();
  htmp1->SetName("energysum");
  htmp1->SetTitle("energysum");
  printf("got histogram, %s\n", htmp1->GetName());

  sprintf( hname,"qHist2D_sig_vslice_%i", 0);
  dir1->GetObject( hname, h2tmp);
  htmp2 = h2tmp->ProjectionX("htmp2",iCaloMin,iCaloMax,"");
  htmp2->Clear();
  htmp2->SetName("xcoordenergysum");
  htmp2->SetTitle("xcoordenergysum");
  printf("got histogram, %s\n", htmp2->GetName());

  for (int islice = 0; islice <= 8; islice++) {

    sprintf( hname,"qHist2D_sig_vslice_%i", islice);
    dir1->GetObject( hname, h2tmp);
    htmp3 = h2tmp->ProjectionX("htmp3",iCaloMin,iCaloMax,"");
    printf("got slice vertical histogram, %s\n", htmp3->GetName());

    for (int ibin = 1; ibin <= htmp3->GetNbinsX(); ibin++) {
      double_t energy = htmp3->GetBinContent(ibin);
      //printf(" bin %i, energy %f, energy slice %f\n", ibin, energy, energy*(islice+1));
      htmp1->AddBinContent( ibin, energy );
      htmp2->AddBinContent( ibin, energy*(islice+1) );
    }
  }

  htmp2->Divide(htmp1);
  htmp2->Rebin(2);
  htmp2->Scale(0.5);


  hhorizontalprofile = new TH1D("horizontalProfile","horizontalProfile", 9, 0, 9);

  for (int islice = 0; islice <= 8; islice++) {
    
    sprintf( hname,"qHist2D_sig_vslice_%i", islice);
    dir1->GetObject( hname, h2tmp);
    htmp3 = h2tmp->ProjectionY( "htmp3", timeMin, timeMax, "");
    printf("got slice vertical histogram, %s\n", htmp3->GetName());
    double energy = htmp3->Integral( iCaloMin, iCaloMax);
    printf(" bin %i, energy %f\n", islice, energy);
    hhorizontalprofile->Fill( islice+0.5, energy );
  }

  hhorizontalprofile->Draw("hist");


  return 1;

}

TF1 *vbof;
int vohorizontal(int iCaloSlice, int ivslice);
int vohorizontal(int iCaloSlice = 0,int ivslice = 0){

  double startT = 55000., endT = 200000.;

  TFile *file1 = TFile::Open( fname);
  printf("got file, %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got directory, %s\n", dname);

  if (iCaloSlice == 0){ 
    sprintf( hname,"qHist1D_sig_hslice_%i", ivslice); // this is sum of all calos
    dir1->GetObject( hname, htmp1);
  } else {
    sprintf( hname,"qHist2D_sig_hslice_%i", ivslice); // this is for single calo
    dir1->GetObject( hname, h2tmp);
    htmp1 = h2tmp->ProjectionX("htmp1",iCaloSlice,iCaloSlice,"");
  }
  htmp1->SetName("row0");
  htmp1->SetTitle("row0");
  printf("got histogram, %s\n", htmp1->GetName());

  if (iCaloSlice == 0){ 
    sprintf( hname,"qHist1D_sig_hslice_%i", 5-ivslice); // this is sum of all calos
    dir1->GetObject( hname, htmp2);
  } else {
    sprintf( hname,"qHist2D_sig_hslice_%i", 5-ivslice); // this is for single calo
    dir1->GetObject( hname, h2tmp2);
    htmp2 = h2tmp2->ProjectionX("htmp2",iCaloSlice,iCaloSlice,"");
  }
  htmp2->SetName("row5");
  htmp2->SetTitle("row5");
  printf("got histogram, %s\n", htmp2->GetName());

  htmp2->Scale( htmp1->Integral( htmp1->FindBin(startT+25000), htmp1->FindBin(endT) ) / htmp2->Integral( htmp1->FindBin(startT+25000), htmp1->FindBin(endT) ) );
  htmp1->Add( htmp2, -1.0);
  htmp1->Draw();

  double A_vbo = 1.6e5, omega_vbo = 1.7707e-2, phi_vbo = 1.5, tau_vbo = 1.8363e4; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
  double A_vbo2 = 1.3e7, omega_vbo2 = 1.3339e-2, phi_vbo2 = 2.8, tau_vbo2 = 0.682e4; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
  double p0 = 1.0e4, p1 = 1.e5;
  double asym = 0.2, omega_a = 1.79991e-3, phi = 3.9;

  vbof = new TF1("vbof", fvo, 49000.0, 200000.0, 13);
  vbof->SetNpx(10000);
  vbof->SetParameters( A_vbo, tau_vbo, omega_vbo, phi_vbo, A_vbo2, tau_vbo2, omega_vbo2, phi_vbo2, p0, p1 );
  vbof->SetParNames( "A_vbo", "tau_vbo", "omega_vbo", "phi_vbo", "A_vbo2", "tau_vbo2", "omega_vbo2", "phi_vbo2", "p0", "p1" );
  vbof->SetParameter(10,asym);
  vbof->SetParName(10,"asym");
  vbof->SetParameter(11,omega_a);
  vbof->SetParName(11,"omega_a");
  vbof->SetParameter(12,phi);
  vbof->SetParName(12,"phi");

  vbof->FixParameter(11, omega_a);
  //vbof->FixParameter(8, 0.0);
  //vbof->FixParameter(12);

  htmp1->Fit( vbof, "", "VRIE", startT, endT); // 5 par fit
  htmp1->Draw();

  for (int ib = 1; ib <= htmp1->GetNbinsX(); ib++){
    double xc = htmp1->GetBinCenter(ib);
    if ( xc < 45000. || xc > 140000 ) htmp1->SetBinContent(ib,0.0);
  }
  
  TVirtualFFT::SetTransform(0);
  hm2 = htmp1->FFT(hm2, "MAG");
  hm2->SetTitle("discrete fourier transform");
  hm2->GetXaxis()->SetRange(0.,hm2->GetNbinsX()/2.); // plot is symmetric, just plot the lower half
  printf("draw residuals FFT\n");
  hm2->Draw("HIST");

  //file1->Close();
  return 0;
}

int plothorizontal();
int plothorizontal(){

  TFile *file1 = TFile::Open( fname);
  printf("got file, %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got directory, %s\n", dname);

  sprintf( hname,"qHist1D_sig_hslice_%i", 0);
  dir1->GetObject( hname, htmp1);
  htmp1->Clear();
  htmp1->SetName("energysum");
  htmp1->SetTitle("energysum");
  printf("got histogram, %s\n", htmp1->GetName());

  sprintf( hname,"qHist1D_sig_hslice_%i", 0);
  dir1->GetObject( hname, htmp2);
  htmp2->Clear();
  htmp2->SetName("xcoordenergysum");
  htmp2->SetTitle("xcoordenergysum");
  printf("got histogram, %s\n", htmp2->GetName());

  for (int islice = 0; islice <= 5; islice++) {

    sprintf( hname,"qHist1D_sig_hslice_%i", islice);
    dir1->GetObject( hname, htmp3);
    printf("got slice horizontal histogram, %s\n", htmp3->GetName());

    for (int ibin = 1; ibin <= htmp3->GetNbinsX(); ibin++) {
      double energy = htmp3->GetBinContent(ibin);
      //printf(" bin %i, energy %f, energy slice %f\n", ibin, energy, energy*(islice+1));
      htmp1->AddBinContent( ibin, energy );
      htmp2->AddBinContent( ibin, energy*(islice+1) );
    }
  }

  htmp2->Divide(htmp1);
  htmp2->Rebin(2);
  htmp2->Scale(0.5);
  htmp2->GetXaxis()->SetRangeUser(25.e3,200.e3);
  htmp2->SetLineColor(kRed);
  htmp2->SetLineWidth(2.0);
  htmp2->Draw("HIST");

  return 1;

}


//timeMin = 482.; timeMax = 615.; //30-40us
//timeMin = 2082.; timeMax = 2748.; //150-200us
//timeMin = 1415.; timeMax = 2081.; //100-150us
//timeMin = 748.; timeMax = 1414.; //50-100us
TH1D *hverticalprofile;

int plothorizontal2D(int iCaloMin, int iCaloMax);
int plothorizontal2D(int iCaloMin, int iCaloMax){

  TFile *file1 = TFile::Open( fname);
  printf("got file, %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got directory, %s\n", dname);

  sprintf( hname,"qHist2D_sig_hslice_%i", 0);
  dir1->GetObject( hname, h2tmp);
  htmp1 = h2tmp->ProjectionX("htmp1",iCaloMin,iCaloMax,"");
  htmp1->Clear();
  htmp1->SetName("energysum");
  htmp1->SetTitle("energysum");
  printf("got histogram, %s\n", h2tmp->GetName());

  sprintf( hname,"qHist2D_sig_hslice_%i", 0);
  dir1->GetObject( hname, h2tmp);
  htmp2 = h2tmp->ProjectionX("htmp2",iCaloMin,iCaloMax,"");
  htmp2->Clear();
  htmp2->SetName("xcoordenergysum");
  htmp2->SetTitle("xcoordenergysum");
  printf("got histogram, %s\n", htmp2->GetName());

  for (int islice = 0; islice <= 5; islice++) {

    sprintf( hname,"qHist2D_sig_hslice_%i", islice);
    dir1->GetObject( hname, h2tmp);
    htmp3 = h2tmp->ProjectionX("htmp3",iCaloMin,iCaloMax,"");
    printf("got slice horizontal histogram, %s\n", htmp3->GetName());

    for (int ibin = 1; ibin <= htmp3->GetNbinsX(); ibin++) {
      double energy = htmp3->GetBinContent(ibin);
      //printf(" bin %i, energy %f, energy slice %f\n", ibin, energy, energy*(islice+1));
      htmp1->AddBinContent( ibin, energy );
      htmp2->AddBinContent( ibin, energy*(islice+1) );
    }
  }

  htmp2->Divide(htmp1);
  htmp2->Rebin(2);
  htmp2->Scale(0.5);
  htmp2->Draw("hist");

  sprintf( hname,"qHist2D_sig_hslice_%i", 0);
  dir1->GetObject( hname, h2tmp);
  htmp4 = h2tmp->ProjectionY( "htmp4", timeMin, timeMax, "");
  htmp4->Reset();
  htmp4->SetName("verticalProfileTimeSlice");
  htmp4->SetTitle("verticalProfileTimSlice");
  printf("got histogram, %s\n", htmp4->GetName());

  hverticalprofile = new TH1D("verticalProfile","verticalProfile", 6, 0, 6);

  for (int islice = 0; islice <= 5; islice++) {
    
    sprintf( hname,"qHist2D_sig_hslice_%i", islice);
    dir1->GetObject( hname, h2tmp);
    htmp3 = h2tmp->ProjectionY( "htmp3", timeMin, timeMax, "");
    printf("got slice horizontal histogram, %s\n", htmp3->GetName());
    double energy = htmp3->Integral( iCaloMin, iCaloMax);
    printf(" bin %i, energy %f\n", islice, energy);
    hverticalprofile->Fill( islice+0.5, energy );
  }

  hverticalprofile->Draw("hist");

  return 1;

}

int plotRebin3New(int cal, int rb){
  
  char fname[64], hname[64];
  //sprintf( fname,"1912_xtal26_1f_%i.root", rb);
  sprintf( fname,"noxtaldata/7860_00_f1_w4_t16.root");
  TFile *file3 = TFile::Open( fname);

  TDirectoryFile *dir3;
  file3->GetObject( "CQBankAnalyzerTim", dir3);

  //dir3->GetObject( "qHist1D_12_12", htmp1);
  //dir3->GetObject( "qHist1D_sig_12_12", htmp2);
  //dir3->GetObject( "qHist1D_ped_all_12_12", htmp3);
  //dir3->GetObject( "qHist1D_ped_usd_12_12", htmp4);
  sprintf( hname,"qHist1D_%i_%i", cal, rb);
  dir3->GetObject( hname, htmp1);
  sprintf( hname,"qHist1D_sig_%i_%i", cal, rb);
  dir3->GetObject( hname, htmp2);
  sprintf( hname,"qHist1D_ped_all_%i_%i", cal, rb);
  dir3->GetObject( hname, htmp3);
  sprintf( hname,"qHist1D_ped_usd_%i_%i", cal, rb);
  dir3->GetObject( hname, htmp4);

  htmp1->GetXaxis()->SetRangeUser(50000.,70000.);
  //htmp1->SetMaximum(20.);
  //htmp1->SetMinimum(-80.);
  htmp1->SetMarkerColor(kGray+2);
  htmp1->SetMarkerStyle(8);
  htmp1->SetMarkerSize(0.7);
  htmp1->Draw("P");

  htmp2->SetLineColor(kRed);
  htmp2->SetLineWidth(2);
  htmp2->Draw("sameHIST");
  htmp3->SetLineColor(kGreen);
  htmp3->SetLineWidth(2);
  htmp3->Draw("sameHIST");
  htmp4->SetLineColor(kViolet);
  htmp4->SetLineWidth(2);
  //htmp4->Draw("sameHIST");

  sprintf( fname, "786x_rawpedsig_rebin%i.png", rb);
  c1->SaveAs( fname);

  return 1;
}

int plotRebin3(int rb){
  
  char fname[64];
  sprintf( fname,"1912_xtal26_1f_%i.root", rb);
  TFile *file3 = TFile::Open( fname);

  TDirectoryFile *dir3;
  file3->GetObject( "CQBankAnalyzerTim", dir3);

  dir3->GetObject( "qHist1D_12_12", htmp1);
  dir3->GetObject( "qHist1D_sig_12_12", htmp2);
  dir3->GetObject( "qHist1D_ped_all_12_12", htmp3);
  dir3->GetObject( "qHist1D_ped_usd_12_12", htmp4);

  htmp1->GetXaxis()->SetRangeUser(360000.,380000.);
  htmp1->SetMaximum(20.);
  htmp1->SetMinimum(-80.);
  htmp1->SetMarkerColor(kGray+2);
  htmp1->SetMarkerStyle(8);
  htmp1->SetMarkerSize(0.4);
  htmp1->Draw("P");

  htmp2->SetLineColor(kRed);
  htmp2->SetLineWidth(2);
  htmp2->Draw("sameHIST");
  htmp3->SetLineColor(kGreen);
  htmp3->SetLineWidth(2);
  htmp3->Draw("sameHIST");
  htmp4->SetLineColor(kViolet);
  htmp4->SetLineWidth(2);
  htmp4->Draw("sameHIST");

  sprintf( fname, "rawpedsig_rebin%i.png", rb);
  c1->SaveAs( fname);

  return 1;
}
int Nrb = 4; // no. of software rebinned histogram to process.

// function goCaloSum() calls summing of  designated calo histograms for a run group
// for 24 calo, no xtal-by-xtal histogram package
int goCaloSum(char* cRun){

  
  caloSumRb(cRun,"qHist1D", Nrb);
  caloSumRb(cRun,"qHist1D_ped_all", Nrb);
  caloSumRb(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRb(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRb(cRun,"qHist1D_sig", Nrb);
  caloSumRb(cRun,"qHist1D_puA", Nrb);
  caloSumRb(cRun,"qHist1D_puL", Nrb);
  caloSumRb(cRun,"qHist1D_puH", Nrb);
  caloSum(cRun,"qHist1D_rms");
  caloSum(cRun,"qHist1D_noise");
  
  caloSumRbHot(cRun,"qHist1D", Nrb);
  caloSumRbHot(cRun,"qHist1D_ped_all", Nrb);
  caloSumRbHot(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRbHot(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRbHot(cRun,"qHist1D_sig", Nrb);
  caloSumRbHot(cRun,"qHist1D_puA", Nrb);
  caloSumRbHot(cRun,"qHist1D_puL", Nrb);
  caloSumRbHot(cRun,"qHist1D_puH", Nrb);
  caloSumHot(cRun,"qHist1D_rms");
  caloSumHot(cRun,"qHist1D_noise");

  caloSumRbCold(cRun,"qHist1D", Nrb);
  caloSumRbCold(cRun,"qHist1D_ped_all", Nrb);
  caloSumRbCold(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRbCold(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRbCold(cRun,"qHist1D_sig", Nrb);
  caloSumRbCold(cRun,"qHist1D_puA", Nrb);
  caloSumRbCold(cRun,"qHist1D_puL", Nrb);
  caloSumRbCold(cRun,"qHist1D_puH", Nrb);
  caloSumCold(cRun,"qHist1D_rms");
  caloSumCold(cRun,"qHist1D_noise");

  //caloSum2D(cRun,"qHist2D_ydiff"); // needs 2D version
  //caloSum2D(cRun,"qHist2D_Coarse"); // needs 2D version

  return 1;
}

// function goCaloSum() calls summing of  designated calo histograms for a run group
// for 24 calo, no xtal-by-xtal histogram package
int goCaloSumPairs(char* cRun){

  caloSumRbPairs(cRun,"qHist1D", Nrb);
  caloSumRbPairs(cRun,"qHist1D_ped_all", Nrb);
  caloSumRbPairs(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRbPairs(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRbPairs(cRun,"qHist1D_sig", Nrb);
  caloSumRbPairs(cRun,"qHist1D_puA", Nrb);
  caloSumRbPairs(cRun,"qHist1D_puL", Nrb);
  caloSumRbPairs(cRun,"qHist1D_puH", Nrb);

  return 1;
}

// function goCaloSum() calls summing of  designated calo histograms for a run group
// for 24 calo, no xtal-by-xtal histogram package
int goCaloSumNear(char* cRun){

  caloSumRbNear(cRun,"qHist1D", Nrb);
  caloSumRbNear(cRun,"qHist1D_ped_all", Nrb);
  caloSumRbNear(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRbNear(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRbNear(cRun,"qHist1D_sig", Nrb);
  caloSumRbNear(cRun,"qHist1D_puA", Nrb);
  caloSumRbNear(cRun,"qHist1D_puL", Nrb);
  caloSumRbNear(cRun,"qHist1D_puH", Nrb);

  return 1;
}

// function goCaloSum() calls summing of  designated calo histograms for a run group
// for 24 calo, no xtal-by-xtal histogram package
int goCaloSumFar(char* cRun){

  caloSumRbFar(cRun,"qHist1D", Nrb);
  caloSumRbFar(cRun,"qHist1D_ped_all", Nrb);
  caloSumRbFar(cRun,"qHist1D_ped_usd", Nrb);
  caloSumRbFar(cRun,"qHist1D_ped_ent", Nrb);
  caloSumRbFar(cRun,"qHist1D_sig", Nrb);
  caloSumRbFar(cRun,"qHist1D_puA", Nrb);
  caloSumRbFar(cRun,"qHist1D_puL", Nrb);
  caloSumRbFar(cRun,"qHist1D_puH", Nrb);

  return 1;
}

// sums designated calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSum(char* cRun, char* cHist){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;  
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"%s_%i_%i", cHist, 1, 1);
  sprintf( hname,"%s_%i", cHist, 1);
  dir1->GetObject( hname, hSum);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"%s_sum", cHist);
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  hSum->Reset();

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if (iCal == 22) continue; // different TQ setting

    sprintf( hname,"%s_%i_%i", cHist, iCal, iCal);
    sprintf( hname,"%s_%i", cHist, iCal);
    dir1->GetObject( hname, h1[iCal]);
    printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
    hSum->Add( h1[iCal], 1.0);
  }

  hSum->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);

  file1->Close();
  return 1;
}

// sums designated calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumHot(char* cRun, char* cHist){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"%s_%i_%i", cHist, 1, 1);
  sprintf( hname,"%s_%i", cHist, 1);
  dir1->GetObject( hname, hSum);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"%s_sum_hot", cHist);
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  hSum->Reset();

  for (int iCal = 1; iCal <= 24; iCal++) {

    if (iCal >= 7 && iCal <= 18) continue; // different TQ setting
    // keep calo's 1-8, 24, reject 9-23

    sprintf( hname,"%s_%i_%i", cHist, iCal, iCal);
    sprintf( hname,"%s_%i", cHist, iCal);
    dir1->GetObject( hname, h1[iCal]);
    printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
    hSum->Add( h1[iCal], 1.0);
  }

  hSum->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);

  file1->Close();
  return 1;
}

// sums designated calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumSlice(char* cRun){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist1D_sig_vslice_0", 0);
  dir1->GetObject( hname, hSum);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"sum_vslice");
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  hSum->Reset();

  for (int iCal = 0; iCal <= 8; iCal++) {

    sprintf( hname,"qHist1D_sig_vslice_%i", iCal);
    dir1->GetObject( hname, h1[iCal]);
    printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
    hSum->Add( h1[iCal], 1.0);
  }

  hSum->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);

  file1->Close();
  return 1;
}

// sums designated calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumCold(char* cRun, char* cHist){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"%s_%i_%i", cHist, 1, 1);
  sprintf( hname,"%s_%i", cHist, 1);
  dir1->GetObject( hname, hSum);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"%s_sum_cold", cHist);
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  hSum->Reset();

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if (iCal < 14 || iCal == 22 || iCal == 24) continue; 
   if (iCal < 7 || iCal > 18) continue; // different TQ setting
    // keep 15-23, reject 1-14, 22, 24

    sprintf( hname,"%s_%i_%i", cHist, iCal, iCal);
    sprintf( hname,"%s_%i", cHist, iCal);
    dir1->GetObject( hname, h1[iCal]);
    printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
    hSum->Add( h1[iCal], 1.0);
  }

  hSum->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);

  file1->Close();
  return 1;
}

// sums designated calo 2D histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSum2D(char* cRun, char* cHist){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"%s_%i_%i", cHist, 1, 1);
  sprintf( hname,"%s_%i", cHist, 1);
  dir1->GetObject( hname, sSum);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"%s_sum", cHist);
  sSum->SetName(hsname);
  sSum->SetTitle(hsname);
  sSum->Reset();

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if (iCal == 22) continue; // different TQ setting - removed Aug 12, 2019

    sprintf( hname,"%s_%i_%i", cHist, iCal, iCal);
    sprintf( hname,"%s_%i", cHist, iCal);
    dir1->GetObject( hname, h2tmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
    sSum->Add( h2tmp, 1.0);
  }

  sSum->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);

  file1->Close();
  return 1;
}

int caloSumDoPU(char* cRun, int rbmax){

  char hname[64];
  TH1D *hSumPUA, *hSumPULH;

  // INPUT / OUTPUT FILE

  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");

  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);
  
  for (int irb = 0; irb < rbmax; irb++) {

    // all calos
    
    sprintf( hname,"qHist1D_sig_sum_%i", irb);
    file1->GetObject( hname, hSumPUA);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puA_sum_%i", irb);
    hSumPUA->SetName(hsname);
    hSumPUA->SetTitle(hsname);

    sprintf( hname,"qHist1D_sig_sum_%i", irb);
    file1->GetObject( hname, hSumPULH);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puLH_sum_%i", irb);
    hSumPULH->SetName(hsname);
    hSumPULH->SetTitle(hsname);
    
    sprintf( hname,"qHist1D_puA_sum_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPUA->Add( htmp, 1.0);

    hSumPUA->Write(); 
    printf("wrote A pileup-corrected histogram %s into root file %s\n", hsname, fname);

    sprintf( hname,"qHist1D_puL_sum_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupL histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    sprintf( hname,"qHist1D_puH_sum_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    hSumPULH->Write(); 
    printf("wrote LH pileup-corrected histogram %s into root file %s\n", hsname, fname);

    // hot calos's

    sprintf( hname,"qHist1D_sig_sum_hot_%i", irb);
    file1->GetObject( hname, hSumPUA);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puA_sum_hot_%i", irb);
    hSumPUA->SetName(hsname);
    hSumPUA->SetTitle(hsname);

    sprintf( hname,"qHist1D_sig_sum_hot_%i", irb);
    file1->GetObject( hname, hSumPULH);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puLH_sum_hot_%i", irb);
    hSumPULH->SetName(hsname);
    hSumPULH->SetTitle(hsname);
    
    sprintf( hname,"qHist1D_puA_sum_hot_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPUA->Add( htmp, 1.0);

    hSumPUA->Write(); 
    printf("wrote A pileup-corrected histogram %s into root file %s\n", hsname, fname);

    sprintf( hname,"qHist1D_puL_sum_hot_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupL histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    sprintf( hname,"qHist1D_puH_sum_hot_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    hSumPULH->Write(); 
    printf("wrote LH pileup-corrected histogram %s into root file %s\n", hsname, fname);

    // cold calos's

    sprintf( hname,"qHist1D_sig_sum_cold_%i", irb);
    file1->GetObject( hname, hSumPUA);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puA_sum_cold_%i", irb);
    hSumPUA->SetName(hsname);
    hSumPUA->SetTitle(hsname);

    sprintf( hname,"qHist1D_sig_sum_cold_%i", irb);
    file1->GetObject( hname, hSumPULH);
    printf("got signal histogram %s from root file %s\n", hname, fname);
    sprintf( hsname,"qHist1D_sig_puLH_sum_cold_%i", irb);
    hSumPULH->SetName(hsname);
    hSumPULH->SetTitle(hsname);
    
    sprintf( hname,"qHist1D_puA_sum_cold_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPUA->Add( htmp, 1.0);

    hSumPUA->Write(); 
    printf("wrote A pileup-corrected histogram %s into root file %s\n", hsname, fname);

    sprintf( hname,"qHist1D_puL_sum_cold_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupL histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    sprintf( hname,"qHist1D_puH_sum_cold_%i", irb);
    file1->GetObject( hname, htmp);
    printf("got pileupA histogram %s from root file %s\n", hname, fname);
    hSumPULH->Add( htmp, 0.5);

    hSumPULH->Write(); 
    printf("wrote LH pileup-corrected histogram %s into root file %s\n", hsname, fname);

  }
  
  
  
  file1->Close();
  
  return 1;
}

// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRb(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    sprintf( hname,"%s_%i_%i", cHist, 1, irb);
    printf("get histogram %s from file %s\n", hname, fname);
    dir1->GetObject( hname, hSum);
    printf("got histogram %s from file %s\n", hname, fname);
    sprintf( hsname,"%s_sum_%i", cHist, irb);
    hSum->SetName(hsname);
    hSum->SetTitle(hsname);
    hSum->Reset();

    for (int iCal = 1; iCal <= 24; iCal++) {

      //if (iCal == 22) continue; // different TQ setting - removed Aug 12, 2019!

      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, h1[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hSum->Add( h1[iCal], 1.0);
    }

    hSum->Write(); 
    printf("wrote histogram %s into file %s\n", hsname, foutname);
  }

  file1->Close();
  
  return 1;
}

// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRbPairs(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    for (int iCal = 1; iCal <= 12; iCal++) {
      
      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, hpair[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      sprintf( hname,"%s_%i_%i", cHist, iCal+12, irb);
      dir1->GetObject( hname, htmp);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hpair[iCal]->Add( htmp, 1.0);
      sprintf( hsname,"%s_pair_%i_%i", cHist, iCal, irb);
      hpair[iCal]->SetName(hsname);
      hpair[iCal]->SetTitle(hsname);
      hpair[iCal]->Write();
      printf("wrote histogram %s into file %s\n", hsname, foutname);
    }

  }

  file1->Close();
  
  return 1;
}

// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRbNear(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    sprintf( hname,"%s_%i_%i", cHist, 1, irb);
    dir1->GetObject( hname, hSum);
    printf("got histogram %s from file %s\n", hname, fname);
    sprintf( hsname,"%s_sum_near_%i", cHist, irb);
    hSum->SetName(hsname);
    hSum->SetTitle(hsname);
    hSum->Reset();

    for (int iCal = 1; iCal <= 12; iCal++) {
      

      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, h1[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hSum->Add( h1[iCal], 1.0);
    }

    hSum->Write(); 
    printf("wrote histogram %s into file %s\n", hsname, foutname);
  }

  file1->Close();
  
  return 1;
}


// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRbFar(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    sprintf( hname,"%s_%i_%i", cHist, 1, irb);
    dir1->GetObject( hname, hSum);
    printf("got histogram %s from file %s\n", hname, fname);
    sprintf( hsname,"%s_sum_far_%i", cHist, irb);
    hSum->SetName(hsname);
    hSum->SetTitle(hsname);
    hSum->Reset();

    for (int iCal = 13; iCal <= 24; iCal++) {
      

      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, h1[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hSum->Add( h1[iCal], 1.0);
    }

    hSum->Write(); 
    printf("wrote histogram %s into file %s\n", hsname, foutname);
  }

  file1->Close();
  
  return 1;
}


// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRbHot(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    sprintf( hname,"%s_%i_%i", cHist, 1, irb);
    dir1->GetObject( hname, hSum);
    printf("got histogram %s from file %s\n", hname, fname);
    sprintf( hsname,"%s_sum_hot_%i", cHist, irb);
    hSum->SetName(hsname);
    hSum->SetTitle(hsname);
    hSum->Reset();

    for (int iCal = 1; iCal <= 24; iCal++) {
      
      if (iCal > 8 && iCal < 24) continue; // different TQ setting
      // keep calo's 1-8, 24, reject 9-23

      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, h1[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hSum->Add( h1[iCal], 1.0);
    }

    hSum->Write(); 
    printf("wrote histogram %s into file %s\n", hsname, foutname);
  }

  file1->Close();
  
  return 1;
}

// sums top3 three horizontal slice for half calo
int caloSumBot(char* cRun){

  char hsname[64];
 
  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist2D_sig_hslice_%01i", 0);
  dir1->GetObject( hname, hbot);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist2D_sig_bot");
  hbot->SetName(hsname);
  hbot->SetTitle(hsname);
  hbot->Reset();

  for (int iSlice = 0; iSlice <= 2; iSlice++) {

    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist2D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, h2tmp);
    printf("got run %s, slice %i, histogram %s\n", cRun, iSlice, hname);
    hbot->Add( h2tmp, 1.0);
  }

  hbot->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// sums top3 three horizontal slice for half calo
int caloSumBot1D(char* cRun);
int caloSumBot1D(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist1D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, hbot1D);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist1D_sig_bot");
  hbot1D->SetName(hsname);
  hbot1D->SetTitle(hsname);
  hbot1D->Reset();

  for (int iSlice = 0; iSlice <= 2; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist1D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, htmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    hbot1D->Add( htmp, 1.0);
  }

  hbot1D->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// middle two rows
int caloSumCore2(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist2D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist2D_sig_core2");
  htop->SetName(hsname);
  htop->SetTitle(hsname);
  htop->Reset();

  for (int iSlice = 2; iSlice <= 3; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist2D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, h2tmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop->Add( h2tmp, 1.0);
  }

  htop->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// middle four rows
int caloSumCore4(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist2D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist2D_sig_core4");
  htop->SetName(hsname);
  htop->SetTitle(hsname);
  htop->Reset();

  for (int iSlice = 1; iSlice <= 4; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist2D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, h2tmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop->Add( h2tmp, 1.0);
  }

  htop->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}


// sums top3 three horizontal slice for half calo
int caloSumTop(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist2D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist2D_sig_top");
  htop->SetName(hsname);
  htop->SetTitle(hsname);
  htop->Reset();

  for (int iSlice = 3; iSlice <= 5; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist2D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, h2tmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop->Add( h2tmp, 1.0);
  }

  htop->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// sums top3 three horizontal slice for half calo
int caloSumTop1D(char* cRun);
int caloSumTop1D(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist1D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop1D);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist1D_sig_top");
  htop1D->SetName(hsname);
  htop1D->SetTitle(hsname);
  htop1D->Reset();

  for (int iSlice = 3; iSlice <= 5; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist1D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, htmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop1D->Add( htmp, 1.0);
  }

  htop1D->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// sums top3 three horizontal slice for half calo
int caloSumCore21D(char* cRun);
int caloSumCore21D(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist1D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop1D);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist1D_sig_core2");
  htop1D->SetName(hsname);
  htop1D->SetTitle(hsname);
  htop1D->Reset();

  for (int iSlice = 2; iSlice <= 3; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist1D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, htmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop1D->Add( htmp, 1.0);
  }

  htop1D->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// sums top3 three horizontal slice for half calo
int caloSumCore41D(char* cRun);
int caloSumCore41D(char* cRun){

  char hsname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname,"qHist1D_sig_hslice_%01i", 3);
  dir1->GetObject( hname, htop1D);
  printf("got histogram %s from file %s\n", hname, fname);
  sprintf( hsname,"qHist1D_sig_core4");
  htop1D->SetName(hsname);
  htop1D->SetTitle(hsname);
  htop1D->Reset();

  for (int iSlice = 1; iSlice <= 4; iSlice++) {
    
    printf("slice %i\n",iSlice);

    sprintf( hname,"qHist1D_sig_hslice_%01i", iSlice);
    dir1->GetObject( hname, htmp);
    printf("got run %s, calo %i, histogram %s\n", cRun, iSlice, hname);
    htop1D->Add( htmp, 1.0);
  }

  htop1D->Write(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file1->Close();
  
  return 1;
}

// sums designated software rebinned calo histogram cHist for a run group
// for 24 calo, no xtal-by-xtal histogram package
int caloSumRbCold(char* cRun, char* cHist, int rbmax){

  char hname[64];

  // INPUT / OUTPUT FILE
  sprintf( fname, "noxtaldata/%s.root", cRun);
  file1 = TFile::Open( fname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", fname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  for (int irb = 0; irb < rbmax; irb++) {

    sprintf( hname,"%s_%i_%i", cHist, 1, irb);
    dir1->GetObject( hname, hSum);
    printf("got histogram %s from file %s\n", hname, fname);
    sprintf( hsname,"%s_sum_cold_%i", cHist, irb);
    hSum->SetName(hsname);
    hSum->SetTitle(hsname);
    hSum->Reset();

    for (int iCal = 1; iCal <= 24; iCal++) {
      
      if (iCal < 14 || iCal == 24) continue; 
      // keep 15-23, reject 1-14, 24

      sprintf( hname,"%s_%i_%i", cHist, iCal, irb);
      dir1->GetObject( hname, h1[iCal]);
      printf("got run %s, calo %i, histogram %s\n", cRun, iCal, hname);
      hSum->Add( h1[iCal], 1.0);
    }

    hSum->Write(); 
    printf("wrote histogram %s into file %s\n", hsname, foutname);
  }

  file1->Close();
  
  return 1;
}


// function call goHist() to sum the xtals of a calo and goSumHist() to
// sum the calos of a run / run group
int goAll(){

  for (int iRun = 1839; iRun<=1839; iRun++){
    //goScat(iRun);
    //goSumScat(iRun);
  }
  goHist(0);
  goSumHist(0);

  return 1;
}

// function goHist() calls summing of  designated histograms for the xtals in each calorimeter
int goHist(int iRun){

    goXtal(iRun,"qHist1D");
    goXtal(iRun,"qHist1D_ped");
    goXtal(iRun,"qHist1D_sig");
    //goXtal(iRun,"qHist1D_thr");
    
    return 1;
}


// function goSumHist() calls summing of  designated histograms for the calos in a run group
int goSumHist(int iRun){

    sumCalo(iRun,"qHist1D");
    sumCalo(iRun,"qHist1D_ped");
    sumCalo(iRun,"qHist1D_sig");
    //sumCalo(iRun,"qHist1D_thr");

    //sumCalo(iRun,"qHist2D");
    //sumCalo(iRun,"qHist2D_Coarse");
    
    return 1;
}

int goSumXtalHist(int iRun){
  
  //sumCaloXtal(iRun,"qHist1D");
  sumCaloXtal(iRun,"qHist1D_thr");
  //sumCaloXtal(iRun,"qHist1D_ped");
  //sumCaloXtal(iRun,"qHist1D_sig");
  //sumCalo(iRun,"qHist2D");
  //sumCalo(iRun,"qHist2D_Coarse");
  
  return 1;
}

int goSumScat(int iRun){

    sumCalo2D(iRun,"qHist2D");
    sumCalo2D(iRun,"qHist2D_Coarse");
    
    return 1;
}


int goScat(int iRun){

   goXtal2D(iRun,"qHist2D");
   goXtal2D(iRun,"qHist2D_Coarse");
   
   return 1;
}

int goXtal(int iRun, char* cHist){

  /*
call sumXtal to sum the histograms for all xtals for particular histogram
of every calo of particular run
  */

  for (int iCal = 1; iCal <= 24; iCal++) {

    // dont need to wed out bad histo's with xtal-dependent threshold from rms
    //if ( iCal == 9 ) continue; // skip calo 9 - ??? 24 Nov 2017
    //if ( iCal == 12 ) continue; // skip calo 12 - this calo is down 6/12/17
    //if ( iCal == 13 ) continue; // skip calo 13 - ??? 24 Nov 2017
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    sumXtal( iRun, iCal, cHist);
  }

  file1->Close();

  return 1;
}

int goXtal2D(int iRun, char* cHist){

  /*
call sumXtal to sum the histograms for all xtals for particular histogram
of every calo of particular run
  */

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 14 - this calo is down 6/12/17
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    sumXtal2D( iRun, iCal, cHist);
  }

  file1->Close();

  return 1;
}

int plotPU(){
  
  char cname[64], finname[64];

  sprintf( finname, "swflsh1.root");
  file1 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_1_26");
  dir1->GetObject( hname, htmp1);
  printf("got histogram %s\n", hname);
  htmp1->SetName("pile-up1");
  htmp1->SetTitle("pile-up 1");
  htmp1->SetLineColor(kViolet);

  sprintf( finname, "swflsh2.root");
  file2 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir2;
  file2->GetObject( dname, dir2);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_sig_1_26");
  dir2->GetObject( hname, htmp2);
  printf("got histogram %s\n", hname);
  htmp2->SetName("pile-up2");
  htmp2->SetTitle("pile-up 2");
  htmp2->SetLineColor(kGreen);
  
  sprintf( finname, "swflsh4.root");
  file3 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir3;
  file3->GetObject( dname, dir3);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_1_26");
  dir3->GetObject( hname, htmp3);
  printf("got histogram %s\n", hname);
  htmp3->SetName("pile-up3");
  htmp3->SetTitle("pile-up 3");
  htmp3->SetLineColor(kRed);

  htmp3->Draw("");
  htmp2->Draw("same");
  htmp1->Draw("same");
  
  printf("done plotting pile-up 1, 2, 4, 8\n");

  return 1;
}

int plotThres( char* iRun){
  
  char cname[64], finname[64];

  sprintf( finname, "%04s/run0%04s_sum_rmsthres_nocalibrate.root", iRun, iRun);
  file1 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  sprintf( hname, "qHist1D_sig_sum_noBadCalo");
  file1->GetObject( hname, htmp1);
  htmp1->SetName("threshold1");
  htmp1->SetTitle("threshold1");
  htmp1->SetLineColor(kViolet);

  sprintf( finname, "%04s/run0%04s_sum_rmsthres2_nocalibrate.root", iRun, iRun);
  file2 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  sprintf( hname, "qHist1D_sig_sum_noBadCalo");
  file2->GetObject( hname, htmp2);
  htmp2->SetName("threshold2");
  htmp2->SetTitle("threshold2");
  htmp2->SetLineColor(kGreen);
  
  sprintf( finname, "%04s/run0%04s_sum_rmsthres3_nocalibrate.root", iRun, iRun);
  file3 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  sprintf( hname, "qHist1D_sig_sum_noBadCalo");
  file3->GetObject( hname, htmp3);
  htmp3->SetName("threshold3");
  htmp3->SetTitle("threshold3");
  htmp3->SetLineColor(kRed);

  htmp1->Rebin(16);
  htmp2->Rebin(16);
  htmp3->Rebin(16);

  htmp3->Draw();
  htmp2->Draw("same");
  htmp1->Draw("same");

  sprintf( hname, "thresPLOT_%04s.png", iRun);
  c1->SaveAs( hname);
  printf("done plotting thresholds 1, 2, 3\n");

  return 1;
}

int plotEntries(int thresRMS){
  
  char cname[64], finname[64];

  int rb = 1;
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file1 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_ped_ent_sum");
  file1->GetObject( hname, htmp1b);
  sprintf( hname, "qHist1D_ped_usd_sum");
  file1->GetObject( hname, htmp1);
  htmp1->Divide(htmp1b);
  sprintf( hname, "1912 rebin %i", rb);
  htmp1->SetName(hname);
  htmp1->SetTitle(hname);
  htmp1->Rebin(8/rb);
  htmp1->SetLineColor(kViolet);

  int minT = 530000, maxT =550000;
  pol0f = new TF1("pol0f","[0]",0,560000);
  pol0f->SetParNames( "pol0");
  pol0f->SetParameters( 0, 1.0);

  //htmp1->Fit( pol0f, "", "WR", minT, maxT);
  //noise1 = pol0f->GetParameter(0);
  //htmp1->Add(pol0f,-1.0,"");

  rb = 2;
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file2 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_ped_ent_sum");
  file2->GetObject( hname, htmp2b);
  sprintf( hname, "qHist1D_ped_usd_sum");
  file2->GetObject( hname, htmp2);
  htmp2->Divide(htmp2b);
  sprintf( hname, "1912 rebin %i", rb);
  htmp2->SetName(hname);
  htmp2->SetTitle(hname);
  htmp2->Rebin(8/rb);
  htmp2->SetLineColor(kOrange);

  //htmp2->Fit( pol0f, "", "WR", minT, maxT);
  //noise2 = pol0f->GetParameter(0);
  //htmp2->Add(pol0f,-1.0,"");

  rb = 4;
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file3 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_ped_ent_sum");
  file3->GetObject( hname, htmp3b);
  sprintf( hname, "qHist1D_ped_usd_sum");
  file3->GetObject( hname, htmp3);
  htmp3->Divide(htmp3b);
  sprintf( hname, "1912 rebin %i", rb);
  htmp3->SetName(hname);
  htmp3->SetTitle(hname);
  htmp3->Rebin(8/rb);
  htmp3->SetLineColor(kRed);

  //htmp3->Fit( pol0f, "", "WR", minT, maxT);
  //noise3 = pol0f->GetParameter(0);
  //htmp3->Add(pol0f,-1.0,"");

  rb = 8;
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file4 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  sprintf( hname, "qHist1D_ped_ent_sum");
  file4->GetObject( hname, htmp4b);
  sprintf( hname, "qHist1D_ped_usd_sum");
  file4->GetObject( hname, htmp4);
  htmp4->Divide(htmp4b);
  sprintf( hname, "1912 rebin %i", rb);
  htmp4->SetName(hname);
  htmp4->SetTitle(hname);
  htmp4->Rebin(8/rb);
  htmp4->SetLineColor(kCyan);

  //htmp4->Fit( pol0f, "", "WR", minT, maxT);
  //noise4 = pol0f->GetParameter(0);
  //htmp4->Add(pol0f,-1.0,"");

  htmp1->Rebin(2);
  htmp2->Rebin(2);
  htmp3->Rebin(2);
  htmp4->Rebin(2);

  //htmp4->SetMaximum(1.e5);
  //htmp4->SetMinimum(-2.5e6);

  htmp4->Draw();
  htmp3->Draw("same");
  htmp2->Draw("same");
  htmp1->Draw("same");

  sprintf( hname, "entriesPLOT_%ithresRMS.png",thresRMS);
  c1->SaveAs( hname);
  printf("done plotting entries for rebins 1, 2, 4, 8\n");

  return 1;
}

int plotRebin(int thresRMS){
  
  char cname[64], finname[64];

  int rb = 1;
  sprintf( finname, "noxtaldata/1912rb%i_%i.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw4.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file1 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_sum");
  //sprintf( hname, "qHist1D_sig_sum");
  //sprintf( hname, "qHist1D_ped_sum");
  //sprintf( hname, "qHist1D_ped_all_sum");
  file1->GetObject( hname, htmp1);
  sprintf( hname, "1912 rebin %i", rb);
  htmp1->SetName(hname);
  htmp1->SetTitle(hname);
  htmp1->Rebin(8/rb);
  htmp1->SetLineColor(kViolet);

  int minT = 530000, maxT =550000;
  pol0f = new TF1("pol0f","[0]",0,560000);
  pol0f->SetParNames( "pol0");
  pol0f->SetParameters( 0, 1.0);

  htmp1->Fit( pol0f, "", "WR", minT, maxT);
  noise1 = pol0f->GetParameter(0);
  //htmp1->Add(pol0f,-1.0,"");

  rb = 2;
  sprintf( finname, "noxtaldata/1912rb%i_%i.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw4.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file2 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_sum");
  //sprintf( hname, "qHist1D_sig_sum");
  //sprintf( hname, "qHist1D_ped_sum");
  //sprintf( hname, "qHist1D_ped_all_sum");
  file2->GetObject( hname, htmp2);
  sprintf( hname, "1912 rebin %i", rb);
  htmp2->SetName(hname);
  htmp2->SetTitle(hname);
  htmp2->Rebin(8/rb);
  htmp2->SetLineColor(kOrange);

  htmp2->Fit( pol0f, "", "WR", minT, maxT);
  noise2 = pol0f->GetParameter(0);
  //htmp2->Add(pol0f,-1.0,"");

  rb = 4;
  sprintf( finname, "noxtaldata/1912rb%i_%i.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw4.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file3 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  sprintf( hname, "qHist1D_sum");
  //sprintf( hname, "qHist1D_sig_sum");
  //sprintf( hname, "qHist1D_ped_sum");
  //sprintf( hname, "qHist1D_ped_all_sum");
  file3->GetObject( hname, htmp3);
  sprintf( hname, "1912 rebin %i", rb);
  htmp3->SetName(hname);
  htmp3->SetTitle(hname);
  htmp3->Rebin(8/rb);
  htmp3->SetLineColor(kRed);

  htmp3->Fit( pol0f, "", "WR", minT, maxT);
  double noise3 = pol0f->GetParameter(0);
  //htmp3->Add(pol0f,-1.0,"");

  rb = 8;
  sprintf( finname, "noxtaldata/1912rb%i_%i.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw4.root", rb, thresRMS);
  sprintf( finname, "noxtaldata/1912rb%i_%i_wndw16.root", rb, thresRMS);
  file4 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  sprintf( hname, "qHist1D_sum");
  //sprintf( hname, "qHist1D_sig_sum");
  //sprintf( hname, "qHist1D_ped_sum");
  //sprintf( hname, "qHist1D_ped_all_sum");
  file4->GetObject( hname, htmp4);
  sprintf( hname, "1912 rebin %i", rb);
  htmp4->SetName(hname);
  htmp4->SetTitle(hname);
  htmp4->Rebin(8/rb);
  htmp4->SetLineColor(kCyan);

  htmp4->Fit( pol0f, "", "WR", minT, maxT);
  noise4 = pol0f->GetParameter(0);
  //htmp4->Add(pol0f,-1.0,"");

  htmp1->Rebin(2);
  htmp2->Rebin(2);
  htmp3->Rebin(2);
  htmp4->Rebin(2);

  htmp4->SetMaximum(1.e5);
  htmp4->SetMinimum(-2.5e6);


  htmp4->Draw();
  htmp3->Draw("same");
  htmp2->Draw("same");
  htmp1->Draw("same");

  sprintf( hname, "rebinPLOT_%ithresRMS.png",thresRMS);
  c1->SaveAs( hname);
  printf("done plotting rebins 1, 2, 4, 8\n");

  return 1;
}

int plotRebin2( int iCal){
  
  char cname[64], finname[64];

  int rb = 1;
  sprintf( finname, "noxtaldata/1912rb%i_16.root", rb);
  sprintf( finname, "noxtaldata/test_1f_%i_16.root", rb);
  file1 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir1;
  file1->GetObject( dname, dir1);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_%i_%i", iCal, iCal);
  sprintf( hname, "qHist1D_ped_%i_%i", iCal, iCal);
  dir1->GetObject( hname, htmp1);
  sprintf( hname, "1912 rebin %i Calo %i", rb, iCal);
  htmp1->SetName(hname);
  htmp1->SetTitle(hname);
  htmp1->SetLineColor(kViolet);
  htmp1->SetLineWidth(rb);
  htmp1->SetLineWidth(2);
  htmp1->Rebin(8./rb);

  rb = 2;
  sprintf( finname, "noxtaldata/1912rb%i_16.root", rb);
  sprintf( finname, "noxtaldata/test_1f_%i_16.root", rb);
  file2 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir2;
  file2->GetObject( dname, dir2);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_%i_%i", iCal, iCal);
  sprintf( hname, "qHist1D_ped_%i_%i", iCal, iCal);
  dir2->GetObject( hname, htmp2);
  sprintf( hname, "1912 rebin %i Calo %i", rb, iCal);
  htmp2->SetName(hname);
  htmp2->SetTitle(hname);
  htmp2->SetLineColor(kOrange);
  htmp2->SetLineWidth(rb);
  htmp2->SetLineWidth(2);
  htmp2->Rebin(8./rb);

  rb = 4;
  sprintf( finname, "noxtaldata/1912rb%i_16.root", rb);
  sprintf( finname, "noxtaldata/test_1f_%i_16.root", rb);
  file3 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir3;
  file3->GetObject( dname, dir3);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_%i_%i", iCal, iCal);
  sprintf( hname, "qHist1D_ped_%i_%i", iCal, iCal);
  dir3->GetObject( hname, htmp3);
  sprintf( hname, "1912 rebin %i Calo %i", rb, iCal);
  htmp3->SetName(hname);
  htmp3->SetTitle(hname);
  htmp3->SetLineColor(kRed);
  htmp3->SetLineWidth(rb);
  htmp3->SetLineWidth(2);
  htmp3->Rebin(8./rb);

  rb = 8;
  sprintf( finname, "noxtaldata/1912rb%i_16.root", rb);
  sprintf( finname, "noxtaldata/test_1f_%i_16.root", rb);
  file4 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir4;
  file4->GetObject( dname, dir4);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "qHist1D_sig_%i_%i", iCal, iCal);
  sprintf( hname, "qHist1D_ped_%i_%i", iCal, iCal);
  dir4->GetObject( hname, htmp4);
  sprintf( hname, "1912 rebin %i Calo %i", rb, iCal);
  htmp4->SetName(hname);
  htmp4->SetTitle(hname);
  htmp4->SetLineColor(kCyan);
  htmp4->SetLineWidth(rb);
  htmp4->SetLineWidth(2);
  htmp4->Rebin(8./rb);

  htmp4->Draw("");
  htmp3->Draw("same");
  htmp2->Draw("same");
  htmp1->Draw("same");

  sprintf( hname, "rebinPLOT.png");
  c1->SaveAs( hname);
  printf("done plotting rebins 1, 2, 4, 8\n");

  return 1;
}

TH1D  *hfr1 = new TH1D("hfr1","hfr1",5000,0,5000);
TH1D  *hfr2 = new TH1D("hfr2","hfr2",2500,0,5000);
TH1D  *hfr3 = new TH1D("hfr3","hfr3",2500,0,5000);

int testFastRot(){
  
  double f = 150./149.2;
  double sumX, lastsumX = 0.0, deltasumX, lastdeltasumX;
  for (int iBin = 1; iBin <= 5000; iBin++) {

    double xBinHi = f*iBin;
    
    if ( ((int)xBinHi)%2 ) {
      sumX = ((int)xBinHi)/2+1;
    } else {
      sumX = ((int)(xBinHi-1))/2+1+(xBinHi-int(xBinHi));
    }
    deltasumX = sumX - lastsumX;

    if (iBin%2) {
      printf("xBinHi %f, sumX %f, lastsumX %f, deltasumX %f\n", xBinHi, sumX, lastsumX, deltasumX);
      hfr1->Fill(iBin,deltasumX);
    } else {
      printf("xBinHi %f, sumX %f, lastsumX %f, deltasumX %f, deltasumX+lastdeltasumX %f\n", xBinHi, sumX, lastsumX, deltasumX, deltasumX+lastdeltasumX );
      hfr1->Fill(iBin,deltasumX);
      hfr2->Fill(iBin,deltasumX+lastdeltasumX);
      hfr3->Fill(iBin,deltasumX-lastdeltasumX);
    }
    lastdeltasumX = deltasumX;
    lastsumX = sumX;
  }
  hfr1->SetLineColor(kRed);
  hfr2->SetLineColor(kMagenta);
  //hfr1->Draw("hist");
  hfr3->Draw("hist");
  hfr2->Draw("samehist");

  return 1;
}

TGraphErrors *gcyc;
TF1 *fastf, *ftf[3]; // array is for piece wise fit function
TH1D *fasth;
TH1D *hhalfspline[24];
TH1D *hhalf[24], *hfull[24];
TH1D *hhalfmean[24], *hfullmean[24];

int plotFastRot3(){
    
  int nCal = 24;

  
  for (int iCal = 1; iCal <= nCal; iCal++) {
    
    icalo = iCal;
    fitCalo( "full", 0, 52000, 200000);

    hhalf[iCal] = (TH1D*)hres->Clone();
    printf("cloned histogram\n");

    sprintf( hname, "hhalf%i", iCal);
    hhalf[iCal]->SetName(hname);
    hhalf[iCal]->SetTitle(hname);
    printf("set histogram name, title \n");

    hhalf[iCal]->GetXaxis()->SetRangeUser(52000,120000);
    hhalf[iCal]->SetLineColor(kBlue);
    hhalf[iCal]->SetMarkerColor(kBlue);
    hhalf[iCal]->SetMaximum(0.05);
    hhalf[iCal]->SetMinimum(-0.05);
  }

  c1->Clear();
  c1->Divide(4,6);
  for (int iCal = 1; iCal <= nCal; iCal++) {

    c1->cd(iCal);
    printf("draw histogram\n");
    hhalf[iCal]->Draw("p");
  }

  return 0;
}

double t[24*3168], r[24*3168], et[24*3168], er[24*3168];
double ti[3168], ri[3168];
TGraph *gcyci[24],*gfrphase;
double dt[24] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
TH1D *hcyc[24];
double calo_index[24], calo_index_err[24] = {0.0},  phase_fr_err[24] = {0.0};
//double phase_fr[24] = { 3.889434, 3.883284, 3.878959, 3.823770, 3.852049, 3.687728, 3.713019, 3.636553, 3.595910, 3.453713, 3.383658, 3.325935, 3.254716, 3.190543, 3.126193, 3.159731, 3.561239, 3.523594, 3.431466, 3.321373, 3.298190, 3.200291, 3.192550, 3.074217};
//double phase_fr[24] = {  0.886423, 0.949525, 1.012864, 1.022316, 1.109354, 1.006378, 1.090126, 1.077057, 1.100187, 1.024070, 1.025143, 1.036220, 1.032247, 1.036091, 1.037612, 1.133586, 1.588337, 1.612988, 1.585450, 1.542004, 1.585416, 1.555550, 1.617694, 1.568187};

//double phase_fr[24] = {  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1 }; // semi-empirical phase correction


int npw = 4;
double pwoverlap = 0.0;
double pwlo[4] = { 30000., 50000., 60000., 75000.};
double pwhi[4] = { 50000.+pwoverlap, 60000.+pwoverlap, 75000.+pwoverlap, 100000.+pwoverlap};

TF1 *ffrpw;

Double_t fnfrpw(Double_t *x, Double_t *par);

Double_t fnfrpw(Double_t *x, Double_t *par){
  
  Double_t xx = x[0], f = 0.0;
 
  // par[0] cyclotron angular freq - common term
  // par[1] cyclotron phase - common term

  // pw1
  // par[2] - amp1
  // par[3] - tau1
  // par[4] - amp2
  // par[5] - tau2
  // par[6] - abeat
  // par[7] beat freq
  // par[8] beat phase

  // par[30] cyclotron ang freq linear change - common term

  // + 7 for pw2 + 7 pw3 = 23 total pars

  int iparoffset = 0;
  if ( xx >= pwlo[0] && xx < pwhi[0] ) iparoffset = 0; 
  if ( xx >= pwlo[1] && xx < pwhi[1] ) iparoffset = 7; 
  if ( xx >= pwlo[2] && xx < pwhi[2] ) iparoffset = 14; 
  if ( xx >= pwlo[3] && xx < pwhi[3] ) iparoffset = 21;
  double phi_cyc = par[1];

  if (doFRfitWithInts) {

    // point at xx is integral over data from xx - 75ns/2  to xx - 75ns/2
    // 
    int np = 10; // points integrated
    double Tb = 75.; // ns bin width
    double xxp; // point coordinate
    double df = 0.0; // function incremental contribution  
    for (int ip = 0; ip < np; ip++){

      xxp = xx - Tb/2. + Tb/2./np + ip*Tb/np;
      // envelope
      df = abs(par[2+iparoffset])*exp( -(xxp-pwlo[0+iparoffset/7])/par[3+iparoffset] ) + abs(par[4+iparoffset])*exp( -(xxp-pwlo[0+iparoffset/7])/par[5+iparoffset] );
      // Tcyc * Tbeat
      df *= sin( par[0]*(xxp-pwlo[0]) + phi_cyc ) * ( 1.0+par[6+iparoffset]*cos( par[7+iparoffset]*(xxp-pwlo[0])+par[8+iparoffset] ) ); 
      f += df/np;
      if (dbg) printf(" ip %i, xxp %f, df %f, f %f\n", ip, xxp, df, f);
    }

  } else { 

    // envelope
    f = abs(par[2+iparoffset])*exp( -(xx-pwlo[0+iparoffset/7])/par[3+iparoffset] ) + abs(par[4+iparoffset])*exp( -(xx-pwlo[0+iparoffset/7])/par[5+iparoffset] );
    // Tcyc * Tbeat
    f *= sin( par[0]*(xx-pwlo[0]) + phi_cyc ) * ( 1.0+par[6+iparoffset]*cos( par[7+iparoffset]*(xx-pwlo[0])+par[8+iparoffset] ) ); 
  }

  return f;
}

bool fixamp1[4] = {1,1,1,1};
bool fixtau1[4] = {1,1,1,1};
bool fixamp2[4] = {1,1,1,1};
bool fixtau2[4] = {1,1,1,1};
bool fixasym[4] = {1,1,1,1};
bool fixomgb[4] = {1,1,1,1};
bool fixphib[4] = {1,1,1,1};
bool fixdcyc = 1;
bool fixomgcyc = 0;
bool fixphicyc = 0;

double omega_cyc = 5.19222e-2/rawBinToNs;
double delta_phi_cyc = 0.0;
double phi_cyc = 1.25700e+01;
double Ampl1[4] = { 0.08143, 0.021866, 0.021866, 0.02};
double Tau1[4] = { 1.5639e4*rawBinToNs, 6.41e9*rawBinToNs, 6.41e9*rawBinToNs, 4.09e9*rawBinToNs};
double Ampl2[4] = { 0.0, 0.0, 0.0, 0.0};
double Tau2[4] = { 7.868e8*rawBinToNs, 1.e4*rawBinToNs, 1.e4*rawBinToNs, 1.e4*rawBinToNs};
double Asym[4] = { -5.7426e-1, 1.01, 1.01, 0.5};
double omega_beat[4] = { 2.25151e-04*rawBinToNs, 2.345e-4*rawBinToNs, 2.345e-4*rawBinToNs, 2.34e-4*rawBinToNs};
double phi_beat[4] = { 1.04527, -2.580, -2.580, -5.42};
int nCalfr = 24;

TF1 *fnoint,*fysint;

bool makeHalfFullHistograms, plotHalfFullHistograms, storeHalfFullHistograms;
TSpline5 *fastRotationSpline;


int plotFastRot2( double fr_start = 52000, double fr_end = 110000, double tcyc = 149.3, bool do_wiggle_fits = 0){

  /* 
     change three things when changing dataset
  1) double phase_fr[24] ={};
  2) file0 = TFile::Open( "hcyc-9day-vbo-pw.root","UPDATE");
  3) sprintf(fname,"noxtaldata/9day-DQC-infillgain-w4t15g11.root");
     sprintf(dname,"QFillByFillAnalyzer");
     sprintf(dataname,"9day");
  */

  // phase shift for 16 ns, 16./149.3*2.*TMath::Pi() =  0.673349

  double phase_fr[24] = {  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0.673, 0.673, 0.673, 0.673,  0.673, 0.673, 0.673, 0.673  }; // all zeros

  //plotFastRot2(50000,110000,149.3) 9day
  //double phase_fr[24] = {  0.000000, 0.032979, 0.045617, 0.033968, 0.065034, -0.057735, -0.000979, 0.016345, 0.075998, 0.045700, 0.090816, 0.116815, 0.117084, 0.116388, 0.059479, 0.092446, 0.192441, 0.193030, 0.197172, 0.193801, 0.302960, 0.329343, 0.394514, 0.336325}; // 9day  semi-empirical phase correction

  //plotFastRot2(50000,110000,149.3) 60hr
  //double phase_fr[24] = { 0.000000, 0.095321, 0.184176, 0.214801, 0.313316, 0.199856, 0.260110, 0.224931, 0.221038, 0.132630, 0.125861, 0.136132, 0.146581, 0.183248, 0.204459, 0.322491, 0.763460, 0.760060, 0.705913, 0.659232, 0.694047, 0.666292, 0.746030, 0.711188}; 

// plotFastRot2(58000,74000,149.3)
//double phase_fr[24] =  {1.421461, 1.443693, 1.459756, 1.472414, 1.528128, 1.405293, 1.473588, 1.478229, 1.525046, 1.473727, 1.505693, 1.526080, 1.533982, 1.525539, 1.480863, 1.532563, 1.666208, 1.674485, 1.682119, 1.659883, 1.726757, 1.740266, 1.797277, 1.740886};  // 9day  semi-empirical phase correction
    
  int n = nCalfr*3168, ni = 3168;  

  if (!do_wiggle_fits){
    // if already generated a correction file with hcyc's from wiggle fits
    file1 = TFile::Open( fnamefastrotation, "UPDATE");
    for (int iCal = 1; iCal <= nCalfr; iCal++) {
      sprintf(hname,"hcyc%i",iCal);
      file1->GetObject( hname, hcyc[iCal-1]);
    }
  }

  // initiate fast rotation analysis
  fastRotationCorrection = 0;
  fastRotationAnalysis = 1;
  // allow analyzer to set
  //sprintf(fitname,"vbo"); 
  //sprintf(fitname,"13par"); 
  //fix_pedestaldrift_amp = 1;
  //set_pedestaldrift_amp = 0;
  //fix_muonloss_amp = 1;
  //set_muonloss_amp = 0;
  //fix_gain_amp = 1;
  //set_gain_amp = 0;

  for (int iCal = 1; iCal <= nCalfr; iCal++) {
    
    if (do_wiggle_fits){
      // if not already generated a correction file of hcyc histograms of residuals from wiggle fits in 
      icalo = iCal;
      fitCalo( "singlecalo", icalo, fr_start, fr_end); // do fits to each calo and get the residuals from fast rotation
      hcyc[iCal-1] = (TH1D*)hres->Clone();
      sprintf(hname,"hcyc%i",iCal);
      hcyc[iCal-1]->SetName(hname);
      hcyc[iCal-1]->SetTitle(hname);
      printf("hcyc iCal %i, bins %i\n", iCal, hcyc[iCal-1]->GetNbinsX());
    }
    
    for (int iBin = 1; iBin <= hcyc[iCal-1]->GetNbinsX(); iBin++) { // superimpose the fast rotation of each calo in time-shifted graph  
      t[nCalfr*(iBin-1)+iCal-1] = hcyc[iCal-1]->GetBinCenter(iBin) + tcyc/24.*((double)(iCal-1)) + tcyc * phase_fr[iCal-1]/(2.*TMath::Pi());
      r[nCalfr*(iBin-1)+iCal-1] = hcyc[iCal-1]->GetBinContent(iBin); // 24 calo's for given iBin are distributed over 2*iBin range
      et[nCalfr*(iBin-1)+iCal-1] = 0.0; // errors
      er[nCalfr*(iBin-1)+iCal-1] = 0.0; // errors
      ti[iBin-1] = hcyc[iCal-1]->GetBinCenter(iBin) + tcyc/24.*((double)(iCal-1)) + tcyc*phase_fr[iCal-1] /(2.*TMath::Pi()); // time coordinate
      ri[iBin-1] = hcyc[iCal-1]->GetBinContent(iBin); // fit residual coordinate
    } 
    
    // individual calo graphs of residuals the time units of gcyci x-axis are 75ns bins
    gcyci[iCal-1] = new TGraph( ni, ti, ri); 
    gcyci[iCal-1]->SetMaximum(0.1);
    gcyci[iCal-1]->SetMinimum(-0.1);
    gcyci[iCal-1]->SetMarkerStyle(20);
    gcyci[iCal-1]->SetMarkerSize(0.5);
    sprintf( hname, "gcyci%i", iCal);
    gcyci[iCal-1]->SetName(hname);
    gcyci[iCal-1]->SetTitle(hname);

    if ((iCal-1)/4 == 0) gcyci[iCal-1]->SetMarkerColor(kBlue+iCal/4e9); // 1-4
    if ((iCal-1)/4 == 1) gcyci[iCal-1]->SetMarkerColor(kMagenta+iCal/4e9);  // 5-8
    if ((iCal-1)/4 == 2) gcyci[iCal-1]->SetMarkerColor(kGreen+iCal/4e9); // 9-12
    if ((iCal-1)/4 == 3) gcyci[iCal-1]->SetMarkerColor(kRed+iCal/4e9); // 13-16
    if ((iCal-1)/4 == 4) gcyci[iCal-1]->SetMarkerColor(kBlack+iCal/4e9); // 17-20
    if ((iCal-1)/4 == 5) gcyci[iCal-1]->SetMarkerColor(kCyan+iCal/4e9); // 21-24
    if ((iCal-1)/4 == 0) gcyci[iCal-1]->SetLineColor(kBlack+iCal/4e9);
    if ((iCal-1)/4 == 1) gcyci[iCal-1]->SetLineColor(kBlue+iCal/4e9);
    if ((iCal-1)/4 == 2) gcyci[iCal-1]->SetLineColor(kGreen+iCal/4e9);
    if ((iCal-1)/4 == 3) gcyci[iCal-1]->SetLineColor(kMagenta+iCal/4e9);
    if ((iCal-1)/4 == 4) gcyci[iCal-1]->SetLineColor(kRed+iCal/4e9);
    if ((iCal-1)/4 == 5) gcyci[iCal-1]->SetLineColor(kGray+iCal/4e9);
    
    printf("draw gcyc[%i] graph\n",iCal-1);
    c1->Clear();
    gcyci[iCal-1]->Draw("ap");
  }

  printf("fill gcyc graph\n");

  // the time units of gcyc x-axis are 75ns bins
  gcyc = new TGraphErrors( n, t, r, et, er);
  gcyc->SetMaximum(1.0);
  gcyc->SetMinimum(-1.0);
  gcyc->SetMarkerStyle(20);
  gcyc->SetMarkerSize(0.5);
  gcyc->SetName("gcyc");
  gcyc->SetTitle("gcyc");
  
  c1->Clear();
  gcyc->Draw("AP");
  
  if (do_wiggle_fits){
    // if not already generated a correction file of hcyc histograms of residuals from wiggle fits in 
    file1 = TFile::Open( fnamefastrotation, "UPDATE");
    for (int iCal = 1; iCal <= nCalfr; iCal++) hcyc[iCal-1]->Write();
    printf("wrote hcyc histograms to file");
    gcyc->Write();
    printf("wrote gcyc");
    file1->Close();
    return 0;
  }
  
  /* old?
  // the time units of ftf x-axis are 75ns bins
  double N_1[3] = { 1.5, 1.5, 1.5},  tau_1[3] = { 100.*60, 100.*60, 100.*60}, N_2[3] = { 0.1, 0.1, 0.1},  tau_2[3] = { 500.*60, 500.*60, 500.*60};
  double omega_f[3] = { 2.*TMath::Pi()/120., 2.*TMath::Pi()/120., 2.*TMath::Pi()/120.},  omega_s[3] = { 2.*TMath::Pi()/500./60., 2.*TMath::Pi()/500./60., 2.*TMath::Pi()/500./60.};
  double phi_f[3] = { 0.0, 0.0, 0.0}, phi_s[3] = { 5.0, 5.0, 5.0}, A_1[3] = { 0.45, 0.45, 0.45}, t0[3] = { 30000., 30000., 30000.};
  */

  /*
  npw = 1; // for npw=1
  pwlo[0] = 50000.; // for npw=1
  pwhi[0] = 110000.; // for npw=1  
  */

  // adjust errors in overlap regions
  /* 
  for (int ip = 1; ip <= gcyc->GetN(); ip++){
    Double_t xi, yi;
    gcyc->GetPoint( ip, xi, yi);
    if (xi >= pwlo[1] && xi < pwlo[1]+pwoverlap) gcyc->SetPointError( ip, 0.0, 1.0e-2/3.);
    if (xi >= pwlo[2] && xi < pwlo[2]+pwoverlap) gcyc->SetPointError( ip, 0.0, 1.0e-2/3.);
  }
  */

  // fit function time in ns
  char pwname[64];
  
  // single, coupled piece-wise fit
  ffrpw = new TF1( "ffrpw", fnfrpw, pwlo[0], pwhi[3], 31);
  ffrpw->SetNpx(10000);
  ffrpw->SetParameters( omega_cyc, phi_cyc, Ampl1[0], Tau1[0], Ampl2[0], Tau2[0], Asym[0], omega_beat[0], phi_beat[0]);
  ffrpw->SetParNames( "omega_cyc", "phi_cyc","A1_0", "tau1_0", "A2_0", "tau2_0", "asym_0", "omega_beat0", "phi_beat0");
  
  ffrpw->SetParameter(30, delta_phi_cyc);
  ffrpw->SetParName(30, "delta_phi_cyc");

  ffrpw->SetParameter(30,delta_phi_cyc);
  ffrpw->SetParameter(2, Ampl1[0]);
  ffrpw->SetParameter(3, Tau1[0]);
  ffrpw->SetParameter(4, Ampl2[0]);
  ffrpw->SetParameter(5, Tau2[0]);
  ffrpw->SetParameter(6, Asym[0]);
  ffrpw->SetParameter(7, omega_beat[0]);
  ffrpw->SetParameter(8, phi_beat[0]);

  if (fixdcyc) ffrpw->FixParameter(30,delta_phi_cyc);
  if (fixomgcyc) ffrpw->FixParameter(0,ffrpw->GetParameter(0));
  if (fixphicyc) ffrpw->FixParameter(1,ffrpw->GetParameter(1));

  if (fixamp1[0]) ffrpw->FixParameter(2, Ampl1[0]);
  if (fixtau1[0]) ffrpw->FixParameter(3, Tau1[0]);
  if (fixamp2[0]) ffrpw->FixParameter(4, Ampl2[0]);
  if (fixtau2[0]) ffrpw->FixParameter(5, Tau2[0]);
  if (fixasym[0]) ffrpw->FixParameter(6, Asym[0]);
  if (fixomgb[0]) ffrpw->FixParameter(7, omega_beat[0]);
  if (fixphib[0]) ffrpw->FixParameter(8, phi_beat[0]);

  if (!fixamp1[0]) ffrpw->ReleaseParameter(2);
  if (!fixtau1[0]) ffrpw->ReleaseParameter(3);
  if (!fixamp2[0]) ffrpw->ReleaseParameter(4);
  if (!fixtau2[0]) ffrpw->ReleaseParameter(5);
  if (!fixasym[0]) ffrpw->ReleaseParameter(6);
  if (!fixomgb[0]) ffrpw->ReleaseParameter(7);
  if (!fixphib[0]) ffrpw->ReleaseParameter(8);
  
  // set pw2 parameters
  ffrpw->SetParameter(9, Ampl1[1]);
  ffrpw->SetParName(9, "A1_1");
  ffrpw->SetParameter(10, Tau1[1]);
  ffrpw->SetParName(10, "Tau1_1");
  ffrpw->SetParameter(11, Ampl2[1]);
  ffrpw->SetParName(11, "A2_1");
  ffrpw->SetParameter(12, Tau2[1]);
  ffrpw->SetParName(12, "Tau2_1");
  ffrpw->SetParameter(13, Asym[1]);
  ffrpw->SetParName(13, "asym_1");
  ffrpw->SetParameter(14, omega_beat[1]);
  ffrpw->SetParName(14, "omega_beat1");
  ffrpw->SetParameter(15, phi_beat[1]);
  ffrpw->SetParName(15, "phi_beat1");
  
  if (fixamp1[1]) ffrpw->FixParameter(9, Ampl1[1]);
  if (fixtau1[1]) ffrpw->FixParameter(10, Tau1[1]);
  if (fixamp2[1]) ffrpw->FixParameter(11, Ampl2[1]);
  if (fixtau2[1]) ffrpw->FixParameter(12, Tau2[1]);
  if (fixasym[1]) ffrpw->FixParameter(13, Asym[1]);
  if (fixomgb[1]) ffrpw->FixParameter(14, omega_beat[1]);
  if (fixphib[1]) ffrpw->FixParameter(15, phi_beat[1]);

  if (!fixamp1[1]) ffrpw->ReleaseParameter(9);
  if (!fixtau1[1]) ffrpw->ReleaseParameter(10);
  if (!fixamp2[1]) ffrpw->ReleaseParameter(11);
  if (!fixtau2[1]) ffrpw->ReleaseParameter(12);
  if (!fixasym[1]) ffrpw->ReleaseParameter(13);
  if (!fixomgb[1]) ffrpw->ReleaseParameter(14);
  if (!fixphib[1]) ffrpw->ReleaseParameter(15);

  // set pw3 parameters
  ffrpw->SetParameter(16, Ampl1[2]);
  ffrpw->SetParName(16, "A1_2");
  ffrpw->SetParameter(17, Tau1[2]);
  ffrpw->SetParName(17, "Tau1_2");
  ffrpw->SetParameter(18, Ampl2[2]);
  ffrpw->SetParName(18, "A2_2");
  ffrpw->SetParameter(19, Tau2[2]);
  ffrpw->SetParName(19, "Tau2_2");
  ffrpw->SetParameter(20, Asym[2]);
  ffrpw->SetParName(20, "asym_2");
  ffrpw->SetParameter(21, omega_beat[2]);
  ffrpw->SetParName(21, "omega_beat2");
  ffrpw->SetParameter(22, phi_beat[2]);
  ffrpw->SetParName(22, "phi_beat2");

  // fix pw3 if out of range of fit
  if (fixamp1[2]) ffrpw->FixParameter(16, Ampl1[2]);
  if (fixtau1[2]) ffrpw->FixParameter(17, Tau1[2]);
  if (fixamp2[2]) ffrpw->FixParameter(18, Ampl2[2]);
  if (fixtau2[2]) ffrpw->FixParameter(19, Tau2[2]);
  if (fixasym[2]) ffrpw->FixParameter(20, Asym[2]);
  if (fixomgb[2]) ffrpw->FixParameter(21, omega_beat[2]);
  if (fixphib[2]) ffrpw->FixParameter(22, phi_beat[2]);

  if (!fixamp1[2]) ffrpw->ReleaseParameter(16);
  if (!fixtau1[2]) ffrpw->ReleaseParameter(17);
  if (!fixamp2[2]) ffrpw->ReleaseParameter(18);
  if (!fixtau2[2]) ffrpw->ReleaseParameter(19);
  if (!fixasym[2]) ffrpw->ReleaseParameter(20);
  if (!fixomgb[2]) ffrpw->ReleaseParameter(21);
  if (!fixphib[2]) ffrpw->ReleaseParameter(22);

  // set pw4 parameters
  ffrpw->SetParameter(23, Ampl1[3]);
  ffrpw->SetParName(23, "A1_3");
  ffrpw->SetParameter(24, Tau1[3]);
  ffrpw->SetParName(24, "Tau1_3");
  ffrpw->SetParameter(25, Ampl2[3]);
  ffrpw->SetParName(25, "A2_3");
  ffrpw->SetParameter(26, Tau2[3]);
  ffrpw->SetParName(26, "Tau2_3");
  ffrpw->SetParameter(27, Asym[3]);
  ffrpw->SetParName(27, "asym_3");
  ffrpw->SetParameter(28, omega_beat[3]);
  ffrpw->SetParName(28, "omega_beat3");
  ffrpw->SetParameter(29, phi_beat[3]);
  ffrpw->SetParName(29, "phi_beat3");

  // fix pw4 if out of range of fit
  if (fixamp1[3]) ffrpw->FixParameter(23, Ampl1[3]);
  if (fixtau1[3]) ffrpw->FixParameter(24, Tau1[3]);
  if (fixamp2[3]) ffrpw->FixParameter(25, Ampl2[3]);
  if (fixtau2[3]) ffrpw->FixParameter(26, Tau2[3]);
  if (fixasym[3]) ffrpw->FixParameter(27, Asym[3]);
  if (fixomgb[3]) ffrpw->FixParameter(28, omega_beat[3]);
  if (fixphib[3]) ffrpw->FixParameter(29, phi_beat[3]);

  if (!fixamp1[3]) ffrpw->ReleaseParameter(23);
  if (!fixtau1[3]) ffrpw->ReleaseParameter(24);
  if (!fixamp2[3]) ffrpw->ReleaseParameter(25);
  if (!fixtau2[3]) ffrpw->ReleaseParameter(26);
  if (!fixasym[3]) ffrpw->ReleaseParameter(27);
  if (!fixomgb[3]) ffrpw->ReleaseParameter(28);
  if (!fixphib[3]) ffrpw->ReleaseParameter(29);

  printf("fit with range %f, %f and bin-integral %i\n", pwlo[0], pwhi[0], doFRfitWithInts);
  gcyc->Fit( ffrpw, "NR", "",  pwlo[0], pwhi[0]);
  ffrpw->Draw("hist");

  if (!makeHalfFullHistograms) return 0;

  /*
  // plotting piecewise fit for individual calos
  gcyc->GetYaxis()->SetRangeUser(-0.1,0.1);
  TGraphErrors* gcyc1 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc2 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc3 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc4 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc5 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc6 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc7 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc8 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc9 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc10 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc11 = (TGraphErrors*) gcyc->Clone();
  TGraphErrors* gcyc12 = (TGraphErrors*) gcyc->Clone();
  c1->Clear();
  c1->Divide(3,4);
  c1->cd(1);
  gcyc1->GetXaxis()->SetRangeUser(46000,51000);
  gcyc1->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(2);
  gcyc2->GetXaxis()->SetRangeUser(51000,56000);
  gcyc2->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(3);
  gcyc3->GetXaxis()->SetRangeUser(56000,61000);
  gcyc3->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(4);
  gcyc4->GetXaxis()->SetRangeUser(61000,66000);
  gcyc4->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(5);
  gcyc5->GetXaxis()->SetRangeUser(66000,71000);
  gcyc5->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(6);
  gcyc6->GetXaxis()->SetRangeUser(71000,76000);
  gcyc6->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(7);
  gcyc7->GetXaxis()->SetRangeUser(76000,81000);
  gcyc7->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(8);
  gcyc8->GetXaxis()->SetRangeUser(81000,86000);
  gcyc8->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(9);
  gcyc9->GetXaxis()->SetRangeUser(86000,91000);
  gcyc9->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(10);
  gcyc10->GetXaxis()->SetRangeUser(91000,96000);
  gcyc10->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(11);
  gcyc11->GetXaxis()->SetRangeUser(96000,101000);
  gcyc11->Draw("ap");
  ffrpw->Draw("histsame");
  c1->cd(12);
  gcyc12->GetXaxis()->SetRangeUser(101000,106000);
  gcyc12->Draw("ap");
  ffrpw->Draw("histsame");
  */

#if 0 // if adjusting phases of individual calo's

  int ipw = 0; // use single fit function version
  double bestphasefr = ftf[0]->GetParameter(5);

  ftf[ipw]->FixParameter( 0, ftf[ipw]->GetParameter(0));
  ftf[ipw]->FixParameter( 1, ftf[ipw]->GetParameter(1));
  ftf[ipw]->FixParameter( 2, ftf[ipw]->GetParameter(2));
  ftf[ipw]->FixParameter( 3, ftf[ipw]->GetParameter(3));
  ftf[ipw]->FixParameter( 4, ftf[ipw]->GetParameter(4));
  ftf[ipw]->FixParameter( 5, ftf[ipw]->GetParameter(5));
 
  // adjust fast rotation phase  
  // this code is old and would need new code for new function
  for (int iCal = 1; iCal <= nCalfr; iCal++) {
    ftf[ipw]->SetParameter(5,bestphasefr);
    ftf[ipw]->ReleaseParameter(5);
    ftf[ipw]->SetParLimits( 5, 0.0, 2.*TMath::Pi());
    gcyci[iCal-1]->Fit( ftf[ipw], "", "RIE", pwlo[ipw], pwhi[ipw]);
    gcyci[iCal-1]->Draw("ap");
    phase_fr[iCal-1] = ftf[ipw]->GetParameter(5);
    phase_fr_err[iCal-1] = ftf[ipw]->GetParError(5);
    calo_index[iCal-1] = iCal;
    for (int iBin = 1; iBin <= hcyc[0]->GetNbinsX(); iBin++ ){  
      double time = hcyc[iCal-1]->GetBinCenter(iBin) + tcyc/24.*((double)(iCal-1)) + tcyc*phase_fr[iCal-1]/(2.*TMath::Pi()); 
      hhalf[iCal-1]->SetBinContent( iBin, ftf[ipw]->Eval(time));
    }
  }
  
  // make plot of fast rotation phase for individual calorimeters
  gfrphase = new TGraphErrors( nCalfr, calo_index, phase_fr, calo_index_err, phase_fr_err); 
  gfrphase->SetMarkerStyle(20);
  gfrphase->SetMarkerSize(0.5);
  gfrphase->SetMarkerColor(kRed);
  gfrphase->SetLineColor(kRed);
  gfrphase->SetName("gfrphase");
  gfrphase->SetTitle("fast rotation phase versus calo index after 2pi/24 shift applied");
  gfrphase->GetXaxis()->SetTitle("calorimeter index");
  gfrphase->GetYaxis()->SetTitle("phase (radians)");
  gfrphase->Draw("ap");
  
  printf("{ ");
  for (int iCal = 1; iCal <= nCalfr; iCal++) printf(" %f,",fmod( phase_fr[iCal-1] - phase_fr[0], 2.*TMath::Pi()));
  printf("} \n");

#endif

  // make correction histograms
  for (int iCal = 1; iCal <= nCalfr; iCal++){ // define 60-tick, 75ns binned histograms
    
    sprintf( hname, "hhalf%i", iCal);
    hhalf[iCal-1] = (TH1D*)hcyc[0]->Clone();
    hhalf[iCal-1]->Reset();
    hhalf[iCal-1]->SetName(hname);
    hhalf[iCal-1]->SetTitle(hname); 
    hhalf[iCal-1]->SetMarkerSize(1.0); 
    hhalf[iCal-1]->GetXaxis()->SetRangeUser(fr_start,fr_end); 
    hhalf[iCal-1]->SetLineColor(kBlue);
    hhalf[iCal-1]->SetMarkerColor(kBlue);

    sprintf( hname, "hhalfmean%i", iCal);
    hhalfmean[iCal-1] = (TH1D*)hcyc[0]->Clone();
    hhalfmean[iCal-1]->Reset();
    hhalfmean[iCal-1]->SetName(hname);
    hhalfmean[iCal-1]->SetTitle(hname); 
    hhalfmean[iCal-1]->SetMarkerSize(1.0); 
    hhalfmean[iCal-1]->GetXaxis()->SetRangeUser(fr_start,fr_end); 
    hhalfmean[iCal-1]->SetLineColor(kMagenta);
    hhalfmean[iCal-1]->SetMarkerColor(kMagenta);
  }

  doFRfitWithInts = 0; // no bin integral so ffrpw->Eval(time) gives unconvoluted fast rotation correction

  for (int iCal = 1; iCal <= nCalfr; iCal++) {

    for (int iBin = 1; iBin <= hcyc[0]->GetNbinsX(); iBin++ ) { 

      double time = hcyc[iCal-1]->GetBinCenter(iBin) + tcyc/24.*((double)(iCal-1)) + tcyc*phase_fr[iCal-1]/(2.*TMath::Pi()); 

      int np = 10; // points integrated
      double Tb = 75.; // ns bin width
      double dtime; // point coordinate
      double sum_yt = 0.0, sum_y = 0.0;

      for (int ip = 0; ip < np; ip++){

      double dtime = time - Tb/2. + Tb/2./np + ip*Tb/np;
      double dyfr = 1. + ffrpw->Eval(dtime); // +1 as correction is histogram + 1 multiplier 
      sum_yt += dyfr*dtime;
      sum_y += dyfr;
      if (dbg) printf(" ip %i, xxp %f, df %f, sum y %f, sum y*t %f\n", ip, dtime, dyfr, sum_y, sum_yt);
    }
      double dtmn = sum_yt/sum_y - time;

      if ( time - tcyc/24.*((double)(iCal-1)) > 30500 && time -tcyc/24.*((double)(iCal-1)) < 69500) {
	hhalf[iCal-1]->SetBinContent( iBin, sum_y/np - 1. );
	hhalfmean[iCal-1]->SetBinContent( iBin, dtmn );
      } else {
	hhalf[iCal-1]->SetBinContent( iBin, 0.0 );   
	hhalfmean[iCal-1]->SetBinContent( iBin, 0.0 );   
      }
    }  
  }
  printf("made hhalf correction histograms\n");

  c1->Clear();
  c1->Divide(4,6);
  
  for (int iCal = 1; iCal <= nCalfr; iCal++){
    
    hfull[iCal-1] = (TH1D*)hhalf[iCal-1]->Clone();
    sprintf( hname, "hfull%i", iCal);
    hfull[iCal-1]->SetName(hname);
    hfull[iCal-1]->SetTitle(hname); 
    hfull[iCal-1]->SetLineColor(kRed);
    hfull[iCal-1]->SetLineWidth(2.0);
    hfull[iCal-1]->SetMarkerColor(kRed);
    hfull[iCal-1]->SetMarkerStyle(20);
    hfull[iCal-1]->SetMarkerStyle(20);
    hfull[iCal-1]->SetMarkerSize(0.5);
    hfull[iCal-1]->SetMarkerSize(0.5);
    hfull[iCal-1]->Rebin(2);
    
    hfullmean[iCal-1] = (TH1D*)hhalfmean[iCal-1]->Clone();
    sprintf( hname, "hfullmean%i", iCal);
    hfullmean[iCal-1]->SetName(hname);
    hfullmean[iCal-1]->SetTitle(hname); 
    hfullmean[iCal-1]->SetLineColor(kOrange);
    hfullmean[iCal-1]->SetLineWidth(2.0);
    hfullmean[iCal-1]->SetMarkerColor(kOrange);
    hfullmean[iCal-1]->SetMarkerStyle(20);
    hfullmean[iCal-1]->SetMarkerStyle(20);
    hhalfmean[iCal-1]->SetMarkerSize(0.5);
    hfullmean[iCal-1]->SetMarkerSize(0.5);
    hfullmean[iCal-1]->Rebin(2);
    hfullmean[iCal-1]->Reset(); // cant calc with rebin(2)
    
    for (int iBin = 1; iBin <= hfullmean[0]->GetNbinsX(); iBin++ ) { 
      
      double time = hfullmean[iCal-1]->GetBinCenter(iBin) + tcyc/24.*((double)(iCal-1)) + tcyc*phase_fr[iCal-1]/(2.*TMath::Pi()); 
      int np = 10; // points integrated
      double Tb = 150.; // ns bin width
      double dtime; // point coordinate
      double sum_yt = 0.0, sum_y = 0.0;
      
      for (int ip = 0; ip < np; ip++) {
	
	double dtime = time - Tb/2. + Tb/2./np + ip*Tb/np;
	double dyfr = 1 + ffrpw->Eval(dtime); // +1 as correction is histogram + 1 multiplier 
	sum_yt += dyfr*dtime;
	sum_y += dyfr;
	if (dbg) printf(" ip %i, xxp %f, df %f, sum y %f, sum y*t %f\n", ip, dtime, dyfr, sum_y, sum_yt);
      }
      double dtmn = sum_yt/sum_y - time;
      
      if ( time - tcyc/24.*((double)(iCal-1)) > 30500 && time -tcyc/24.*((double)(iCal-1)) < 69500) {
	hfullmean[iCal-1]->SetBinContent( iBin, dtmn ); // units are time
      } else { 
	hfullmean[iCal-1]->SetBinContent( iBin, 0.0 );  // units are time
      }
    }
  
    //fix boundary bins of piece function - very awkward, needs better function 
    /*
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[0]), 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[0])-1, 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[0])+1, 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[1]), 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[1])-1, 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[1])+1, 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[2]), 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[2])-1, 0.0 );
    hfull[iCal-1]->SetBinContent( hfull[iCal-1]->FindBin(pwhi[2])+1, 0.0 );
    */
  
    c1->cd(iCal);
    hfullmean[iCal-1]->Draw("hist");
    hfull[iCal-1]->Draw("histsame");
  }
  printf("made hfull correction histograms\n");

  if (!storeHalfFullHistograms) return 0;

  c1->Update();  
  c1->Clear();  

  for (int iCal = 1; iCal <= nCalfr; iCal++){

    if (plotHalfFullHistograms) {
      hcyc[iCal-1]->Draw("hist");
      hhalf[iCal-1]->Draw("same");
      hhalfmean[iCal-1]->Draw("same");
      hfull[iCal-1]->Draw("same");
      hfullmean[iCal-1]->Draw("same");
    }
    c1->Update();
    sprintf(hname,"calo-fr-%i.png",iCal);
    c1->SaveAs(hname);

  }

  printf("write half, full histos\n");
  if (!do_wiggle_fits){
    for (int iCal = 1; iCal <= nCalfr; iCal++) hhalf[iCal-1]->Write();
    for (int iCal = 1; iCal <= nCalfr; iCal++) hfull[iCal-1]->Write();
    for (int iCal = 1; iCal <= nCalfr; iCal++) hhalfmean[iCal-1]->Write();
    for (int iCal = 1; iCal <= nCalfr; iCal++) hfullmean[iCal-1]->Write();
    ffrpw->Write();
    //file1->Close();
  }
  
  return 0;  
}


int goFastRot2();
int goFastRot2(){

  makeHalfFullHistograms=0;
  if (doShiftBins) phi_cyc = phi_cyc + 3.14;
  fixamp1[0]=0;
  fixtau1[0]=0;
  fixasym[0]=1;
  fixomgb[0]=1;
  fixphib[0]=1;
  pwlo[0]=29800;
  pwhi[0]=36000;
  pwlo[1]=36000;
  pwhi[1]=36000;
  pwlo[2]=36000;
  pwhi[2]=36000;
  pwlo[3]=36000;
  pwhi[3]=36000;
  fixdcyc=1;
  doFRfitWithInts=0;
  plotFastRot2( 29800, 36000, 149.3, 0);
  gcyc->Draw("ap");
  ffrpw->Draw("histsame");

  Ampl1[0]=ffrpw->GetParameter(2);
  Tau1[0]=ffrpw->GetParameter(3);
  fixamp1[0]=0;
  fixtau1[0]=0;
  fixasym[0]=0;
  fixomgb[0]=1;
  fixphib[0]=1;
  pwlo[0]=29800;
  pwhi[0]=50000;
  pwlo[1]=50000;
  pwhi[1]=50000;
  pwlo[2]=50000;
  pwhi[2]=50000;
  pwlo[3]=50000;
  pwhi[3]=50000;
  fixdcyc=1;
  doFRfitWithInts=0;
  plotFastRot2( 29800, 50000, 149.3, 0);
  gcyc->Draw("ap");
  ffrpw->Draw("histsame");

  Tau1[0]=ffrpw->GetParameter(3);
  fixamp1[0]=0;
  fixtau1[0]=1;
  fixasym[0]=0;
  fixomgb[0]=0;
  fixphib[0]=0;;
  pwlo[0]=29800;
  pwhi[0]=50000;
  pwlo[1]=50000;
  pwhi[1]=50000;
  pwlo[2]=50000;
  pwhi[2]=50000;
  pwlo[3]=50000;
  pwhi[3]=50000;
  fixdcyc=1;
  doFRfitWithInts=0;
  plotFastRot2( 29800, 69000, 149.3, 0);
  gcyc->Draw("ap");
  ffrpw->Draw("histsame");

  Asym[0]=ffrpw->GetParameter(6);
  omega_beat[0]=ffrpw->GetParameter(7);
  phi_beat[0]=ffrpw->GetParameter(8);
  fixamp1[0]=0;
  fixtau1[0]=1;
  fixasym[0]=1;
  fixomgb[0]=0;
  fixphib[0]=0;
  pwlo[0]=29800;
  pwhi[0]=69000;
  pwlo[1]=69000;
  pwhi[1]=69000;
  pwlo[2]=69000;
  pwhi[2]=69000;
  pwlo[3]=69000;
  pwhi[3]=69000;
  doFRfitWithInts=1;
  plotFastRot2( 29800, 69000, 149.3, 0);
  gcyc->Draw("ap");
  ffrpw->Draw("histsame");

  Ampl1[0]=ffrpw->GetParameter(2);
  Tau1[0]=ffrpw->GetParameter(3);
  Asym[0]=ffrpw->GetParameter(6);
  omega_beat[0]=ffrpw->GetParameter(7);
  phi_beat[0]=ffrpw->GetParameter(8);
  fixamp1[0]=0;
  fixtau1[0]=1;
  fixasym[0]=1;
  fixomgb[0]=0;
  fixphib[0]=0;
  pwlo[0]=29800;
  pwhi[0]=69000;
  pwlo[1]=69000;
  pwhi[1]=69000;
  pwlo[2]=69000;
  pwhi[2]=69000;
  pwlo[3]=69000;
  pwhi[3]=69000;
  makeHalfFullHistograms = 1;
  storeHalfFullHistograms = 1;
  doFRfitWithInts=1;
  plotFastRot2( 29800, 69000, 149.3, 0);
  gcyc->Draw("ap");
  ffrpw->Draw("histsame");

  printf("doFRfitWithInts %i\n", doFRfitWithInts);

  return 0;
  }


int plotFastRot(){
  
  char cname[64], finname[64];

  //sprintf( finname, "noxtaldata/xmas_w4_t16.root");
  sprintf( finname, "noxtaldata/60hr-newerrors-withcalib-w2.root");

  file0 = TFile::Open( finname,"UPDATE");
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_sig_%i_0", 1);
  dir0->GetObject( hname, htmp);
  printf("got input histogram object %s\n", hname);

  hfastrotation = new TH1D("hfastrotation","hfastrotation", 24*htmp->GetNbinsX()/2, 0, 24*htmp->GetNbinsX()/2);
  
  for (int iCal = 1; iCal <= 24; iCal++) {
   
    sprintf( hname, "qHist1D_sig_%i_0", iCal);
    dir0->GetObject( hname, htmp);
   
    for (int ib = 1; ib <= htmp->GetNbinsX()/2; ib++){

      double r;
      if  ( htmp->GetBinContent(2*ib) + htmp->GetBinContent(2*ib-1) != 0) {
	r = ( htmp->GetBinContent(2*ib) - htmp->GetBinContent(2*ib-1) ) / ( htmp->GetBinContent(2*ib) + htmp->GetBinContent(2*ib-1) );
      } else {
	r = 0;
      }

      int ibin  = iCal + 24*(ib-1);
      printf(" ib %i, iCal %i, ibin %i, r %f \n", ib, iCal, ibin, r);
      hfastrotation->SetBinContent( ibin, r);
    }

    /*
    h1[iCal-1] = (TH1D*)htmp->Clone();
    h2[iCal-1] = (TH1D*)htmp->Clone();
    sprintf( hname, "qHist1D_fastrot1_%i_0", iCal);
    h1[iCal-1]->SetName(hname);
    h1[iCal-1]->SetTitle(hname);
    sprintf( hname, "qHist1D_fastrot2_%i_0", iCal);
    h2[iCal-1]->SetName(hname);
    h2[iCal-1]->SetTitle(hname);
    */

    /*
    sprintf( hname, "qHist1D_%i_0", iCal);
    dir0->GetObject( hname, htmp);
    double a1 = htmp->Integral(200,220);
    double a2 = htmp->Integral(2900,2920);
    printf("a1 %f, a2 %f, a1/a2 %f", a1, a2, a1/a2);
    */

    /*
    printf("h1, iCal %i, %s\n",iCal,h1[iCal-1]->GetName());
    printf("h2, iCal %i, %s\n",iCal,h2[iCal-1]->GetName());
    h2[iCal-1]->Rebin(2);

    for (int iBin = 1; iBin <= h2[iCal-1]->GetNbinsX(); iBin++) {
      double y1 = h1[iCal-1]->GetBinContent(2*iBin-1);
      double y2 = h1[iCal-1]->GetBinContent(2*iBin);
      double dy=0.0;
      if (y1+y2 > 0) dy = (y1-y2)/(y1+y2);
      h2[iCal-1]->SetBinContent(iBin,dy);
    }
    */

    //printf("h2, iCal %i, %s\n",iCal,h2[iCal-1]->GetName());
    //h2[iCal-1]->GetXaxis()->SetRangeUser(30000.,100000.);
    //h2[iCal-1]->Draw("hist");
    //c1->Update();
    //sleep(1);

    //h2[iCal-1]->Write("",TObject::kOverwrite);
    //h2[iCal-1]->Write();
    //printf("wrote histogram %s into file %s\n", hsname, foutname);

    //htmp->SetMinimum(0.);
    //htmp->SetMaximum(2*2000000.);

    //c1->Update();
    //sleep(0.3);
    //sprintf( hname, "qHist1D_sig_%i_0_2.png", iCal);
    //c1->SaveAs(hname);
  }

  hfastrotation->Draw("hist");
  c1->Update();
  sleep(1);

  //file0->Close();
  return 1;
}


int plotPedSize(){
  
  char cname[64], finname[64];

  sprintf( finname, "noxtaldata/xmas_w4_t16.root");
  file0 = TFile::Open( finname,"UPDATE");
  printf("got input root file %s\n", finname);

  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  double Ni = 24, xi[24], yi[24];
  
  for (int iCal = 1; iCal <= 24; iCal++) {
   

    sprintf( hname, "qHist1D_%i_0", iCal);
    dir0->GetObject( hname, htmp);
    yi[iCal-1] = htmp->Integral(200,220);
    xi[iCal-1] = (float)iCal;
    printf("xi %f, yi %f\n", xi[iCal-1], yi[iCal-1]);
  }

  file0->Close();

  gPed = new TGraph(Ni, xi, yi);
  gPed->SetName("gPed");
  gPed->SetTitle("pedestal amplitude versus calo index");
  gPed->GetXaxis()->SetTitle("calo index");
  gPed->GetYaxis()->SetTitle("ped amplitude (arb. units)");
  gPed->GetYaxis()->SetTitleOffset(0.5);
  gPed->GetYaxis()->SetTitleSize(0.07);
  gPed->GetXaxis()->SetTitleOffset(0.5);
  gPed->GetXaxis()->SetTitleSize(0.07);
  gPed->SetMarkerStyle(8);
  gPed->SetMarkerSize(1.0);
  gPed->SetMarkerColor(kRed+1);
  gPed->SetLineWidth(2);
  gPed->SetLineColor(kRed+1);
  gPed->Draw("AP");

  return 1;
}

int plotCBO(){
  
  char cname[64], finname[64];
    
  sprintf( finname, "noxtaldata/8129xx.root");
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);

  c1->Clear();
  gPad->SetLogy();
  
  for (int iCal = 12; iCal <= 24; iCal++) {
    
    sprintf( hname, "qHist2D_ydiff_%i", iCal);
    dir0->GetObject( hname, h2tmp);
    sprintf( hname, "qHist2D_ydiff_%i_x", iCal);
    h1[iCal-1] = h2tmp->ProjectionX( hname, 500, 590);

    /*
    h2tmp->Draw("colz");
    c1->Update();
    h1[iCal-1]->Draw();
    c1->Update();
    */

    h1[iCal-1]->Rebin(2);
    h1[iCal-1]->SetMinimum(2000.);
    h1[iCal-1]->SetMaximum(40000.);
    h1[iCal-1]->GetXaxis()->SetRangeUser(30000.,60000.);
    if (iCal >= 1 && iCal <= 4) h1[iCal-1]->SetLineColor(kRed);
    if (iCal >= 5 && iCal <= 8) h1[iCal-1]->SetLineColor(kBlue);
    if (iCal >= 9 && iCal <= 12) h1[iCal-1]->SetLineColor(kGreen);
    if (iCal >= 13 && iCal <= 16) h1[iCal-1]->SetLineColor(kMagenta);
    if (iCal >= 17 && iCal <= 20) h1[iCal-1]->SetLineColor(kYellow);
    if (iCal >= 21 && iCal <= 24) h1[iCal-1]->SetLineColor(kAzure);
    double wid = ((double) ((iCal-1)%4));
    h1[iCal-1]->SetLineWidth(0.5*wid + 1.);
    if (iCal==12){
      h1[iCal-1]->Draw("");
      c1->Update();
    } else {
      h1[iCal-1]->Draw("same");
      c1->Update();
 }
    if (iCal==24) c1->SaveAs("muonlossintermediate2.png");
  }

  
  /*
  h1[0]->SetLineColor(25);
  h1[0]->Draw();
  printf("calo %i plotted\n", 1);

  for (int iCal = 2; iCal <= 24; iCal++) {
    h1[iCal-1]->SetLineColor(24+iCal);
    h1[iCal-1]->Draw("same");
    printf("calo %i plotted\n", iCal);
  }
  */

  file0->Close();
  return 1;
}


int plotRMS( char* iRun){

  hRMS = new TH1D("hRMS","calo / xtal distribution of noise RMS", 54*24, 0.5, 0.5+54.*24. );
  sRMS = new TH2D("sRMS","calo / xtal distribution of noise RMS", 54, 0.5, 54.5, 24, 0.5, 24.5 );
  
  char cname[64], finname[64];
    
  for (int iCal = 1; iCal <= 24; iCal++) {

    sprintf( finname, "%04s/run0%04s_calo%02i_rmsthres_nocalibrate.root", iRun, iRun, iCal);
    file0 = TFile::Open( finname);
    printf("got input root file %s\n", finname);
  
    TDirectoryFile *dir0;
    file0->GetObject( dname, dir0);
    printf("got input root directory %s\n", dname);

    for (int iSeg = 0; iSeg < 54; iSeg++) {
 
      sprintf( hname, "qHist1D_rms_%i_%i", iCal, iSeg);
      dir0->GetObject( hname, htmp);
      double rms = htmp->GetRMS();
      int iBin = 54*(iCal-1)+iSeg+1;
      hRMS->SetBinContent( iBin, rms);
      sRMS->SetBinContent( iSeg+1, iCal, rms);
      printf("got calo %i, xtal %i, ibin %i, rms %f, histogram %s\n", iCal, iSeg, iBin, rms, hname);
      htmp->Clear(); 
    }

    file0->Close();
  }

  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  hRMS->Draw();
  c1->cd(2);
  sRMS->Draw("colz");

  sprintf( hname, "rmsPLOT_%04s.png", iRun);
  c1->SaveAs( hname);
  printf("done plotting RMS run %04s\n", iRun);

  return 1;
}

int plotCnts( char* iRun){
  
  hCnts = new TH1D("hCnts","calo / xtal distribution of above threshold counts", 54*24, 0.5, 0.5+54.*24. );
  sCnts = new TH2D("sCnts","calo / xtal distribution of above threshold counts", 54, 0.5, 54.5, 24, 0.5, 24.5 );

  char cname[64], finname[64];
  TDirectoryFile *dir0;
    
  for (int iCal = 1; iCal <= 24; iCal++) {

    sprintf( finname, "%04s/run0%04s_calo%02i_rmsthres_nocalibrate.root", iRun, iRun, iCal);
    file0 = TFile::Open( finname);
    printf("got input root file %s\n", finname);
  
    file0->GetObject( dname, dir0);
    printf("got input root directory %s\n", dname);

    for (int iSeg = 0; iSeg < 54; iSeg++) {
 
      sprintf( hname, "qHist1D_sig_%i_%i", iCal, iSeg);
      dir0->GetObject( hname, htmp);
      double Cnts = htmp->Integral( 5000, 35000);
      int iBin = 54*(iCal-1)+iSeg+1;
      hCnts->SetBinContent( iBin, Cnts);
      sCnts->SetBinContent( iSeg+1, iCal, Cnts);
      printf("got calo %i, xtal %i, ibin %i, sum %f, histogram %s\n", iCal, iSeg, iBin, Cnts, hname);
      htmp->Clear(); 
    }

    file0->Close();
  }

  // make signals positive
  hCnts->Scale(-1.0);
  sCnts->Scale(-1.0);

  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  hCnts->Draw();
  c1->cd(2);
  sCnts->Draw("colz");

  printf("done plotting Cnts run %04s\n", iRun);

  return 1;
}

// looking at spike at ~81us in run group 18xx
int plotRF( char* iRun){

  char cname[64], finname[64];

  sprintf( finname, "%04s/run0%04s_sum_rmsthres_nocalibrate.root", iRun, iRun);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  c1->Clear();
  c1->Divide(4,6);
    
  for (int iCal = 1; iCal <= 24; iCal++) {

      sprintf( hname, "qHist1D_sig_%i_sum", iCal);
      file0->GetObject( hname, htmp);
      printf("got calo %i, histogram %s\n", iCal, hname);
      sprintf( hname, "hsig_%i_sum", iCal);
      htmp->SetName(hname);
      htmp->SetTitle(hname);
      htmp->SetTitleSize(0.05);
      htmp->SetLineWidth(2.5);
      htmp->Rebin(8);
      htmp->GetXaxis()->SetRange(1200.,1300.);
      c1->cd(iCal);
      htmp->Draw();
      c1->Update();
  }

  //file0->Close();

  printf("done plotting RF peak run %04s\n", iRun);

  return 1;
}

int goDrawSync(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawSync( iRun, iCal);
  }

  return 1;
}

int drawSync(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "1912/run%05i_calo%02i_rmsthres_nocalibrate.root", iRun, iCal);
  //sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->GetXaxis()->SetRange( 610, 650 );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int chanPerMod = 5;
  c1->Clear();
  c1->Divide(3,4);
  c1->cd(1);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg/chanPerMod + 1;
    int iColor = iSeg%chanPerMod + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->GetXaxis()->SetRange( 610, 650 );
    
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }

  sprintf( cname, "png/run%05i_%02i_SyncPulse.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_SyncPulse.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawFlash(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {
    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawFlash( iRun, iCal);
  }

  return 1;
}

int drawFlashNew(){
  
  char cname[64], finname[64];
  
  //sprintf( finname, "noxtaldata/786x_w4_t16.root");
  sprintf( finname, "8129.root");
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  c1->Clear();
  c1->Divide(6,4);

  for (int iCal = 1; iCal <= 24; iCal++) {

    /*
    sprintf( hname, "qHist1D_%i_0", iCal);
    dir0->GetObject( hname, htmp1);
    htmp1->SetName(hsname);
    htmp1->SetTitle(hsname);
    htmp1->SetLineColor( 1 );
    htmp1->SetLineWidth( 2 );
    htmp1->GetXaxis()->SetRangeUser( 24000, 26000 );
    */

    sprintf( hname, "qHist2D_ydiff_%i", iCal);
    dir0->GetObject( hname, h2tmp);
    h2tmp->SetName(hname);
    h2tmp->SetTitle(hname);
    h2tmp->RebinX(2);
    h2tmp->RebinY(8);
    h2tmp->SetMaximum(40);
    h2tmp->GetXaxis()->SetRangeUser( 24000, 35000 );

    c1->cd(iCal);
    
    //h1[0]->SetMinimum(-2.e6);
    //h1[0]->SetMaximum(2.e5);
    //h1[0]->SetMinimum(-8000);
    //h1[0]->SetMaximum(1000);
    //h1[0]->SetMinimum( 2.*h1[0]->GetMinimum());
    //h1[0]->SetMaximum( 2.*h1[0]->GetMaximum());

    //htmp1->Draw("HIST");
    h2tmp->Draw("colz");
    
  }
  sprintf( cname, "8129_flash.png"); 
  c1->SaveAs( cname );
  printf("done plotting run %s", cname);

  return 1;
}

TH1D *hflash = new TH1D("hflash","flash versus calorimeter", 24, 0., 24. );
TH1D *hpositron = new TH1D("hpositron","positron versus calorimeter", 24, 0., 24. );
int drawFlashDistbn(){
  
  char cname[64], finname[64];
  
  //sprintf( finname, "noxtaldata/lego8129.root");
  sprintf( finname, "noxtaldata/r13680.root");
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  c1->Clear();

  for (int iCal = 1; iCal <= 24; iCal++) {

    sprintf( hname, "qHist2D_sig_flashcalo_%i", iCal);
    dir0->GetObject( hname, h2tmp);
    double flash = h2tmp->Integral();
    sprintf( hname, "qHist2D_sig_legocalo_%i", iCal);
    dir0->GetObject( hname, h2tmp);
    double positron = h2tmp->Integral();

    hflash->SetBinContent( iCal, flash);     
    hpositron->SetBinContent( iCal, positron);     
    printf("position, flash %i, %f %f\n", iCal, flash, positron);
   
  }

  hflash->SetLineWidth(4.);
  hpositron->SetLineWidth(4.);
  hflash->SetLineColor(kRed);
  hpositron->SetLineColor(kBlue);

  hflash->Draw("HIST");
  hpositron->Draw("HISTsame");
  sprintf( cname, "8129flashvcalo.png"); 
  c1->SaveAs( cname );
  printf("done plotting run %s", cname);

  return 1;
}

int drawFlash(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "1912/run%05i_calo%02i_rmsthres_nocalibrate.root", iRun, iCal);
  //sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->GetXaxis()->SetRange( 3160, 3320 );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);
  //h1[0]->SetMinimum(-2.e6);
  //h1[0]->SetMaximum(2.e5);
  //h1[0]->SetMinimum(-8000);
  //h1[0]->SetMaximum(1000);
  h1[0]->SetMinimum( 2.*h1[0]->GetMinimum());
  h1[0]->SetMaximum( 2.*h1[0]->GetMaximum());
  sprintf( hname,"Calo %i Xtal %i row", iCal, 0);
  h1[0]->SetTitle(hname);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg%xtalPerRow + 1;
    int iColor = iSeg/xtalPerRow + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
    h1[0]->SetTitle(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->GetXaxis()->SetRange( 3160, 3320 );
    
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      //h1[iSeg]->SetMinimum(-2.e6);
      //h1[iSeg]->SetMaximum(2.e5);
      //h1[iSeg]->SetMinimum(-8000);
      //h1[iSeg]->SetMaximum(2000);
      h1[iSeg]->SetMinimum( 2.*h1[iSeg]->GetMinimum());
      h1[iSeg]->SetMaximum( 2.*h1[iSeg]->GetMaximum());
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }

  sprintf( cname, "png/run%05i_%02i_Flash_flush1.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_Flash_flush1.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawThrSplash(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawThrSplash( iRun, iCal);
  }

  return 1;
}

int drawThrSplash(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "1912/run%05i_calo%02i_rmsthres_nocalibrate.root", iRun, iCal);
  //sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_sig_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_sig_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->Rebin( 16 );
  h1[0]->GetXaxis()->SetRange( 12000/8, 22000/8  );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);
  h1[0]->SetMinimum(-8.e2);
  h1[0]->SetMaximum(1.e2);
  sprintf( hname,"Calo %i Xtal %i row", iCal, 0);
  h1[0]->SetTitle(hname);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg%xtalPerRow + 1;
    int iColor = iSeg/xtalPerRow + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_sig_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
    h1[0]->SetTitle(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->Rebin( 16 );
    h1[iSeg]->GetXaxis()->SetRange( 12000/8, 22000/8 );
     
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      h1[iSeg]->SetMinimum(-8.e2);
      h1[iSeg]->SetMaximum(1.e2);
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }
  

  sprintf( cname, "png/run%05i_%02i_SihSplash.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_SigSplash.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawThrStore(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawThrStore( iRun, iCal);
  }

  return 1;
}

int drawThrStore(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_thr_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_thr_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->Rebin( 16 );
  h1[0]->GetXaxis()->SetRange( 1800/16, 8000/16  );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);
  h1[0]->SetMinimum(-8.e2);
  h1[0]->SetMaximum(1.e2);
  sprintf( hname,"Calo %i Xtal %i row", iCal, 0);
  h1[0]->SetTitle(hname);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg%xtalPerRow + 1;
    int iColor = iSeg/xtalPerRow + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_thr_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
    h1[0]->SetTitle(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->Rebin( 16 );
    h1[iSeg]->GetXaxis()->SetRange( 1800/16, 8000/16 );
     
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      h1[iSeg]->SetMinimum(-8.e2);
      h1[iSeg]->SetMaximum(1.e2);
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }
  

  sprintf( cname, "png/run%05i_%02i_ThrStore.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_ThrStore.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawPed1(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawPed1( iRun, iCal);
  }

  return 1;
}

int drawPed1(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  //sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  sprintf( finname, "Data/run%05i_%i_flush1.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->GetXaxis()->SetRange(  1500, 3500 );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);
  //h1[0]->SetMinimum(-5.e3);
  //h1[0]->SetMaximum(2.e4);
  h1[0]->SetMinimum(-2.e1);
  h1[0]->SetMaximum(1.e2);
  sprintf( hname,"Calo %i Xtal %i row", iCal, 0);
  h1[0]->SetTitle(hname);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg%xtalPerRow + 1;
    int iColor = iSeg/xtalPerRow + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
    h1[0]->SetTitle(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->GetXaxis()->SetRange( 1500, 3500 );
    
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      //h1[iSeg]->SetMinimum(-5.e3);
      //h1[iSeg]->SetMaximum(2.e4);
      h1[iSeg]->SetMinimum(-2.e1);
      h1[iSeg]->SetMaximum(1.e2);
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }
  

  sprintf( cname, "png/run%05i_%02i_Ped1_flush1.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_Ped1_flush1.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawPed2(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawPed2( iRun, iCal);
  }

  return 1;
}

int drawPed2(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  sprintf( hname, "qHist1D_%i_0", iCal);
  dir0->GetObject( hname, h1[0]);
  sprintf( hsname,"qHist1D_%i_%i", iCal, 0);
  h1[0]->SetName(hsname);
  h1[0]->SetTitle(hsname);
  h1[0]->SetLineColor( 1 );
  h1[0]->SetLineWidth( 2 );
  h1[0]->GetXaxis()->SetRange(  1500, 1900 );
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);
  h1[0]->SetMinimum(-5.e3);
  h1[0]->SetMaximum(2.e4);
  sprintf( hname,"Calo %i Xtal %i row", iCal, 0);
  h1[0]->SetTitle(hname);
  h1[0]->Draw();
  printf("iSeg %i, iPanel %i, iColor %i\n", 0, 1, 1);    

  for (int iSeg = 1; iSeg < 54; iSeg++) {
 
    int iPanel = iSeg%xtalPerRow + 1;
    int iColor = iSeg/xtalPerRow + 1; 
    printf("iSeg %i, iPanel %i, iColor %i\n", iSeg, iPanel, iColor);    

    sprintf( hname, "qHist1D_%i_%i", iCal, iSeg);    
    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
    h1[0]->SetTitle(hname);
    h1[iSeg]->SetLineColor( iColor );
    h1[iSeg]->SetLineWidth( 2 );
    h1[iSeg]->GetXaxis()->SetRange( 1500, 1900 );
    
    c1->cd( iPanel );
    if ( !(iColor-1) ) {
      h1[iSeg]->SetMinimum(-5.e3);
      h1[iSeg]->SetMaximum(2.e4);
      h1[iSeg]->Draw();
    } else {
      h1[iSeg]->Draw("same");
    }
  }
  

  sprintf( cname, "png/run%05i_%02i_Ped2.png", iRun, iCal); 
  c1->SaveAs( cname );
  sprintf( cname, "pdf/run%05i_%02i_Ped2.pdf", iRun, iCal); 
  c1->SaveAs( cname );

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int goDrawPedvSample(int iRun){

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 12 ) continue; // skip calo 11 - this calo is noisy
    //if ( iCal == 14 ) continue; // skip calo 14 - this calo is down 6/12/17
    drawPedvSample( iRun, iCal);
  }

  return 1;
}

int drawPedvSample(int iRun, int iCal){
  
  char cname[64], finname[64];
  
  sprintf( finname, "Data/run%05i_%i_Coarse.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);   
  
  int xtalPerRow = 9;
  c1->Clear();
  c1->Divide(3,3);
  c1->cd(1);

  int iPanel;
  for (int iRow = 1; iRow <= 6; iRow++) { 

    for (int iSeg = 9*(iRow-1); iSeg < 9*iRow; iSeg++) { 
      
      iPanel = iSeg%xtalPerRow + 1;
      iRow = iSeg/xtalPerRow + 1; 
      printf("iSeg %i, iPanel %i, iRow %i\n", iSeg, iPanel, iRow);    
      
      sprintf( hname, "qHist2D_%i_%i", iCal, iSeg);    
      dir0->GetObject( hname, s1[iSeg]);
      s1[iSeg]->SetName(hname);
      sprintf( hname,"Calo %i Xtal %i row", iCal, iSeg);
      s1[0]->SetTitle(hname);
      c1->cd( iPanel );
      s1[iSeg]->Draw("colz");
      
    }

    sprintf( cname, "png/run%05i_%02i_Row%01i_PedvSample.png", iRun, iCal, iRow); 
    c1->SaveAs( cname );
    //sprintf( cname, "pdf/run%05i_%02i_Row%01i_PedvSample.pdf", iRun, iCal, iRow); 
    //c1->SaveAs( cname );

  }

  printf("done plotting run %i, calo %i\n", iRun, iCal);

  return 1;
}

int sumXtal(int iRun, int iCal, char* cHist){
 
  /*
sums the histograms for all xtals for a particular histogram, calorimeter, 
and run from individual histograms in finame and writes summed histogram in 
foutname 
  */

  char finname[64];
  TDirectoryFile *dir0;

  // INPUT FILE
  sprintf( finname, "18xx/run018xx_calo%02i_rmsthres3_nocalibrate.root", iCal);
  //sprintf( finname, "1874/run01874_calo%02i_rmsthres3_nocalibrate.root", iCal);
  //sprintf( finname, "1912/run01912_calo%02i_rmsthres3_nocalibrate.root", iCal);
  //sprintf( finname, "18xx/run018xx_calo%02i_rmsthres2_nocalibrate.root", iCal);
  //sprintf( finname, "1874/run01874_calo%02i_rmsthres2_nocalibrate.root", iCal);
  //sprintf( finname, "1912/run01912_calo%02i_rmsthres2_nocalibrate.root", iCal);
  //sprintf( finname, "1874/run01874_calo%02i_rmsthres_nocalibrate.root", iCal);
  //sprintf( finname, "18xx/run018xx_calo%02i_rmsthres_nocalibrate.root", iCal);
  //sprintf( finname, "1912/run01912_calo%02i_rmsthres_nocalibrate.root", iCal);
  //sprintf( finname, "1912/run01912_calo%02i_calibrate.root", iCal);
  //sprintf( finname, "18xx/run018xx_calo%02i_calibrate.root", iCal);
  //sprintf( finname, "1874/run01874_calo%02i_calibrate.root", iCal);
  //sprintf( finname, "1912/run01912_calo%02i_nocalibrate.root", iCal);
  //sprintf( finname, "18xx/run018xx_calo%02i_nocalibrate.root", iCal);
  //sprintf( finname, "1874/run01874_calo%02i_nocalibrate.root", iCal);

  //sprintf( finname, "18xx/calo%i_calibrate.root", iCal);
  //sprintf( finname, "DataErrorBarUpdate/calo%irun%04i.root", iCal, iRun);
  //sprintf( finname, "Data/run%05i_%i_Coarse_testped.root", iRun, iCal);

  // OUTPUT FILE
  sprintf( foutname, "18xx/run018xx_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_calibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_calibrate.root");
  //sprintf( foutname, "1874/run01874_sum_calibrate.root");
  //sprintf( foutname, "1912/run01912_sum_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_nocalibrate.root");

  //sprintf( foutname, "18xx/sum_calibrate.root");
  //sprintf( foutname, "Data/sum%05i_coarse_testped.root", iRun);

  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);

  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "%s_%i_0", cHist, iCal);
  dir0->GetObject( hname, hSum);
  sprintf( hsname,"%s_%i_sum", cHist, iCal);
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);

  // hot, inner xtals are 8,17,26,35,44,53
  // cold, inner xtals are 0,9,18,27,36,45

  for (int iSeg = 1; iSeg < 54; iSeg++) {

    sprintf( hname, "%s_%i_%i", cHist, iCal, iSeg);

    dir0->GetObject( hname, h1[iSeg]);
    h1[iSeg]->SetName(hname);
    printf("got calo %i, xtal %i, histogram %s\n", iCal, iSeg, hname);

    // hot, inner xtals are 8,17,26,35,44,53
    // cold, inner xtals are 0,9,18,27,36,45
    if ( (iSeg%9 != 0) && (iSeg%9 != 1) ) hSum->Add( h1[iSeg], 1.0);
  }

  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("created output root file %s\n", foutname);

  hSum->Write();
  file1->Print();
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file0->Close();
  file1->Close();

  return 1;
}

int sumXtal2D(int iRun, int iCal, char* cHist){
 
  /*
sums the histograms for all xtals for a particular histogram, calorimeter, 
and run from individual histograms in finame and writes summed histogram in 
foutname 
  */

  char finname[64];

  sprintf( finname, "Data/run%05i_%i_Coarse_testped.root", iRun, iCal);
  file0 = TFile::Open( finname);
  printf("got input root file %s\n", finname);
 
  TDirectoryFile *dir0;
  file0->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);

  sprintf( hname, "%s_%i_0", cHist, iCal);
  dir0->GetObject( hname, sSum);
  sprintf( hsname,"%s_%i_sum", cHist, iCal);
  sSum->SetName(hsname);
  sSum->SetTitle(hsname);
  printf("got calo %i, xtal %i, histogram %s\n", iCal, 0, hname);

  for (int iSeg = 1; iSeg < 54; iSeg++) {

    sprintf( hname, "%s_%i_%i", cHist, iCal, iSeg);

    dir0->GetObject( hname, s1[iSeg]);
    s1[iSeg]->SetName(hname);
    printf("got calo %i, xtal %i, histogram %s\n", iCal, iSeg, hname);

    sSum->Add( s1[iSeg], 1.0);
  }

  sprintf( foutname, "Data/sum%05i_coarse_testped.root", iRun);
  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("created output root file %s\n", foutname);

  sSum->Write();
  file1->Print();
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  
  file0->Close();
  file1->Close();

  return 1;
}


int sumCalo(int iRun, char* cHist){

  printf("sum Calo\n");

  // INPUT / OUTPUT FILE
  sprintf( foutname, "18xx/run018xx_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres3_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres2_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_rmsthres_nocalibrate.root");
  //sprintf( foutname, "1912/run01912_sum_calibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_calibrate.root");
  //sprintf( foutname, "1874/run01874_sum_calibrate.root");
  //sprintf( foutname, "1912/run01912_sum_nocalibrate.root");
  //sprintf( foutname, "18xx/run018xx_sum_nocalibrate.root");
  //sprintf( foutname, "1874/run01874_sum_nocalibrate.root");

  //sprintf( foutname, "18xx/sum_calibrate.root");
  //sprintf( foutname, "DataErrorBars/sum%04i.root", iRun);
  //sprintf( foutname, "DataErrorBars/sum%04i.root", iRun);
  //sprintf( foutname, "Data/sum%05i_coarse_testped.root", iRun);

  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("got sum root file %s\n", foutname);

  //sprintf( hname, "qHist1D_sig_%i_sum", 1);
  sprintf( hname,"%s_%i_sum", cHist, 1);
  file1->GetObject( hname, hSum);
  //sprintf( hsname, "qHist1D_sig_calo_sum"); // includes bad Calos
  sprintf( hsname,"%s_sum_noBadCalo", cHist); // excludes bad Calos
  hSum->SetName(hsname);
  hSum->SetTitle(hsname);
  printf("got run %i, calo %i, histogram %s\n", iRun, 1, hname);

  for (int iCal = 1; iCal <= 24; iCal++) {

    //if ( iCal == 4 ) continue; // skip - noise Xtal 39
    //if ( iCal == 12 ) continue; // skip - bad noise
    //if ( iCal == 14 ) continue; // skip - calo off
    //if ( iCal == 17 ) continue; // skip - noise Xtal 2,47,39,31,52,53
    //if ( iCal == 19 ) continue; // skip - noise Xtal 46,20,31,32
    //if ( iCal == 23 ) continue; // skip - noise Xtal 1,29,31,51,52
    //if ( iCal == 24 ) continue; // skip - noise Xtal 43

    sprintf( hname, "%s_%i_sum", cHist, iCal);
    file1->GetObject( hname, h2[iCal]);
    hSum->Add( h2[iCal], 1.0);
    printf("got run %i, calo %i, histogram %s\n", iRun, iCal, hname);

  }

  hSum->Write(); 
  file1->Print(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  file1->Close();

  return 1;
}

int sumCaloXtal(int iRun, char* cHist){

  char finname[64];
  
  // open calo 1 file
  sprintf( finname, "Data/run%05i_%i_Coarse_testped.root", iRun, 1);
  file1 = TFile::Open( finname, "READ");
  if ( file1->IsOpen() ) printf("read input root file %s\n", finname);
  
  TDirectoryFile *dir0;
  file1->GetObject( dname, dir0);
  printf("got input root directory %s\n", dname);
  
  for (int iSeg = 0; iSeg < 53; iSeg++) {
    
    sprintf( hname,"%s_%i_%i", cHist, 1, iSeg);
    dir0->GetObject( hname, h1[iSeg]);
    printf("got run %i, calo %i, histogram %s\n", iRun, 1, hname);
    
    sprintf( hsname,"%s_CaloSum_%i", cHist, iSeg);
    h1[iSeg]->SetName(hsname);
    h1[iSeg]->SetTitle(hsname);
    printf("renamed histogram %s\n", hsname);

  } // loop over xtals


  for (int iCal = 2; iCal <= 24; iCal++) {

    sprintf( finname, "Data/run%05i_%i_Coarse_testped.root", iRun, iCal);
    file2 = TFile::Open( finname, "READ");
    if ( file2->IsOpen() ) printf("read input root file %s\n", finname);
    
    file2->GetObject( dname, dir0);
    printf("got input root directory %s\n", dname);
    
    //if ( iCal == 4 ) continue; // skip - noise Xtal 39
    //if ( iCal == 12 ) continue; // skip - bad noise
    //if ( iCal == 14 ) continue; // skip - calo off
    //if ( iCal == 17 ) continue; // skip - noise Xtal 2,47,39,31,52,53
    //if ( iCal == 19 ) continue; // skip - noise Xtal 46,20,31,32
    //if ( iCal == 23 ) continue; // skip - noise Xtal 1,29,31,51,52
    //if ( iCal == 24 ) continue; // skip - noise Xtal 43

    for (int iSeg = 0; iSeg < 53; iSeg++) {
      
      sprintf( hname,"%s_%i_%i", cHist, iCal, iSeg);
      dir0->GetObject( hname, htmp);

      printf("got run %i, calo %i, xtal %i, histogram %s ...n", iRun, iCal, iSeg, hname);
      
      h1[iSeg]->Add( htmp, 1.0);
      printf(" done add\n");
    } // loop over xtals

    printf(" done calo %i\n", iCal);
    file2->Close();    
  }
  
  // open sum calo output file
  sprintf( foutname, "Data/sum%05i_coarse_testped.root", iRun);
  file0 = TFile::Open( foutname, "UPDATE");
  if ( file0->IsOpen() ) printf("update output root file %s\n", foutname);

  for (int iSeg = 0; iSeg < 53; iSeg++) {
    h1[iSeg]->Write(); 
    file0->Print();
    printf("wrote histogram %s into file %s\n", h1[iSeg]->GetTitle(), foutname);
  }
  file0->Close();
  
  return 1;
}

int sumCalo2D(int iRun, char* cHist){

  sprintf( foutname, "Data/sum%05i_coarse_testped.root", iRun);
  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("created output root file %s\n", foutname);

  //sprintf( hname, "qHist1D_sig_%i_sum", 1);
  sprintf( hname,"%s_%i_sum", cHist, 1);
  file1->GetObject( hname, sSum);
  //sprintf( hsname, "qHist1D_sig_calo_sum");
  //sprintf( hsname,"%s_sum", cHist); // includes bad Calos
  sprintf( hsname,"%s_sum_noBadCalo", cHist); // excludes bad Calos
  sSum->SetName(hsname);
  sSum->SetTitle(hsname);
  printf("got run %i, xtal %i, histogram %s\n", iRun, 0, hname);

  for (int iCal = 2; iCal <= 24; iCal++) {

    //if ( iCal == 4 ) continue; // skip - noise Xtal 39
    //if ( iCal == 12 ) continue; // skip - bad noise
    //if ( iCal == 14 ) continue; // skip - calo off
    //if ( iCal == 17 ) continue; // skip - noise Xtal 2,47,39,31,52,53
    //if ( iCal == 19 ) continue; // skip - noise Xtal 46,20,31,32
    //if ( iCal == 23 ) continue; // skip - noise Xtal 1,29,31,51,52
    //if ( iCal == 24 ) continue; // skip - noise Xtal 43

    sprintf( hname, "%s_%i_sum", cHist, iCal);
    file1->GetObject( hname, s2[iCal]);
    sSum->Add( s2[iCal], 1.0);
    printf("got run %i, xtal %i, histogram %s\n", iRun, iCal, hname);

  }

  sSum->Write(); 
  file1->Print(); 
  printf("wrote histogram %s into file %s\n", hsname, foutname);
  file1->Close();

  return 1;
}

int NCalo = 24;
double GCalonum[24], GCBOmeanX[24], dGCalonum[24], dGCBOmeanX[24];

int goFitCalo( int iS, int iF){

  for (int ic = 1; ic <= 24; ic++){
    ICalo = ic;
    fitCalo( "", 0, iS, iF); 
  }

  gACBOmeanX = new TGraphErrors(NCalo, GCalonum, GCBOmeanX, dGCalonum, dGCBOmeanX);
  gACBOmeanX->SetName("gACBOmeanX");
  gACBOmeanX->SetTitle("Q-Method CBO ampl. mean Y versus calorimeter index");
  gACBOmeanX->GetYaxis()->SetTitle("Q-Method CBO ampl. mean Y");
  gACBOmeanX->GetXaxis()->SetTitle("calorimeter index");
  gACBOmeanX->GetYaxis()->SetTitleOffset(0.5);
  gACBOmeanX->GetYaxis()->SetTitleSize(0.07);
  gACBOmeanX->GetXaxis()->SetTitleOffset(0.5);
  gACBOmeanX->GetXaxis()->SetTitleSize(0.07);
  gACBOmeanX->SetMarkerStyle(8);
  gACBOmeanX->SetMarkerSize(1.1);
  gACBOmeanX->SetMarkerColor(kMagenta);
  gACBOmeanX->SetLineWidth(2);
  gACBOmeanX->SetLineColor(kMagenta);

  c1->Clear();
  gACBOmeanX->Draw("AP");

  return 1;
}

int gorunsum();

int gorunsum(){
  
  char fname[128];
  int runnos[62] =
    {15921, 15922, 15923, 15924, 15925, 15926, 15927, 15928, 15929, 
     15930, 15931, 15932, 15933, 15934, 15935, 15936, 15937, 15938, 15939,
     15940, 15941, 15942, 15943, 15944, 15945, 15946, 15947, 15948, 15949,
     15950, 15951, 15952, 15953, 15954, 15955, 15957, 15958, 15959,
     15960, 15961,        15963, 15964,               15967, 15968, 15969,
     15970,        15972, 15973, 15974, 15975,        15977, 15978,
     15980, 15981, 15982,        15984, 15985, 15986, 15987, 15988, 15989,
     15991};
  
  for (I = 0; I < 62; I++){
    sprintf( fname, "v4w%it%02iprd%s_%05i", iwind, ithreshold, aname, runnos[I] );
    goCaloSum(fname);
  }

  return 0;
}

int fitCalo( char* cFit, int iC, int iStart, int iEnd){

  double startTime = iStart;
  double endTime = iEnd;

  int runnos[62] =
    {15921, 15922, 15923, 15924, 15925, 15926, 15927, 15928, 15929, 
     15930, 15931, 15932, 15933, 15934, 15935, 15936, 15937, 15938, 15939,
     15940, 15941, 15942, 15943, 15944, 15945, 15946, 15947, 15948, 15949,
     15950, 15951, 15952, 15953, 15954, 15955, 15957, 15958, 15959,
     15960, 15961,        15963, 15964,               15967, 15968, 15969,
     15970,        15972, 15973, 15974, 15975,        15977, 15978,
     15980, 15981, 15982,        15984, 15985, 15986, 15987, 15988, 15989,
     15991};

  fpedestaldrift = TFile::Open( fnamepedestal);
  fmuonloss = TFile::Open( fnamemuonloss);

  // single calo histogram fit
  if ( strcmp( "singlecalo", cFit) == 0 ) {
    icalo = iC-1;
    sprintf( hname, "qHist1D_sig_%i_%i", iC, iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
  }

  // opposite pair calos histogram fit
  if ( strcmp( "paircalo", cFit) == 0 ) {
    icalo = iC-1;
    sprintf( hname, "qHist1D_sig_pair_%i_%i", iC, iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
  }

  // hot calosum histogram fit
  if ( strcmp( "hot", cFit) == 0 ) {
    Ifit = 0;
    icalo = iC;
    GPar_QM[Ifit] = Ifit;
    dGPar_QM[Ifit] = 0;
    if (iC == 0) sprintf( hname, "qHist1D_sig_sum_hot_%i", iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }

  // hot calosum histogram fit
  if ( strcmp( "cold", cFit) == 0 ) {
    Ifit = 0;
    icalo = iC;
    GPar_QM[Ifit] = Ifit;
    dGPar_QM[Ifit] = 0;
    if (iC == 0) sprintf( hname, "qHist1D_sig_sum_cold_%i", iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }

  // calosum histogram fit
  if ( strcmp( "full", cFit) == 0 ) {
    Ifit = 0;
    icalo = iC;
    GPar_QM[Ifit] = Ifit;
    dGPar_QM[Ifit] = 0;

    if (iC == 0) sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2); // histogram name 
    if ( strcmp( "core4", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_core4"); // test, drops the top row and bottom row 
    if ( strcmp( "core2", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_core2"); // test, drops the top two rows and bottom two rows 
    if ( strcmp( "puA", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_puA_sum_%i", iswrb/2);
    if ( strcmp( "puLH", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_puLH_sum_%i", iswrb/2);
    if ( strcmp( "top", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_top"); // calo top half
    if ( strcmp( "bot", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_bot"); // calo botttom half

    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }
  
  
  if ( strcmp( "near", cFit) == 0 ) {
    Ifit = 0;
    icalo = iC;
    GPar_QM[Ifit] = Ifit;
    dGPar_QM[Ifit] = 0;

    if (iC == 0) sprintf( hname, "qHist1D_sig_sum_near_%i", iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }
  
  if ( strcmp( "far", cFit) == 0 ) {
    Ifit = 0;
    icalo = iC;
    GPar_QM[Ifit] = Ifit;
    dGPar_QM[Ifit] = 0;

    if (iC == 0) sprintf( hname, "qHist1D_sig_sum_far_%i", iswrb/2); // histogram name 
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }
  
  // do a calorimeter scan, slice scan, etc
  int Nf;
  if (iC > 0) {
    Nf = iC+1;
  } else {
    Nf = Nfit; // full scan
  }

  if ( strcmp( "runscan", cFit) == 0 ) {
      for (Ifit = 0; Ifit < Nfit; Ifit++){
	
	// do a start time scan
	double startTime = iStart;
	double endTime = iEnd;
	GPar_QM[Ifit] = runnos[Ifit];
	dGPar_QM[Ifit] = 0.0;
	sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
	sprintf( fname, "noxtaldata/v4w%it%02iprd%s_%05i.root", iwind, ithreshold, aname, runnos[Ifit]); 
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      }
    }
  
  if ( strcmp( "starttimescan", cFit) == 0 ) {
    for (Ifit = 0; Ifit < Nfstarttimescan; Ifit++){
      
	// do a start time scan
	double startTime = 30000. + Ifit*4370./4.; // in tenths of T_a
	double endTime = 2.155e5;
	GPar_QM[Ifit] = startTime;
	dGPar_QM[Ifit] = 0.0;
	if ( strcmp( "all", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
	if ( strcmp( "near", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_near_%i", iswrb/2);
	if ( strcmp( "far", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_far_%i", iswrb/2);
	if ( strcmp( "hot", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_hot_%i", iswrb/2);
	if ( strcmp( "cold", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_cold_%i", iswrb/2);
	if ( strcmp( "top", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_top");
	if ( strcmp( "bot", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_bot");
	//sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd1e9shadowtest9dayp2.root", iwind, inoise, ithreshold); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%sshadowtestendgamefraction.root", iwind, inoise, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%swithgaincorrection.root", iwind, inoise, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%swithgaincorrection-gap%02i.root", iwind, inoise, ithreshold, aname, igap);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      }
    }

  if ( strcmp( "tauscan", cFit) == 0 ) {
    for (Ifit = 0; Ifit <= 10; Ifit++){
      
	// do a start time scan
        fix_gammatau = 1;
        set_gammatau = 64200+40.*Ifit;
	GPar_QM[Ifit] = set_gammatau;
	dGPar_QM[Ifit] = 0.0;
	if ( strcmp( "all", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
	if ( strcmp( "near", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_near_%i", iswrb/2);
	if ( strcmp( "far", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_far_%i", iswrb/2);
	if ( strcmp( "hot", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_hot_%i", iswrb/2);
	if ( strcmp( "cold", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_sum_cold_%i", iswrb/2);
	if ( strcmp( "top", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_top");
	if ( strcmp( "bot", starttimescantype ) == 0 ) sprintf( hname, "qHist1D_sig_bot");
	//sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd1e9shadowtest9dayp2.root", iwind, inoise, ithreshold); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%sshadowtestendgamefraction.root", iwind, inoise, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%swithgaincorrection.root", iwind, inoise, ithreshold, aname); 
	//sprintf( fname, "noxtaldata/v4w%in%it%02iprd%swithgaincorrection-gap%02i.root", iwind, inoise, ithreshold, aname, igap);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      }
    }

  // no correction _0, _2 and puLH _0 have too small error bars?
  if ( strcmp( "puscan", cFit) == 0 ) {
    int Ithr[3] = {15, 15, 15}; 
    for (int Ith = 0; Ith <= 0; Ith++){ // loop over threshold settings
      ithreshold = Ithr[Ith];
      iswrb = 1; 
      for (int Irb = 0; Irb <= 2; Irb++){ // loop over rebinning factors

	Ifit = 0 + Irb*4 + Ith*12;
	GPar_QM[Ifit] = 1 + Irb*4 + Ith*12;
	dGPar_QM[Ifit] = 0;

        // get file of chosen dataset, threshold setting, pedestal window, noise cut w/o pulse rejection corrections
	sprintf( fname, "noxtaldata/v9.16-w%in%it%02iprd%swithgaincorrection-gap%02i-60hr.root", iwind, inoise, ithreshold, "1e9", igap); // Feb 8, 2019 

        // get no pu correction  histogram 
	sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime); 

	Ifit = 1 + Irb*4 + Ith*12;
	GPar_QM[Ifit] = 2 + Irb*4 + Ith*12;
	dGPar_QM[Ifit] = 0;

        // get adjacent shadow pileup method histogram 
	sprintf( hname, "qHist1D_sig_puA_sum_%i", iswrb/2);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);

	Ifit = 2 + Irb*4 + Ith*12;
	GPar_QM[Ifit] = 3 + Irb*4  + Ith*12;
	dGPar_QM[Ifit] = 0;

        // get left-right shadow pileup method histogram 
	sprintf( hname, "qHist1D_sig_puLH_sum_%i", iswrb/2);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);

	Ifit = 3 + Irb*4 + Ith*12;
	GPar_QM[Ifit] = 4 + Irb*4 + Ith*12;
	dGPar_QM[Ifit] = 0;

        // get file of chosen dataset, threshold setting, pedestal window, noise cut with pulse rejection corrections
	sprintf( fname, "noxtaldata/v9.16-w%in%it%02iprd%swithgaincorrection-gap%02i-60hr.root", iwind, inoise, ithreshold, "1e0", igap); // Feb 8, 2019 

        // get pulse rejection pileup method histogram 
	sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);

	iswrb *= 2; // iswrb goes 1,2,4, iswrb/2 0,1,2

      }
    }
  }
  
  if ( strcmp( "noisescan", cFit) == 0 ) {
    for (Ifit = iC; Ifit < Nf; Ifit++){
      
      int in = 2*Ifit + 4; // noise rejection multiplier
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( hname, "qHist1D_sig_%i_%i", 23, iswrb/2);
      sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2);
      sprintf( fname, "noxtaldata/v4w%in%it%02iprd%s-60hr.root", iwind, in, ithreshold, aname); 
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    }
  }

  if ( strcmp( "top", cFit) == 0 ) {

    icalo = iC;
    sprintf( hname, "qHist2D_sig_top");
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }

  if ( strcmp( "bot", cFit) == 0 ) {

    icalo = iC;
    sprintf( hname, "qHist2D_sig_bot");
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
    return 0;
  }

  if ( strcmp( "horizontalslice", cFit) == 0 ) {

    icalo = iC;
    sprintf( hname, "qHist2D_sig_hslice_%i", iC);
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);    
  }

  if ( strcmp( "verticalslice", cFit) == 0 ) {

    icalo = iC;
    sprintf( hname, "qHist2D_sig_vslice_%i", iC);
    viewCalo( cFit, fname, hname, iswrb, startTime, endTime);    
  }
  

  if ( ( strcmp( "sequencescan", cFit) == 0 ) || ( strcmp( "pairscan", cFit) == 0 ) || ( strcmp( "caloscan", cFit) == 0 ) || ( strcmp( "verticalscan", cFit) == 0 ) || ( strcmp( "horizontalscan", cFit) == 0 ) ) {
    
    if ( strcmp( "pairscan", cFit) == 0 ) Nf = 12;
    if ( strcmp( "caloscan", cFit) == 0 ) Nf = 24;
    if ( strcmp( "verticalscan", cFit) == 0 ) Nf = 9;
    if ( strcmp( "horizontalscan", cFit) == 0 ) Nf = 6;
    if ( strcmp( "sequencescan", cFit) == 0 ) Nf = 8;

    for (Ifit = 0; Ifit < Nf; Ifit++){
      
      GPar_QM[Ifit] = Ifit+1;
      dGPar_QM[Ifit] = 0;
      if ( strcmp( "pairscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_pair_%i_%i", Ifit+1, iswrb/2);
      if ( strcmp( "caloscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_%i_%i", Ifit+1, iswrb/2);
      if ( strcmp( "verticalscan", cFit) == 0 ) sprintf( hname, "qHist2D_sig_vslice_%i", Ifit);
      if ( strcmp( "horizontalscan", cFit) == 0 ) sprintf( hname, "qHist2D_sig_hslice_%i", Ifit);
      if ( strcmp( "sequencescan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_seqnum_%i", Ifit);
      
      icalo = Ifit;
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);    
      if ( strcmp( "caloscan", cFit) == 0 )  caloslowtermsf[Ifit] = (TF1*)slowtermsf->Clone();
      
    }
  }
   
    if ( strcmp( "windowscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2); // histogram name 
    if ( strcmp( "thresholdscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2); // histogram name 
    if ( strcmp( "coldthresholdscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_sum_cold_%i", iswrb/2); // histogram name 
    if ( strcmp( "hotthresholdscan", cFit) == 0 ) sprintf( hname, "qHist1D_sig_sum_hot_%i", iswrb/2); // histogram name 

    if ( strcmp( "thresholdscan", cFit) == 0 ) {

      int Nthr, thresarray[11] = {12,15,18,24,30,35,40,45,50,60,70};
      int thresarray_hk[11] = {12,15,18,24,30,40,50,60,1000,1000,1000};
 
     if ( strcmp( "endgame", dataname) == 0 ) Nthr = 11;
      int thresenrgy_eg[11] = {271,339,407,543,708,791,942,1017,1130,1356,1582};
      if ( strcmp( "9day", dataname) == 0 ) Nthr = 11;
      int thresenrgy_9d[11] = {212,265,417,556,708,619,707,795,884,1060,1237};
      if ( strcmp( "60hr", dataname) == 0 ) Nthr = 5;
      int thresenrgy_60[5] = {283,354,425,567,709};
      if ( strcmp( "hk", dataname) == 0 ) Nthr = 8;
      int thresenrgy_hk[11] = {212,265,318,424,515,687,859,1031,10000,10000,10000};
      if ( strcmp( "lk", dataname) == 0 ) Nthr = 5;
      int thresenrgy_lk[11] = {212,265,417,556,708,619,707,795,884,1060,1237};
      

      for (int ithr = 0; ithr < Nthr; ithr++){
	ithreshold = thresarray[ithr];
	if ( strcmp( "hk", dataname) == 0 ) ithreshold = thresarray_hk[ithr];
	Ifit = ithr;
	GPar_QM[Ifit] = ithreshold;
	if ( strcmp( "endgame", dataname) == 0 ) GPar_QM[Ifit] = thresenrgy_eg[ithr];
	if ( strcmp( "9day", dataname) == 0 ) GPar_QM[Ifit] = thresenrgy_9d[ithr];
	if ( strcmp( "60hr", dataname) == 0 ) GPar_QM[Ifit] = thresenrgy_60[ithr];
	if ( strcmp( "hk", dataname) == 0 ) GPar_QM[Ifit] = thresenrgy_hk[ithr];
	if ( strcmp( "lk", dataname) == 0 ) GPar_QM[Ifit] = thresenrgy_lk[ithr];
	dGPar_QM[Ifit] = 0;

        // different windows are in different files
	if ( strcmp( "endgame", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "endgame", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
       
	if ( strcmp( "9day", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "9day", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "60hr", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-60hr-pufix-gold-withpederr-gainbugfix-newenergyscale.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "60hr", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-60hr-pufix-gold-withpederr-gainbugfix-newenergyscale.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "hk", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.21-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "hk", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.21-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "lk", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-lk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "lk", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-lk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2); 
	if ( strcmp( "core4", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_core4"); // test, drops the top row and bottom row 
	if ( strcmp( "core2", thresholdscantype ) == 0 ) sprintf( hname, "qHist1D_sig_core2"); // test, drops the top two rows and bottom two rows 
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      }
    }

    if ( strcmp( "windowscan", cFit) == 0 ) {
      
      int nwnd, windwarray[5] = {2,4,6,8,10};
      if ( strcmp( "endgame", dataname) == 0 ) nwnd = 5;
      if ( strcmp( "9day", dataname) == 0 ) nwnd = 5;
      if ( strcmp( "hk", dataname) == 0 ) nwnd = 5;
      if ( strcmp( "lk", dataname) == 0 ) nwnd = 5;
      if ( strcmp( "60hr", dataname) == 0 ) nwnd = 5;
      for (int iwnd = 0; iwnd < nwnd; iwnd++){
	iwind = windwarray[iwnd];
	Ifit = iwnd;
	GPar_QM[Ifit] = iwind;
	dGPar_QM[Ifit] = 0;

        // different windows are in different files
        // ithreshold=18; // only threshold with window scan
	if ( strcmp( "endgame", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "endgame", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
        // ithreshold=15; // only threshold with window scan
	if ( strcmp( "9day", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "9day", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "hk", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "hk", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "lk", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-lk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "lk", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-lk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	if ( strcmp( "60hr", dataname) == 0 ) sprintf( fname, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-60hr-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);
	if ( strcmp( "60hr", dataname) == 0 ) sprintf( fnamepedestal, "noxtaldata/v9.17-w%in%it%02iprd%sinfillwithstdp-OOFyes-gap%02i-60hr-pufix-withpederr-gainbugfix-newenergyscale-SRC.root", iwind, inoise, ithreshold, aname, igap);

	sprintf( hname, "qHist1D_sig_sum_%i", iswrb/2); 
	viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      }
    }

    if ( ( strcmp( "coldthresholdscan", cFit) == 0 ) || ( strcmp( "hotthresholdscan", cFit) == 0 ) ) {
      Ifit = 0;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 8, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);    
      Ifit = 1;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 12, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      Ifit = 2;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 16, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      Ifit = 3;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 20, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      Ifit = 4;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 24, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      Ifit = 5;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 28, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      Ifit = 6;
      GPar_QM[Ifit] = Ifit;
      dGPar_QM[Ifit] = 0;
      sprintf( fname, "noxtaldata/v4w%it%02iprd%s-60hr.root", iwind, 32, aname); // no point reject
      viewCalo( cFit, fname, hname, iswrb, startTime, endTime);
      
    }
    //*/
    
    printf("iC %i\n", iC);
    if (iC >= 1) return 1;
    if ( ( strcmp( "horizontalslice", cFit) == 0 ) || ( strcmp( "verticalslice", cFit) == 0 ) ) return 1;
    
  gOmega_afit = new TGraphErrors(Nfit, GPar_QM, omega_aF_QM, dGPar_QM, domega_aF_QM);
  gOmega_afit->SetName("gOmega_afit");
  gOmega_afit->SetTitle("Q-Method blind R[ppm]");
  gOmega_afit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gOmega_afit->GetXaxis()->SetTitle("run group");
  //gOmega_afit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gOmega_afit->GetYaxis()->SetTitle("blind R[ppm]");
  gOmega_afit->GetYaxis()->SetTitleOffset(0.5);
  gOmega_afit->GetYaxis()->SetTitleSize(0.07);
  gOmega_afit->GetXaxis()->SetTitleOffset(0.5);
  gOmega_afit->GetXaxis()->SetTitleSize(0.07);
  gOmega_afit->SetMarkerStyle(8);
  gOmega_afit->SetMarkerSize(1.0);
  gOmega_afit->SetMarkerColor(kOrange+1);
  gOmega_afit->SetLineWidth(2);
  gOmega_afit->SetLineColor(kOrange+1);

  gTaufit = new TGraphErrors(Nfit, GPar_QM, tauF_QM, dGPar_QM, dtauF_QM);
  gTaufit->SetName("gTaufit");
  gTaufit->SetTitle("Q-Method tau");
  gTaufit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gTaufit->GetXaxis()->SetTitle("run group");
  //gTaufit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gTaufit->GetYaxis()->SetTitle("tau ( ns )");
  gTaufit->GetYaxis()->SetTitleOffset(0.5);
  gTaufit->GetYaxis()->SetTitleSize(0.07);
  gTaufit->GetXaxis()->SetTitleOffset(0.5);
  gTaufit->GetXaxis()->SetTitleSize(0.07);
  gTaufit->SetMarkerStyle(8);
  gTaufit->SetMarkerSize(1.0);
  gTaufit->SetMarkerColor(kRed);
  gTaufit->SetLineWidth(2);
  gTaufit->SetLineColor(kRed);

  gAfit = new TGraphErrors(Nfit, GPar_QM, AF_QM, dGPar_QM, dAF_QM);
  gAfit->SetName("gAfit");
  gAfit->SetTitle("Q-Method amplitude");
  gAfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gAfit->GetXaxis()->SetTitle("run group");
  //gAfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gAfit->GetYaxis()->SetTitle("amplitude");
  gAfit->GetYaxis()->SetTitleOffset(0.5);
  gAfit->GetYaxis()->SetTitleSize(0.07);
  gAfit->GetXaxis()->SetTitleOffset(0.5);
  gAfit->GetXaxis()->SetTitleSize(0.07);
  gAfit->SetMarkerStyle(8);
  gAfit->SetMarkerSize(1.0);
  gAfit->SetMarkerColor(kAzure);
  gAfit->SetLineWidth(2);
  gAfit->SetLineColor(kAzure);

  gArelaxfit = new TGraphErrors(Nfit, GPar_QM, AFrelax_QM, dGPar_QM, dAFrelax_QM);
  gArelaxfit->SetName("gArelaxfit");
  gArelaxfit->SetTitle("Asymmetry relaxation");
  gArelaxfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gArelaxfit->GetXaxis()->SetTitle("run group");
  //gArelax->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gArelaxfit->GetYaxis()->SetTitle("asymmetry relaxation");
  gArelaxfit->GetYaxis()->SetTitleOffset(0.5);
  gArelaxfit->GetYaxis()->SetTitleSize(0.07);
  gArelaxfit->GetXaxis()->SetTitleOffset(0.5);
  gArelaxfit->GetXaxis()->SetTitleSize(0.07);
  gArelaxfit->SetMarkerStyle(8);
  gArelaxfit->SetMarkerSize(1.0);
  gArelaxfit->SetMarkerColor(kAzure);
  gArelaxfit->SetLineWidth(2);
  gArelaxfit->SetLineColor(kAzure);

  gPhifit = new TGraphErrors(Nfit, GPar_QM, phaseF_QM, dGPar_QM, dphaseF_QM);
  gPhifit->SetName("gPhifit");
  gPhifit->SetTitle("Q-Method phase");
  gPhifit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gPhifit->GetXaxis()->SetTitle("run group");
  //gPhifit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gPhifit->GetYaxis()->SetTitle("phase (rad)");
  gPhifit->GetYaxis()->SetTitleOffset(0.5);
  gPhifit->GetYaxis()->SetTitleSize(0.07);
  gPhifit->GetXaxis()->SetTitleOffset(0.5);
  gPhifit->GetXaxis()->SetTitleSize(0.07);
  gPhifit->SetMarkerStyle(8);
  gPhifit->SetMarkerSize(1.0);
  gPhifit->SetMarkerColor(kViolet);
  gPhifit->SetLineWidth(2);
  gPhifit->SetLineColor(kViolet);

  gNfit = new TGraphErrors(Nfit, GPar_QM, NF_QM, dGPar_QM, dNF_QM);
  gNfit->SetName("gNfit");
  gNfit->SetTitle("Q-Method normalization");
  gNfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gNfit->GetXaxis()->SetTitle("run group");
  //gNfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gNfit->GetYaxis()->SetTitle("normalization");
  gNfit->GetYaxis()->SetTitleOffset(0.5);
  gNfit->GetYaxis()->SetTitleSize(0.07);
  gNfit->GetXaxis()->SetTitleOffset(0.5);
  gNfit->GetXaxis()->SetTitleSize(0.07);
  gNfit->SetMarkerStyle(8);
  gNfit->SetMarkerSize(1.0);
  gNfit->SetMarkerColor(kGray+2);
  gNfit->SetLineWidth(2);
  gNfit->SetLineColor(kGray+2);

  gBfit = new TGraphErrors(Nfit, GPar_QM, BF_QM, dGPar_QM, dBF_QM);
  gBfit->SetName("gBfit");
  gBfit->SetTitle("Q-Method background");
  gBfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gBfit->GetXaxis()->SetTitle("run group");
  //gBfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gBfit->GetYaxis()->SetTitle("background");
  gBfit->GetYaxis()->SetTitleOffset(0.5);
  gBfit->GetYaxis()->SetTitleSize(0.07);
  gBfit->GetXaxis()->SetTitleOffset(0.5);
  gBfit->GetXaxis()->SetTitleSize(0.07);
  gBfit->SetMarkerStyle(8);
  gBfit->SetMarkerSize(1.0);
  gBfit->SetMarkerColor(kTeal);
  gBfit->SetLineWidth(2);
  gBfit->SetLineColor(kTeal);

  gACBOfit = new TGraphErrors(Nfit, GPar_QM, ACBOF_QM, dGPar_QM, dACBOF_QM);
  gACBOfit->SetName("gACBOfit");
  gACBOfit->SetTitle("Q-Method CBO ampl");
  gACBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gACBOfit->GetXaxis()->SetTitle("run group");
  //gACBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gACBOfit->GetYaxis()->SetTitle("CBO ampl");
  gACBOfit->GetYaxis()->SetTitleOffset(0.5);
  gACBOfit->GetYaxis()->SetTitleSize(0.07);
  gACBOfit->GetXaxis()->SetTitleOffset(0.5);
  gACBOfit->GetXaxis()->SetTitleSize(0.07);
  gACBOfit->SetMarkerStyle(8);
  gACBOfit->SetMarkerSize(1.1);
  gACBOfit->SetMarkerColor(kMagenta);
  gACBOfit->SetLineWidth(2);
  gACBOfit->SetLineColor(kMagenta);

  gAVWfit = new TGraphErrors(Nfit, GPar_QM, AVWF_QM, dGPar_QM, dAVWF_QM);
  gAVWfit->SetName("gAVWfit");
  gAVWfit->SetTitle("Q-Method VW ampl");
  gAVWfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gAVWfit->GetXaxis()->SetTitle("run group");
  //gAVWfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gAVWfit->GetYaxis()->SetTitle("VW ampl");
  gAVWfit->GetYaxis()->SetTitleOffset(0.5);
  gAVWfit->GetYaxis()->SetTitleSize(0.07);
  gAVWfit->GetXaxis()->SetTitleOffset(0.5);
  gAVWfit->GetXaxis()->SetTitleSize(0.07);
  gAVWfit->SetMarkerStyle(8);
  gAVWfit->SetMarkerSize(1.1);
  gAVWfit->SetMarkerColor(kMagenta);
  gAVWfit->SetLineWidth(2);
  gAVWfit->SetLineColor(kMagenta);

  gomegaVWfit = new TGraphErrors(Nfit, GPar_QM, omegaVWF_QM, dGPar_QM, domegaVWF_QM);
  gomegaVWfit->SetName("gomegaVWfit");
  gomegaVWfit->SetTitle("Q-Method VW omega");
  gomegaVWfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gomegaVWfit->GetXaxis()->SetTitle("run group");
  //gomegaVWfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gomegaVWfit->GetYaxis()->SetTitle("VW ampl");
  gomegaVWfit->GetYaxis()->SetTitleOffset(0.5);
  gomegaVWfit->GetYaxis()->SetTitleSize(0.07);
  gomegaVWfit->GetXaxis()->SetTitleOffset(0.5);
  gomegaVWfit->GetXaxis()->SetTitleSize(0.07);
  gomegaVWfit->SetMarkerStyle(8);
  gomegaVWfit->SetMarkerSize(1.1);
  gomegaVWfit->SetMarkerColor(kMagenta);
  gomegaVWfit->SetLineWidth(2);
  gomegaVWfit->SetLineColor(kMagenta);

  gtauVWfit = new TGraphErrors(Nfit, GPar_QM, tauVWF_QM, dGPar_QM, dtauVWF_QM);
  gtauVWfit->SetName("gtauVWfit");
  gtauVWfit->SetTitle("Q-Method VW tau");
  gtauVWfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gtauVWfit->GetXaxis()->SetTitle("run group");
  //gtauVWfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gtauVWfit->GetYaxis()->SetTitle("VW ampl");
  gtauVWfit->GetYaxis()->SetTitleOffset(0.5);
  gtauVWfit->GetYaxis()->SetTitleSize(0.07);
  gtauVWfit->GetXaxis()->SetTitleOffset(0.5);
  gtauVWfit->GetXaxis()->SetTitleSize(0.07);
  gtauVWfit->SetMarkerStyle(8);
  gtauVWfit->SetMarkerSize(1.1);
  gtauVWfit->SetMarkerColor(kMagenta);
  gtauVWfit->SetLineWidth(2);
  gtauVWfit->SetLineColor(kMagenta);

  gphaseVWfit = new TGraphErrors(Nfit, GPar_QM, phaseVWF_QM, dGPar_QM, dphaseVWF_QM);
  gphaseVWfit->SetName("gphaseVWfit");
  gphaseVWfit->SetTitle("Q-Method VW phase");
  gphaseVWfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gphaseVWfit->GetXaxis()->SetTitle("run group");
  //gphaseVWfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gphaseVWfit->GetYaxis()->SetTitle("VW phase");
  gphaseVWfit->GetYaxis()->SetTitleOffset(0.5);
  gphaseVWfit->GetYaxis()->SetTitleSize(0.07);
  gphaseVWfit->GetXaxis()->SetTitleOffset(0.5);
  gphaseVWfit->GetXaxis()->SetTitleSize(0.07);
  gphaseVWfit->SetMarkerStyle(8);
  gphaseVWfit->SetMarkerSize(1.1);
  gphaseVWfit->SetMarkerColor(kMagenta);
  gphaseVWfit->SetLineWidth(2);
  gphaseVWfit->SetLineColor(kMagenta);

  gAVW2fit = new TGraphErrors(Nfit, GPar_QM, AVW2F_QM, dGPar_QM, dAVW2F_QM);
  gAVW2fit->SetName("gAVW2fit");
  gAVW2fit->SetTitle("Q-Method VW2 ampl");
  gAVW2fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gAVW2fit->GetXaxis()->SetTitle("run group");
  //gAVW2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gAVW2fit->GetYaxis()->SetTitle("VW2 ampl");
  gAVW2fit->GetYaxis()->SetTitleOffset(0.5);
  gAVW2fit->GetYaxis()->SetTitleSize(0.07);
  gAVW2fit->GetXaxis()->SetTitleOffset(0.5);
  gAVW2fit->GetXaxis()->SetTitleSize(0.07);
  gAVW2fit->SetMarkerStyle(8);
  gAVW2fit->SetMarkerSize(1.1);
  gAVW2fit->SetMarkerColor(kMagenta);
  gAVW2fit->SetLineWidth(2);
  gAVW2fit->SetLineColor(kMagenta);

  gomegaVW2fit = new TGraphErrors(Nfit, GPar_QM, omegaVW2F_QM, dGPar_QM, domegaVW2F_QM);
  gomegaVW2fit->SetName("gomegaVW2fit");
  gomegaVW2fit->SetTitle("Q-Method VW2 omega");
  gomegaVW2fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gomegaVW2fit->GetXaxis()->SetTitle("run group");
  //gomegaVW2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gomegaVW2fit->GetYaxis()->SetTitle("VW2 ampl");
  gomegaVW2fit->GetYaxis()->SetTitleOffset(0.5);
  gomegaVW2fit->GetYaxis()->SetTitleSize(0.07);
  gomegaVW2fit->GetXaxis()->SetTitleOffset(0.5);
  gomegaVW2fit->GetXaxis()->SetTitleSize(0.07);
  gomegaVW2fit->SetMarkerStyle(8);
  gomegaVW2fit->SetMarkerSize(1.1);
  gomegaVW2fit->SetMarkerColor(kMagenta);
  gomegaVW2fit->SetLineWidth(2);
  gomegaVW2fit->SetLineColor(kMagenta);

  gtauVW2fit = new TGraphErrors(Nfit, GPar_QM, tauVW2F_QM, dGPar_QM, dtauVW2F_QM);
  gtauVW2fit->SetName("gtauVW2fit");
  gtauVW2fit->SetTitle("Q-Method VW2 tau");
  gtauVW2fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gtauVW2fit->GetXaxis()->SetTitle("run group");
  //gtauVW2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gtauVW2fit->GetYaxis()->SetTitle("VW2 ampl");
  gtauVW2fit->GetYaxis()->SetTitleOffset(0.5);
  gtauVW2fit->GetYaxis()->SetTitleSize(0.07);
  gtauVW2fit->GetXaxis()->SetTitleOffset(0.5);
  gtauVW2fit->GetXaxis()->SetTitleSize(0.07);
  gtauVW2fit->SetMarkerStyle(8);
  gtauVW2fit->SetMarkerSize(1.1);
  gtauVW2fit->SetMarkerColor(kMagenta);
  gtauVW2fit->SetLineWidth(2);
  gtauVW2fit->SetLineColor(kMagenta);

  gphaseVW2fit = new TGraphErrors(Nfit, GPar_QM, phaseVW2F_QM, dGPar_QM, dphaseVW2F_QM);
  gphaseVW2fit->SetName("gphaseVW2fit");
  gphaseVW2fit->SetTitle("Q-Method VW2 phase");
  gphaseVW2fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gphaseVW2fit->GetXaxis()->SetTitle("run group");
  //gphaseVW2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gphaseVW2fit->GetYaxis()->SetTitle("VW2 phase");
  gphaseVW2fit->GetYaxis()->SetTitleOffset(0.5);
  gphaseVW2fit->GetYaxis()->SetTitleSize(0.07);
  gphaseVW2fit->GetXaxis()->SetTitleOffset(0.5);
  gphaseVW2fit->GetXaxis()->SetTitleSize(0.07);
  gphaseVW2fit->SetMarkerStyle(8);
  gphaseVW2fit->SetMarkerSize(1.1);
  gphaseVW2fit->SetMarkerColor(kMagenta);
  gphaseVW2fit->SetLineWidth(2);
  gphaseVW2fit->SetLineColor(kMagenta);

  gAGAINfit = new TGraphErrors(Nfit, GPar_QM, AGAINF_QM, dGPar_QM, dAGAINF_QM);
  gAGAINfit->SetName("gAGAINfit");
  gAGAINfit->SetTitle("Q-Method GAIN ampl");
  gAGAINfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gAGAINfit->GetXaxis()->SetTitle("run group");
  //gAGAINfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gAGAINfit->GetYaxis()->SetTitle("GAIN ampl");
  gAGAINfit->GetYaxis()->SetTitleOffset(0.5);
  gAGAINfit->GetYaxis()->SetTitleSize(0.07);
  gAGAINfit->GetXaxis()->SetTitleOffset(0.5);
  gAGAINfit->GetXaxis()->SetTitleSize(0.07);
  gAGAINfit->SetMarkerStyle(8);
  gAGAINfit->SetMarkerSize(1.0);
  gAGAINfit->SetMarkerColor(kTeal);
  gAGAINfit->SetLineWidth(2);
  gAGAINfit->SetLineColor(kTeal);

  gtauGAINfit = new TGraphErrors(Nfit, GPar_QM, tauGAINF_QM, dGPar_QM, dtauGAINF_QM);
  gtauGAINfit->SetName("gtauGAINfit");
  gtauGAINfit->SetTitle("Q-Method GAIN tau");
  gtauGAINfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gtauGAINfit->GetXaxis()->SetTitle("run group");
  //gtauGAINfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gtauGAINfit->GetYaxis()->SetTitle("GAIN ampl");
  gtauGAINfit->GetYaxis()->SetTitleOffset(0.5);
  gtauGAINfit->GetYaxis()->SetTitleSize(0.07);
  gtauGAINfit->GetXaxis()->SetTitleOffset(0.5);
  gtauGAINfit->GetXaxis()->SetTitleSize(0.07);
  gtauGAINfit->SetMarkerStyle(8);
  gtauGAINfit->SetMarkerSize(1.0);
  gtauGAINfit->SetMarkerColor(kTeal);
  gtauGAINfit->SetLineWidth(2);
  gtauGAINfit->SetLineColor(kTeal);

  gmuonlossfit = new TGraphErrors(Nfit, GPar_QM, muonlossF_QM, dGPar_QM, dmuonlossF_QM);
  gmuonlossfit->SetName("gmuonlossfit");
  gmuonlossfit->SetTitle("Q-Method muon loss ampl.");
  gmuonlossfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gmuonlossfit->GetXaxis()->SetTitle("run group");
  //gmuonlossfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gmuonlossfit->GetYaxis()->SetTitle("muon loss ampl.");
  gmuonlossfit->GetYaxis()->SetTitleOffset(0.5);
  gmuonlossfit->GetYaxis()->SetTitleSize(0.07);
  gmuonlossfit->GetXaxis()->SetTitleOffset(0.5);
  gmuonlossfit->GetXaxis()->SetTitleSize(0.07);
  gmuonlossfit->SetMarkerStyle(8);
  gmuonlossfit->SetMarkerSize(1.0);
  gmuonlossfit->SetMarkerColor(kBlack);
  gmuonlossfit->SetLineWidth(2);
  gmuonlossfit->SetLineColor(kBlack);

  gpeddriftfit = new TGraphErrors(Nfit, GPar_QM, peddriftF_QM, dGPar_QM, dpeddriftF_QM);
  gpeddriftfit->SetName("gpeddriftfit");
  gpeddriftfit->SetTitle("Q-Method ped drift ampl.");
  gpeddriftfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gpeddriftfit->GetXaxis()->SetTitle("run group");
  //gpeddriftfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gpeddriftfit->GetYaxis()->SetTitle("muon loss ampl.");
  gpeddriftfit->GetYaxis()->SetTitleOffset(0.5);
  gpeddriftfit->GetYaxis()->SetTitleSize(0.07);
  gpeddriftfit->GetXaxis()->SetTitleOffset(0.5);
  gpeddriftfit->GetXaxis()->SetTitleSize(0.07);
  gpeddriftfit->SetMarkerStyle(8);
  gpeddriftfit->SetMarkerSize(1.0);
  gpeddriftfit->SetMarkerColor(kBlack);
  gpeddriftfit->SetLineWidth(2);
  gpeddriftfit->SetLineColor(kBlack);

  gRPILEUPfit = new TGraphErrors(Nfit, GPar_QM, RPILEUPF_QM, dGPar_QM, dRPILEUPF_QM);
  gRPILEUPfit->SetName("gRPILEUPfit");
  gRPILEUPfit->SetTitle("Q-Method PILEUP term");
  gRPILEUPfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gRPILEUPfit->GetXaxis()->SetTitle("run group");
  //gRPILEUPfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gRPILEUPfit->GetYaxis()->SetTitle("PILEUP term");
  gRPILEUPfit->GetYaxis()->SetTitleOffset(0.5);
  gRPILEUPfit->GetYaxis()->SetTitleSize(0.07);
  gRPILEUPfit->GetXaxis()->SetTitleOffset(0.5);
  gRPILEUPfit->GetXaxis()->SetTitleSize(0.07);
  gRPILEUPfit->SetMarkerStyle(8);
  gRPILEUPfit->SetMarkerSize(1.0);
  gRPILEUPfit->SetMarkerColor(kTeal);
  gRPILEUPfit->SetLineWidth(2);
  gRPILEUPfit->SetLineColor(kTeal);

  gtauCBOfit = new TGraphErrors(Nfit, GPar_QM, tauCBOF_QM, dGPar_QM, dtauCBOF_QM);
  gtauCBOfit->SetName("gtauCBOfit");
  gtauCBOfit->SetTitle("Q-Method CBO tau");
  gtauCBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gtauCBOfit->GetXaxis()->SetTitle("run group");
  //gtauCBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gtauCBOfit->GetYaxis()->SetTitle("CBO tau");
  gtauCBOfit->GetYaxis()->SetTitleOffset(0.5);
  gtauCBOfit->GetYaxis()->SetTitleSize(0.07);
  gtauCBOfit->GetXaxis()->SetTitleOffset(0.5);
  gtauCBOfit->GetXaxis()->SetTitleSize(0.07);
  gtauCBOfit->SetMarkerStyle(8);
  gtauCBOfit->SetMarkerSize(1.0);
  gtauCBOfit->SetMarkerColor(38);
  gtauCBOfit->SetLineWidth(2);
  gtauCBOfit->SetLineColor(38);

  gphaseCBOfit = new TGraphErrors(Nfit, GPar_QM, phaseCBOF_QM, dGPar_QM, dphaseCBOF_QM);
  gphaseCBOfit->SetName("gphaseCBOfit");
  gphaseCBOfit->SetTitle("Q-Method CBO phase");
  gphaseCBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gphaseCBOfit->GetXaxis()->SetTitle("run group");
 //gphaseCBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gphaseCBOfit->GetYaxis()->SetTitle("CBO phase");
  gphaseCBOfit->GetYaxis()->SetTitleOffset(0.5);
  gphaseCBOfit->GetYaxis()->SetTitleSize(0.07);
  gphaseCBOfit->GetXaxis()->SetTitleOffset(0.5);
  gphaseCBOfit->GetXaxis()->SetTitleSize(0.07);
  gphaseCBOfit->SetMarkerStyle(8);
  gphaseCBOfit->SetMarkerSize(1.1);
  gphaseCBOfit->SetMarkerColor(kRed);
  gphaseCBOfit->SetLineWidth(2);
  gphaseCBOfit->SetLineColor(kRed);

  gfreqCBOfit = new TGraphErrors(Nfit, GPar_QM, freqCBOF_QM, dGPar_QM, dfreqCBOF_QM);
  gfreqCBOfit->SetName("gfreqCBOfit");
  gfreqCBOfit->SetTitle("Q-Method CBO freq");
  gfreqCBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gfreqCBOfit->GetXaxis()->SetTitle("run group");
  //gfreqCBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gfreqCBOfit->GetYaxis()->SetTitle("CBO freq");
  gfreqCBOfit->GetYaxis()->SetTitleOffset(0.5);
  gfreqCBOfit->GetYaxis()->SetTitleSize(0.07);
  gfreqCBOfit->GetXaxis()->SetTitleOffset(0.5);
  gfreqCBOfit->GetXaxis()->SetTitleSize(0.07);
  gfreqCBOfit->SetMarkerStyle(8);
  gfreqCBOfit->SetMarkerSize(1.0);
  gfreqCBOfit->SetMarkerColor(kGreen);
  gfreqCBOfit->SetLineWidth(2);
  gfreqCBOfit->SetLineColor(kGreen);

  gdfreqCBOfit = new TGraphErrors(Nfit, GPar_QM, deltafreqCBOF_QM, dGPar_QM, ddeltafreqCBOF_QM);
  gdfreqCBOfit->SetName("gdfreqCBOfit");
  gdfreqCBOfit->SetTitle("Q-Method CBO delta freq");
  gdfreqCBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gdfreqCBOfit->GetXaxis()->SetTitle("run group");
  //gdfreqCBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gdfreqCBOfit->GetYaxis()->SetTitle("CBO delta freq");
  gdfreqCBOfit->GetYaxis()->SetTitleOffset(0.5);
  gdfreqCBOfit->GetYaxis()->SetTitleSize(0.07);
  gdfreqCBOfit->GetXaxis()->SetTitleOffset(0.5);
  gdfreqCBOfit->GetXaxis()->SetTitleSize(0.07);
  gdfreqCBOfit->SetMarkerStyle(8);
  gdfreqCBOfit->SetMarkerSize(1.0);
  gdfreqCBOfit->SetMarkerColor(kGreen);
  gdfreqCBOfit->SetLineWidth(2);
  gdfreqCBOfit->SetLineColor(kGreen);

  gACBO2fit = new TGraphErrors(Nfit, GPar_QM, ACBO2F_QM, dGPar_QM, dACBO2F_QM);
  gACBO2fit->SetName("gACBO2fit");
  gACBO2fit->SetTitle("Q-Method CBO2 ampl (for A_a)");
  gACBO2fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gACBO2fit->GetXaxis()->SetTitle("run group");
  //gACBO2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gACBO2fit->GetYaxis()->SetTitle("CBO ampl");
  gACBO2fit->GetYaxis()->SetTitleOffset(0.5);
  gACBO2fit->GetYaxis()->SetTitleSize(0.07);
  gACBO2fit->GetXaxis()->SetTitleOffset(0.5);
  gACBO2fit->GetXaxis()->SetTitleSize(0.07);
  gACBO2fit->SetMarkerStyle(8);
  gACBO2fit->SetMarkerSize(1.0);
  gACBO2fit->SetMarkerColor(28);
  gACBO2fit->SetLineWidth(2);
  gACBO2fit->SetLineColor(28);

  gACBO3fit = new TGraphErrors(Nfit, GPar_QM, ACBO3F_QM, dGPar_QM, dACBO3F_QM);
  gACBO3fit->SetName("gACBO3fit");
  gACBO3fit->SetTitle("Q-Method CBO3 ampl (for phi_a)");
  gACBO3fit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gACBO3fit->GetXaxis()->SetTitle("run group");
  //gACBO3fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gACBO3fit->GetYaxis()->SetTitle("CBO ampl");
  gACBO3fit->GetYaxis()->SetTitleOffset(0.5);
  gACBO3fit->GetYaxis()->SetTitleSize(0.07);
  gACBO3fit->GetXaxis()->SetTitleOffset(0.5);
  gACBO3fit->GetXaxis()->SetTitleSize(0.07);
  gACBO3fit->SetMarkerStyle(8);
  gACBO3fit->SetMarkerSize(1.0);
  gACBO3fit->SetMarkerColor(28);
  gACBO3fit->SetLineWidth(2);
  gACBO3fit->SetLineColor(28);

  g2ACBOfit = new TGraphErrors(Nfit, GPar_QM, A2CBOF_QM, dGPar_QM, dA2CBOF_QM);
  g2ACBOfit->SetName("g2ACBOfit");
  g2ACBOfit->SetTitle("Q-Method 2xCBO ampl");
  gACBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  g2ACBOfit->GetXaxis()->SetTitle("run group");
  //g2ACBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  g2ACBOfit->GetYaxis()->SetTitle("2xCBO ampl");
  g2ACBOfit->GetYaxis()->SetTitleOffset(0.5);
  g2ACBOfit->GetYaxis()->SetTitleSize(0.07);
  g2ACBOfit->GetXaxis()->SetTitleOffset(0.5);
  g2ACBOfit->GetXaxis()->SetTitleSize(0.07);
  g2ACBOfit->SetMarkerStyle(8);
  g2ACBOfit->SetMarkerSize(1.0);
  g2ACBOfit->SetMarkerColor(28);
  g2ACBOfit->SetLineWidth(2);
  g2ACBOfit->SetLineColor(28);

  gphase2CBOfit = new TGraphErrors(Nfit, GPar_QM, phaseCBO2F_QM, dGPar_QM, dphaseCBO2F_QM);
  gphase2CBOfit->SetName("gphase2CBOfit");
  gphase2CBOfit->SetTitle("Q-Method CBO2 phase (for A_a)");
  gphase2CBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gphase2CBOfit->GetXaxis()->SetTitle("run group");
 //gphase2CBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gphase2CBOfit->GetYaxis()->SetTitle("CBO phase");
  gphase2CBOfit->GetYaxis()->SetTitleOffset(0.5);
  gphase2CBOfit->GetYaxis()->SetTitleSize(0.07);
  gphase2CBOfit->GetXaxis()->SetTitleOffset(0.5);
  gphase2CBOfit->GetXaxis()->SetTitleSize(0.07);
  gphase2CBOfit->SetMarkerStyle(8);
  gphase2CBOfit->SetMarkerSize(1.0);
  gphase2CBOfit->SetMarkerColor(7);
  gphase2CBOfit->SetLineWidth(2);
  gphase2CBOfit->SetLineColor(7);

  gphase3CBOfit = new TGraphErrors(Nfit, GPar_QM, phaseCBO3F_QM, dGPar_QM, dphaseCBO3F_QM);
  gphase3CBOfit->SetName("gphase3CBOfit");
  gphase3CBOfit->SetTitle("Q-Method CBO3 phase (for phi_a)");
  gphase3CBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gphase3CBOfit->GetXaxis()->SetTitle("run group");
 //gphase3CBOfit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gphase3CBOfit->GetYaxis()->SetTitle("CBO phase");
  gphase3CBOfit->GetYaxis()->SetTitleOffset(0.5);
  gphase3CBOfit->GetYaxis()->SetTitleSize(0.07);
  gphase3CBOfit->GetXaxis()->SetTitleOffset(0.5);
  gphase3CBOfit->GetXaxis()->SetTitleSize(0.07);
  gphase3CBOfit->SetMarkerStyle(8);
  gphase3CBOfit->SetMarkerSize(1.0);
  gphase3CBOfit->SetMarkerColor(7);
  gphase3CBOfit->SetLineWidth(2);
  gphase3CBOfit->SetLineColor(7);

  g2phaseCBOfit = new TGraphErrors(Nfit, GPar_QM, phase2CBOF_QM, dGPar_QM, dphase2CBOF_QM);
  g2phaseCBOfit->SetName("g2phaseCBOfit");
  g2phaseCBOfit->SetTitle("Q-Method 2xCBO phase");
  g2phaseCBOfit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  g2phaseCBOfit->GetXaxis()->SetTitle("run group");
 //g2phaseCBO2fit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  g2phaseCBOfit->GetYaxis()->SetTitle("2xCBO phase");
  g2phaseCBOfit->GetYaxis()->SetTitleOffset(0.5);
  g2phaseCBOfit->GetYaxis()->SetTitleSize(0.07);
  g2phaseCBOfit->GetXaxis()->SetTitleOffset(0.5);
  g2phaseCBOfit->GetXaxis()->SetTitleSize(0.07);
  g2phaseCBOfit->SetMarkerStyle(8);
  g2phaseCBOfit->SetMarkerSize(1.0);
  g2phaseCBOfit->SetMarkerColor(7);
  g2phaseCBOfit->SetLineWidth(2);
  g2phaseCBOfit->SetLineColor(7);

  /*

  //GCBOmeanX[ICalo] = gACBOfit->GetRMS(1);
  double xtmp, ytmp, sxy = 0.0, sy = 0.0, meanx;
  for (int imean = 0; imean < Nfit; imean++){
    gACBOfit->GetPoint( imean, xtmp, ytmp);
    printf("xtmp, ytmp %f, %f\n", xtmp, ytmp);
    sxy += xtmp*ytmp;
    sy += ytmp;
  }
  meanx = sxy/sy;
  printf("sy, sxy, meanx %f, %f, %f\n", sy, sxy, meanx);
  
  GCBOmeanX[ICalo-1] = meanx;
  dGCBOmeanX[ICalo-1] = 0;
  GCalonum[ICalo-1] = ICalo;
  dGCalonum[ICalo-1] = 0;
  printf("calo number %f, meanX %f\n", GCalonum[ICalo], GCBOmeanX[ICalo]);

  */

  gChifit = new TGraphErrors(Nfit, GPar_QM, chisq_QM, dGPar_QM, dchisq_QM);
  gChifit->SetName("gChifit");
  gChifit->SetTitle("Q-Method chi-squared/pdf");
  gChifit->GetXaxis()->SetTitle("start time (800 MSPS units)");
  gChifit->GetXaxis()->SetTitle("run group");
  //gChifit->GetXaxis()->SetTitle("end time (800 MSPS units)");
  gChifit->GetYaxis()->SetTitle("Chi_squared/pdf");
  gChifit->GetYaxis()->SetTitleOffset(0.5);
  gChifit->GetYaxis()->SetTitleSize(0.07);
  gChifit->GetXaxis()->SetTitleOffset(0.5);
  gChifit->GetXaxis()->SetTitleSize(0.07);
  gChifit->SetMarkerStyle(8);
  gChifit->SetMarkerSize(1.0);
  gChifit->SetMarkerColor(9);
  gChifit->SetLineWidth(2);
  gChifit->SetLineColor(9);

  /*
  c1->Clear();
  c1->Divide(1,1);
  gNfit->Draw("AP");
  c1->SaveAs("param1scan.png");
  */

  c1->Clear();
  c1->Divide(1,1);
  gTaufit->SetMaximum(66000);
  gTaufit->SetMinimum(63000);
  gTaufit->Draw("AP");
  c1->SaveAs("param2scan.png");

  /*
  c1->Clear();
  c1->Divide(1,1);
  gAfit->Draw("AP");
  c1->SaveAs("param3scan.png");
  */

  c1->Clear();
  c1->Divide(1,1);
  gOmega_afit->Draw("AP");
  c1->SaveAs("param4scan.png");

  c1->Clear();
  c1->Divide(1,1);
  gPhifit->Draw("AP");
  c1->SaveAs("param5scan.png");

  /*
  c1->Clear();
  c1->Divide(1,1);
  //gBfit->Draw("AP");
  gAGAINfit->Draw("AP");
  c1->SaveAs("param6scan.png");

  c1->Clear();
  c1->Divide(1,1);
  //gBfit->Draw("AP");
  gRPILEUPfit->Draw("AP");
  c1->SaveAs("RPILUEPscan.png");

  c1->Clear();
  c1->Divide(1,1);
  gACBOfit->Draw("AP");
  c1->SaveAs("ACBOscan.png");

  c1->Clear();
  c1->Divide(1,1);
  gACBO2fit->Draw("AP");
  c1->SaveAs("ACBO2scan.png");

  c1->Clear();
  c1->Divide(1,1);
  gACBO3fit->Draw("AP");
  c1->SaveAs("ACBO3scan.png");


  c1->Clear();
  c1->Divide(1,1);
  gtauCBOfit->Draw("AP");
  c1->SaveAs("tauCBOscan.png");

  c1->Clear();
  c1->Divide(1,1);
  gphaseCBOfit->Draw("AP");
  c1->SaveAs("phaseCBOscan.png");

  c1->Clear();
  c1->Divide(1,1);
  gphase2CBOfit->Draw("AP");
  c1->SaveAs("phaseCBO2scan.png");

  c1->Clear();
  c1->Divide(1,1);
  gphase3CBOfit->Draw("AP");
  c1->SaveAs("phaseCBO3scan.png");

  c1->Clear();
  c1->Divide(1,1);
  gfreqCBOfit->Draw("AP");
  c1->SaveAs("freqCBOscan.png");

  c1->Clear();
  c1->Divide(1,1);
  //gBfit->Draw("AP");
  gmuonlossfit->Draw("AP");
  c1->SaveAs("muonlossscan.png");

  c1->Clear();
  c1->Divide(1,1);
  gChifit->Draw("AP");
  c1->SaveAs("Chisqscan.png");
  */

  gChifit->GetYaxis()->SetLabelSize(0.08);
  gNfit->GetYaxis()->SetLabelSize(0.08);
  gTaufit->GetYaxis()->SetLabelSize(0.08);
  gAfit->GetYaxis()->SetLabelSize(0.08);
  gArelaxfit->GetYaxis()->SetLabelSize(0.08);
  gPhifit->GetYaxis()->SetLabelSize(0.08);
  gOmega_afit->GetYaxis()->SetLabelSize(0.08);
  gAGAINfit->GetYaxis()->SetLabelSize(0.08);
  gtauGAINfit->GetYaxis()->SetLabelSize(0.08);
  gmuonlossfit->GetYaxis()->SetLabelSize(0.08);
  gpeddriftfit->GetYaxis()->SetLabelSize(0.08);
  gRPILEUPfit->GetYaxis()->SetLabelSize(0.08);
  gACBOfit->GetYaxis()->SetLabelSize(0.08);
  gACBO2fit->GetYaxis()->SetLabelSize(0.08);
  gACBO3fit->GetYaxis()->SetLabelSize(0.08);
  gtauCBOfit->GetYaxis()->SetLabelSize(0.08);
  gfreqCBOfit->GetYaxis()->SetLabelSize(0.08);
  gdfreqCBOfit->GetYaxis()->SetLabelSize(0.08);
  gphaseCBOfit->GetYaxis()->SetLabelSize(0.08);
  gphase2CBOfit->GetYaxis()->SetLabelSize(0.08);
  gphase3CBOfit->GetYaxis()->SetLabelSize(0.08);
  gAVWfit->GetYaxis()->SetLabelSize(0.08);
  gomegaVWfit->GetYaxis()->SetLabelSize(0.08);
  gtauVWfit->GetYaxis()->SetLabelSize(0.08);
  gphaseVWfit->GetYaxis()->SetLabelSize(0.08);
  gAVW2fit->GetYaxis()->SetLabelSize(0.08);
  gomegaVW2fit->GetYaxis()->SetLabelSize(0.08);
  gtauVW2fit->GetYaxis()->SetLabelSize(0.08);
  gphaseVW2fit->GetYaxis()->SetLabelSize(0.08);
  gChifit->GetYaxis()->SetTitleSize(0.12);
  gNfit->GetYaxis()->SetTitleSize(0.12);
  gTaufit->GetYaxis()->SetTitleSize(0.12);
  gAfit->GetYaxis()->SetTitleSize(0.12);
  gArelaxfit->GetYaxis()->SetTitleSize(0.12);
  gPhifit->GetYaxis()->SetTitleSize(0.12);
  gOmega_afit->GetYaxis()->SetTitleSize(0.12);
  gAGAINfit->GetYaxis()->SetTitleSize(0.12);
  gtauGAINfit->GetYaxis()->SetTitleSize(0.12);
  gmuonlossfit->GetYaxis()->SetTitleSize(0.12);
  gpeddriftfit->GetYaxis()->SetTitleSize(0.12);
  gRPILEUPfit->GetYaxis()->SetTitleSize(0.12);
  gACBOfit->GetYaxis()->SetTitleSize(0.12);
  gACBO2fit->GetYaxis()->SetTitleSize(0.12);
  gACBO3fit->GetYaxis()->SetTitleSize(0.12);
  gtauCBOfit->GetYaxis()->SetTitleSize(0.12);
  gfreqCBOfit->GetYaxis()->SetTitleSize(0.12);
  gdfreqCBOfit->GetYaxis()->SetTitleSize(0.12);
  gphaseCBOfit->GetYaxis()->SetTitleSize(0.12);
  gphase2CBOfit->GetYaxis()->SetTitleSize(0.12);
  gphase3CBOfit->GetYaxis()->SetTitleSize(0.12);
  gChifit->GetYaxis()->SetTitleSize(0.12);
  gNfit->GetYaxis()->SetTitleSize(0.12);
  gTaufit->GetYaxis()->SetTitleSize(0.12);
  gAfit->GetYaxis()->SetTitleSize(0.12);
  gArelaxfit->GetYaxis()->SetTitleSize(0.12);
  gPhifit->GetYaxis()->SetTitleSize(0.12);
  gOmega_afit->GetYaxis()->SetTitleSize(0.12);
  gAGAINfit->GetYaxis()->SetTitleSize(0.12);
  gtauGAINfit->GetYaxis()->SetTitleSize(0.12);
  gmuonlossfit->GetYaxis()->SetTitleSize(0.12);
  gpeddriftfit->GetYaxis()->SetTitleSize(0.12);
  gRPILEUPfit->GetYaxis()->SetTitleSize(0.12);
  gACBOfit->GetYaxis()->SetTitleSize(0.12);
  gACBO2fit->GetYaxis()->SetTitleSize(0.12);
  gACBO3fit->GetYaxis()->SetTitleSize(0.12);
  gtauCBOfit->GetYaxis()->SetTitleSize(0.12);
  gfreqCBOfit->GetYaxis()->SetTitleSize(0.12);
  gdfreqCBOfit->GetYaxis()->SetTitleSize(0.12);
  gphaseCBOfit->GetYaxis()->SetTitleSize(0.12);
  gphase2CBOfit->GetYaxis()->SetTitleSize(0.12);
  gphase3CBOfit->GetYaxis()->SetTitleSize(0.12);
  gAVWfit->GetYaxis()->SetTitleSize(0.12);
  gomegaVWfit->GetYaxis()->SetTitleSize(0.12);
  gtauVWfit->GetYaxis()->SetTitleSize(0.12);
  gphaseVWfit->GetYaxis()->SetTitleSize(0.12);
  gAVW2fit->GetYaxis()->SetTitleSize(0.12);
  gomegaVW2fit->GetYaxis()->SetTitleSize(0.12);
  gtauVW2fit->GetYaxis()->SetTitleSize(0.12);
  gphaseVW2fit->GetYaxis()->SetTitleSize(0.12);
  gChifit->GetYaxis()->SetTitleOffset(-0.22);
  gNfit->GetYaxis()->SetTitleOffset(-0.22);
  gTaufit->GetYaxis()->SetTitleOffset(-0.22);
  gAfit->GetYaxis()->SetTitleOffset(-0.22);
  gArelaxfit->GetYaxis()->SetTitleOffset(-0.22);
  gPhifit->GetYaxis()->SetTitleOffset(-0.22);
  gOmega_afit->GetYaxis()->SetTitleOffset(-0.22);
  gAGAINfit->GetYaxis()->SetTitleOffset(-0.22);
  gtauGAINfit->GetYaxis()->SetTitleOffset(-0.22);
  gmuonlossfit->GetYaxis()->SetTitleOffset(-0.22);
  gpeddriftfit->GetYaxis()->SetTitleOffset(-0.22);
  gRPILEUPfit->GetYaxis()->SetTitleOffset(-0.22);
  gACBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gACBO2fit->GetYaxis()->SetTitleOffset(-0.22);
  gACBO3fit->GetYaxis()->SetTitleOffset(-0.22);
  gtauCBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gfreqCBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gphaseCBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gphase2CBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gphase3CBOfit->GetYaxis()->SetTitleOffset(-0.22);
  gAVWfit->GetYaxis()->SetTitleOffset(-0.22);
  gomegaVWfit->GetYaxis()->SetTitleOffset(-0.22);
  gtauVWfit->GetYaxis()->SetTitleOffset(-0.22);
  gphaseVWfit->GetYaxis()->SetTitleOffset(-0.22);
  gAVW2fit->GetYaxis()->SetTitleOffset(-0.22);
  gomegaVW2fit->GetYaxis()->SetTitleOffset(-0.22);
  gtauVW2fit->GetYaxis()->SetTitleOffset(-0.22);
  gphaseVW2fit->GetYaxis()->SetTitleOffset(-0.22);
  

  c1->Clear();
  c1->Divide(2,8);
  c1->cd(1);
  gNfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gNfit->Draw("AP");
  c1->cd(2);
  gTaufit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gTaufit->Draw("AP");
  c1->cd(3);
  gAfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gAfit->Draw("AP");
  c1->cd(4);
  gOmega_afit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gOmega_afit->Draw("AP");
  c1->cd(5);
  gPhifit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gPhifit->Draw("AP");
  c1->cd(6);
  gACBOfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gACBOfit->Draw("AP");
  c1->cd(7);  
  gfreqCBOfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gfreqCBOfit->Draw("AP");
  //gdfreqCBOfit->Draw("AP");
  c1->cd(8);  
  gphaseCBOfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gphaseCBOfit->Draw("AP");
  c1->cd(9);  
  //gAGAINfit->Draw("AP");
  //gpeddriftfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  //gpeddriftfit->Draw("AP");
  gmuonlossfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gmuonlossfit->Draw("AP");
  c1->cd(10); 
  gAVWfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gAVWfit->Draw("AP");
  c1->cd(11);
  gphaseVWfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gphaseVWfit->Draw("AP");
  c1->cd(12);
  gAVW2fit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gAVW2fit->Draw("AP");
  c1->cd(13);
  gAGAINfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gAGAINfit->Draw("AP");
  c1->cd(14);
  gtauGAINfit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gtauGAINfit->Draw("AP");
  c1->cd(15);
  gphaseVW2fit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gphaseVW2fit->Draw("AP");
  c1->cd(16);
  gChifit->GetXaxis()->SetRangeUser(0.5,Ifit+0.5);
  gChifit->Draw("AP");

  /*
  c1->Clear();
  c1->Divide(2,3);
  c1->cd(1);
  gACBOfit->Draw("AP");
  c1->cd(2);
  gtauCBOfit->Draw("AP");
  c1->cd(3);
  gfreqCBOfit->Draw("AP");
  c1->cd(4);
  gphaseCBOfit->Draw("AP");
  c1->cd(5);
  gChifit->Draw("AP");
  */

  sprintf( hname, "cboVersusCalo%02i.png", ICalo);
  c1->SaveAs(hname);


  // save a graph
  sprintf( foutname, "graphs-w%i-t%i-g%i-swrb%i-prd%s-%s-%s-%s-%s.root", iwind, ithreshold, igap, iswrb, aname, cFit, fitname, dataname, datasubset);
  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("created output root file %s\n", foutname);

  gNfit->Write();
  gBfit->Write();
  gTaufit->Write();
  gAfit->Write();
  gArelaxfit->Write();
  gOmega_afit->Write();
  gPhifit->Write();
  gACBOfit->Write();
  gfreqCBOfit->Write();
  gdfreqCBOfit->Write();
  gphaseCBOfit->Write();
  gtauCBOfit->Write();
  gACBO2fit->Write();
  gphase2CBOfit->Write();
  g2ACBOfit->Write();
  g2phaseCBOfit->Write();
  gAGAINfit->Write();
  gtauGAINfit->Write();
  gmuonlossfit->Write();
  gpeddriftfit->Write();
  gChifit->Write();
  gAVWfit->Write();
  gomegaVWfit->Write();
  gtauVWfit->Write();
  gphaseVWfit->Write();
  gAVW2fit->Write();
  gomegaVW2fit->Write();
  gtauVW2fit->Write();
  gphaseVW2fit->Write();
  gRPILEUPfit->Write();

  file1->Print();
  printf("wrote graph into file %s\n", foutname);  
  file1->Close();

  //return 0;

  // fitting a sine curve to calo dependence

  p0f = new TF1("p0f","[0] +[1]*sin(2*3.14*x/24+[2])",0,24);
  p0f->SetParameters( 0, -35.0);
  p0f->SetParameter( 1, 0.0);
  p0f->SetParameter( 2, 3.14);
  p0f->SetParNames( "R", "dR", "phase");
  p0f->FixParameter(1, 0.0);
  p0f->FixParameter(2, 3.14);
  gOmega_afit->Fit( p0f, "", "VR", 0, 24);
  printf("p0f %f +/- %f, chsq %f, ndf %i, cdsq/ndf %f\n", p0f->GetParameter(0),p0f->GetParError(0), p0f->GetChisquare(), p0f->GetNDF(), p0f->GetChisquare()/p0f->GetNDF());
  printf("probability, %f", p0f->GetProb());

  //gOmega_afit->Fit( p0f, "", "VR", 0, 8);
  //printf("hot p0f %f +/- %f, chsq %f, ndf %i, cdsq/ndf %f\n", p0f->GetParameter(0),p0f->GetParError(0), p0f->GetChisquare(), p0f->GetNDF(), p0f->GetChisquare()/p0f->GetNDF());
  //gOmega_afit->Fit( p0f, "", "VR", 14, 22);
  //printf("cold p0f %f +/- %f, chsq %f, ndf %i, cdsq/ndf %f\n", p0f->GetParameter(0),p0f->GetParError(0), p0f->GetChisquare(), p0f->GetNDF(), p0f->GetChisquare()/p0f->GetNDF());
 
  //

  c1->Clear();
  gStyle->SetOptFit(1111);
  gOmega_afit->Draw("AP");

  return 1;
}


double A_vbo = 0.0, omega_vbo, phi_vbo, tau_vbo; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
double A_vbo2 = 0.0, omega_vbo2, phi_vbo2, tau_vbo2; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
double A_vbo2asym, phi_vbo2asym; // asymmetry term in vbo as seen in FFT

int viewCalo( char* cFit, char* cRun, char* cHis, int swrb, double startTime, double endTime){  // in 800 MSPS units
  
  // open root file
  sprintf( foutname, "%s", cRun);
  file1 = TFile::Open( foutname, "READ");
  printf("got input root file %s\n", foutname);

  //find energy fraction on individual calos in fit region (used in old fast rotation correction)
  TH1D *hix; 
  char hn[60];
  TDirectoryFile *dr;
  double calosm = 0.0, calofrac[24]={0.0};
  file1->GetObject(dname,dr);
  for (int ix = 1; ix <= 24; ix++){
    sprintf( hn, "qHist1D_sig_%i_0", ix);
    dr->GetObject( hn, hix);
    calofrac[ix-1] = hix->Integral(hix->FindBin(startTime*rawBinToNs), hix->FindBin(endTime*rawBinToNs));
    calosm += calofrac[ix-1];
  }
  for (int ix = 1; ix <= 24; ix++) {
    calofrac[ix-1] /= calosm; 
    printf("ix %i, calofrac %f\n", ix, calofrac[ix-1]);
  }

  // ***** histogram isn't in folder CQBankAnalyzerTim
  // for 1D histograms
  if ( ( strcmp( "full", cFit) == 0 ) || ( strcmp( "near", cFit) == 0 )  || ( strcmp( "far", cFit) == 0 ) || ( strcmp( "hot", cFit) == 0 )  || ( strcmp( "cold", cFit) == 0 )  || ( strcmp( "paircalo", cFit) == 0 ) || ( strcmp( "pairscan", cFit) == 0 ) || ( strcmp( "runscan", cFit) == 0 ) || ( strcmp( "noisescan", cFit) == 0 ) || ( strcmp( "puscan", cFit) == 0 ) || ( strcmp( "starttimescan", cFit) == 0 ) || ( strcmp( "tauscan", cFit) == 0 ) || ( strcmp( "windowscan", cFit) == 0 ) || ( strcmp( "thresholdscan", cFit) == 0 ) || ( strcmp( "coldthresholdscan", cFit) == 0 ) || ( strcmp( "hotthresholdscan", cFit) == 0 ) ) {
    sprintf( hname, "%s", cHis);
    printf( "get 1D, root folder %s\n", hname);
    file1->GetObject( hname, htmp); 
    printf( "got 1D, root folder %s\n", hname);
  }
  // for 2D histograms
  if ( ( strcmp( "top", cFit) == 0 ) || ( strcmp( "bot", cFit) == 0 ) ) { 
      file1->GetObject( hname, h2tmp);
      htmp = h2tmp->ProjectionX("htmp",firstcaloinslice,lastcaloinslice,"");
      printf("got 2D histogram\n");
    }

  // ***** histogram is in folder CQBankAnalyzerTim
  // for 1D histograms
  if ( ( strcmp( "singlecalo", cFit) == 0 ) || ( strcmp( "sequencescan", cFit) == 0 ) || ( strcmp( "caloscan", cFit) == 0 ) ) {
    TDirectoryFile *dir;
    file1->GetObject(dname,dir);
    printf("got directory\n");
    // if 1D
    dir->GetObject( hname, htmp);
    printf("got 1D histogram\n");
  }
  // ***** histogram is in sub-folder
  // for 1D histograms
  if ( ( strcmp( "verticalscan", cFit) == 0 ) || ( strcmp( "horizontalscan", cFit) == 0 ) || ( strcmp( "verticalslice", cFit) == 0 ) || ( strcmp( "horizontalslice", cFit) == 0 ) ) {
    TDirectoryFile *dir;
    file1->GetObject( dname, dir);
    printf("got directory %s for histogram %s\n", dname, hname);
    // if 2D
    dir->GetObject( hname, h2tmp);
    htmp = h2tmp->ProjectionX("htmp",firstcaloinslice,lastcaloinslice,"");
    printf("got 2D histogram\n");
  }
  
  // changes for run 2
  // rebin to 60 ct bins 
  if (run_2) {    
    printf("do changes for run 2\n");
    htmp->Rebin(4);
    htmp->Draw("hist");
  }

  printf( "do axes, etc %s\n", hname);
  htmp->GetXaxis()->SetTitle("time (ns)");
  htmp->GetYaxis()->SetTitle("ADC counts");
  htmp->GetYaxis()->SetTitleOffset(0.5);
  htmp->GetYaxis()->SetTitleSize(0.07);
  htmp->GetXaxis()->SetTitleOffset(0.5);
  htmp->GetXaxis()->SetTitleSize(0.07);
  htmp->SetMarkerStyle(8);
  htmp->SetMarkerSize(0.1);
  htmp->SetMarkerColor(kGray+2);
  htmp->SetLineColor(kGray+2);
  printf( "done axes, etc %s\n", hname);

  if (fixbug){
    // FIX for bug in raw bin - decimated bin correspondance, Oct21, 2018
    //printf("htmp before adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp->GetNbinsX(), htmp->GetBinWidth(1), htmp->GetBinLowEdge(1), htmp->GetBinLowEdge( htmp->GetNbinsX() + htmp->GetBinWidth(htmp->GetNbinsX()) ) );
    htmp->SetBins(3168/swrb,20001,210021+60);
      //printf("htmp after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp->GetNbinsX(), htmp->GetBinWidth(1), htmp->GetBinLowEdge(1), htmp->GetBinLowEdge( htmp->GetNbinsX() + htmp->GetBinWidth(htmp->GetNbinsX()) ) );
  }
    
  // run group 786x, ... have time decimation 60, flush factor 1
  printf( "rebin %s\n", hname);
  int extrarb;

  if (swrb == 1 && fastRotationAnalysis ) extrarb = 1; // fast rotation analysis 
  if (swrb == 1 && !fastRotationAnalysis ) extrarb = 2; // regular wiggle analysis
  if (swrb == 2) extrarb = 1; // for swrb 4 cant make 150ns bins
  if (swrb == 4) extrarb = 1; // for swrb 4 cant make 150ns bins

// do a shift of -1 bins to test ideas on fr correction
  htmpshift = (TH1D*)htmp->Clone();
  printf("unshifted example bin N(1000) = %f\n", htmp->GetBinContent( 1000) );
  for (int ib = 2; ib < htmp->GetNbinsX()-2; ib++ ) {

    // 1+2+2+3 for blending and shifting // 
    htmpshift->SetBinContent( ib, htmp->GetBinContent( ib+1 ) ); 
    htmpshift->SetBinError( ib, htmp->GetBinError( ib+1 ) );

    // 2+3+3+4 for blending not shifting //
    //htmpshift->SetBinContent( ib, htmp->GetBinContent( ib+2 ) ); 
    //htmpshift->SetBinError( ib, htmp->GetBinError( ib+2 ) );
    //htmp->SetBinContent( ib, htmp->GetBinContent( ib+1 ) ); 
    //htmp->SetBinError( ib, htmp->GetBinError( ib+1 ) );
  }

  if (doShiftBins) {
    printf("shifted example bins, example bin N(1000) = %f\n", htmpshift->GetBinContent( 1000) );
    for (int ib = 1+1; ib < htmp->GetNbinsX()-1; ib++ ) { // approx error bar fix-up
      htmp->SetBinContent( ib, htmpshift->GetBinContent( ib ) );
      htmp->SetBinError( ib, htmpshift->GetBinError( ib ) );
    }
  }

  htmpshift->Rebin(extrarb);    
  if (doBlendBins) {
    for (int ib = 1; ib < htmp->GetNbinsX(); ib ++ ) { // approx error bar fix-up                                                                                              
       htmpshift->SetBinContent( ib, htmp->GetBinContent( 2*ib) + 0.5*( htmp->GetBinContent( 2*ib-1) + htmp->GetBinContent( 2*ib+1) ) );
       htmpshift->SetBinError( ib, sqrt( htmp->GetBinError( 2*ib)*htmp->GetBinError( 2*ib) + 0.25*( htmp->GetBinError( 2*ib-1)*htmp->GetBinError( 2*ib-1) + htmp->GetBinError( 2*ib+1)*htmp->GetBinError( 2*ib+1) ) ) );
    }
    htmp = (TH1D*)htmpshift->Clone();
      //htmp->Scale(0.5);
      //htmp->Add(htmpshift,0.5);
  } else {
    htmp->Rebin(extrarb);
  }

  printf("do analysis rebinning\n");
  int rb=extrarb*60; // gpu rebin * analysis rebin
  rb = swrb*rb; // multiply by software rebin in gm2 analysis

  // "after the fact" infill gain correction (so doesnt handle the threshold effect)
#if 0
  // if serious about this need to fix for swrb = 2, 4
  if (icalo == 0) {
  printf("do after-the-fact infill gain correction\n");
  sprintf(hname,"Hinfillgain");
  filefr->GetObject( hname, htmp4);

  if (fixbug){
    // FIX for bug in raw bin - decimated bin correspondance, Oct21, 2018
    //printf("htmp4 before adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp4->GetNbinsX(), htmp4->GetBinWidth(1), htmp4->GetBinLowEdge(1), htmp4->GetBinLowEdge( htmp4->GetNbinsX() + htmp4->GetBinWidth(htmp->GetNbinsX()) ) );
    htmp4->SetBins(3168,20001,210021+60);
    //printf("htmp4 after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp4->GetNbinsX(), htmp4->GetBinWidth(1), htmp4->GetBinLowEdge(1), htmp4->GetBinLowEdge( htmp4->GetNbinsX() + htmp4->GetBinWidth(htmp->GetNbinsX()) ) );
  }

  htmp4->Rebin(2);
  htmp4->Scale(0.5);

  //adjust magnitude of gain correction
  for (int ib = 1; ib <= htmp4->GetNbinsX(); ib++ ) htmp4->SetBinContent( ib, htmp4->GetBinContent(ib)-1); // subtract 1
  htmp4->Scale(infillgainscaling); // scaling factor for magnitude of gain correction
  for (int ib = 1; ib <= htmp4->GetNbinsX(); ib++ ) htmp4->SetBinContent( ib, htmp4->GetBinContent(ib)+1); // add 1

  htmp->Divide( htmp4); // do gain correction.
}
#endif

  if (doPUErrorBarCorr) {
  // test error bar correction, divide error bars by overestimate from simulation
// do gain correction for pileup effects on error bar, see docdb
 for (int ib = 1; ib <= htmp->GetNbinsX(); ib++ ) 
   // see page 5, https://gm2-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=14983
   // threshold 12
   htmp->SetBinError( ib, ( 0.9888 + 0.0283*exp( -(htmp->GetBinCenter( ib))/(6.45e4) ) )*htmp->GetBinError(ib) ); 
   // threshold 12
   //htmp->SetBinError( ib, ( 0.995 + 0.0274*exp( -(htmp->GetBinCenter( ib))/(6.29e4) ) )*htmp->GetBinError(ib) ); 
  }
 
 //  fit start time, end time in time decimated, in 800MSPS units
 double minT = startTime, maxT = endTime;
 double minNs = minT, maxNs = maxT;

  // fit start time, end time in time decimated, rebinned units
  int minB = (startTime - htmp->GetBinLowEdge(1))/rb, maxB = (endTime - htmp->GetBinLowEdge(1))/rb;

  for (int ib =  minB; ib <= maxB; ib++){
    // define bin error as sqrt of bin content
    //htmp->SetBinError(  ib, sqrt( htmp->GetBinContent( ib)));
  }
                                                                              
  printf(" minB, maxB %i, %i, minT, maxT %f, %f, minNs, maxNs %f, %f\n", minB, maxB, minT, maxT, minNs, maxNs);
  
  // time dilated muon lifetime tau_mu*gamma
  //2.1969811*sqrt(3.094*3.094+0.10566*0.10566)/0.10566 =  64.3708 us

  // wiggle plot fit function
  //double N = swrb*3.8e8, tau = 6.4371e4/1.25, ampl = 0.21, omega_a = 1.8e-3, phi = 3.02, B = 0.; // all calos
  //double A_cbo = 8.e5, tau_cbo = 1.261e5, omega_cbo = 3.279e-3, phi_cbo = 4.80; // all calos

  // low rate, individual calos  
  //double N = swrb*5.e6, tau = 6.4371e4/1.25, ampl = 0.22, omega_a = 1.79987e-3, phi = 3.02, B = 0.; // calo 24
  //double A_cbo = 10.0e4, tau_cbo = 1.261e5, omega_cbo = 3.279e-3, phi_cbo = 4.80;

  // high rate (60hr), sum calos  
  double N = 24*1.985e7, tau = set_gammatau, ampl = 0.217, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; // calo 24
  double A_cbo = 24*6.15e5, tau_cbo = 1.261e5*rawBinToNs, omega_cbo = 2.90171e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.004099/rawBinToNs;

  double A2_cbo = set_anomasym_amp, phi2_cbo = 0.0, A3_cbo = set_anomphi_amp, phi3_cbo = set_anomphi_phi, A_2cbo = set_2cbo_amp/N, phi_2cbo = set_2cbo_phi; // calo 24
  double A_gain = set_gain_amp, tau_gain = set_gain_tau, t0_gain = 2.8e4; // empirical gain correction parameters
  omega_vbo = 1.7707e-2/rawBinToNs, phi_vbo = 1.5, tau_vbo = 1.8363e4*rawBinToNs; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)

  // this improved chi-squared at early times in 9day, this is VW term

  if ( ( strcmp( "full", cFit) == 0 ) || ( strcmp( "top", cFit) == 0 ) || ( strcmp( "bot", cFit) == 0 ) || ( strcmp( "near", cFit) == 0 )  || ( strcmp( "far", cFit) == 0 )   || ( strcmp( "hot", cFit) == 0 )  || ( strcmp( "cold", cFit) == 0 ) || ( strcmp( "tauscan", cFit) == 0 )  || ( strcmp( "starttimescan", cFit) == 0 ) || ( strcmp( "noisescan", cFit) == 0 ) || ( strcmp( "puscan", cFit) == 0 ) || ( strcmp( "windowscan", cFit) == 0 ) || ( strcmp( "thresholdscan", cFit) == 0 ) || ( strcmp( "hotthresholdscan", cFit) == 0 ) || ( strcmp( "coldthresholdscan", cFit) == 0 ) ) { 

    // best guesses for calorimeter sum, 60hr dataset
    if ( strcmp( "60hr", dataname) == 0 ) {
      printf("60hr dataset initial guesses\n");
      N = 5.8e8, tau = set_gammatau, ampl = 0.226, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; 
      A_cbo = 2.52e6/N, tau_cbo = 1.371e5*rawBinToNs, omega_cbo = 2.90131e-03/rawBinToNs, phi_cbo = 1.7369, delta_omega_cbo = -0.004099/rawBinToNs;
    }

    // high rate sum calos  (note important difference in time and cbo_freq when fitting 9 day dataset)
    if ( strcmp( "9day", dataname) == 0 ) {
      printf("9day dataset initial guesses\n");
      N = 50.e9/24., tau = set_gammatau, ampl = 0.229, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.;
      A_cbo = 2.062e6/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.25088e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
    }    

    // high rate sum calos  (note important difference in time and cbo_freq when fitting hk day dataset)
    if ( strcmp( "hk", dataname) == 0 ) {
      printf("high kick dataset initial guesses\n");
      N = 50.e9/24., tau = set_gammatau, ampl = 0.229, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.;
      A_cbo = 2.062e6/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.25088e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
    }    

    // high rate sum calos  (note important difference in time and cbo_freq when fitting lk day dataset)
    if ( strcmp( "lk", dataname) == 0 ) {
      printf("high kick dataset initial guesses\n");
      N = 50.e9/24., tau = set_gammatau, ampl = 0.229, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.;
      A_cbo = 2.062e6/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.25088e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
    }    

    // best guesses for calorimeter sum, partial endgame
    // tg != -> == feb 13, 2020
    if ( strcmp( "endgame", dataname) == 0 ) {
      printf("endgame dataset initial guesses\n");
      N = 25.e8, tau = set_gammatau, ampl = 0.226, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; 
      A_cbo = 2.52e6/N, tau_cbo = 1.371e5*rawBinToNs, omega_cbo = 2.90131e-03/rawBinToNs, phi_cbo = 1.7369, delta_omega_cbo = -0.004099/rawBinToNs;
    }
  }
  
  if ( strcmp( "runscan", cFit) == 0 ) { // best guesses for individual calos, 60hr dataset
    N = 24./62.*1.985e7, tau = set_gammatau, ampl = 0.217, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; 
    A_cbo = 6.15e5/N, tau_cbo = 1.261e5*rawBinToNs, omega_cbo = 2.890e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
  }

  if ( ( strcmp( "pairscan", cFit) == 0 ) || ( strcmp( "paircalo", cFit) == 0 ) || ( strcmp( "caloscan", cFit) == 0 ) || ( strcmp( "singlecalo", cFit) == 0 ) ) { 
    
    // best guesses for individual calos, 60hr dataset
    if ( strcmp( "60hr", dataname) ==  0 ) {
      printf("60hr dataset initial guesses\n");
      N = 5.17e8/24, tau = set_gammatau, ampl = 0.266, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; 
      A_cbo = 5*1.85e6/24/N, tau_cbo = 1.361e5*rawBinToNs, omega_cbo = 2.90131e-03/rawBinToNs, phi_cbo = 1.7369, delta_omega_cbo = -0.004099/rawBinToNs;
    }
    
    // high rate sum calos  (note important difference in time and cbo_freq when fitting 9 day dataset)
    if ( strcmp( "9day", dataname) == 0 ) {
      printf("9day dataset initial guesses\n");
      N = 1.e9/24., tau = set_gammatau, ampl = 0.269, omega_a = 1.79991e-3/rawBinToNs, phi = 2.23, B = 0.; 
      A_cbo = 6.062e5/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.24188e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00389/rawBinToNs;
    }

    if ( strcmp( "hk", dataname) == 0 ) {
      printf("hk dataset initial guesses\n");
      N = 1.e9/24., tau = set_gammatau, ampl = 0.269, omega_a = 1.79991e-3/rawBinToNs, phi = 2.23, B = 0.; 
      A_cbo = 6.062e5/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.24188e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00389/rawBinToNs;
    }

    if ( strcmp( "lk", dataname) == 0 ) {
      printf("hk dataset initial guesses\n");
      N = 1.e9/24., tau = set_gammatau, ampl = 0.269, omega_a = 1.79991e-3/rawBinToNs, phi = 2.23, B = 0.; 
      A_cbo = 6.062e5/N, tau_cbo = 1.55e5*rawBinToNs, omega_cbo = 3.24188e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00389/rawBinToNs;
    }

    if ( strcmp( "endgame", dataname) ==  0 ) {
      printf("endgame dataset initial guesses\n");
      N =8.e8, tau = set_gammatau, ampl = 0.267, omega_a = 1.79991e-3/rawBinToNs, phi = 2.32, B = 0.; 
      A_cbo = 0.01, tau_cbo = 1.351e5*rawBinToNs, omega_cbo = 2.868e-03/rawBinToNs, phi_cbo = 1.7369, delta_omega_cbo = -0.004399/rawBinToNs;
    }
  }

  if ( strcmp( "sequencescan", cFit) == 0 ) { // best guesses for individual calos, 60hr dataset
    N = 24./8.*1.985e7, tau = set_gammatau, ampl = 0.217, omega_a = 1.79991e-3/rawBinToNs, phi = 3.036, B = 0.; 
    A_cbo = 6.15e5/N, tau_cbo = 1.261e5*rawBinToNs, omega_cbo = 2.890e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
  }

  if ( strcmp( "verticalscan", cFit) == 0 || strcmp( "verticalslice", cFit) == 0 ) { // best guesses for individual calos, 60hr dataset
    N = 24./9.*1.985e7, tau = set_gammatau, ampl = 0.217, omega_a = 1.79991e-3/rawBinToNs, phi = 3.036, B = 0.; 
    A_cbo = 6.15e5/N, tau_cbo = 1.261e5*rawBinToNs, omega_cbo = 2.890e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
  }

  if ( strcmp( "horizontalscan", cFit) == 0 || strcmp( "horizontalslice", cFit) == 0 ) { // best guesses for individual calos, 60hr dataset
    N = 24./6.*1.985e7, tau = set_gammatau, ampl = 0.217, omega_a = 1.79991e-3/rawBinToNs, phi = 3.036, B = 0.; 
    A_cbo = 6.15e5/N, tau_cbo = 1.261e5*rawBinToNs, omega_cbo = 2.890e-03/rawBinToNs, phi_cbo = 1.740, delta_omega_cbo = -0.00457/rawBinToNs;
  }

  double blindR = 0.0;

  /*
From Denver talk, the histograms on page 11, 12:
The upper left plot is the loss function (L(t)) area normalized. In your fit function this should be included. The upper right plot  is the loss function normalized to the number of decay electrons, for this distribution, I plot N_loss/N0*Exp(t/64.4), where N0 is the number of decay positrons in the time bin of 30 us and t is the decay time.
And the lower right plot is an example of what you'll see if you take the upper left plot (L(t)) without any special normalization and multiply the exp(t'/tau') and integrate over t0  to t. 
So the correction for loss function Lamba (t) = 1 - C exp (-t0/tau) integral [my loss histogram * exp(t'/tau') dt']
So the lower right plot actually shows the expression in red above.
  */

  // empirical handling of muon loss 
  //double p0 =-8.5670e5, p1=1.74573e6, p2=-6.2772e4, p3=4.8396e2;
  double p0 = 250.0e6, p1=1.74573e6, p2 =0.0, p3 = 0.0;
  
  //muonlossf = new TF1( "muonlossf", "[0] + [1]*sqrt(x*1.25e-3 - 30.) + [2]*(x*1.25e-3 - 30.) + [3]*(x*1.25e-3 - 30.)*(x*1.25e-3 -30.)", 45000, 210000. ); // *1.25e-3 converts from Sudeshna use'd  to raw bins (800MHz)
  muonlossf = new TF1( "muonlossf", "[0] + [1]*sqrt(x) + [2]*(x) + [3]*(x^2)", 20., 210. ); // *1.25e-3 converts from Sudeshna usec's to raw bins (800MHz)
  muonlossf->SetParNames( "p0", "p1", "p2");
  muonlossf->SetParameters( p0, p1, p2, p3);
  muonlossf->SetLineColor(kBlue);

  // hmuonloss is the correction to the positron time distribution for lossing parent muons
  printf("get muonloss histogram\n");
  fmuonloss->GetObject( "hI", hmuonloss);
  hmuonloss->SetName("hmuonloss");
  hmuonloss->Fit( muonlossf, "", "VRIE", 20., 210. );
  hmuonloss->Draw("hist");
  printf("got muonloss histogram\n");

  // hmuonloss is the correction to the positron time distribution for lossing parent muons
  // only do for full, single calo, calo scan due to availability of qHist1D_*
  char fpname[64];
  sprintf(fpname, "qHist1D_%i_%i", icalo, iswrb/2);
  if ( ( strcmp( "singlecalo", cFit) == 0 ) ) fpedestaldrift->GetObject( fpname, hpedestaldrift);
  sprintf(fpname, "qHist1D_%i_%i", Ifit, iswrb/2);
  if ( ( strcmp( "caloscan", cFit) == 0 ) ) fpedestaldrift->GetObject( fpname, hpedestaldrift);
  if ( ( strcmp( "full", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "horizontalscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); // approx handling of ringing
  if ( ( strcmp( "verticalscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); // approx handling of ringing
  if ( ( strcmp( "horizontalslice", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); // approx handling of ringing
  if ( ( strcmp( "verticalslice", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); // approx handling of ringing
  if ( ( strcmp( "tauscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "starttimescan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "thresholdscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "windowscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "puscan", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift);
  if ( ( strcmp( "top", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); //fix me
  if ( ( strcmp( "bot", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_0", hpedestaldrift); //fix me
  if ( ( strcmp( "near", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_near_0", hpedestaldrift);
  if ( ( strcmp( "far", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_far_0", hpedestaldrift);
  if ( ( strcmp( "hot", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_hot_0", hpedestaldrift);
  if ( ( strcmp( "cold", cFit) == 0 ) ) fpedestaldrift->GetObject( "qHist1D_sum_cold_0", hpedestaldrift);
  if ( ( strcmp( "pairscan", cFit) == 0 ) || ( strcmp( "paircalo", cFit) == 0 ) ) {
    sprintf(hpname,"qHist1D_pair_%i_0",icalo+1);
    fpedestaldrift->GetObject( hpname, hpedestaldrift);
  }
  if ( ( strcmp( "caloscan", cFit) == 0 ) || ( strcmp( "singlecalo", cFit) == 0 ) ) { 
    TDirectoryFile *dirp;
    fpedestaldrift->GetObject( dname, dirp);
    sprintf(hpname,"qHist1D_%i_0",icalo+1);
    dirp->GetObject( hpname, hpedestaldrift);
  }
  hpedestaldrift->SetName("hpedestaldrift");
  printf("got pedestal histogram \n");

  htmp5 = (TH1D*) hpedestaldrift->Clone();
  htmp5->Reset();  
  htmp5b = (TH1D*) hpedestaldrift->Clone();
  htmp5b->Reset();  
  int pwind = 4, pgap = 1, ipb, jpb;
  for (ipb = 1; ipb <= hpedestaldrift->GetNbinsX(); ipb++){
  int pcnt = 0;
  double pedsampleaverage = 0.0;
    for (jpb = ipb - pwind - pgap; jpb <= ipb + pwind + pgap; jpb++){
      if (jpb >= 1 && jpb <= hpedestaldrift->GetNbinsX()) {// bookend histogram bins
	if (jpb < ipb - pgap)  { // exclude trigger / gap samples
	  pedsampleaverage += hpedestaldrift->GetBinContent(jpb);
          pcnt++;
	}
	if (jpb > ipb + pgap)  { // exclude trigger / gap samples
	  pedsampleaverage += hpedestaldrift->GetBinContent(jpb);
          pcnt++;
	}
      }
    }
    pedsampleaverage /= 2*pwind;
    htmp5->SetBinContent( ipb, pedsampleaverage);  
    htmp5b->SetBinContent( ipb, pedsampleaverage);  
  }
  htmp5b->Add( hpedestaldrift, -1.0);  

  hpedestaldrift->Draw("hist");
  htmp5->SetLineColor(kRed);
  htmp5->Draw("histsame");
  htmp5b->SetLineColor(kGreen);
  htmp5b->Draw("histsame");
  hpedestaldrift->Reset();
  hpedestaldrift->Add( htmp5b, 1.0);
  // basic idea is that we caculate the leakage of pedestal drift/ringing through trigger - pedestal calc
   
  // hlostmuontimedist is the time distribution for the lost muons themselves
  fmuonloss->GetObject( "time_1", hlostmuontimedist);
  hlostmuontimedist->SetName("hlostmuontimedist");

  double A_muonloss = set_muonloss_amp; // for muon loss effect on positron time distribution
  double A_detectmuon = set_detectmuon_amp; // for muon loss effect of detecting lost muons
  double A_pedestaldrift = set_pedestaldrift_amp; // for pedestal drift correction in fit


  printf("define hmlcorrection");
  hmlcorrection = (TH1F*)hmuonloss->Clone(); // histogram of lost muon correction
  hmlcorrection->SetName("hmlcorrection");
  hmlcorrection->Reset();
  hdmcorrection = (TH1F*)hmuonloss->Clone(); // histogram of detected muon correction
  hdmcorrection->SetName("hdmcorrection"); 
  hdmcorrection->Reset();
  higcorrection = (TH1F*)hmuonloss->Clone(); // histogram of infill gain correction
  higcorrection->SetName("higcorrection");
  higcorrection->Reset();
  printf("defined hmlcorrection");

  // define swrbglobal in order to inspect its value in root
  if ( fastRotationAnalysis )   swrbglobal = 2*1600;
  if ( (swrb == 1 || swrb == 2) && !fastRotationAnalysis )   swrbglobal = 1600;
  if (swrb == 4 && !fastRotationAnalysis ) swrbglobal = 1600/2;

  // make bins in ns with t=0 corresponding to injection
  resett0 = 2.49e4; // injection time in 1.25 ns bins - run 1
  if (run_2) {    
    resett0 = 1.048e5; // injection time in 1.25 ns bins - run 2
  }

  resetbins = htmp->GetNbinsX();  // bins in wiggle plot (can by 60ct, 120ct, 240 ct)
  resetloedge = htmp->GetBinLowEdge(1); // lower range of histogram
  resethiedge = htmp->GetBinLowEdge(resetbins)+htmp->GetBinWidth(resetbins); // upper range of histogram
  printf("hsig: resetbins %i, resetloedge %f, resethiedge %f\n", resetbins, resetloedge, resethiedge);
  htmp->SetBins(resetbins,rawBinToNs*(resetloedge-resett0), rawBinToNs*(resethiedge-resett0));
  //if (fastRotationCorrection) htmp2->SetBins(resetbins,rawBinToNs*(resetloedge-resett0), rawBinToNs*(resethiedge-resett0));

  resetbins = hpedestaldrift->GetNbinsX(); 
  resetloedge = hpedestaldrift->GetBinLowEdge(1);
  resethiedge = hpedestaldrift->GetBinLowEdge(resetbins)+hpedestaldrift->GetBinWidth(resetbins);
  printf("hped: resetbins %i, resetloedge %f, resethiedge %f\n", resetbins, resetloedge, resethiedge);
  hpedestaldrift->SetBins(resetbins,rawBinToNs*(resetloedge-resett0),rawBinToNs*(resethiedge-resett0));
  htmp5->SetBins(resetbins,rawBinToNs*(resetloedge-resett0),rawBinToNs*(resethiedge-resett0));
  htmp5b->SetBins(resetbins,rawBinToNs*(resetloedge-resett0),rawBinToNs*(resethiedge-resett0));
  minB = htmp->FindBin(startTime);
  maxB = htmp->FindBin(endTime);


 printf("begin fast rotation correction %i, analysis %i\n", fastRotationCorrection, fastRotationAnalysis);

 /* 

for correct fast rotation correction to 120ct data use (this avoids the mistake below)

fastRotationCorrection=1
fastRotationAnalysis=1

this does correction to 60ct data then rebins (rather than rebin then correcction)

for no fast rotation correction with 120ct data use 

fastRotationCorrection=0
fastRotationAnalysis=0

this avoids correction and uses 120ct rebinnied histogram

for obtaining residuals for calculating the FR correction do

fastRotationCorrection=0
fastRotationAnalysis=1

 */

 filefr = TFile::Open(fnamefastrotation);

 // fast rotation correction for individual calo's
 if ( fastRotationCorrection ) {
   printf("begin fast rotation correction for individual calos %i %i, file %s\n", fastRotationCorrection, fastRotationAnalysis, fnamefastrotation);
   
   if ( strcmp( "caloscan", cFit) == 0 || strcmp( "singlecalo", cFit) == 0 ) {
     
     if ( !fastRotationAnalysis ) sprintf(hname,"hfull%i",icalo+1);
     if ( fastRotationAnalysis ) sprintf(hname,"hhalf%i",icalo+1); // fast rotation analysis, testing correction on 60ct data
     printf("do fast rotation correction - individual calo's, calo %i, histogram %s\n", icalo+1, hname);
     filefr->GetObject( hname, htmp2);
     printf("time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );
     
     if (fixbug){
       // FIX for bug in raw bin - decimated bin correspondance, Oct21, 2018
       //printf("htmp2 before adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
       
       htmp2->SetBins(3168/swrb,20001,210021+60);
       if ( (swrb == 1 || swrb == 2) && !fastRotationAnalysis ) htmp2->SetBins(1584,20001,210021+60);
       if (swrb == 4 ) htmp2->SetBins(1584/2,20001,210021+60);
       //printf("htmp2 after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fast rotation histogram for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) htmp2->AddBinContent(ib,1.0);
     htmp->Multiply(htmp2);
     
     if ( fastRotationAnalysis ) { // only need to do for a fast rotation analysis
       htmp->Rebin(2); // test
       htmp2->Rebin(2); //test
       rb *= 2;
     }
   } 
 }
 printf("after fast rotation correction for individual calos %i %i\n", fastRotationCorrection, fastRotationAnalysis);

      
 if ( fastRotationCorrection ) { // for calo sum
   if ( strcmp( "full", cFit) == 0) {
     printf("begin fast rotation correction for summed calos %i %i, file %s\n", fastRotationCorrection, fastRotationAnalysis, fnamefastrotation);
     
     if ( !fastRotationAnalysis ) sprintf(hname,"hfull%i",1);
     if ( fastRotationAnalysis ) sprintf(hname,"hhalf%i",icalo+1); // fast rotation analysis, testing correction on 60ct data
     printf("do fast rotation correction - sum calo's, calo %i, histogram %s\n", icalo+1, hname);
     filefr->GetObject( hname, htmp2);
     printf("time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );
     
     if (fixbug){
       // FIX for bug in raw bin - decimated bin correspondance, Oct21, 2018
       //printf("htmp2 before adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
       if ( (swrb == 1 || swrb == 2) && !fastRotationAnalysis ) htmp2->SetBins(1584,20001,210021+60);
       if (swrb == 4 ) htmp2->SetBins(1584/2,20001,210021+60);
       //printf("htmp2 after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
     }
     
     htmp2->Scale(calofrac[0]); // is proprtion of individual calo contribution to fast roation correction 
     for (int icfr = 2; icfr <= 24; icfr++){
       sprintf(hname,"hfull%i",icfr);
       if ( fastRotationAnalysis ) sprintf(hname,"hhalf%i",icfr); // fast rotation analysis, testing correction on 60ct data
       filefr->GetObject( hname, htmp3);
       
       if (fixbug){
	 // FIX for bug in raw bin - decimated bin correspondance, Oct21, 2018
	 //printf("htmp3 before adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
	 if ( (swrb == 1 || swrb == 2) && !fastRotationAnalysis ) htmp3->SetBins(1584,20001,210021+60);
	 if (swrb == 4 ) htmp3->SetBins(1584/2,20001,210021+60);
	 //printf("htmp3 after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
       }
       
       htmp2->Add( htmp3, calofrac[icfr-1]);
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fast rotation histogram for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) {
       htmp2->AddBinContent(ib,1.0);
       htmp2->SetBinError(ib,0.0);
     } 

     printf("htmp after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp->GetNbinsX(), htmp->GetBinWidth(1), htmp->GetBinLowEdge(1), htmp->GetBinLowEdge( htmp->GetNbinsX() + htmp->GetBinWidth(htmp->GetNbinsX()) ) );
     printf("htmp2 after adjust: nbins, bin width, lower edge, upper edge %i %f, %f, %f\n", htmp2->GetNbinsX(), htmp2->GetBinWidth(1), htmp2->GetBinLowEdge(1), htmp2->GetBinLowEdge( htmp2->GetNbinsX() + htmp2->GetBinWidth(htmp->GetNbinsX()) ) );
     htmp->Multiply(htmp2);     
     printf("done multiply");
     
     if ( fastRotationAnalysis ) { // only need to do for fast rotation analysis 
       htmp->Rebin(2); // test
       htmp2->Rebin(2); //test
       rb *= 2;
     }

   }
   
   if ( strcmp( "hot", cFit) == 0 && !fastRotationAnalysis ) {

     
     sprintf(hname,"hfull%i",1);
     filefr->GetObject( hname, htmp2);
     htmp2->Reset();
     
     printf("hot calos, time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );

     for (int icfr = 1; icfr <= 24; icfr++){
       if (icfr > 8 && icfr < 24) continue; 
       sprintf(hname,"hfull%i",icfr);
       filefr->GetObject( hname, htmp3);       
       htmp2->Add( htmp3, calofrac[icfr-1]);
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fas rotation correction for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) {
       htmp2->AddBinContent(ib,1.0);
       htmp2->SetBinError(ib,0.0);
     }

     htmp->Multiply(htmp2);    
   }

   if ( strcmp( "near", cFit) == 0 && !fastRotationAnalysis ) {
     
     sprintf(hname,"hfull%i",1);
     filefr->GetObject( hname, htmp2);
     htmp2->Reset();
     
     printf("near calos, time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );

     for (int icfr = 1; icfr <= 12; icfr++){
       sprintf(hname,"hfull%i",icfr);
       filefr->GetObject( hname, htmp3);       
       htmp2->Add( htmp3, calofrac[icfr-1]);
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fas rotation correction for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) {
       htmp2->AddBinContent(ib,1.0);
       htmp2->SetBinError(ib,0.0);
     }

     htmp->Multiply(htmp2);    
   }

   if ( strcmp( "cold", cFit) == 0 && !fastRotationAnalysis ) {
     printf("begin fast rotation correction for cold calos %i %i, file %s\n", fastRotationCorrection, fastRotationAnalysis, fnamefastrotation);
     
     sprintf(hname,"hfull%i",1);
     filefr->GetObject( hname, htmp2);
     htmp2->Reset();
     printf("cold calos, time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );
     
     for (int icfr = 1; icfr <= 24; icfr++){
       if (icfr < 7 || icfr > 18) continue; 
       sprintf(hname,"hfull%i",icfr);
       filefr->GetObject( hname, htmp3);       
       htmp2->Add( htmp3, calofrac[icfr-1]);
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fas rotation correction for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) {
       htmp2->AddBinContent(ib,1.0);
       htmp2->SetBinError(ib,0.0);
     }
     htmp->Multiply(htmp2);     
   }

   if ( strcmp( "far", cFit) == 0 && !fastRotationAnalysis ) {
     printf("begin fast rotation correction for far calos %i %i, file %s\n", fastRotationCorrection, fastRotationAnalysis, fnamefastrotation);
     
     sprintf(hname,"hfull%i",1);
     filefr->GetObject( hname, htmp2);
     htmp2->Reset();
     printf("far calos, time histogram bins %i, corr histogram bins %i\n", htmp->GetNbinsX(), htmp2->GetNbinsX() );
     
     for (int icfr = 13; icfr <= 24; icfr++){

       sprintf(hname,"hfull%i",icfr);
       filefr->GetObject( hname, htmp3);       
       htmp2->Add( htmp3, calofrac[icfr-1]);
     }
     
     htmp2->Scale(-1.0); // undo fast rotation effects
     htmp2b = (TH1D*) htmp2->Clone(); // clone fas rotation correction for testing
     htmp2->Scale(scalefastrotation); // scale fast rotation effects by energy in fit range
     for (int ib = 1; ib <= htmp2->GetNbinsX(); ib++) {
       htmp2->AddBinContent(ib,1.0);
       htmp2->SetBinError(ib,0.0);
     }
     htmp->Multiply(htmp2);     
   }

 }
 
 printf("after fast rotation correction for summed calos %i %i\n", fastRotationCorrection, fastRotationAnalysis);

 //char funcname[32];
 //sprintf(funcname,"ffrpw");
 //filefr->GetObject( funcname, ffastrotation);
 //BuildFrInt();

 if (doConvoluted) {
   precessf = new TF1("precessf", fprecconv, 0.0, 250000.0, 47); // convoluted 0.5 + 1.0 + 0.5
 } else {
   precessf = new TF1("precessf", fprec, 0.0, 250000.0, 47); // unconvoluted
 }
 precessf->SetParameters( N, tau, ampl, blindR, phi, B, A_cbo, tau_cbo, omega_cbo, phi_cbo, A_gain);
 precessf->SetParNames( "norm", "tau", "amplitude", "blindR", "phi", "B", "A_cbo", "tau_cbo", "omega_cbo", "phi_cbo", "A_gain");

 minTglobal = minT;
 maxTglobal = maxT;
  
 //ROOT::Math::WrappedMultiTF1 fitFunction(*precessf, precessf->GetNdim() );
 //fitter.SetFunction( fitFunction, false);
 
  // chisquare function
  //fitter.SetFCN(chisquare);

  // wiggle parameters
  //precessf->FixParameter(0,N); // time-dilated muon lifetime
  //precessf->FixParameter(1,tau); // time-dilated muon lifetime
  //precessf->FixParameter(2,ampl); // wiggle amplitude
  //precessf->FixParameter(3,omega_a); // anomalous precession frequency
  //precessf->FixParameter(4,phi); // anomalous precession phase
  //precessf->SetParLimits( 4, -2.*3.141593, 2.*3.141593); // limit range of wiggle phase

  double delta_omega_a = set_delta_omega_a; // just a test - 60 hr didnt want a delta_omega_a 
  precessf->SetParName(29,"delta_omega_a");
  if (fix_delta_omega_a) precessf->FixParameter(29,delta_omega_a); // anomalous precession frequency change


  // empirical background term for tests
  if (fix_R) precessf->FixParameter(3,set_R);

  // empirical background term for tests
  B = set_bkd;
  precessf->FixParameter(5,B);

  // empirical gain correction term for gain sag
  precessf->SetParameter(10, A_gain);
  precessf->SetParName(10, "A_gain");
  precessf->SetParameter(11, tau_gain);
  precessf->SetParName(11, "tau_gain");
  precessf->FixParameter(10, A_gain); // switch off gain sag term
  precessf->FixParameter(11, tau_gain); // fix gain sag time constant
 
  // muon loss parameter
  precessf->SetParName( 21, "A_muonloss");
  precessf->SetParameter( 21, A_muonloss); 
  //precessf->FixParameter( 21, A_muonloss);

  // additional muon loss parameter for detected muons
  precessf->SetParName( 30, "A_detectmuon");
  precessf->SetParameter( 30, A_detectmuon); 
  //precessf->FixParameter( 30, A_detectmuon);

  // pedestal drift / ringing parameter
  precessf->SetParName( 35, "A_pedestaldrift");
  precessf->SetParameter( 35, A_pedestaldrift); 
  //precessf->FixParameter( 35, A_pedestaldrift);


  // CBO modulation of normalization

  // tracker data
  // < 15970
  //2.3051 rad/musec = 2.3051 x 10^-3 rad/nsec = 1.25*2.3051 x 10^-3 = 0.00288138 
  // 
  // > 15970
  //2.3046 rad/musec = 2.3046 x 10^-3 rad/nsec = 1.25*2.3046 x 10^-3 = 0.00288075 

  // 60hr dataset
  if ( strcmp( "60hr", dataname) == 0 ) omega_cbo = 0.00288138/rawBinToNs; // tracker data
  //if ( strcmp( "60hr", dataname) == 0 ) omega_cbo = 2.88670e-03;   // best fit calo sum

  //9day dataset
  if ( strcmp( "9day", dataname) == 0 ) omega_cbo = 3.25088e-03/rawBinToNs;

  //hk dataset
  if ( strcmp( "hk", dataname) == 0 ) omega_cbo = 3.25088e-03/rawBinToNs;

  //lk dataset
  if ( strcmp( "lk", dataname) == 0 ) omega_cbo = 3.25088e-03/rawBinToNs;

  // endgame dataset
  if ( strcmp( "60hr", dataname) == 0 ) omega_cbo = 0.00288138/rawBinToNs; // tracker data

  if (fix_omega_cbo) omega_cbo = set_omega_cbo; 
  precessf->FixParameter(8, omega_cbo); // fix  cbo freq

  // CBO freq change for norm, asymmetry, and phase terms

  // tracker data
  // < 15970
  // 1.25*1.86x10^-8 = 2.32500e-08
  //
  // > 15970
  // 1.25*1.91x10^-8 = 2.38750e-08

  // 60hr dataset
  //if ( strcmp( "60hr", dataname) == 0 ) delta_omega_cbo = 2.32500e-08/rawBinToNs; // tracker value
  //if ( strcmp( "9day", dataname) == 0 ) delta_omega_cbo =  1.4822e-8/rawBinToNs; // best fit to calo sum
  //if ( strcmp( "endgame", dataname) == 0 ) delta_omega_cbo = 2.32500e-08/rawBinToNs; // tracker value
  if ( strcmp( "60hr", dataname) == 0 ) delta_omega_cbo = set_delta_omega_cbo; // tracker value
  if ( strcmp( "9day", dataname) == 0 ) delta_omega_cbo =  set_delta_omega_cbo; // best fit to calo sum
  if ( strcmp( "hk", dataname) == 0 ) delta_omega_cbo =  set_delta_omega_cbo; // best fit to calo sum
  if ( strcmp( "lk", dataname) == 0 ) delta_omega_cbo =  set_delta_omega_cbo; // best fit to calo sum
  if ( strcmp( "endgame", dataname) == 0 ) delta_omega_cbo = set_delta_omega_cbo; // tracker value

  precessf->SetParameter(20, delta_omega_cbo);
  precessf->SetParName(20, "delta_omega_cbo");
  if (!fix_delta_omega_cbo) precessf->ReleaseParameter(20); // release  cbo freq

  //fitter.Config().ParSettings(20).SetStepSize(0.00001*delta_omega_cbo); // set step size, this parameter is tough to fit
  //printf("print step size %e\n",fitter.Config().ParSettings(20).StepSize());
  precessf->FixParameter(20, delta_omega_cbo); // fix  cbo linear freq change

  // A*exp(-t/tau) term of CBO freq in tracker parameterization
  double AExpTermCBO = 0.0;
  precessf->SetParName( 23, "AExpTermCBO");
  precessf->SetParameter( 23, AExpTermCBO); 
  precessf->FixParameter( 23, AExpTermCBO);
  //precessf->FixParameter( 23, 0.0); // set exp term to zero 

  // B*exp(-t/tau) term of CBO freq in tracker parameterization
  double BExpTermCBO = 0.0;
  precessf->SetParName( 24, "BExpTermCBO");
  precessf->SetParameter( 24, BExpTermCBO); 
  precessf->FixParameter( 24, BExpTermCBO); 
  //precessf->FixParameter( 24, 0.0); // set exp term to zero 

  // CBO envelope time constant for norm, asymmetry, and phase terms
  if (fix_tau_cbo) tau_cbo = set_tau_cbo;
  if (!fix_tau_cbo) precessf->ReleaseParameter(7); // SWITCH ON cbo tau envelope variation 
  precessf->FixParameter(7, tau_cbo); // fix cbo tau envelope

  // amplitude of CBO modulation of wiggle asymmetry 
  precessf->SetParameter(12,A2_cbo);
  precessf->SetParName(12,"A2_cbo (for A_a)");
  //precessf->FixParameter(12,0.0); // SWITCH OFF CBO asymmetry term

  // phase of CBO modulation of wiggle asymmetry 
  precessf->SetParameter(13, phi2_cbo);
  precessf->SetParName(13, "phi2_cbo (for A_a)");
  precessf->SetParLimits( 13, -2*3.141593, 2.*3.141593); // limit range of phase // TG TEST 9/3/19
  //precessf->FixParameter(13,0.0); // SWITCH OFF CBO asymmetry term

  // amplitude of CBO modulation of wiggle phase (not needed in 60-hr dataset) 
  precessf->SetParameter(14,A3_cbo);
  precessf->SetParName(14,"A3_cbo (for phi_a)");
  precessf->FixParameter(14,0.0); // SWITCH OFF CBO phase

  // phase of CBO modulation of wiggle phase (not needed in 60-hr dataset) 
  precessf->SetParameter(15, phi3_cbo);
  precessf->SetParName(15, "phi3_cbo (for phi_a)");
  precessf->SetParLimits( 15, -2*3.141593, 2.*3.141593); // limit range of phase 
  precessf->FixParameter(15,0.0); // SWITCH OFF CBO phase

  // amplitude of 2*CBO modulation of normalization (helpful for individual calos) 
  precessf->SetParameter(36,A_2cbo);
  precessf->SetParName(36,"A_2cbo (for norm)");
  precessf->FixParameter(36,0.0); // SWITCH OFF CBO phase

  // phase of CBO modulation of wiggle phase (not needed in 60-hr dataset) 
  precessf->SetParameter(37, phi_2cbo);
  precessf->SetParName(37, "phi_2cbo (for norm)");
  precessf->SetParLimits( 37, -2*3.141593, 2.*3.141593); // limit range of phase 
  precessf->FixParameter(37, 0.0); // SWITCH OFF CBO phase

  // pedestal ringing term
  //f *= (1.0 - par[31]/par[0] * exp( -( xx - t0 ) / par[33] ) * sin( par[32] * ( xx - t0 ) + par[34] ) );

  double A_ring = set_fitring_amp, omega_ring = set_fitring_omega, phi_ring = 2.8, tau_ring = set_fitring_tau; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
  precessf->SetParameter(31,A_ring);
  precessf->SetParName(31, "A_ring");
  precessf->FixParameter(31, 0); // SWITCH OFF ringing
  precessf->SetParName(32, "omega_ring");
  precessf->FixParameter(32, omega_ring); // SWITCH OFF ringing
  precessf->SetParameter(33, tau_ring);
  precessf->SetParName(33, "tau_ring");
  precessf->FixParameter(33, tau_ring); // SWITCH OFF ringing
  precessf->SetParameter(34, phi_ring);
  precessf->SetParName(34, "phi_ring");
  precessf->FixParameter(34, phi_ring); // SWITCH OFF ringing

  double A_ring2 = set_fitring2_amp, omega_ring2 = set_fitring2_omega, phi_ring2 = 2.8, tau_ring2 = set_fitring2_tau; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
  precessf->SetParameter(49,A_ring2);
  precessf->SetParName(40, "A_ring2");
  precessf->FixParameter(40, 0); // SWITCH OFF ringing2
  precessf->SetParName(41, "omega_ring2");
  precessf->FixParameter(41, omega_ring2); // SWITCH OFF ringing2
  precessf->SetParameter(42, tau_ring2);
  precessf->SetParName(42, "tau_ring2");
  precessf->FixParameter(42, tau_ring2); // SWITCH OFF ringing2
  precessf->SetParameter(43, phi_ring2);
  precessf->SetParName(43, "phi_ring2");
  precessf->FixParameter(43, phi_ring2); // SWITCH OFF ringing2

  // relaxation parameter
  double asymtau = set_asymtau;
  precessf->SetParameter(44, asymtau);
  precessf->SetParName(44, "asymtau");
  precessf->FixParameter(44, asymtau); // SWITCH OFF ringing2

  // 2*omega_a parameter
  precessf->SetParameter(45, 0.0);
  precessf->SetParName(45, "A_2omega_a");
  precessf->FixParameter(45, 0.0); 
  precessf->SetParameter(46, 0.0);
  precessf->SetParName(46, "phi_2omega_a");
  precessf->FixParameter(46, 0.0); 
   
  // amplitude of vertical BO modulation of wiggle normalization, VW term???
  precessf->SetParameter(16,A_vbo);
  precessf->SetParName(16, "A_vbo");
  precessf->FixParameter(16, 0.0); // SWITCH OFF vertical BO

  // vertical waist term
  // frequency of horizontal BO modulation of wiggle normalization
  precessf->SetParameter(17, omega_vbo);
  precessf->SetParName(17, "omega_vbo");
  precessf->FixParameter(17, omega_vbo); // SWITCH OFF vertical BO

  // time constant of vertical BO modulation of wiggle normalization 
  precessf->SetParameter(22, tau_vbo);
  precessf->SetParName(22, "tau_vbo");
  precessf->FixParameter(22, tau_vbo); // SWITCH OFF vertical BO

  // phase of horizontal BO modulation of wiggle normalization 
  precessf->SetParameter(18, phi_vbo);
  precessf->SetParName(18, "phi_vbo");
  //precessf->SetParLimits( 18, -2*3.141593, 2.*3.141593);
  precessf->FixParameter(18, 2*1.57); // SWITCH OFF vertical BOfix_gain_amp = 0

  omega_vbo2 = 1.3339e-2, phi_vbo2 = 2.8, tau_vbo2 = 1.83563e+05; // vertical BO ang freq 1.25*2.*3.14159/(449.ns)
  A_vbo2asym = 0.0, phi_vbo2asym = 0.0; // asymmetry term in vbo as seen in FFT

  // obtain tau's for vbo1, vbo2 from difference in top, bottom vertical slices of calos using vohorizontal() function
  //if ( strcmp( "60hr", dataname) == 0 ) tau_vbo2 = 2.122e4;
  //if ( strcmp( "60hr", dataname) == 0 ) omega_vbo2 = 1.30968e-2;

  // summary for vertical frequencies
  // "60hr"   omega_vbo = 1.77073e-2; //fc-2fy 
  // "60hr"   omega_vbo2 = 1.3339e-2; //fy
  /// 60hr ratio 1.33
  // "9day"   omega_vbo = 1.5745e-2;  //fc-2fy
  // "9day"   omega_vbo2 = 1.8332e-2; //fy
  //  9day ratio 0.86

  // amplitude of vertical BO modulation 2 of wiggle normalization, not VW term???
  precessf->SetParameter(25,A_vbo2);
  precessf->SetParName(25, "A_vbo2");
  precessf->FixParameter(25, 0.0); // SWITCH OFF vertical BO 2

  // frequency of horizontal BO modulation 2 of wiggle normalization, not VW term???
  precessf->SetParameter(26, omega_vbo2);
  precessf->SetParName(26, "omega_vbo2");
  precessf->FixParameter(26, omega_vbo2); // SWITCH OFF vertical BO 2

  // time constant of vertical BO modulation 2 of wiggle normalization, not VW term???
  precessf->SetParameter(27, tau_vbo2);
  precessf->SetParName(27, "tau_vbo2");
  precessf->FixParameter(27, tau_vbo2); // SWITCH OFF vertical BO 2

  // phase of horizontal BO modulation of 2 wiggle normalization, not VW term???
  precessf->SetParameter(28, phi_vbo2);
  precessf->SetParName(28, "phi_vbo2");
  precessf->FixParameter(28, phi_vbo2); // SWITCH OFF vertical BO 2

  // asymmetry term - amplitude of vertical BO modulation 2 of wiggle normalization, not VW term???
  precessf->SetParameter(38,A_vbo2asym);
  precessf->SetParName(38, "A_vbo2asym");
  precessf->FixParameter(38, 0.0); // SWITCH OFF vertical BO 2

  // asymmetry - phase of horizontal BO modulation of 2 wiggle normalization, not VW term???
  precessf->SetParameter(39, phi_vbo2asym);
  precessf->SetParName(39, "phi_vbo2asym");
  precessf->FixParameter(39, phi_vbo2asym); // SWITCH OFF vertical BO 2

  // empirical rate effect parameter (in place of using a pileup correction term)
  double R_pileup = 0.0;
  precessf->SetParameter(19, R_pileup);
  precessf->SetParName(19, "R_pileup");
  precessf->FixParameter(19, 0.0); //  SWITCH off pileup term

  // fit function draw options
  precessf->SetLineColor(kRed);
  precessf->SetLineWidth(2);
  precessf->SetNpx(10000);

  // gain term helpers  
  //precessf->FixParameter(10,-0.007); // switch off gain change term correlates too strongly with muon loss term
  //precessf->FixParameter(10,0.0); // switch off gain change term correlates too strongly with muon loss term (small compared to pu?)
  //precessf->ReleaseParameter(10); // switch on gain term 
  
  // muon loss term helpers
  //precessf->FixParameter(21,1.6e-9); // pu corrected value for muon loss term (typical -ileup corrected fit value
  //precessf->FixParameter(21, precessf->GetParameter(21)); // muon loss
  //precessf->ReleaseParameter(21); // switch on muon loss term 
  //precessf->FixParameter(21,0.0); // switch off muon loss term 

  // cbo helpers
  //precessf->FixParameter(18, precessf->GetParameter(18)); // Asym cbo
  //precessf->ReleaseParameter(20); // switch on cbo frequency time dependence variation 

   if ( strcmp( "60hr", dataname) == 0 ) omega_vbo = set_omega_vbo_60hr; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
   if ( strcmp( "60hr", dataname) == 0 ) tau_vbo = set_tau_vbo_60hr;
   if ( strcmp( "60hr", dataname) == 0 ) A_vbo = set_amp_vbo_60hr;
   if ( strcmp( "60hr", dataname) == 0 ) phi_vbo = set_phi_vbo_60hr;

   if ( strcmp( "9day", dataname) == 0 ) omega_vbo = set_omega_vbo_9day; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "9day", dataname) == 0 ) tau_vbo = set_tau_vbo_9day;
   if ( strcmp( "9day", dataname) == 0 ) A_vbo = set_amp_vbo_9day;
   if ( strcmp( "9day", dataname) == 0 ) phi_vbo = set_phi_vbo_9day;

   if ( strcmp( "hk", dataname) == 0 ) omega_vbo = set_omega_vbo_hk; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "hk", dataname) == 0 ) tau_vbo = set_tau_vbo_hk;
   if ( strcmp( "hk", dataname) == 0 ) A_vbo = set_amp_vbo_hk;
   if ( strcmp( "hk", dataname) == 0 ) phi_vbo = set_phi_vbo_hk;

   if ( strcmp( "lk", dataname) == 0 ) omega_vbo = set_omega_vbo_lk; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "lk", dataname) == 0 ) tau_vbo = set_tau_vbo_lk;
   if ( strcmp( "lk", dataname) == 0 ) A_vbo = set_amp_vbo_lk;
   if ( strcmp( "lk", dataname) == 0 ) phi_vbo = set_phi_vbo_lk;

   if ( strcmp( "endgame", dataname) == 0 ) omega_vbo = set_omega_vbo_endgame; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
   if ( strcmp( "endgame", dataname) == 0 ) tau_vbo = set_tau_vbo_endgame;
   if ( strcmp( "endgame", dataname) == 0 ) A_vbo = set_amp_vbo_endgame;
   if ( strcmp( "endgame", dataname) == 0 ) phi_vbo = set_phi_vbo_endgame;

   if ( strcmp( "60hr", dataname) == 0 ) omega_vbo2 = set_omega_vbo2_60hr; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
   if ( strcmp( "60hr", dataname) == 0 ) tau_vbo2 = set_tau_vbo2_60hr;
   if ( strcmp( "60hr", dataname) == 0 ) A_vbo2 = set_amp_vbo2_60hr;
   if ( strcmp( "60hr", dataname) == 0 ) phi_vbo2 = set_phi_vbo2_60hr;

   if ( strcmp( "9day", dataname) == 0 ) omega_vbo2 = set_omega_vbo2_9day; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "9day", dataname) == 0 ) tau_vbo2 = set_tau_vbo2_9day;
   if ( strcmp( "9day", dataname) == 0 ) A_vbo2 = set_amp_vbo2_9day;
   if ( strcmp( "9day", dataname) == 0 ) phi_vbo2 = set_phi_vbo2_9day;

   if ( strcmp( "hk", dataname) == 0 ) omega_vbo2 = set_omega_vbo2_hk; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "hk", dataname) == 0 ) tau_vbo2 = set_tau_vbo2_hk;
   if ( strcmp( "hk", dataname) == 0 ) A_vbo2 = set_amp_vbo2_hk;
   if ( strcmp( "hk", dataname) == 0 ) phi_vbo2 = set_phi_vbo2_hk;

   if ( strcmp( "lk", dataname) == 0 ) omega_vbo2 = set_omega_vbo2_lk; // based on FFT of absolute residuals and peak at 477 -> T_vbo = 190020./477. = 398.365 ct, omega_cbo = 2pi/T_vbo = 0.0157724
   if ( strcmp( "lk", dataname) == 0 ) tau_vbo2 = set_tau_vbo2_lk;
   if ( strcmp( "lk", dataname) == 0 ) A_vbo2 = set_amp_vbo2_lk;
   if ( strcmp( "lk", dataname) == 0 ) phi_vbo2 = set_phi_vbo2_lk;

   if ( strcmp( "endgame", dataname) == 0 ) omega_vbo2 = set_omega_vbo2_endgame; // based on FFT of absolute residuals and peak at 534 -> T_vbo = 190020./534. = 355.843 ct, omega_cbo = 2pi/T_vbo = 0.0176573
   if ( strcmp( "endgame", dataname) == 0 ) tau_vbo2 = set_tau_vbo2_endgame;
   if ( strcmp( "endgame", dataname) == 0 ) A_vbo2 = set_amp_vbo2_endgame;
   if ( strcmp( "endgame", dataname) == 0 ) phi_vbo2 = set_phi_vbo2_endgame;

  precessf->FixParameter(7, tau_cbo);  // no tau variation for cbo terms
  precessf->FixParameter(6, 0.0);  // no cbo normalization term
  precessf->FixParameter(9, phi_cbo);  // no phase variation for cbo normalization term
  precessf->FixParameter(12, 0.0); // no cbo asymmetry term
  precessf->FixParameter(13, phi2_cbo);  // no phase variation for cbo asymmetry term
  precessf->FixParameter(14, 0.0); // no cbo phase term
  precessf->FixParameter(15, phi3_cbo);  // no phase variation for cbo asymmetry term

  precessf->FixParameter(16, 0.0); // vertical bo term
  precessf->FixParameter(17, omega_vbo);  //no freq variation for vertical bo term
  precessf->FixParameter(18, phi_vbo);  // no phase variation for vertical bo term
  precessf->FixParameter(22, tau_vbo);  // no tau variation for vertical bo term

  precessf->FixParameter(25, 0.0); // vertical bo term
  precessf->FixParameter(26, omega_vbo2);  //no freq variation for vertical bo term
  precessf->FixParameter(28, phi_vbo2);  // no phase variation for vertical bo term
  precessf->FixParameter(27, tau_vbo2);  // no tau variation for vertical bo term

  precessf->FixParameter(10, 0.0); // no empircal gain term
  precessf->FixParameter(19, 0.0); // no empircal pileup term
  precessf->FixParameter(21, 0.0); // no muon loss term
  precessf->FixParameter(30, 0.0); // no detected muon term
  precessf->FixParameter(35, 0.0); // no pedestal drift
  precessf->SetParLimits( 4, -2.*3.141593, 2.*3.141593); // limit range of wiggle phase 

  // update for new energy scale / energy calibration in gm2 analysis (23.62 is scale change to make energies, thresholds in MeV)
  N *= 23.63;
  //A_cbo *= 23.63; 
  //A_vbo *= 23.63; 
  //A_vbo2 *= 23.63; 

  precessf->SetParameter(0,htmp->Integral(htmp->FindBin(30000),htmp->FindBin(215000))/254); // empirical estimate N for 150ns binning
  precessf->FixParameter(1, tau);
  if (fix_asym) ampl = set_asym;
  precessf->FixParameter(2, ampl); // muon lifetime
  if (!fix_R) precessf->SetParameter(3, set_R);

  // for drift correction
  if (doDriftCorrection) htmp=(TH1D*)htmpdriftcorrection->Clone();

  printf("**** 5par, non-bin-integral fit ****\n");
  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit

  precessf->SetParameter(0,htmp->Integral(htmp->FindBin(30000),htmp->FindBin(215000))/254); // empirical estimate N for 150ns binning
  precessf->FixParameter(1, tau);
  precessf->FixParameter(2, ampl); // muon lifetime
  if (!fix_R) precessf->SetParameter(3, set_R);

  printf("**** 5par, bin-integral fit ****\n");
  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit

  if (!fix_gammatau) precessf->ReleaseParameter(1); // muon lifetime
  if (!fix_asym) precessf->ReleaseParameter(2); // wiggle asym 

  printf("**** 5par, non-bin-integral fit ****\n");
  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit

  if (!fix_gammatau) precessf->ReleaseParameter(1); // muon lifetime
  if (!fix_asym) precessf->ReleaseParameter(2); // wiggle asym 
 
  printf("**** 5par, bin-integral fit ****\n");
  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit

  precessf->SetParameter(6, A_cbo);
  if (fix_cbo_amp) precessf->SetParameter(6, set_cbo_amp);

  if ( strcmp( "5par", fitname) != 0 ) {

    precessf->FixParameter( 0, precessf->GetParameter(0));
    precessf->FixParameter( 1, precessf->GetParameter(1));
    precessf->FixParameter( 2, precessf->GetParameter(2));

    // old model
    //omega_cbo = 2.89875e-03/rawBinToNs; // no A,B exp terms
    //delta_omega_cbo = 2.60823e-8/rawBinToNs; // no A,B exp terms
    AExpTermCBO = 0.0; // no A,B exp terms
    BExpTermCBO = 0.0; // no A,B exp terms

    if (use_new_cbomodel) {   // old model from James Mott email
      if ( strcmp( "60hr", dataname) == 0 ) AExpTermCBO = set_AExpTermCBO_60hr; // 81.8 us lifetime
      if ( strcmp( "60hr", dataname) == 0 ) BExpTermCBO = set_BExpTermCBO_60hr; //  7.7 us lifetime
      if ( strcmp( "9day", dataname) == 0 ) AExpTermCBO = set_AExpTermCBO_9day; // 81.8 us lifetime
      if ( strcmp( "9day", dataname) == 0 ) BExpTermCBO = set_BExpTermCBO_9day; //  7.7 us lifetime
      if ( strcmp( "hk", dataname) == 0 ) AExpTermCBO = set_AExpTermCBO_hk; // 81.8 us lifetime
      if ( strcmp( "hk", dataname) == 0 ) BExpTermCBO = set_BExpTermCBO_hk; //  7.7 us lifetime
      if ( strcmp( "lk", dataname) == 0 ) AExpTermCBO = set_AExpTermCBO_lk; // 81.8 us lifetime
      if ( strcmp( "lk", dataname) == 0 ) BExpTermCBO = set_BExpTermCBO_lk; //  7.7 us lifetime
      if ( strcmp( "endgame", dataname) == 0 ) AExpTermCBO = set_AExpTermCBO_eg; // from fit to calo
      if ( strcmp( "endgame", dataname) == 0 ) BExpTermCBO = set_BExpTermCBO_eg; //  7.7 us lifetime
    } else { // old model from Joe
      AExpTermCBO = -5.04e-2; // no A,B exp terms
      BExpTermCBO = -13.10e-2; // no A,B exp terms
    }

    //precessf->FixParameter(7, tau_cbo);  // no tau variation for cbo terms
    //precessf->FixParameter(8, omega_cbo); // set cbo freq
    //precessf->FixParameter(20, delta_omega_cbo); // set  cbo linear freq chage
    precessf->FixParameter(23, AExpTermCBO); // set  cbo freq exp A change
    precessf->FixParameter(24, BExpTermCBO); // set  cbo freq exp B change
    
    precessf->SetParameter(6, A_cbo);  // no cbo normalization term
    precessf->ReleaseParameter(6);  // cbo normalization term
    if (fix_cbo_amp) precessf->FixParameter(6, set_cbo_amp);
    precessf->SetParameter(9, phi_cbo);  // no phase variation for cbo normalization term
    precessf->ReleaseParameter(9);  // phase for cbo normalization term
    
    if (doBinInts) {
      htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
    } else {
      htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
    }

    precessf->ReleaseParameter(0);
    if (!fix_gammatau) precessf->ReleaseParameter(1);
    if (!fix_asym) precessf->ReleaseParameter(2);
    if (!fix_tau_cbo) precessf->ReleaseParameter(7);
    if (!fix_omega_cbo) precessf->ReleaseParameter(8);  // freq for cbo normalization term

    if (doBinInts) {
      htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
    } else {
      htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
    }

    precessf->SetParameter(16, 0.0); //vbo off unless vbo fit
    precessf->SetParameter(25, 0.0); //vbo2 off unless vbo2 fit
    
    if ( strcmp( "5par", fitname) != 0 && strcmp( "9par", fitname) != 0 ) { 
      
      if (!fix_delta_omega_cbo)  precessf->ReleaseParameter(20); // fix cbo dfreq term

      if (doBinInts) {
	htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
      } else {
	htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
      }
      
      if ( strcmp( "5par", fitname) != 0 && strcmp( "9par", fitname) != 0 && strcmp( "10par", fitname) != 0 ) { 
	
	//precessf->FixParameter(7, precessf->GetParameter(7));  // fix tau for cbo normalization term
	//precessf->FixParameter(8, precessf->GetParameter(8));  // fix freq for cbo normalization term
	precessf->SetParameter( 21, A_muonloss); 
	if (!fix_muonloss_amp) precessf->ReleaseParameter( 21);
	precessf->SetParameter( 30, A_detectmuon); 
	if (!fix_detectmuon_amp) precessf->ReleaseParameter( 30);
	precessf->SetParameter( 35, A_pedestaldrift); 
	if (!fix_pedestaldrift_amp) precessf->ReleaseParameter( 35);
	
	if (doBinInts) {
	  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	} else {
	  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	}

	if (!fix_gain_amp) precessf->SetParameter( 10, 0.0); // along empirical gain amplitude to vary
	if (!fix_gain_amp) precessf->ReleaseParameter( 10); // along empirical gain amplitude to vary
	if (!fix_gain_tau) precessf->SetParameter( 11, tau_gain); // along empirical gain amplitude to vary
	if (!fix_gain_tau) precessf->ReleaseParameter( 11); // along empirical gain amplitude to vary
	
	if (doBinInts) {
	  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	} else {
	  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	}
	
	if ( strcmp( "5par", fitname) != 0 && strcmp( "9par", fitname) != 0 && strcmp( "10par", fitname) != 0 && strcmp( "11par", fitname) != 0 ) { 
      
	  //precessf->ReleaseParameter(8); // release  cbo freq
	  //precessf->ReleaseParameter(20); // release  cbo dfreq
	  //precessf->SetParLimits(20,  0.5*1.86e-8, 2.0*1.86e-8); // limit range of wiggle phase
	  
	  //htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (released freq, dfreq) + muonloss fit
	  
          if (!fix_anomasym_amp) {
	    precessf->SetParameter(12, A2_cbo); // cbo asymmetry term
	    precessf->SetParameter(13, phi2_cbo);  // phase variation for cbo asymmetry term
	    precessf->ReleaseParameter(12);  // cbo asymmetry term
	    precessf->ReleaseParameter(13);  // phase variation for cbo asymmetry term
	  }

	  if (doBinInts) {
	    htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  } else {
	    htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  }

          if (!fix_anomphi_amp) {
	    precessf->SetParameter(14, A3_cbo); // cbo asymmetry term
	    precessf->SetParameter(15, phi3_cbo);  // phase variation for cbo asymmetry term
	    precessf->ReleaseParameter(14);  // cbo asymmetry term
	    precessf->ReleaseParameter(15);  // phase variation for cbo asymmetry term
	  }

	  
	  if (doBinInts) {
	    htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  } else {
	    htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  }
	  
          if (!fix_2cbo_amp) {
	    precessf->SetParameter(36, A_2cbo); // cbo asymmetry term
	    precessf->SetParameter(37, phi_2cbo);  // phase variation for cbo asymmetry term
	    precessf->ReleaseParameter(36);  // cbo asymmetry term
	    precessf->ReleaseParameter(37);  // phase variation for cbo asymmetry term
	  }

	  if (doBinInts) {
	    htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  } else {
	    htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  }

	  if (!fix_pileup) precessf->ReleaseParameter(19); // add-in pileup term
	  
	  if (doBinInts) {
	    htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  } else {
	    htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	  }

	  if ( strcmp( "5par", fitname) != 0 && strcmp( "9par", fitname) != 0 && strcmp( "10par", fitname) != 0 && strcmp( "11par", fitname) != 0 && strcmp( "13par", fitname) != 0 ) { 
	    
	    if ( strcmp( "vbo", fitname) == 0 || strcmp( "vbo2", fitname) == 0 ) { // for early start time add vw term

	      //release ringing term    
              if (!fix_fitring_amp) {
		precessf->FixParameter(31, set_fitring_amp);
		precessf->ReleaseParameter(31); // ringing ampl
		precessf->ReleaseParameter(34); //ringing phase	    
		if (!fix_fitring_omega) precessf->ReleaseParameter(32); // vbo freq
		if (!fix_fitring_tau) precessf->ReleaseParameter(33); // vbo tau
		if (doBinInts) {
		  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		} else {
		  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		}
	      }

	      //release ringing2 term    
              if (!fix_fitring2_amp) {
		precessf->FixParameter(40, set_fitring2_amp);
		precessf->ReleaseParameter(40); // ringing ampl
		precessf->ReleaseParameter(43); //ringing phase	    
		if (!fix_fitring2_omega) precessf->ReleaseParameter(41); // vbo freq
		if (!fix_fitring2_tau) precessf->ReleaseParameter(42); // vbo tau
		if (doBinInts) {
		  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		} else {
      htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		}
	      }

              //fix vbo amplitudes   	      
	      precessf->FixParameter(16, A_vbo); // vertical bo term
	      precessf->FixParameter(25, A_vbo2); // vertical bo term

	      //fix wiggle
	      precessf->FixParameter(0, precessf->GetParameter(0)); 
	      precessf->FixParameter(1, precessf->GetParameter(1)); 
	      precessf->FixParameter(2, precessf->GetParameter(2)); 
	      precessf->FixParameter(3, precessf->GetParameter(3)); 
	      precessf->FixParameter(4, precessf->GetParameter(4)); 
	      //fix norm cbo
	      precessf->FixParameter(6, precessf->GetParameter(6)); 
	      if (!fix_tau_cbo) precessf->FixParameter(7, precessf->GetParameter(7)); 
	      precessf->FixParameter(8, precessf->GetParameter(8)); 
	      precessf->FixParameter(9, precessf->GetParameter(9)); 
	      precessf->FixParameter(20, precessf->GetParameter(20)); 
	      //fix asym cbo
	      precessf->FixParameter(12, precessf->GetParameter(12)); 
	      precessf->FixParameter(13, precessf->GetParameter(13));
              // pileup
              set_pileup = precessf->GetParameter(19);
	  
	      //release vw term    
	      if (!fix_vbo_amp) precessf->ReleaseParameter(16); //vbo ampl
	      if (!fix_vbo_phi) precessf->ReleaseParameter(18); //vbo phase
	      //precessf->SetParLimits( 18, 0, 2.*3.141593); // limit range of wiggle phase
	       
	      if (doBinInts) {
		htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      } else {
		htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      }
	      
	      if (!fix_vbo_tau) precessf->ReleaseParameter(22); // vbo tau
	      if (!fix_vbo_omega) precessf->ReleaseParameter(17); // vbo freq
	       
	      if (doBinInts) {
		htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      } else {
		htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      }

	      //early time additional VBO2
	      
	      if ( strcmp( "vbo2", fitname) == 0 ) { // very early times has additional freq (whatis is it?)
				
		// release vbo2 term parameters
		if (!fix_vbo_amp2) precessf->ReleaseParameter(25); //vbo2 ampl
		if (!fix_vbo_phi2) precessf->ReleaseParameter(28); //vbo2 phase
		 
		if (doBinInts) {
		  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		} else {
		  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		}
		
		if (!fix_vbo_omega2) precessf->ReleaseParameter(26); //vbo2 freq
		if (!fix_vbo_tau2) precessf->ReleaseParameter(27); //vbo2 envelope time constant
		 
		if (doBinInts) {
		  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		} else {
		  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
		}

		precessf->FixParameter(25, precessf->GetParameter(25)); 
		precessf->FixParameter(26, precessf->GetParameter(26)); 
		precessf->FixParameter(27, precessf->GetParameter(27)); 
		precessf->FixParameter(28, precessf->GetParameter(28)); 
	      }
	      
	      //unfix wiggle
	      precessf->ReleaseParameter(0); 
	      if (!fix_gammatau) precessf->ReleaseParameter(1); 
	      if (!fix_asym) precessf->ReleaseParameter(2); 
	      if (!fix_R) precessf->ReleaseParameter(3); 
	      precessf->ReleaseParameter(4); 
	      //unfix norm cbo
	      precessf->ReleaseParameter(6); 
	      if (fix_cbo_amp) precessf->FixParameter(6, set_cbo_amp);
	      if (!fix_tau_cbo) precessf->ReleaseParameter(7); 
	      if (!fix_omega_cbo) precessf->ReleaseParameter(8); 
	      precessf->ReleaseParameter(9); 
	      if (!fix_delta_omega_cbo) precessf->ReleaseParameter(20); // release  cbo freq variation
	      if (!fix_AExpTermCBO) precessf->ReleaseParameter(23);
	      if (!fix_BExpTermCBO) precessf->ReleaseParameter(24);
	      //unfix asym cbo
	      if (!fix_anomasym_amp) precessf->ReleaseParameter(12); 
	      if (!fix_anomasym_amp) precessf->ReleaseParameter(13); 
              // unfix pileup
	      if (!fix_pileup) precessf->ReleaseParameter(19);

	    }
	    
	    if (doBinInts) {
	      htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	    } else {
	      htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	    }

	    // free vbo2 parameters and fit if vbo2 fit
	    if (strcmp( "vbo2", fitname) == 0) {
	      
	      if (!fix_vbo_amp2) precessf->ReleaseParameter(25);
	      if (!fix_vbo_phi2) precessf->ReleaseParameter(28);
	      if (!fix_vbo_omega2) precessf->ReleaseParameter(26);
	      if (!fix_vbo_tau2) precessf->ReleaseParameter(27);
	      
	      if (!fix_vbo_Aasym2) precessf->ReleaseParameter(38); // asymmetry term in vbo2
	      if (!fix_vbo_Aasym2) precessf->ReleaseParameter(39); // asymmetry term in vbo2
	      
	      if (doBinInts) {
		htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      } else {
		htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	      }
	    }
      
      sprintf( foutname, "graphs-w%i-swrb%i-prd%s-icalo%i-Ncbo1-Acbo1-df1.root",  iwind, iswrb, aname, icalo);
            
      // dont yet know how to set the errors on certain datasets
      if ( ( strcmp( "sequencescan", cFit) == 0 ) || ( strcmp( "hotthresholdscan", cFit) == 0 ) || ( strcmp( "coldthresholdscan", cFit) == 0 )) {
	for (int ib =  minB; ib <= maxB; ib++){
	  // define bin error using chi-squared and fit function
	  htmp->SetBinError(  ib, htmp->GetBinError( ib) * sqrt(precessf->GetChisquare() / precessf->GetNDF()) );
	}
      }
      
      //htmp->Fit( precessf, "", "NRI", minT, maxT);  
      
      // iterate as dont yet know how to set the errors on certain datasets
      if ( ( strcmp( "sequencescan", cFit) == 0 ) || ( strcmp( "hotthresholdscan", cFit) == 0 ) || ( strcmp( "coldthresholdscan", cFit) == 0 )) {
	for (int ib =  minB; ib <= maxB; ib++){
	  // define bin error using chi-squared and fit function  for improperly analyzed sequnce number histograms     
	  htmp->SetBinError(  ib, htmp->GetBinError( ib) * sqrt(precessf->GetChisquare() / precessf->GetNDF()) );
	}
      }
      
      // helper for fitting horizontal slices with large VBO on top, bottom slices
      if ( strcmp( "horizontalscan", cFit) == 0 || strcmp( "horizontalslice", cFit) == 0 ) {
	
	precessf->FixParameter(0, precessf->GetParameter(0)); // wiggle 
	precessf->FixParameter(1, precessf->GetParameter(1)); // wiggle
	precessf->FixParameter(2, precessf->GetParameter(2)); // wiggle
	precessf->FixParameter(3, precessf->GetParameter(3)); // wiggle
	precessf->FixParameter(4, precessf->GetParameter(4)); // wiggle
	if (!fix_tau_cbo) precessf->FixParameter(7, precessf->GetParameter(7)); // tau cbo
	precessf->FixParameter(10, precessf->GetParameter(10)); // Norm cbo  
	precessf->FixParameter(12, precessf->GetParameter(12)); // Asym cbo
	precessf->FixParameter(13, precessf->GetParameter(13)); // Asym cbo
	precessf->ReleaseParameter(16); // SWITCH on vertical BO amp
	precessf->ReleaseParameter(18); // SWITCH on vertical BO phase
	
	if (doBinInts) {
	  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	} else {
	  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	}

	precessf->ReleaseParameter(0); // wiggle
	if (!fix_gammatau) precessf->ReleaseParameter(1); // wiggle
	if (!fix_asym) precessf->ReleaseParameter(2); // wiggle
	if (!fix_R) precessf->ReleaseParameter(3); // wiggle
	precessf->ReleaseParameter(4); // wiggle
	if (!fix_tau_cbo) precessf->ReleaseParameter(7); // tau cbo
	if (!fix_gain_amp) precessf->ReleaseParameter(10); // norm cbo
	if (!fix_anomasym_amp) precessf->ReleaseParameter(12); // asym cbo
	if (!fix_anomasym_amp) precessf->ReleaseParameter(13); // asym cbo
	precessf->FixParameter(16, precessf->GetParameter(16)); //  SWITCH on vertical BO ampl
	precessf->FixParameter(18, precessf->GetParameter(18)); //  SWITCH on vertical BO hase
	
	if (doBinInts) {
	  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	} else {
	  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	}
	
	precessf->ReleaseParameter(16); // SWITCH on vertical BO ampl
	if (!fix_vbo_omega) precessf->ReleaseParameter(17); // SWITCH on vertical BO freq
	precessf->ReleaseParameter(18); // SWITCH on vertical BO phase
	
	if (doBinInts) {
	  htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	} else {
	  htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
	}
      }
 
      if (!fix_asymtau) precessf->ReleaseParameter(44); // release spin relaxation term
      if (!fix_bkd) precessf->ReleaseParameter(5); // release empirical background term

      if (!fix_2omega_a) { // free 2omega_a term
	precessf->ReleaseParameter(45);
	precessf->ReleaseParameter(46); 
      }
      
      printf("fit w/o using covariance matrix, integrate fn over bin-width %i\n", doBinInts);
      if (doBinInts) {
	htmp->Fit( precessf, "", "NRIE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
      } else {
	htmp->Fit( precessf, "", "NRE", minT, maxT); // 5par + cbo norm (fixed freq, dfreq) fit
      }
      
      // option U uses user fcn
      if (doUserCovariance) {
	gFitter = TVirtualFitter::Fitter( htmp);
	gFitter->SetFCN( chsq);  
	if (doMinos) {
	  if (!doBinInts) {
	    printf("fit w/ using covariance matrix, w/ minos, dont integrate fn over bin-width %i\n", doBinInts);
	    htmp->Fit(precessf,"RUE","", minT, maxT);
	  }
	  if (doBinInts) {
	    printf("fit w/ using covariance matrix, w/ minos, integrate fn over bin-width %i\n", doBinInts);
	    htmp->Fit(precessf,"RUIE","", minT, maxT);
	  }
	} else {
	  if (!doBinInts) {
	    printf("fit w/ using covariance matrix, w/o minos, dont integrate fn over bin-width %i\n", doBinInts);
	    htmp->Fit(precessf,"RUE","", minT, maxT);
	  }
	  if (doBinInts) {
	    printf("fit w/ using covariance matrix, w/o minos, integrate fn over bin-width %i\n", doBinInts);
	    htmp->Fit(precessf,"RUIE","", minT, maxT);
	  }
	}
	ndf_QM[Ifit] = fcn_ndf;
	chisq_QM[Ifit] = fcn_chi2;
        printf("covariance matrix based fit ndf %f, fcn %f from chsq()\n", ndf_QM[Ifit], chisq_QM[Ifit]);
      }

      if (doTimeShift) {
	printf("build graph time shifts\n");
	sprintf(hname,"hfullmean%i",icalo+1);
	filefr->GetObject( hname, htmp4);
	double gx[1600], gdx[1600], gy[1600], gdy[1600]; 
	for (int ib =  1; ib <= htmp->GetNbinsX(); ib++) {
	  gx[ib-1] = htmp->GetBinCenter( ib);
	  if (doTimeShift) gx[ib-1] += timeShiftMultiplier*htmp4->GetBinContent( ib);
	  gdx[ib-1] = 0.0;
	  gy[ib-1] = htmp->GetBinContent( ib);
	  gdy[ib-1] = htmp->GetBinError( ib);
	}
	gWiggle = new TGraphErrors(1600, gx, gy, gdx, gdy);
	gWiggle->Fit( precessf, "R", "", minT, maxT);
      }

      /*
      if (!doFRStudy) {
	filefrstudy = TFile::Open( "frstudy.root");
	sprintf(hname,"ds%i",icalo+1); 
	filefrstudy->GetObject( hname, htmp2);
	if (doShiftBins) htmp->Add(htmp2, -0.5);
	if (!doShiftBins) htmp->Add(htmp2, +0.5);
	
	printf("next fit after fr adjustment\n");
	htmp->Fit( precessf, "R", "", minT, maxT);
      }
      */      
      
      //printf("last fit w/ bin integration\n");
      //htmp->Fit( precessf, "IER", "", minT, maxT);

      // try improve fit results
      //htmp->Fit( precessf, "", "NRM", minT, maxT);

      // helper for horizontal slice fits

      // all fit, also includes vbo, vbo2 terms depending of fit start time      
      } // 13 parameter fit, 11 par  + ampl, phase for cbo asymmetry term
     } // 11 parameter fit, 10 par + muon loss term for cbo normalization term
    } // 10 parameter fit, 9 par with cbo freq change for cbo normalization term
   } // 9 parameter fit, wiggle + amp, tau, freq, phase of cbo normalization term
  } // 5 parameter fit, theoretical wiggle

  // ringing fit function
  ringf = new TF1("ringf", fring, 0.0, 220000.0, 4);
  ringf->SetParameters( precessf->GetParameter(31), precessf->GetParameter(33), precessf->GetParameter(32), precessf->GetParameter(34));
  ringf->SetParNames( "A_ring", "omega_ring", "tau_ring", "phi_ring");
  ringf->SetNpx(10000);

  // ringing fit function
  ring2f = new TF1("ring2f", fring, 0.0, 220000.0, 4);
  ring2f->SetParameters( precessf->GetParameter(40), precessf->GetParameter(42), precessf->GetParameter(41), precessf->GetParameter(43));
  ring2f->SetParNames( "A_ring2",  "omega_ring2", "tau_ring2", "phi_ring2");
  ring2f->SetNpx(10000);

  // slow terms function
  slowtermsf = new TF1("slowtermsf", fslowterms, 0.0, 220000.0, 5);
  slowtermsf->SetParameters( 1.0, precessf->GetParameter(1), precessf->GetParameter(10), precessf->GetParameter(11), precessf->GetParameter(21) );
  slowtermsf->SetParNames( "N0", "tau_mu", "A_gain", "tau_gain", "muonloss");
  slowtermsf->SetNpx(10000);

  NF_QM[Ifit] = precessf->GetParameter(0);
  dNF_QM[Ifit] = precessf->GetParError(0);
  tauF_QM[Ifit] = precessf->GetParameter(1);
  dtauF_QM[Ifit] = precessf->GetParError(1);
  AF_QM[Ifit] = precessf->GetParameter(2); 
  AF_QM[Ifit] = sqrt(AF_QM[Ifit]*AF_QM[Ifit]); // handle +/- sign
  dAF_QM[Ifit] = precessf->GetParError(2);
  omega_aF_QM[Ifit] = precessf->GetParameter(3); // parameter change to blind R
  domega_aF_QM[Ifit] = precessf->GetParError(3); // parameter change to blind R
  if ( strcmp( "60hr", dataname) == 0 ) omega_aF_QM[Ifit] -= 5.09020; // fix to independent blinding of 60 hr dataset

  phaseF_QM[Ifit] = precessf->GetParameter(4);
  //if (phaseF_QM[Ifit] <= 0.0) phaseF_QM[Ifit] += 3.141593; // deal with phase ambiquity
  dphaseF_QM[Ifit] = precessf->GetParError(4);
  BF_QM[Ifit] = precessf->GetParameter(5);
  dBF_QM[Ifit] = precessf->GetParError(5);

     // time dependent cbo leading-order nomalization term - exponential envelope
     //f *= (1.0 - par[6]/par[0] * exp( -( xx - t0 ) / par[7] ) * sin( omega_cbo * ( xx - t0 )   +  par[9] ) );

  ACBOF_QM[Ifit] = abs(precessf->GetParameter(6));
  dACBOF_QM[Ifit] = precessf->GetParError(6);

     // vertical BO modifier, VW?
     //f *= (1.0 - par[16]/par[0] * exp( -( xx - t0 ) / par[22] ) * sin( par[17] * ( xx - t0 ) + par[18] ) );
     // vertical BO modifier 2, not VW?
     //f *= (1.0 - par[25]/par[0] * exp( -( xx - t0 ) / par[27] ) * sin( par[26] * ( xx - t0 ) + par[28] ) );

  AVWF_QM[Ifit] = precessf->GetParameter(16);
  dAVWF_QM[Ifit] = precessf->GetParError(16);
  omegaVWF_QM[Ifit] = precessf->GetParameter(17);
  domegaVWF_QM[Ifit] = precessf->GetParError(17);
  tauVWF_QM[Ifit] = precessf->GetParameter(22);
  dtauVWF_QM[Ifit] = precessf->GetParError(22);
  phaseVWF_QM[Ifit] = fmod( precessf->GetParameter(18), TMath::Pi() );
  dphaseVWF_QM[Ifit] = precessf->GetParError(18);
  AVW2F_QM[Ifit] = precessf->GetParameter(25);
  dAVW2F_QM[Ifit] = precessf->GetParError(25);
  omegaVW2F_QM[Ifit] = precessf->GetParameter(26);
  domegaVW2F_QM[Ifit] = precessf->GetParError(26);
  tauVW2F_QM[Ifit] = precessf->GetParameter(27);
  dtauVW2F_QM[Ifit] = precessf->GetParError(27);
  phaseVW2F_QM[Ifit] = fmod( precessf->GetParameter(28), TMath::Pi() );
  dphaseVW2F_QM[Ifit] = precessf->GetParError(28);

  ACBO2F_QM[Ifit] = precessf->GetParameter(12);
  dACBO2F_QM[Ifit] = precessf->GetParError(12);
  A2CBOF_QM[Ifit] = abs(precessf->GetParameter(36)); // stop it flipping sign
  dA2CBOF_QM[Ifit] = precessf->GetParError(36);
  ACBO3F_QM[Ifit] = precessf->GetParameter(14);
  dACBO3F_QM[Ifit] = precessf->GetParError(14);
  tauCBOF_QM[Ifit] = precessf->GetParameter(7);
  dtauCBOF_QM[Ifit] = precessf->GetParError(7);
  freqCBOF_QM[Ifit] = precessf->GetParameter(8);
  dfreqCBOF_QM[Ifit] = precessf->GetParError(8);
  deltafreqCBOF_QM[Ifit] = precessf->GetParameter(20);
  ddeltafreqCBOF_QM[Ifit] = precessf->GetParError(20);
  phaseCBOF_QM[Ifit] = fmod( precessf->GetParameter(9), TMath::Pi() );
  dphaseCBOF_QM[Ifit] = precessf->GetParError(9);
  phaseCBO2F_QM[Ifit] = fmod( precessf->GetParameter(13), TMath::Pi() );
  dphaseCBO2F_QM[Ifit] = precessf->GetParError(13);
  phaseCBO3F_QM[Ifit] = fmod( precessf->GetParameter(15), TMath::Pi() );
  dphaseCBO3F_QM[Ifit] = precessf->GetParError(15);
  phase2CBOF_QM[Ifit] = fmod( precessf->GetParameter(37), TMath::Pi() );
  dphase2CBOF_QM[Ifit] = precessf->GetParError(37);
  AGAINF_QM[Ifit] = precessf->GetParameter(10);
  dAGAINF_QM[Ifit] = precessf->GetParError(10);
  tauGAINF_QM[Ifit] = precessf->GetParameter(11);
  dtauGAINF_QM[Ifit] = precessf->GetParError(11);
  muonlossF_QM[Ifit] = precessf->GetParameter(21);
  dmuonlossF_QM[Ifit] = precessf->GetParError(21);
  peddriftF_QM[Ifit] = precessf->GetParameter(35);
  dpeddriftF_QM[Ifit] = precessf->GetParError(35);
  RPILEUPF_QM[Ifit] = precessf->GetParameter(19);
  dRPILEUPF_QM[Ifit] = precessf->GetParError(19);
  AFrelax_QM[Ifit] = precessf->GetParameter(44);
  dAFrelax_QM[Ifit] = precessf->GetParError(44);
  if ( !( strcmp( "vbo", fitname) == 0 || strcmp( "vbo2", fitname) == 0 ) || !doUserCovariance ) { // otherwise chisq calculated above
    ndf_QM[Ifit] = precessf->GetNDF();
    chisq_QM[Ifit] = precessf->GetChisquare();
    dchisq_QM[Ifit] = chisq_QM[Ifit]/sqrt(ndf_QM[Ifit]);
    printf("non-covariance chi-squared, ndf %f, chi2 %f\n", ndf_QM[Ifit], chisq_QM[Ifit]);
  }
  if (chisq_QM[Ifit] > 1.0) domega_aF_QM[Ifit] = domega_aF_QM[Ifit]*sqrt(chisq_QM[Ifit]); // inflate blind R error

  // inflate error bars appropriately

  if ( inflate_errors && ( chisq_QM[Ifit] > 1.0 ) ) domega_aF_QM[Ifit] = sqrt(chisq_QM[Ifit])*domega_aF_QM[Ifit];

  printf("Norm  %e +/- %e \n", NF_QM[Ifit], dNF_QM[Ifit] );
  printf("tau  %e +/- %e ns \n", tauF_QM[Ifit], dtauF_QM[Ifit] );
  printf("A  %e +/- %e \n", AF_QM[Ifit], dAF_QM[Ifit] );
  printf("omega_a  %12.7e +/- %e rad/ns \n", omega_aF_QM[Ifit], domega_aF_QM[Ifit] );
  printf("phi  %e +/- %e rad \n", phaseF_QM[Ifit], dphaseF_QM[Ifit] );
  printf("B  %e +/- %e \n", BF_QM[Ifit], dBF_QM[Ifit] );
  printf("omega_cbo  %12.7e +/- %e rad/ns \n", freqCBOF_QM[Ifit], dfreqCBOF_QM[Ifit] );
  printf("delta_omega_cbo  %12.7e +/- %e rad/ns \n", dfreqCBOF_QM[Ifit], ddeltafreqCBOF_QM[Ifit] );
  printf("chisq %e, ndf %e, chisq/ndf %e \n \n", chisq_QM[Ifit], ndf_QM[Ifit], chisq_QM[Ifit]/ndf_QM[Ifit] );
  printf("probability, %f", precessf->GetProb());

  hres=(TH1D*)htmp->Clone();
  //hres = new TH1D("hres","hres",maxB-minB+1,minB*htmp->GetBinWidth(1),(maxB+1)*htmp->GetBinWidth(1));
  hres->Reset();
  hres->SetName("hres");
  hres->SetTitle("(data - fit)/sigma residuals hits time distribution");
  hres->SetMarkerColor(kGray+2);
  hres->SetLineColor(kGray+2);
  TH1D *hresProj;
  if ( !fracResiUnits ) hresProj = new TH1D("hresProj","distribution of residuals (data - fit)/sigma", 1000, -2.5e5, +2.5e5 ); // absolute scale
  if ( fracResiUnits ) hresProj = new TH1D("hresProj","distribution of residuals (data - fit)/sigma", 100, -10, +10 ); // relative scale
  double resi;
  for (int ib = minB; ib <= maxB; ib++){
    resi = 0.0;
    //double funcbinaverage = precessf->Eval( htmp->GetBinCenter(ib) + htmp4->GetBinContent( ib) ); // if doing graph fit
    double funcbinaverage = precessf->Eval( htmp->GetBinCenter(ib) ); // of doing histo fit
    if ( !fracResiUnits && htmp->GetBinError(ib) > 0) resi = ( htmp->GetBinContent(ib)  - funcbinaverage ); // absolute units
    if ( fracResiUnits && htmp->GetBinError(ib) > 0) resi = ( htmp->GetBinContent(ib) - funcbinaverage ) / funcbinaverage; //  fraction units
    hres->SetBinContent( ib, resi );
    if ( !fracResiUnits ) hres->SetBinError( ib, htmp->GetBinError(ib) ); // absolute units
    if ( fracResiUnits ) hres->SetBinError( ib, htmp->GetBinError(ib) / funcbinaverage ); // fractional units
    hresProj->Fill( resi );
     }

  printf("finished building residuals\n");

  /*
  int fr_start = 60000, fr_end = 90000;
  Double_t A_fr = 0.1, tau_fr = 5.e4, omega_fr =2.*TMath::Pi()*1.25/149.7, phi_fr = -1.40+2.*TMath::Pi();
  omega_fr=0.0525;
  frf = new TF1("frf", ffr, fr_start, fr_end, 4);
  frf->SetNpx(10000);
  frf->SetParameters( A_fr, tau_fr, omega_fr, phi_fr);
  frf->SetParNames( "A_fr", "tau_fr", "omega_fr", "phi_fr");
  frf->FixParameter( 0, A_fr);
  //frf->FixParameter( 1, tau_fr);
  //frf->FixParameter( 2, omega_fr);
  //hres->Fit( frf, "", "RIE", fr_start, fr_end);
  */

  hres->SetStats(1);
  hresProj->SetStats(1);
  hresProj->SetLineWidth(3);
  hresProj->SetMarkerColor(kViolet+2);
  hresProj->SetLineColor(kViolet+2);
  hresProj->GetYaxis()->SetTitleOffset(0.5);
  hresProj->GetYaxis()->SetTitleSize(0.04);
  hresProj->GetXaxis()->SetTitleOffset(0.5);
  hresProj->GetXaxis()->SetTitleSize(0.04);
  hresProj->GetXaxis()->SetTitle("sigma");
  hresProj->GetYaxis()->SetTitle("counts");
  //htmp->SetMaximum(1.5e5);
  htmp->GetXaxis()->SetRange(minB-100, maxB+100);
  hres->GetXaxis()->SetRange(minB-100, maxB+100);       
  htmp->SetStats(1);

  printf("finished setting scales\n");

  //gStyle->SetOptStat("iourmrn");
  //gStyle->SetOptFit(1111);
  
  if (doDraw) {

    c1->Clear();
    if (FFTon) {
      c1->Divide(1,3);
    } else {
      c1->Divide(1,2);
    }
    
    printf("finished making canvas arrangement\n");
    
    int cinp;
    
    gStyle->SetOptFit(1111);
    c1->cd(1);
    htmp->GetXaxis()->SetRangeUser(30000.,215500.);
    //htmp->SetMaximum(1.25*htmp->GetBinContent(320/swrb));
    //htmp->SetMinimum(0.75*htmp->GetBinContent(1500/swrb));
    printf("draw time distribution\n");
    htmp->Draw("hist");
    htmp->Draw("psame");
    //printf( "continue?");
    //cinp = getchar();
    
    c1->cd(2);
    //hres->SetMaximum(+0.005);
    //hres->SetMinimum(-0.005);
    hres->GetXaxis()->SetRangeUser(30000., 215500.);
    //hres->GetXaxis()->SetRangeUser(50000, 80000); // see strange spike
    hres->SetTitle("(data - 5-par wiggle fit) residuals time distribution");
    printf("draw fit residuals\n");
    hres->Draw("HIST");
    //printf( "continue?");
    //cinp = getchar();
    
    if (FFTon) {
      
      c1->cd(3);
      //TH1 *hm2 =0; TVirtualFFT::SetTransform(0); hm2 = hres->FFT(hm2, "MAG"); hm2->SetTitle("discrete fourier transform");hm2->Draw();
      TVirtualFFT::SetTransform(0);
      hm2 = hres->FFT(hm2, "MAG");
      hm2->SetTitle("discrete fourier transform");
      hm2->SetBins(hm2->GetNbinsX(), 0, 1./75.e-9*1.e-6*TMath::Pi() ); // put in unuts of rad
      hm2->GetXaxis()->SetRange(0.,hm2->GetNbinsX()/2.); // plot is symmetric, just plot the lower half
      printf("draw residuals FFT\n");
      hm2->Draw("HIST");
      //printf( "continue?");
      //cinp = getchar();
      //hresProj->Draw();    
    }
    
    char foutsave[128];
    sprintf( foutsave, "png/graphs-w%i-t%i-g%i-swrb%i-prd%s-calo%i-%s-%s-%s-%s.png", iwind, ithreshold, igap, iswrb, aname, icalo+1, cFit, fitname, dataname, datasubset);
    c1->SaveAs(foutsave);

  } // end doDraw=1

  return 0;

  // save a graph
  file1 = TFile::Open( foutname, "UPDATE");
  if ( file1->IsOpen() ) printf("created output root file %s\n", foutname);

  hm2->Write();
  file1->Print();
  printf("wrote graph into file %s\n", foutname);  
  file1->Close();
  

  //sprintf( hname, "hslice%i-calo%02i-residuals-7par.png", Ifit, ICalo);
  //c1->SaveAs(hname);

  //sprintf( hname, "qHist1D_sig_All_fit%02i.png",Ifit);
  //c1->SaveAs(hname);
 
  /*
  char fname[64];
  sprintf( fname, "cbo-6par-pileupterm-residuals-lowrate.root");
  sprintf( fname, "cbo-6par-pileupterm-residuals-highrate.root");
  //sprintf( fname, "cbo-5par-residuals.root");
  file1 = TFile::Open( fname, "update");
  //sprintf( hname, "hot");
  sprintf( hname, "cbo_calo_%02i", Ifit);
  //sprintf( hname, "cbo_sumi");
  hres->SetName(hname);
  hres->SetTitle(hname);
  hres->Write();
  file1->Close();
  */
  

  /*  
  A_cbo = precessf->GetParameter(6)/precessf->GetParameter(0);
  tau_cbo = precessf->GetParameter(7);
  omega_cbo = precessf->GetParameter(8);
  phi_cbo = precessf->GetParameter(9);
  */
  A_cbo = 0.002;
  tau_cbo = 1.0e9;
  omega_cbo = 3.2285e-03;
  phi_cbo = 0.0;
  float p0_cbo = 0.0;

  cbof = new TF1("cbof", "[0] * exp( -( x - 2.7e4 ) / [1] ) * sin( [2]*x + [3] ) + [4]", minT, maxT);
  cbof->SetParameters( A_cbo, tau_cbo, omega_cbo, phi_cbo, p0_cbo);
  cbof->SetParNames( "A_cbo", "tau_cbo", "omega_cbo", "phi_cbo", "p0_cbo" );
  cbof->SetParLimits(3,0,2*3.14);
  //cbof->SetParLimits(0,0.0,1000000.);
  //cbof->FixParameter(0,A_cbo);
  cbof->FixParameter(1,tau_cbo);
  //cbof->FixParameter(2,omega_cbo);
  cbof->FixParameter(4,p0_cbo);


  cbof->SetLineWidth(2);
  cbof->SetLineColor(kRed);
  cbof->SetNpx(10000);
  hres->SetLineWidth(2);
  hres->Fit( cbof, "", "RIE", minT, maxT);
  int newminT = 58000, newmaxT = 195000;
  hres->Fit( cbof, "", "RIE", newminT, newmaxT);
  hres->GetXaxis()->SetRangeUser(50000, 160000);
  hres->Draw("hist"); 
  cbof->Draw("same"); 
  c1->Update();

  // for 5-par fit followed by cbo fit store CBO amplitude / phase
  ACBOF_QM[Ifit] = sqrt(cbof->GetParameter(0)*cbof->GetParameter(0));
  dACBOF_QM[Ifit] = cbof->GetParError(0);
  phaseCBOF_QM[Ifit] = cbof->GetParameter(3);
  dphaseCBOF_QM[Ifit] = cbof->GetParError(3);

  return 1;

  hres2=(TH1D*)hres->Clone();
  hres2->Reset();
  hres2->SetName("hres");
  hres2->SetTitle("(data - fit)/sigma CBO residuals distribution");
  hres2->SetMarkerColor(kGray+2);
  hres2->SetLineColor(kGray+2);
  hres2->GetXaxis()->SetRangeUser(30000, 205000);
  for (int ib = minB; ib <= maxB; ib++){
    if ( hres->GetBinError(ib) > 0) resi = ( hres->GetBinContent(ib) - cbof->Eval( hres->GetBinCenter(ib) )); // absolute units
    hres2->SetBinContent( ib, resi );
    hres2->SetBinError( ib, hres->GetBinError(ib) ); // absolute units
  }

  /*
  hres3=(TH1D*)hres->Clone();
  hres3->Reset();
  hres3->SetName("hres");
  hres3->SetTitle("(data - fit)/sigma CBO residuals modulu cbo period");
  float cboPeriodInBins = 2.*3.14159265/cbof->GetParameter(2)/hres->GetBinWidth(1);
  printf("cbo period in bins %f\n", cboPeriodInBins);
  for (int ib = minB; ib <= maxB; ib++){
    hres3->SetBinContent( fmod((float)ib,cboPeriodInBins), hres->GetBinContent(ib) );
  }
  hres3->Draw();
  return 1;
  */

  float  xmodulo[4000], ymodulo[4000];
  int NCBOPoints = maxB - minB + 1;
  NCBOPoints = 200;
  printf("make cbo x,y's\n");
  for (int ib = minB; ib <= minB + NCBOPoints - 1; ib++){
    float tInHistBins = ib*hres->GetBinWidth(1);
    float omegaTPlusPhi = tInHistBins*cbof->GetParameter(2) + cbof->GetParameter(3);
    xmodulo[ib-minB] = cbof->Eval( hres->GetBinCenter(ib) );
    ymodulo[ib-minB] = ( hres->GetBinContent(ib) - cbof->Eval( hres->GetBinCenter(ib)) );
    printf("ib, omegaTPlusPhi, sinomegaTPlusPhi), x, y (residual) %i, %f, %f, %f, %f \n",  ib, omegaTPlusPhi, sin(omegaTPlusPhi), xmodulo[ib-minB], ymodulo[ib-minB] );
  }

  printf("make cbo graph\n");
  c1->Clear();
  c1->Divide(1,2);
  g1 = new TGraph( NCBOPoints, xmodulo, ymodulo);
  g1->SetName("gcbomodulo");
  g1->SetMarkerStyle(8);
  g1->SetMarkerColor(kRed);
  c1->cd(1);
  hres->Draw();
  c1->cd(2);
  g1->Draw("AP");
  return 1;

  c1->Clear();
  c1->Divide(1,2);

  c1->cd(1);
  hres->Draw("");
  cbof->Draw("same");

  c1->cd(2);
  hres2->Draw("");

  sprintf( hname, "cboCalo%02i.png",Ifit);
  c1->SaveAs(hname);


  ACBOF_QM[Ifit] = cbof->GetParameter(0)/precessf->GetParameter(0);
  dACBOF_QM[Ifit] = cbof->GetParError(0)/precessf->GetParameter(0);

  tauCBOF_QM[Ifit] = cbof->GetParameter(1);
  dtauCBOF_QM[Ifit] = cbof->GetParError(1);

  freqCBOF_QM[Ifit] = cbof->GetParameter(2)/rawBinToNs;
  dfreqCBOF_QM[Ifit] = cbof->GetParError(2)/rawBinToNs;

  phaseCBOF_QM[Ifit] = cbof->GetParameter(3);
  dphaseCBOF_QM[Ifit] = cbof->GetParError(3);

  ndf_QM[Ifit] = cbof->GetNDF();
  chisq_QM[Ifit] = cbof->GetChisquare()/ndf_QM[Ifit];
  dchisq_QM[Ifit] = chisq_QM[Ifit]/sqrt(ndf_QM[Ifit]);


  return 1;

  // TH1 *hm =0;TVirtualFFT::SetTransform(0);hm = hres->FFT(hm, "MAG"); hm->SetTitle("discrete fourier transform");hm->Draw();
//Compute the transform and look at the magnitude of the output
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = hres->FFT(hm, "MAG");
  hm->SetTitle("discrete fourier transform");
  c1->Clear();
  hm->Draw();

  sprintf( hname, "qHist1D_sig_All_fit%02i.png",Ifit);
  c1->SaveAs(hname);

  return 1;
}

int energyPlots(){

  sprintf( foutname, "Data/sum_coarse_testped.root");
  file1 = TFile::Open( foutname, "READ");
  printf("got input root file %s\n", foutname);

  int xtalPerRow = 9;
  int rb = 14; // rebin factor

  c1->Clear();
  c1->Divide(3,3);

  for (int iCol = 0; iCol<9; iCol++){ 

    int iSeg = iCol;
    sprintf( hname, "qHist1D_thr_CaloSum_%i", iSeg);
    file1->GetObject( hname, h1[iCol]);
    sprintf( hsname,"qHist1D_thr_ColSum_%i", iCol);
    h1[iCol]->SetName(hsname);
    h1[iCol]->SetTitle(hsname);
    printf("got xtal %i, column %i, histogram %s\n", iSeg, iCol, hname);

    for (int iRow = 1; iRow<6; iRow++){ 

      iSeg = iCol + xtalPerRow*iRow;
      sprintf( hname, "qHist1D_thr_CaloSum_%i", iSeg);
      file1->GetObject( hname, htmp);
      h1[iCol]->Add(htmp, 1.0);
      printf("add xtal %i, column %i, histogram %s\n", iSeg, iCol, hname);
    }

  c1->cd(iCol+1);
    
  h1[iCol]->Rebin(rb);
  h1[iCol]->SetMinimum(-2000*rb);
  h1[iCol]->SetMaximum(200*rb);
  h1[iCol]->Draw();
  printf("draw column %i\n", iCol);    
  }

  c1->SaveAs("timedistCaloColumns.png");

  TH1D  *hEnergyColumn_e = new TH1D("hEnergyColumn_e","hEnergyColumn_e",9,0,9);
  hEnergyColumn_e->SetLineColor(kRed);
  hEnergyColumn_e->SetLineWidth(4.0);
  TH1D  *hEnergyColumn_p = new TH1D("hEnergyColumn_p","hEnergyColumn_p",9,0,9);
  hEnergyColumn_p->SetLineColor(kBlue);
  hEnergyColumn_p->SetLineWidth(4.0);

  for (int iCol = 0; iCol<9; iCol++){

    hEnergyColumn_e->SetBinContent( 9-iCol, h1[iCol]->Integral(224,851)); // e+ time range 50-190 us
    hEnergyColumn_p->SetBinContent( 9-iCol, h1[iCol]->Integral(921,1500)); // p time range 50-190 us
  }

  c1->Clear();

  hEnergyColumn_p->Draw();
  hEnergyColumn_e->Draw("same");
  c1->SaveAs("energydistCaloColumns.png");

  c1->Clear();
  c1->Divide(2,3);

  for (int iRow = 0; iRow<6; iRow++){ 

    int iSeg = 9*iRow;
    sprintf( hname, "qHist1D_thr_CaloSum_%i", iSeg);
    file1->GetObject( hname, h2[iRow]);
    sprintf( hsname,"qHist1D_thr_RowSum_%i", iRow);
    h2[iRow]->SetName(hsname);
    h2[iRow]->SetTitle(hsname);
    printf("got xtal %i, row %i, histogram %s\n", iSeg, iRow, hname);

    for (int iC = 1; iC<6; iC++){ 

      iSeg = iC + xtalPerRow*iRow;
      sprintf( hname, "qHist1D_thr_CaloSum_%i", iSeg);
      file1->GetObject( hname, htmp);
      h2[iRow]->Add(htmp, 1.0);
      printf("add xtal %i, row %i, histogram %s\n", iSeg, iRow, hname);
    }

  c1->cd(iRow+1);
    
  h2[iRow]->Rebin(rb);
  h2[iRow]->SetMinimum(-2000*rb);
  h2[iRow]->SetMaximum(200*rb);
  h2[iRow]->Draw();
  printf("draw row %i\n", iRow);    
  }

  c1->SaveAs("timedistCaloRows.png");


  TH1D  *hEnergyRow_e = new TH1D("hEnergyRow_e","hEnergyRow_e",6,0,6);
  hEnergyRow_e->SetLineColor(kRed);
  hEnergyRow_e->SetLineWidth(4.0);
  TH1D  *hEnergyRow_p1 = new TH1D("hEnergyRow_p1","hEnergyRow_p1",6,0,6);
  hEnergyRow_p1->SetLineColor(kBlue);
  hEnergyRow_p1->SetLineWidth(4.0);
  TH1D  *hEnergyRow_p2 = new TH1D("hEnergyRow_p2","hEnergyRow_p2",6,0,6);
  hEnergyRow_p2->SetLineColor(kRed);
  hEnergyRow_p2->SetLineWidth(4.0);
  TH1D  *hEnergyRow_p3 = new TH1D("hEnergyRow_p3","hEnergyRow_p3",6,0,6);
  hEnergyRow_p3->SetLineColor(kGreen);
  hEnergyRow_p3->SetLineWidth(4.0);
  TH1D  *hEnergyRow_p4 = new TH1D("hEnergyRow_p4","hEnergyRow_p4",6,0,6);
  hEnergyRow_p4->SetLineColor(kMagenta);
  hEnergyRow_p4->SetLineWidth(4.0);

  for (int iRow = 0; iRow<6; iRow++){

    hEnergyRow_e->SetBinContent( iRow+1, h2[iRow]->Integral(224,851)); // e+ time range 50-190 us
    //hEnergyRow_p->SetBinContent( iRow+1, h2[iRow]->Integral(921,1500)); // p time range 50-190 us
    hEnergyRow_p1->SetBinContent( iRow+1, h2[iRow]->Integral( 921,1050)); // p time range 50-190 us 
    hEnergyRow_p2->SetBinContent( iRow+1, h2[iRow]->Integral(1051,1200)); // p time range 50-190 us
    hEnergyRow_p3->SetBinContent( iRow+1, h2[iRow]->Integral(1201,1350)); // p time range 50-190 us
    hEnergyRow_p4->SetBinContent( iRow+1, h2[iRow]->Integral(1351,1500)); // p time range 50-190 us
  }

  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1);
  hEnergyRow_p1->Draw();
  c1->cd(2);
  hEnergyRow_p2->Draw();
  c1->cd(3);
  hEnergyRow_p3->Draw();
  c1->cd(4);
  hEnergyRow_p4->Draw();
  //hEnergyRow_e->Draw("same");
  c1->SaveAs("energydistCaloRows-timecuts.png");

  return 1;
}

 TH1D *h_alla, *h_allb; 
 TH1D *h_0a, *h_0b, *h_1a, *h_1b, *h_2a, *h_2b,  *h_3a, *h_3b,  *h_4a, *h_4b, *h_5a, *h_5b;
 TH1D *h_0c, *h_1c, *h_2c, *h_3c, *h_4c, *h_5c;
 TH1D *hfft0,  *hfft1,  *hfft2,  *hfft3,  *hfft4,  *hfft5;
 double R[6], dR[6], phi[6], dphi[6];
 double Re[6], dRe[6], phie[6], dphie[6];
 double Rl[6], dRl[6], phil[6], dphil[6];
 double DR[6], DdR[6], Dphi[6], Ddphi[6];
 double xR[6]={1,2,3,4,5,6};
 double dxR[6]={0};
 double omega[6], domega[6];
 double Slope[6], dSlope[6];
 double Offst[5], dOffst[6];
 double Rc[6], dRc[6], Rp[6], dRp[6]; 
 double flow[6];
 double delta[6]={0};
 double frac[6]={0};
 TGraphErrors *gphi, *gRS, *gphiS, *gSlope;
 TF1 *fn0, *fn1, *fn2, *fn3, *fn4, *fn5; 
 TGraphErrors *gRbefore, *gRafter, *gRdiff; 
 TH1D *h_delta0, *h_delta1, *h_delta2, *h_delta3, *h_delta4, *h_delta5;
 TH1D *h_flow0, *h_flow1, *h_flow2, *h_flow3, *h_flow4, *h_flow5;
 double Norig[6];
 TGraphErrors *gNorig;
 double Nnew[6];
 TGraphErrors *gNnew;
 double dlta[6];
 TGraphErrors *gdlta;
 double fl[6];
 TGraphErrors *gfl;
 double Phic[6], dPhic[6]; 
 double Phip[6], dPhip[6];
 double diffR[6];
 double diffPhi[6];
 TGraphErrors  *gPhibefore,  *gPhiafter,  *gPhidiff;
 TH1D *h_fracflow0, *h_fracflow1, *h_fracflow2, *h_fracflow3, *h_fracflow4, *h_fracflow5;
 TH1D *h_fracdelta0, *h_fracdelta1, *h_fracdelta2, *h_fracdelta3, *h_fracdelta4, *h_fracdelta5;
 TH1D *h_0d, *h_1d, *h_2d, *h_3d, *h_4d, *h_5d; 
 TH1D *h_0e, *h_1e, *h_2e, *h_3e, *h_4e, *h_5e; 
 

void driftcorrection();
void driftcorrection(){

 // free muon lifetime
 fix_gammatau=0;
 set_gammatau=64500;
 
 // no fit-method fr correction
 fastRotationAnalysis = 0;
 fastRotationCorrection = 0;
 scalefastrotation=1.0;
 timeShiftMultiplier=1.0;

 // switch on muon loss correction
 fix_muonloss_amp = 0;
 set_muonloss_amp = 0.0;
 // switch on pedestal ringing correction
 fix_pedestaldrift_amp = 0;
 set_pedestaldrift_amp = 0.0;

 // set cbo parameters, vary vw ,vo
 set_omega_vbo_60hr=9.75196e-01; // from "full" (old 9.75112e-01;)
 set_omega_vbo2_60hr=1.03036e+00; // from "full" (old 1.03043e+00;)
 set_omega_vbo_9day=9.67629e-01; // from "full" (old 9.67436e-01;)
 set_omega_vbo2_9day=1.00726e+00; // from "full" (old 1.00717e+00;)
 set_omega_vbo_hk=9.67629e-01; // from "full" 
 set_omega_vbo2_hk=1.00726e+00; // from "full"
 set_omega_vbo_endgame = 9.74178e-01; // from "full" (old 9.75382e-01;)
 set_omega_vbo2_endgame = 1.00776e+00; // from "full" (old 1.00662e+00;) 
 
 // these are poorly deterimed from "full" fit
 set_tau_vbo_60hr=2.04e+04;  //from "full" (old 8.79011e+04;)
 set_tau_vbo2_60hr=1.74268e+04; //from "full" (old 1.87389e+04;)
 set_tau_vbo_9day=4.28368e+04; // from "full" (old 4.44523e+04;)
 set_tau_vbo2_9day=8.21351e+04; // from "full" (old 7.90928e+04;)
 set_tau_vbo_hk=3.88298e+04; // from "full" 
 set_tau_vbo2_hk=7.36571e+04; // from "full"
 set_tau_vbo_endgame=3.75748e+04; // from "full" (old 2.08418e+04;)
 set_tau_vbo2_endgame=2.e4; //2.84686e+06; // from "full" (old 5.84952e+04;)
 
 // these are estimates from horizontal slice fits
 set_tau_vbo_60hr=5.e+04;  //WHY LONGER? F_VW ~ F_VB ? - from "full" (old 8.79011e+04;)
 set_tau_vbo2_60hr=2.e+04; //WHY SHORTER F_VW ~ F_VB ? - from "full" (old 1.87389e+04;)
 set_tau_vbo_endgame=3.e+04; // from "full" (old 2.08418e+04;)
 set_tau_vbo2_endgame=1.0e5; //2.84686e+06; // from "full" (old 5.84952e+04;)
 set_tau_vbo_9day=5.e4; // from "full" (old 4.44523e+04;)
 set_tau_vbo2_9day=1.1e5; // from "full" (old 7.90928e+04;)
 set_tau_vbo_hk=2.e+04; // from "full" 
 set_tau_vbo2_hk=1.2e+05; // from "full"
 
 // set cbo freq change parameters
 // pre 1/9/2020 values
 set_AExpTermCBO_60hr = 2.90; // from tracker
 set_BExpTermCBO_60hr = 5.12; // from tracker
 set_AExpTermCBO_9day = 2.86;
 set_BExpTermCBO_9day = 5.50;
 set_AExpTermCBO_hk = 2.86;
 set_BExpTermCBO_hk = 5.50;
 set_AExpTermCBO_lk = 2.86;
 set_BExpTermCBO_lk = 5.50;
 set_AExpTermCBO_eg = 6.626; // best fit calo data
 set_BExpTermCBO_eg = 6.821; // best fit calo data
// post 1/9/2020 values
 set_AExpTermCBO_60hr = 2.79; // from tracker
 set_BExpTermCBO_60hr = 5.63; // from tracker
 set_AExpTermCBO_9day = 2.80;
 set_BExpTermCBO_9day = 6.18;
 set_AExpTermCBO_hk = 2.86; // TO DO
 set_BExpTermCBO_hk = 5.50; // TO DO
 set_AExpTermCBO_lk = 2.86; // TO DO
 set_BExpTermCBO_lk = 5.50; // TO DO
 set_AExpTermCBO_eg = 6.82; // from tracker
 set_BExpTermCBO_eg = 5.42; // from tracker
 
 // pre 1/9/2020 values
 set_AExpTauTermCBO_60hr = 81.8e3; // from tracker
 set_BExpTauTermCBO_60hr = 7.7e3; // from tracker
 set_AExpTauTermCBO_9day = 72.8e3;
 set_BExpTauTermCBO_9day = 8.5e3;
 set_AExpTauTermCBO_hk = 72.8e3; // TO DO
 set_BExpTauTermCBO_hk = 8.5e3; // TO DO
 set_AExpTauTermCBO_lk = 72.8e3; // TO DO
 set_BExpTauTermCBO_lk = 8.5e3; // TO DO
 set_AExpTauTermCBO_eg = 81.8e3; // from tracker
 set_BExpTauTermCBO_eg = 7.7e3; // from tracker
 // post 1/9/2020 values
 set_AExpTauTermCBO_60hr = 61.1e3; // from tracker
 set_BExpTauTermCBO_60hr = 6.07e3; // from tracker
 set_AExpTauTermCBO_9day = 56.6e3;
 set_BExpTauTermCBO_9day = 6.32e3;
 set_AExpTauTermCBO_hk = 56.6e3;  // TO DO
 set_BExpTauTermCBO_hk = 6.32e3;  // TO DO
 set_AExpTauTermCBO_lk = 56.6e3;  // TO DO
 set_BExpTauTermCBO_lk = 6.32e3;  // TO DO
 set_AExpTauTermCBO_eg = 78.3e3; // from tracker
 set_BExpTauTermCBO_eg = 6.54e3; // from tracker

 // vary horizontal waist, 2-3 sigma effect (60hr)
 fix_2cbo_amp = 0;
 fix_2cbo_phi = 0;
 // fix horizontal cbo term in precession phase ~1-2 sigma effect (60hr)
 fix_anomphi_amp=0;
 fix_anomphi_phi=0;

 // absolute units for residuals
 fracResiUnits=1;

 // free relax, no pu, no bkd
 fix_asymtau=0;
 set_asymtau=1.e8;
 fix_pileup=1;
 fix_bkd=1;

 // vertical drift
 use_verticaldrift_correction = 1;
 set_gain_amp = 0.e-3;
 fix_gain_amp = 0;
 set_gain_tau = 1.e12;
 fix_gain_tau = 1;

 // nominal caloscan
 fix_AExpTermCBO=1;
 fix_BExpTermCBO=1;
 fix_vbo_tau = 1;
 fix_vbo_tau2 = 1;
 fix_vbo_omega = 1;
 fix_vbo_omega2 = 1;
 fix_2cbo_amp = 0;
 fix_2cbo_phi = 0;
 fix_anomphi_amp=0;
 fix_anomphi_phi=0;
 fix_asymtau=1;
 set_asymtau=7.0e7;

 // do minos, not migrad, errors
 doBlendBins=0;
 doUserCovariance=0;
 doMinos=0;

 // configure fit and run

 sprintf(datasubset,"driftcorrection");

 sprintf(fname,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamepedestal,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamemuonloss,"muonloss/60hours-DQC.root");
 sprintf(dname,"QFillByFillAnalyzer");
 sprintf(dataname,"endgame");
 sprintf(fitname, "vbo2");
 ithreshold=15;
 iwind=4;

 /*
 sprintf(fname,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamepedestal,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-9day-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamemuonloss,"muonloss/9days-DQC.root");
 sprintf(dname,"QFillByFillAnalyzer");
 sprintf(dataname,"9day");
 ithreshold=15;
 iwind=4;

 sprintf(fname,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamepedestal,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-hk-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamemuonloss,"muonloss/9days-DQC.root");
 sprintf(dname,"QFillByFillAnalyzer");
 sprintf(dataname,"hk");
 ithreshold=15;
 iwind=4;

 sprintf(fname,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-60hr-pufix-gold-withpederr-gainbugfix-newenergyscale.root");
 sprintf(fnamepedestal,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-60hr-pufix-gold-withpederr-gainbugfix-newenergyscale.root");
 sprintf(fnamemuonloss,"muonloss/60hours-DQC.root");
 sprintf(dname,"QFillByFillAnalyzer");
 sprintf(dataname,"60hr");
 ithreshold=15;
 iwind=4;
 */

 sprintf(fitname,"vbo2");
 fitCalo("full",0,30000,215500);
 h_alla=(TH1D*)htmp->Clone();
 h_allb=(TH1D*)htmp->Clone();
 h_allc=(TH1D*)htmp->Clone();
 Rall=precessf->GetParameter(3);
 phiall=precessf->GetParameter(4);
 dRall=precessf->GetParError(3);
 dphiall=precessf->GetParError(4);
 fitCalo("horizontalslice",0,30000,215500);
 h_0a=(TH1D*)htmp->Clone();
 h_0b=(TH1D*)htmp->Clone();
 h_0c=(TH1D*)htmp->Clone();
 R[0]=precessf->GetParameter(3);
 phi[0]=precessf->GetParameter(4);
 dR[0]=precessf->GetParError(3);
 dphi[0]=precessf->GetParError(4);
 fitCalo("horizontalslice",1,30000,215500);
 h_1a=(TH1D*)htmp->Clone();
 h_1b=(TH1D*)htmp->Clone();
 h_1c=(TH1D*)htmp->Clone();
 R[1]=precessf->GetParameter(3);
 phi[1]=precessf->GetParameter(4);
 dR[1]=precessf->GetParError(3);
 dphi[1]=precessf->GetParError(4);
 fitCalo("horizontalslice",2,30000,215500);
 h_2a=(TH1D*)htmp->Clone();
 h_2b=(TH1D*)htmp->Clone();
 h_2c=(TH1D*)htmp->Clone();
 R[2]=precessf->GetParameter(3);
 phi[2]=precessf->GetParameter(4);
 dR[2]=precessf->GetParError(3);
 dphi[2]=precessf->GetParError(4);
 fitCalo("horizontalslice",3,30000,215500);
 h_3a=(TH1D*)htmp->Clone();
 h_3b=(TH1D*)htmp->Clone();
 h_3c=(TH1D*)htmp->Clone();
 R[3]=precessf->GetParameter(3);
 phi[3]=precessf->GetParameter(4);
 dR[3]=precessf->GetParError(3);
 dphi[3]=precessf->GetParError(4);
 fitCalo("horizontalslice",4,30000,215500);
 h_4a=(TH1D*)htmp->Clone();
 h_4b=(TH1D*)htmp->Clone();
 h_4c=(TH1D*)htmp->Clone();
 R[4]=precessf->GetParameter(3);
 phi[4]=precessf->GetParameter(4);
 dR[4]=precessf->GetParError(3);
 dphi[4]=precessf->GetParError(4);
 fitCalo("horizontalslice",5,30000,215500);
 h_5a=(TH1D*)htmp->Clone();
 h_5b=(TH1D*)htmp->Clone();
 h_5c=(TH1D*)htmp->Clone();
 R[5]=precessf->GetParameter(3);
 phi[5]=precessf->GetParameter(4);
 dR[5]=precessf->GetParError(3);
 dphi[5]=precessf->GetParError(4);
 for (int ip = 0; ip < 6; ip++) printf("%i %f %f \n", ip, phi[ip], dphi[ip]);
 
 h_0c->Reset();
 h_1c->Reset();
 h_2c->Reset();
 h_3c->Reset();
 h_4c->Reset();
 h_5c->Reset();
 
 h_0a->SetLineColor(kBlack);
 h_0b->SetLineColor(kBlack);
 h_0c->SetLineColor(kBlack);
 h_1a->SetLineColor(kRed);
 h_1b->SetLineColor(kRed);
 h_1c->SetLineColor(kRed);
 h_2a->SetLineColor(kGreen);
 h_2b->SetLineColor(kGreen);
 h_2c->SetLineColor(kGreen);
 h_3a->SetLineColor(kBlue);
 h_3b->SetLineColor(kBlue);
 h_3c->SetLineColor(kBlue);
 h_4a->SetLineColor(kMagenta);
 h_4b->SetLineColor(kMagenta);
 h_4c->SetLineColor(kMagenta);
 h_5a->SetLineColor(kOrange);
 h_5b->SetLineColor(kOrange);
 h_5c->SetLineColor(kOrange);
 
 h_0a->SetLineWidth(2);
 h_1a->SetLineWidth(2);
 h_2a->SetLineWidth(2);
 h_3a->SetLineWidth(2);
 h_4a->SetLineWidth(2);
 h_5a->SetLineWidth(2);
 h_0b->SetLineWidth(2);
 h_1b->SetLineWidth(2);
 h_2b->SetLineWidth(2);
 h_3b->SetLineWidth(2);
 h_4b->SetLineWidth(2);
 h_5b->SetLineWidth(2);
 h_0c->SetLineWidth(1);
 h_1c->SetLineWidth(1);
 h_2c->SetLineWidth(1);
 h_3c->SetLineWidth(1);
 h_4c->SetLineWidth(1);
 h_5c->SetLineWidth(1);
 
 h_0c->Add(h_0a,1.0);
 h_1c->Add(h_1a,1.0);
 h_2c->Add(h_2a,1.0);
 h_3c->Add(h_3a,1.0);
 h_4c->Add(h_4a,1.0);
 h_5c->Add(h_5a,1.0);
 
 h_0b->Scale(h_allb->Integral()/h_0b->Integral());
 h_1b->Scale(h_allb->Integral()/h_1b->Integral());
 h_2b->Scale(h_allb->Integral()/h_2b->Integral());
 h_3b->Scale(h_allb->Integral()/h_3b->Integral());
 h_4b->Scale(h_allb->Integral()/h_4b->Integral());
 h_5b->Scale(h_allb->Integral()/h_5b->Integral());
 
 h_0b->Divide(h_allb);
 h_1b->Divide(h_allb);
 h_2b->Divide(h_allb);
 h_3b->Divide(h_allb);
 h_4b->Divide(h_allb);
 h_5b->Divide(h_allb);
 
 /*
 hfft0 = (TH1D*) h_0b->FFT(hfft0, "MAG");
 hfft1 = (TH1D*) h_1b->FFT(hfft1, "MAG");
 hfft2 = (TH1D*) h_2b->FFT(hfft2, "MAG");
 hfft3 = (TH1D*) h_3b->FFT(hfft3, "MAG");
 hfft4 = (TH1D*) h_4b->FFT(hfft4, "MAG");
 hfft5 = (TH1D*) h_5b->FFT(hfft5, "MAG");
 
 hfft0->SetLineColor(kBlack);
 hfft1->SetLineColor(kRed);
 hfft2->SetLineColor(kGreen);
 hfft3->SetLineColor(kBlue);
 hfft4->SetLineColor(kMagenta);
 hfft5->SetLineColor(kOrange);
 
 hfft0->SetLineWidth(2);
 hfft1->SetLineWidth(2);
 hfft2->SetLineWidth(2);
 hfft3->SetLineWidth(2);
 hfft4->SetLineWidth(2);
 hfft5->SetLineWidth(2);
 
 hfft0->Scale(hfft3->Integral(10,600)/hfft0->Integral(10,600));
 hfft1->Scale(hfft3->Integral(10,600)/hfft1->Integral(10,600));
 hfft2->Scale(hfft3->Integral(10,600)/hfft2->Integral(10,600));
 hfft4->Scale(hfft3->Integral(10,600)/hfft4->Integral(10,600));
 hfft5->Scale(hfft3->Integral(10,600)/hfft5->Integral(10,600));
 
 c1->Clear();
 c1->Divide(1,1);
 c1->cd(1);
 hfft0->GetXaxis()->SetRangeUser(10,100);
 hfft0->Draw("hist");
 hfft1->Draw("histsame");
 hfft2->Draw("histsame");
 hfft3->Draw("histsame");
 hfft4->Draw("histsame");
 hfft5->Draw("histsame");
 */

 // get slopes due to drift
 fn0 = new TF1("fn0","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn0->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn0->SetNpx(10000);
 fn0->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn0->FixParameter(0,0.03);
 fn0->FixParameter(3,1.0);
 h_0b->Fit( fn0, "", "NRE", 30000, 215500);
 //fn0->FixParameter(1,fn2->GetParameter(1));
 //fn0->FixParameter(2,fn2->GetParameter(2));
 fn0->ReleaseParameter(0);
 fn0->ReleaseParameter(1);
 fn0->ReleaseParameter(2);
 fn0->ReleaseParameter(3);
 fn0->ReleaseParameter(4);
 h_0b->Fit( fn0, "", "NRE", 30000, 215500);
 omega[0] = fn0->GetParameter(1);
 domega[0] = fn0->GetParError(1);
 Slope[0] = fn0->GetParameter(4);
 dSlope[0] = fn0->GetParError(4);
 fn1 = new TF1("fn1","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn1->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn1->SetNpx(10000);
 fn1->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn1->FixParameter(0,0.03);
 fn1->FixParameter(3,1.0);
 h_1b->Fit( fn1, "", "NRE", 30000, 215500);
 //fn1->FixParameter(1,fn1->GetParameter(1));
 //fn1->FixParameter(2,fn1->GetParameter(2));
 fn1->ReleaseParameter(0);
 fn1->ReleaseParameter(1);
 fn1->ReleaseParameter(2);
 fn1->ReleaseParameter(3);
 fn1->ReleaseParameter(4);
 h_1b->Fit( fn1, "", "NRE", 30000, 215500);
 omega[1] = fn1->GetParameter(1);
 domega[1] = fn1->GetParError(1);
 Slope[1] = fn1->GetParameter(4);
 dSlope[1] = fn1->GetParError(4);
 fn2 = new TF1("fn2","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn2->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn2->SetNpx(10000);
 fn2->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn2->FixParameter(0,0.03);
 fn2->FixParameter(3,1.0);
 h_2b->Fit( fn2, "", "NRE", 30000, 215500);
 //fn2->FixParameter(1,fn2->GetParameter(1));
 //fn2->FixParameter(2,fn2->GetParameter(2));
 fn2->ReleaseParameter(0);
 fn2->ReleaseParameter(1);
 fn2->ReleaseParameter(2);
 fn2->ReleaseParameter(3);
 fn2->ReleaseParameter(4); 
h_2b->Fit( fn2, "", "NRE", 30000, 215500);
 omega[2] = fn2->GetParameter(1);
 domega[2] = fn2->GetParError(1);
 Slope[2] = fn2->GetParameter(4);
 dSlope[2] = fn2->GetParError(4);
 fn3 = new TF1("fn3","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn3->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn3->SetNpx(10000);
 fn3->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn3->FixParameter(0,0.03);
 fn3->FixParameter(3,1.0);
 h_3b->Fit( fn3, "", "NRE", 30000, 215500);
 //fn3->FixParameter(1,fn3->GetParameter(1));
 //fn3->FixParameter(2,fn3->GetParameter(2));
 fn3->ReleaseParameter(0);
 fn3->ReleaseParameter(1);
 fn3->ReleaseParameter(2);
 fn3->ReleaseParameter(3);
 fn3->ReleaseParameter(4);
 h_3b->Fit( fn3, "", "NRE", 30000, 215500);
 omega[3] = fn3->GetParameter(1);
 domega[3] = fn3->GetParError(1);
 Slope[3] = fn3->GetParameter(4);
 dSlope[3] = fn3->GetParError(4);
 fn4 = new TF1("fn4","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn4->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn4->SetNpx(10000);
 fn4->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn4->FixParameter(0,0.03);
 fn4->FixParameter(3,1.0);
 h_4b->Fit( fn4, "", "NRE", 30000, 215500);
 //fn4->FixParameter(1,fn4->GetParameter(1));
 //fn4->FixParameter(2,fn4->GetParameter(2));
 fn4->ReleaseParameter(0);
 fn4->ReleaseParameter(1);
 fn4->ReleaseParameter(2);
 fn4->ReleaseParameter(3);
 fn4->ReleaseParameter(4);
 h_4b->Fit( fn4, "", "NRE", 30000, 215500);
 omega[4] = fn4->GetParameter(1);
 domega[4] = fn4->GetParError(1);
 Slope[4] = fn4->GetParameter(4);
 dSlope[4] = fn4->GetParError(4);
 fn5 = new TF1("fn5","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
 fn5->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
 fn5->SetNpx(10000);
 fn5->SetParameters(0.07,0.00143,0.0,1.0,0.0);
 fn5->FixParameter(0,0.03);
 fn5->FixParameter(3,1.0);
 h_5b->Fit( fn5, "", "NRE", 30000, 215500);
 //fn5->FixParameter(1,fn5->GetParameter(1));
 //fn5->FixParameter(2,fn5->GetParameter(2));
 fn5->ReleaseParameter(0);
 fn5->ReleaseParameter(1);
 fn5->ReleaseParameter(2);
 fn5->ReleaseParameter(3);
 fn5->ReleaseParameter(4);
 h_5b->Fit( fn5, "", "NRE", 30000, 215500);
 omega[5] = fn5->GetParameter(1);
 domega[5] = fn5->GetParError(1);
 Slope[5] = fn5->GetParameter(4);
 dSlope[5] = fn5->GetParError(4);
 for (int ip = 0; ip < 6; ip++)  printf(" %i %f %f %14.8e %14.8e  %e %e \n", ip, R[ip], dR[ip], omega[ip], domega[ip], Slope[ip], dSlope[ip]);
 
 c1->Clear();
 c1->Divide(1,3);

 c1->cd(1);
 gSlope = new TGraphErrors(6,xR,Slope,dxR,dSlope);
 gSlope->SetTitle("slope parameter v's y-coordinate");
 gSlope->SetMarkerStyle(8);
 gSlope->SetMarkerSize(2);
 gSlope->SetMarkerColor(kBlue);
 gSlope->Draw("AP");
 
 c1->cd(2);
 gphi = new TGraphErrors(6,xR,phi,dxR,dphi);
 gphi->SetTitle("phase v's y-coordinate");
 gphi->SetMarkerStyle(8);
 gphi->SetMarkerSize(2);
 gphi->SetMarkerColor(kRed);
 gphi->Draw("AP");
 
 c1->cd(3);
 gdR = new TGraphErrors(6,xR,R,dxR,dR);
 gdR->SetTitle("R versus y-coordinate");
 gdR->SetMarkerStyle(8);
 gdR->SetMarkerSize(2);
 gdR->SetMarkerColor(kGreen);
 gdR->Draw("AP");
 
 // delta's are histograms of changes in slices v's time
 h_delta0=(TH1D*)h_0c->Clone();
 h_delta1=(TH1D*)h_1c->Clone();
 h_delta2=(TH1D*)h_2c->Clone();
 h_delta3=(TH1D*)h_3c->Clone();
 h_delta4=(TH1D*)h_4c->Clone();
 h_delta5=(TH1D*)h_5c->Clone();
 h_delta0->Reset();
 h_delta1->Reset();
 h_delta2->Reset();
 h_delta3->Reset();
 h_delta4->Reset();
 h_delta5->Reset();
 h_delta0->SetLineColor(kBlack);
 h_delta1->SetLineColor(kRed);
 h_delta2->SetLineColor(kGreen);
 h_delta3->SetLineColor(kBlue);
 h_delta4->SetLineColor(kMagenta);
 h_delta5->SetLineColor(kOrange);
 h_delta0->SetLineWidth(2);
 h_delta1->SetLineWidth(2);
 h_delta2->SetLineWidth(2);
 h_delta3->SetLineWidth(2);
 h_delta4->SetLineWidth(2);
 h_delta5->SetLineWidth(2);

 h_delta0->SetTitle("delta_n versus time, slice 0");
 h_delta1->SetTitle("delta_n versus time, slice 1");
 h_delta2->SetTitle("delta_n versus time, slice 2");
 h_delta3->SetTitle("delta_n versus time, slice 3");
 h_delta4->SetTitle("delta_n versus time, slice 4");
 h_delta5->SetTitle("delta_n versus time, slice 5");

 h_delta0->GetXaxis()->SetTitle("time (ns)");
 h_delta1->GetXaxis()->SetTitle("time (ns)");
 h_delta2->GetXaxis()->SetTitle("time (ns)");
 h_delta3->GetXaxis()->SetTitle("time (ns)");
 h_delta4->GetXaxis()->SetTitle("time (ns)");
 h_delta5->GetXaxis()->SetTitle("time (ns)");

 
 // flow's are histograms of flow into slices v's time
 h_flow0=(TH1D*)h_0c->Clone();
 h_flow1=(TH1D*)h_1c->Clone();
 h_flow2=(TH1D*)h_2c->Clone();
 h_flow3=(TH1D*)h_3c->Clone();
 h_flow4=(TH1D*)h_4c->Clone();
 h_flow5=(TH1D*)h_5c->Clone();
 h_flow0->Reset();
 h_flow1->Reset();
 h_flow2->Reset();
 h_flow3->Reset();
 h_flow4->Reset();
 h_flow5->Reset();
 h_flow0->SetLineColor(kBlack);
 h_flow1->SetLineColor(kRed);
 h_flow2->SetLineColor(kGreen);
 h_flow3->SetLineColor(kBlue);
 h_flow4->SetLineColor(kMagenta);
 h_flow5->SetLineColor(kOrange);
 h_flow0->SetLineWidth(2);
 h_flow1->SetLineWidth(2);
 h_flow2->SetLineWidth(2);
 h_flow3->SetLineWidth(2);
 h_flow4->SetLineWidth(2);
 h_flow5->SetLineWidth(2);

 h_flow0->SetTitle("flow_n versus time, -> slice0");
 h_flow1->SetTitle("flow_n versus time, -> slice1");
 h_flow2->SetTitle("flow_n versus time, -> slice2");
 h_flow3->SetTitle("flow_n versus time, -> slice3");
 h_flow4->SetTitle("flow_n versus time, -> slice4");
 h_flow5->SetTitle("flow_n versus time, -> slice5");

 h_flow0->GetXaxis()->SetTitle("time (ns)");
 h_flow1->GetXaxis()->SetTitle("time (ns)");
 h_flow2->GetXaxis()->SetTitle("time (ns)");
 h_flow3->GetXaxis()->SetTitle("time (ns)");
 h_flow4->GetXaxis()->SetTitle("time (ns)");
 h_flow5->GetXaxis()->SetTitle("time (ns)");
 
 // get change in slices from drift
 for (int ib = 1; ib <= h_0c->GetNbinsX(); ib++) h_delta0->SetBinContent( ib, Slope[0]*h_delta0->GetBinCenter(ib)*h_0a->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_1c->GetNbinsX(); ib++) h_delta1->SetBinContent( ib, Slope[1]*h_delta1->GetBinCenter(ib)*h_1a->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_2c->GetNbinsX(); ib++) h_delta2->SetBinContent( ib, Slope[2]*h_delta2->GetBinCenter(ib)*h_2a->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_3c->GetNbinsX(); ib++) h_delta3->SetBinContent( ib, Slope[3]*h_delta3->GetBinCenter(ib)*h_3a->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_4c->GetNbinsX(); ib++) h_delta4->SetBinContent( ib, Slope[4]*h_delta4->GetBinCenter(ib)*h_4a->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_5c->GetNbinsX(); ib++) h_delta5->SetBinContent( ib, Slope[5]*h_delta5->GetBinCenter(ib)*h_5a->GetBinContent(ib)); 
 printf("%f\n",h_delta5->GetBinContent(500));
 printf("%f\n",h_delta4->GetBinContent(500));
 printf("%f\n",h_delta3->GetBinContent(500));
 printf("%f\n",h_delta2->GetBinContent(500));
 printf("%f\n",h_delta1->GetBinContent(500));
 printf("%f\n",h_delta0->GetBinContent(500));
 printf("%f\n",h_5a->GetBinContent(500));
 printf("%f\n",h_4a->GetBinContent(500));
 printf("%f\n",h_3a->GetBinContent(500));
 printf("%f\n",h_2a->GetBinContent(500));
 printf("%f\n",h_1a->GetBinContent(500));
 printf("%f\n",h_0a->GetBinContent(500));
 
 // get flow i+1 -> i in slices from drift
 for (int ib = 1; ib <= h_5c->GetNbinsX(); ib++) h_flow5->SetBinContent( ib, 0.0); 
 for (int ib = 1; ib <= h_4c->GetNbinsX(); ib++) h_flow4->SetBinContent( ib, -h_delta5->GetBinContent(ib)+h_flow5->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_3c->GetNbinsX(); ib++) h_flow3->SetBinContent( ib, -h_delta4->GetBinContent(ib)+h_flow4->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_2c->GetNbinsX(); ib++) h_flow2->SetBinContent( ib, -h_delta3->GetBinContent(ib)+h_flow3->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_1c->GetNbinsX(); ib++) h_flow1->SetBinContent( ib, -h_delta2->GetBinContent(ib)+h_flow2->GetBinContent(ib)); 
 for (int ib = 1; ib <= h_0c->GetNbinsX(); ib++) h_flow0->SetBinContent( ib, -h_delta1->GetBinContent(ib)+h_flow1->GetBinContent(ib)); 
 printf("%f\n",h_flow5->GetBinContent(500));
 printf("%f\n",h_flow4->GetBinContent(500));
 printf("%f\n",h_flow3->GetBinContent(500));
 printf("%f\n",h_flow2->GetBinContent(500));
 printf("%f\n",h_flow1->GetBinContent(500));
 printf("%f\n",h_flow0->GetBinContent(500));
 
 h_fracdelta0=(TH1D*)h_delta0->Clone();
 h_fracdelta1=(TH1D*)h_delta1->Clone();
 h_fracdelta2=(TH1D*)h_delta2->Clone();
 h_fracdelta3=(TH1D*)h_delta3->Clone();
 h_fracdelta4=(TH1D*)h_delta4->Clone();
 h_fracdelta5=(TH1D*)h_delta5->Clone();
 h_fracdelta0->Divide(h_0a);
 h_fracdelta1->Divide(h_1a);
 h_fracdelta2->Divide(h_2a);
 h_fracdelta3->Divide(h_3a);
 h_fracdelta4->Divide(h_4a);
 h_fracdelta5->Divide(h_5a);
 h_fracdelta0->SetTitle("fractional delta_n versus time, slice 0");
 h_fracdelta1->SetTitle("fractional delta_n versus time, slice 1");
 h_fracdelta2->SetTitle("fractional delta_n versus time, slice 2");
 h_fracdelta3->SetTitle("fractional delta_n versus time, slice 3");
 h_fracdelta4->SetTitle("fractional delta_n versus time, slice 4");
 h_fracdelta5->SetTitle("fractional delta_n versus time, slice 5");
 h_fracdelta0->SetMaximum(0.03);
 h_fracdelta0->SetMinimum(-0.03);
 h_fracdelta1->SetMaximum(0.03);
 h_fracdelta1->SetMinimum(-0.03);
 h_fracdelta2->SetMaximum(0.03);
 h_fracdelta2->SetMinimum(-0.03);
 h_fracdelta3->SetMaximum(0.03);
 h_fracdelta3->SetMinimum(-0.03);
 h_fracdelta4->SetMaximum(0.03);
 h_fracdelta4->SetMinimum(-0.03);
 h_fracdelta5->SetMaximum(0.03);
 h_fracdelta5->SetMinimum(-0.03);

 h_fracflow0=(TH1D*)h_flow0->Clone();
 h_fracflow1=(TH1D*)h_flow1->Clone();
 h_fracflow2=(TH1D*)h_flow2->Clone();
 h_fracflow3=(TH1D*)h_flow3->Clone();
 h_fracflow4=(TH1D*)h_flow4->Clone();
 h_fracflow5=(TH1D*)h_flow5->Clone();
 h_fracflow0->Divide(h_0a);
 h_fracflow1->Divide(h_1a);
 h_fracflow2->Divide(h_2a);
 h_fracflow3->Divide(h_3a);
 h_fracflow4->Divide(h_4a);
 h_fracflow5->Divide(h_5a);
 h_fracflow0->SetTitle("fractional flow_n versus time, slice 0");
 h_fracflow1->SetTitle("fractional flow_n versus time, slice 1");
 h_fracflow2->SetTitle("fractional flow_n versus time, slice 2");
 h_fracflow3->SetTitle("fractional flow_n versus time, slice 3");
 h_fracflow4->SetTitle("fractional flow_n versus time, slice 4");
 h_fracflow5->SetTitle("fractional flow_n versus time, slice 5");
 h_fracflow0->SetMaximum(0.03);
 h_fracflow0->SetMinimum(-0.03);
 h_fracflow1->SetMaximum(0.03);
 h_fracflow1->SetMinimum(-0.03);
 h_fracflow2->SetMaximum(0.03);
 h_fracflow2->SetMinimum(-0.03);
 h_fracflow3->SetMaximum(0.03);
 h_fracflow3->SetMinimum(-0.03);
 h_fracflow4->SetMaximum(0.03);
 h_fracflow4->SetMinimum(-0.03);
 h_fracflow5->SetMaximum(0.03);
 h_fracflow5->SetMinimum(-0.03);

 h_delta0->SetMaximum(1.e7);
 h_delta0->SetMinimum(-1.e7);
 h_delta1->SetMaximum(1.e7);
 h_delta1->SetMinimum(-1.e7);
 h_delta2->SetMaximum(1.e7);
 h_delta2->SetMinimum(-1.e7);
 h_delta3->SetMaximum(1.e7);
 h_delta3->SetMinimum(-1.e7);
 h_delta4->SetMaximum(1.e7);
 h_delta4->SetMinimum(-1.e7);
 h_delta5->SetMaximum(1.e7);
 h_delta5->SetMinimum(-1.e7);

 h_flow0->SetMaximum(2.e7);
 h_flow0->SetMinimum(-0.e7);
 h_flow1->SetMaximum(2.e7);
 h_flow1->SetMinimum(-0.e7);
 h_flow2->SetMaximum(2.e7);
 h_flow2->SetMinimum(-0.e7);
 h_flow3->SetMaximum(2.e7);
 h_flow3->SetMinimum(-0.e7);
 h_flow4->SetMaximum(2.e7);
 h_flow4->SetMinimum(-0.e7);
 h_flow5->SetMaximum(2.e7);
 h_flow5->SetMinimum(-0.e7);

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 c1->cd(1); 
 h_delta0->Draw("hist");
 c1->cd(2);
 h_delta1->Draw("hist");
 c1->cd(3);
 h_delta2->Draw("hist");
 c1->cd(4);
 h_delta3->Draw("hist");
 c1->cd(5);
 h_delta4->Draw("hist");
 c1->cd(6);
 h_delta5->Draw("hist");

 c1->SaveAs("deltaN_v_time_per_slice.png");

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 c1->cd(1); 
 h_fracdelta0->Draw("hist");
 c1->cd(2);
 h_fracdelta1->Draw("hist");
 c1->cd(3);
 h_fracdelta2->Draw("hist");
 c1->cd(4);
 h_fracdelta3->Draw("hist");
 c1->cd(5);
 h_fracdelta4->Draw("hist");
 c1->cd(6);
 h_fracdelta5->Draw("hist");

 c1->SaveAs("frac_deltaN_v_time_per_slice.png");

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 c1->cd(1); 
 h_flow0->Draw("hist");
 c1->cd(2);
 h_flow1->Draw("hist");
 c1->cd(3);
 h_flow2->Draw("hist");
 c1->cd(4);
 h_flow3->Draw("hist");
 c1->cd(5);
 h_flow4->Draw("hist");
 c1->cd(6);
 h_flow5->Draw("hist");

 c1->SaveAs("flowN_v_time_per_slice.png");

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 c1->cd(1); 
 h_fracflow0->Draw("hist");
 c1->cd(2);
 h_fracflow1->Draw("hist");
 c1->cd(3);
 h_fracflow2->Draw("hist");
 c1->cd(4);
 h_fracflow3->Draw("hist");
 c1->cd(5);
 h_fracflow4->Draw("hist");
 c1->cd(6);
 h_fracflow5->Draw("hist");

 c1->SaveAs("frac_flowN_v_time_per_slice.png");

 // top to bottom: orange, magenta, blue, green, red, black
 
 // correct original slices for flow into slice from above

 // cant measure the flow out of bottom slice 0
 // cant measure the flow into top slice 5
 h_0c->Add(h_flow0, -1.0);

 h_1c->Add(h_flow1, -1.0); // remove from from 2 into 1
 h_1c->Add(h_flow0, 1.0); // add flow into 0 from 1

 h_2c->Add(h_flow2, -1.0);
 h_2c->Add(h_flow1, 1.0);

 h_3c->Add(h_flow3, -1.0);
 h_3c->Add(h_flow2, 1.0);

 h_4c->Add(h_flow4, -1.0);
 h_4c->Add(h_flow3, 1.0);

 h_5c->Add(h_flow4, 1.0); 
 //h_5c->Add(h_flow5, -1.0);

 h_0d=(TH1D*)h_0c->Clone();
 h_1d=(TH1D*)h_1c->Clone();
 h_2d=(TH1D*)h_2c->Clone();
 h_3d=(TH1D*)h_3c->Clone();
 h_4d=(TH1D*)h_4c->Clone();
 h_5d=(TH1D*)h_5c->Clone();
 
 h_0d->Add(h_0a,-1.0);
 h_1d->Add(h_1a,-1.0);
 h_2d->Add(h_2a,-1.0);
 h_3d->Add(h_3a,-1.0);
 h_4d->Add(h_4a,-1.0);
 h_5d->Add(h_5a,-1.0);
 
 h_0d->Divide(h_0a);
 h_1d->Divide(h_1a);
 h_2d->Divide(h_2a);
 h_3d->Divide(h_3a);
 h_4d->Divide(h_4a);
 h_5d->Divide(h_5a);

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 h_3a->GetXaxis()->SetRangeUser(30000,220000);
 //c1->cd(1); 
 h_3a->Draw("hist");
 //c1->cd(2); 
 h_1a->Draw("histsame");
 //c1->cd(3); 
 h_2a->Draw("histsame");
 //c1->cd(4); 
 h_0a->Draw("histsame");
 //c1->cd(5); 
 h_4a->Draw("histsame");
 //c1->cd(6); 
 h_5a->Draw("histsame");
 c1->SaveAs("uncorrectedslices-timedist.png"); 

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 h_3c->GetXaxis()->SetRangeUser(30000,220000);
 //c1->cd(1); 
 h_3c->Draw("hist");
 //c1->cd(2); 
 h_1c->Draw("histsame");
 //c1->cd(3); 
 h_2c->Draw("histsame");
 //c1->cd(4); 
 h_0c->Draw("histsame");
 //c1->cd(5); 
 h_4c->Draw("histsame");
 //c1->cd(6); 
 h_5c->Draw("histsame");
 c1->SaveAs("correctedslices-timedist.png"); 

 h_0d->SetLineColor(kBlack);
 h_1d->SetLineColor(kRed);
 h_2d->SetLineColor(kGreen);
 h_3d->SetLineColor(kBlue);
 h_4d->SetLineColor(kMagenta);
 h_5d->SetLineColor(kOrange);
 h_0d->SetLineWidth(2);
 h_1d->SetLineWidth(2);
 h_2d->SetLineWidth(2);
 h_3d->SetLineWidth(2);
 h_4d->SetLineWidth(2);
 h_5d->SetLineWidth(2);

 c1->Clear();
 c1->Divide(2,3);
 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 h_0d->GetXaxis()->SetRangeUser(30000,220000);
 //c1->cd(1); 
 h_0d->Draw("hist");
 //c1->cd(2); 
 h_1d->Draw("histsame");
 //c1->cd(3); 
 h_2d->Draw("histsame");
 //c1->cd(4); 
 h_3d->Draw("histsame");
 //c1->cd(5); 
 h_4d->Draw("histsame");
 //c1->cd(6); 
 h_5d->Draw("histsame");
 c1->SaveAs("fractionalcorrectionslices-timedist.png"); 

 Norig[0]=h_0a->GetBinContent(500);
 Norig[1]=h_1a->GetBinContent(500);
 Norig[2]=h_2a->GetBinContent(500);
 Norig[3]=h_3a->GetBinContent(500);
 Norig[4]=h_4a->GetBinContent(500);
 Norig[5]=h_5a->GetBinContent(500);
 for (int ip = 0; ip < 6; ip++) printf("%i %f \n", ip, Norig[ip]);

 gNorig = new TGraphErrors(6,xR,Norig,dxR,0);
 gNorig->SetTitle("Ne-before versus y-coordinate for bin 500, time 65us");
 gNorig->SetMarkerStyle(8);
 gNorig->SetMarkerSize(2);
 gNorig->SetMarkerColor(kBlack);
 
 Nnew[0]=h_0c->GetBinContent(500);
 Nnew[1]=h_1c->GetBinContent(500);
 Nnew[2]=h_2c->GetBinContent(500);
 Nnew[3]=h_3c->GetBinContent(500);
 Nnew[4]=h_4c->GetBinContent(500);
 Nnew[5]=h_5a->GetBinContent(500);
 for (int ip = 0; ip < 6; ip++) printf("%i %f \n", ip, Nnew[ip]);

 gNnew = new TGraphErrors(6,xR,Nnew,dxR,0);
 gNnew->SetTitle("Ne-after versus y-coordinate for bin 500, time 65us");
 gNnew->SetMarkerStyle(8);
 gNnew->SetMarkerSize(2);
 gNnew->SetMarkerColor(kRed);
 dlta[0]=h_delta0->GetBinContent(500);
 dlta[1]=h_delta1->GetBinContent(500);
 dlta[2]=h_delta2->GetBinContent(500);
 dlta[3]=h_delta3->GetBinContent(500);
 dlta[4]=h_delta4->GetBinContent(500);
 dlta[5]=h_delta5->GetBinContent(500);
 for (int ip = 0; ip < 6; ip++) printf("%i %f \n", ip, dlta[ip]);

 gdlta = new TGraphErrors(6,xR,dlta,dxR,0);
 gdlta->SetTitle("delta versus y-coordinate for bin 500, time 65us");
 gdlta->SetMarkerStyle(8);
 gdlta->SetMarkerSize(2);
 gdlta->SetMarkerColor(kMagenta);
 
 fl[0]=h_flow0->GetBinContent(500);
 fl[1]=h_flow1->GetBinContent(500);
 fl[2]=h_flow2->GetBinContent(500);
 fl[3]=h_flow3->GetBinContent(500);
 fl[4]=h_flow4->GetBinContent(500);
 fl[5]=h_flow5->GetBinContent(500);
 
 gfl = new TGraphErrors(6,xR,fl,dxR,0);
 gfl->SetTitle("flow versus y-coordinate for bin 500, time 65us");
 gfl->SetMarkerStyle(8);
 gfl->SetMarkerSize(2);
 gfl->SetMarkerColor(kGreen);
 
 for (int ip = 5; ip >= 0; ip--) printf("%i %e %e %e %f %f\n", ip, Norig[ip], Nnew[ip], Slope[ip], dlta[ip], fl[ip]);
 
 c1->Clear();
 c1->Divide(1,4);
 c1->cd(1);
 gNorig->Draw("AP");
 gNnew->Draw("Psame");
 c1->cd(2);
 gSlope->Draw("AP");
 c1->cd(3);
 gdlta->Draw("AP");
 c1->cd(4);
 gfl->Draw("AP");
 
 c1->Clear();
 c1->Divide(1,3);
 c1->cd(1);
 h_3a->Draw("hist");
 h_1a->Draw("histsame");
 h_2a->Draw("histsame");
 h_0a->Draw("histsame");
 h_4a->Draw("histsame");
 h_5a->Draw("histsame");
 c1->cd(2);
 h_3b->Draw("hist");
 h_1b->Draw("histsame");
 h_2b->Draw("histsame");
 h_1b->Draw("histsame");
 h_4b->Draw("histsame");
 h_5b->Draw("histsame");
 c1->cd(3);
 h_3c->Draw("hist");
 h_1c->Draw("histsame");
 h_2c->Draw("histsame");
 h_0c->Draw("histsame");
 h_4c->Draw("histsame");
 h_5c->Draw("histsame");
 // top to bottom: orange, magenta, blue, green, red, black
 
 // "c" check fit numbers before correction for drift 
 doDriftCorrection=1;
 htmpdriftcorrection=(TH1D*)h_0a->Clone();
 fitCalo("horizontalslice",0,30000,215500);
 Rc[0]=precessf->GetParameter(3);
 dRc[0]=precessf->GetParError(3);
 Phic[0]=precessf->GetParameter(4);
 dPhic[0]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_0c->Clone();
 fitCalo("horizontalslice",0,30000,215500);
 Rp[0]=precessf->GetParameter(3);
 dRp[0]=precessf->GetParError(3);
 Phip[0]=precessf->GetParameter(4);
 dPhip[0]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_1a->Clone();
 fitCalo("horizontalslice",1,30000,215500);
 Rc[1]=precessf->GetParameter(3);
 dRc[1]=precessf->GetParError(3);
 Phic[1]=precessf->GetParameter(4);
 dPhic[1]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_1c->Clone();
 fitCalo("horizontalslice",1,30000,215500);
 Rp[1]=precessf->GetParameter(3);
 dRp[1]=precessf->GetParError(3);
 Phip[1]=precessf->GetParameter(4);
 dPhip[1]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_2a->Clone();
 fitCalo("horizontalslice",2,30000,215500);
 Rc[2]=precessf->GetParameter(3);
 dRc[2]=precessf->GetParError(3);
 Phic[2]=precessf->GetParameter(4);
 dPhic[2]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_2c->Clone();
 fitCalo("horizontalslice",2,30000,215500);
 Rp[2]=precessf->GetParameter(3);
 dRp[2]=precessf->GetParError(3);
 Phip[2]=precessf->GetParameter(4);
 dPhip[2]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_3a->Clone();
 fitCalo("horizontalslice",3,30000,215500);
 Rc[3]=precessf->GetParameter(3);
 dRc[3]=precessf->GetParError(3);
 Phic[3]=precessf->GetParameter(4);
 dPhic[3]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_3c->Clone();
 fitCalo("horizontalslice",3,30000,215500);
 Rp[3]=precessf->GetParameter(3);
 dRp[3]=precessf->GetParError(3);
 Phip[3]=precessf->GetParameter(4);
 dPhip[3]=precessf->GetParError(4); 
 htmpdriftcorrection=(TH1D*)h_4a->Clone();
 fitCalo("horizontalslice",4,30000,215500);
 Rc[4]=precessf->GetParameter(3);
 dRc[4]=precessf->GetParError(3);
 Phic[4]=precessf->GetParameter(4);
 dPhic[4]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_4c->Clone();
 fitCalo("horizontalslice",4,30000,215500);
 Rp[4]=precessf->GetParameter(3);
 dRp[4]=precessf->GetParError(3);
 Phip[4]=precessf->GetParameter(4);
 dPhip[4]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_5a->Clone();
 fitCalo("horizontalslice",5,30000,215500);
 Rc[5]=precessf->GetParameter(3);
 dRc[5]=precessf->GetParError(3);
 Phic[5]=precessf->GetParameter(4);
 dPhic[5]=precessf->GetParError(4);
 htmpdriftcorrection=(TH1D*)h_5c->Clone();
 fitCalo("horizontalslice",5,30000,215500);
 Rp[5]=precessf->GetParameter(3);
 dRp[5]=precessf->GetParError(3);
 Phip[5]=precessf->GetParameter(4);
 dPhip[5]=precessf->GetParError(4);
 doDriftCorrection=0;
 for (int ip = 0; ip < 6; ip++)  printf(" %i %f %f %f %f %f \n", ip, Rp[ip], dRp[ip], Rc[ip], dRc[ip], Rp[ip]-Rc[ip]);
 for (int ip = 0; ip < 6; ip++)  printf(" %i %f %f %f %f %f \n", ip, Phip[ip], dPhip[ip], Phic[ip], dPhic[ip], Phip[ip]-Phic[ip]);

for (int ip = 0; ip < 6; ip++) diffR[ip] = Rp[ip]-Rc[ip];
 for (int ip = 0; ip < 6; ip++) diffPhi[ip] = Phip[ip]-Phic[ip];
 
 c1->Clear();
 c1->Divide(1,2);
 c1->cd(1);
 gRbefore = new TGraphErrors(6,xR,Rc,dxR,dRc);
 gRbefore->SetTitle("Rbefore vs y-coordinate");
 gRbefore->SetMarkerStyle(8);
 gRbefore->SetMarkerSize(2);
 gRbefore->SetMarkerColor(kRed);
 gRbefore->Draw("AP");
 gRafter = new TGraphErrors(6,xR,Rp,dxR,dRp);
 gRafter->SetTitle("Rafter vs y-coordinate");
 gRafter->SetMarkerStyle(8);
 gRafter->SetMarkerSize(2);
 gRafter->SetMarkerColor(kGreen);
 gRafter->Draw("psame");
 c1->cd(2);
 gRdiff = new TGraphErrors(6,xR,diffR,dxR,dRp);
 gRdiff->SetTitle("Rdiff vs y-coordinate");
 gRdiff->SetMarkerStyle(8);
 gRdiff->SetMarkerSize(2);
 gRdiff->SetMarkerColor(kBlue);
 gRdiff->Draw("ap");
 c1->SaveAs("gRdiff.png");

 c1->Clear();
 c1->Divide(1,2);
 c1->cd(1);
 gPhibefore = new TGraphErrors(6,xR,Phic,dxR,dPhic);
 gPhibefore->SetTitle("Phibefore vs y-coordinate");
 gPhibefore->SetMarkerStyle(8);
 gPhibefore->SetMarkerSize(2);
 gPhibefore->SetMarkerColor(kRed);
 gPhibefore->Draw("AP");
 gPhiafter = new TGraphErrors(6,xR,Phip,dxR,dPhip);
 gPhiafter->SetTitle("Phiafter vs y-coordinate");
 gPhiafter->SetMarkerStyle(8);
 gPhiafter->SetMarkerSize(2);
 gPhiafter->SetMarkerColor(kGreen);
 gPhiafter->Draw("psame");
 c1->cd(2);
 gPhidiff = new TGraphErrors(6,xR,diffPhi,dxR,dPhip);
 gPhidiff->SetTitle("Phidiff vs y-coordinate");
 gPhidiff->SetMarkerStyle(8);
 gPhidiff->SetMarkerSize(2);
 gPhidiff->SetMarkerColor(kBlue);
 gPhidiff->Draw("ap");
c1->SaveAs("gPhidiff.png");

 return;
}

void configurefit();
void configurefit(){

 // free muon lifetime
 fix_gammatau=0;
 set_gammatau=64500;
 
 // no fit-method fr correction
 fastRotationAnalysis = 0;
 fastRotationCorrection = 0;
 scalefastrotation=1.0;
 timeShiftMultiplier=1.0;

 // switch on muon loss correction
 fix_muonloss_amp = 0;
 set_muonloss_amp = 0.0;
 // switch on pedestal ringing correction
 fix_pedestaldrift_amp = 0;
 set_pedestaldrift_amp = 0.0;

 // set cbo parameters, vary vw ,vo
 set_omega_vbo_60hr=9.75196e-01; // from "full" (old 9.75112e-01;)
 set_omega_vbo2_60hr=1.03036e+00; // from "full" (old 1.03043e+00;)
 set_omega_vbo_9day=9.67629e-01; // from "full" (old 9.67436e-01;)
 set_omega_vbo2_9day=1.00726e+00; // from "full" (old 1.00717e+00;)
 set_omega_vbo_hk=9.67629e-01; // from "full" 
 set_omega_vbo2_hk=1.00726e+00; // from "full"
 set_omega_vbo_endgame = 9.74178e-01; // from "full" (old 9.75382e-01;)
 set_omega_vbo2_endgame = 1.00776e+00; // from "full" (old 1.00662e+00;) 
 
 // these are poorly deterimed from "full" fit
 set_tau_vbo_60hr=2.04e+04;  //from "full" (old 8.79011e+04;)
 set_tau_vbo2_60hr=1.74268e+04; //from "full" (old 1.87389e+04;)
 set_tau_vbo_9day=4.28368e+04; // from "full" (old 4.44523e+04;)
 set_tau_vbo2_9day=8.21351e+04; // from "full" (old 7.90928e+04;)
 set_tau_vbo_hk=3.88298e+04; // from "full" 
 set_tau_vbo2_hk=7.36571e+04; // from "full"
 set_tau_vbo_endgame=3.75748e+04; // from "full" (old 2.08418e+04;)
 set_tau_vbo2_endgame=2.e4; //2.84686e+06; // from "full" (old 5.84952e+04;)
 
 // these are estimates from horizontal slice fits
 set_tau_vbo_60hr=5.e+04;  //WHY LONGER? F_VW ~ F_VB ? - from "full" (old 8.79011e+04;)
 set_tau_vbo2_60hr=2.e+04; //WHY SHORTER F_VW ~ F_VB ? - from "full" (old 1.87389e+04;)
 set_tau_vbo_endgame=3.e+04; // from "full" (old 2.08418e+04;)
 set_tau_vbo2_endgame=1.0e5; //2.84686e+06; // from "full" (old 5.84952e+04;)
 set_tau_vbo_9day=5.e4; // from "full" (old 4.44523e+04;)
 set_tau_vbo2_9day=1.1e5; // from "full" (old 7.90928e+04;)
 set_tau_vbo_hk=2.e+04; // from "full" 
 set_tau_vbo2_hk=1.2e+05; // from "full"
 
 // set cbo freq change parameters
 // pre 1/9/2020 values
 set_AExpTermCBO_60hr = 2.90; // from tracker
 set_BExpTermCBO_60hr = 5.12; // from tracker
 set_AExpTermCBO_9day = 2.86;
 set_BExpTermCBO_9day = 5.50;
 set_AExpTermCBO_hk = 2.86;
 set_BExpTermCBO_hk = 5.50;
 set_AExpTermCBO_lk = 2.86;
 set_BExpTermCBO_lk = 5.50;
 set_AExpTermCBO_eg = 6.626; // best fit calo data
 set_BExpTermCBO_eg = 6.821; // best fit calo data
// post 1/9/2020 values
 set_AExpTermCBO_60hr = 2.79; // from tracker
 set_BExpTermCBO_60hr = 5.63; // from tracker
 set_AExpTermCBO_9day = 2.80;
 set_BExpTermCBO_9day = 6.18;
 set_AExpTermCBO_hk = 2.86; // TO DO
 set_BExpTermCBO_hk = 5.50; // TO DO
 set_AExpTermCBO_lk = 2.86; // TO DO
 set_BExpTermCBO_lk = 5.50; // TO DO
 set_AExpTermCBO_eg = 6.82; // from tracker
 set_BExpTermCBO_eg = 5.42; // from tracker
 
 // pre 1/9/2020 values
 set_AExpTauTermCBO_60hr = 81.8e3; // from tracker
 set_BExpTauTermCBO_60hr = 7.7e3; // from tracker
 set_AExpTauTermCBO_9day = 72.8e3;
 set_BExpTauTermCBO_9day = 8.5e3;
 set_AExpTauTermCBO_hk = 72.8e3; // TO DO
 set_BExpTauTermCBO_hk = 8.5e3; // TO DO
 set_AExpTauTermCBO_lk = 72.8e3; // TO DO
 set_BExpTauTermCBO_lk = 8.5e3; // TO DO
 set_AExpTauTermCBO_eg = 81.8e3; // from tracker
 set_BExpTauTermCBO_eg = 7.7e3; // from tracker
 // post 1/9/2020 values
 set_AExpTauTermCBO_60hr = 61.1e3; // from tracker
 set_BExpTauTermCBO_60hr = 6.07e3; // from tracker
 set_AExpTauTermCBO_9day = 56.6e3;
 set_BExpTauTermCBO_9day = 6.32e3;
 set_AExpTauTermCBO_hk = 56.6e3;  // TO DO
 set_BExpTauTermCBO_hk = 6.32e3;  // TO DO
 set_AExpTauTermCBO_lk = 56.6e3;  // TO DO
 set_BExpTauTermCBO_lk = 6.32e3;  // TO DO
 set_AExpTauTermCBO_eg = 78.3e3; // from tracker
 set_BExpTauTermCBO_eg = 6.54e3; // from tracker

 // vary horizontal waist, 2-3 sigma effect (60hr)
 fix_2cbo_amp = 0;
 fix_2cbo_phi = 0;
 // fix horizontal cbo term in precession phase ~1-2 sigma effect (60hr)
 fix_anomphi_amp=0;
 fix_anomphi_phi=0;

 // absolute units for residuals
 fracResiUnits=1;

 // free relax, no pu, no bkd
 fix_asymtau=0;
 set_asymtau=1.e8;
 fix_pileup=1;
 fix_bkd=1;

 // vertical drift
 use_verticaldrift_correction = 1;
 set_gain_amp = 0.e-3;
 fix_gain_amp = 0;
 set_gain_tau = 1.e12;
 fix_gain_tau = 1;

 // nominal caloscan
 fix_AExpTermCBO=1;
 fix_BExpTermCBO=1;
 fix_vbo_tau = 1;
 fix_vbo_tau2 = 1;
 fix_vbo_omega = 1;
 fix_vbo_omega2 = 1;
 fix_2cbo_amp = 0;
 fix_2cbo_phi = 0;
 fix_anomphi_amp=0;
 fix_anomphi_phi=0;
 fix_asymtau=1;
 set_asymtau=7.0e7;

 // do minos, not migrad, errors
 doBlendBins=0;
 doUserCovariance=0;
 doMinos=0;
 sprintf(datasubset,"early-late");

 sprintf(fname,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamepedestal,"noxtaldata/v9.17-w4n10t15prd1e0infillwithstdp-OOFyes-gap11-eg-pufix-withpederr-gainbugfix-newenergyscale-SRC.root");
 sprintf(fnamemuonloss,"muonloss/60hours-DQC.root");
 sprintf(dname,"QFillByFillAnalyzer");
 sprintf(dataname,"endgame");
 sprintf(fitname, "vbo2");
 ithreshold=15;
 iwind=4;

 return 0;
}

double xT[10]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
double BR[10][6], dBR[10][6], Phi[10][6], dPhi[10][6];
double N[10][6], dN[10][6], A[10][6], dA[10][6];
double BRfull[10], dBRfull[10], Phifull[10], dPhifull[10];
double Nfull[10], dNfull[10], Afull[10], dAfull[10];
double BRfit[6], dBRfit[6], Phifit[6], dPhifit[6];
double N_Ampl_fit[10], N_dAmpl_fit[10], N_Cent_fit[10], N_dCent_fit[10], N_Sgma_fit[10], N_dSgma_fit[10];
double A_Ampl_fit[10], A_dAmpl_fit[10], A_Cent_fit[10], A_dCent_fit[10], A_Sgma_fit[10], A_dSgma_fit[10];
double Phi_Ampl_fit[10], Phi_dAmpl_fit[10], Phi_Cent_fit[10], Phi_dCent_fit[10], Phi_Sgma_fit[10], Phi_dSgma_fit[10], Phi_Ofst_fit[10], Phi_dOfst_fit[10];
TGraphErrors *grphN_Ampl_fit, *grphN_Cent_fit, *grphN_Sgma_fit;
TGraphErrors *grphA_Ampl_fit, *grphA_Cent_fit, *grphA_Sgma_fit;
TGraphErrors *grphPhi_Ampl_fit, *grphPhi_Cent_fit, *grphPhi_Sgma_fit, *grphPhi_Ofst_fit;
double Afit[10][6], dAfit[6];
TGraphErrors *grphBR[10], *grphN[10], *grphPhi[10], *grphA[10];
TGraphErrors *grphBRdiff[10], *grphNdiff[10], *grphPhidiff[10], *grphAdiff[10];
TGraphErrors *grphBRfdiff[10], *grphNfdiff[10], *grphPhifdiff[10], *grphAfdiff[10]; 
TGraph *grphTs, *grphTe, *grphTi, *grphTm;
TGraphErrors *grphPhiMean, *grphPhiMean2, *grphPhiMean3, *grphPhifull;
double  tstart[10], tend[10], tincr[10], tmean[10];
char pngname[128];
int Nt;
TF1 *f_N[10], *f_A[10], *f_Phi[10];
TF1 *f_Nspline[10], *f_Aspline[10], *f_Phispline[10];
double Ny[10]={0}, Ny2[10]={0}, Ay[10]={0}, Phiy[10]={0}, meanphi[10]={0}, meanphi2[10]={0}, meanphi3[10]={0};
bool Ngausfit = 0, Phigausfit = 0;
bool FIXSlide = 0, FIXStretch = 0, FIXAmplitude = 0;
bool USEAsymmetry = 1;
TF1 *polf;
char vname[64], fpolname[64];

TH1D  *hN0 = new TH1D( "hN0", "hN0", 8, -0.5, 7.5);
TSpline3 *hSN0;
TH1D  *hA0 = new TH1D( "hA0", "hA0", 8, -0.5, 7.5);
TSpline3 *hSA0;
TH1D  *hPhi0 = new TH1D( "hPhi0", "hPhi0", 8, -0.5, 7.5);
TSpline3 *hSPhi0;

 
TH1 *hslce[6], *hslcesum, *hfll; 
TF1 *fnslce[6];

Double_t N_fspline(Double_t *x, Double_t *par);
Double_t N_fspline(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     // par[1] - stretch
     // par[2] - amplitude
     
    double f;
    double xx = x[0];
 
    hN0->SetBins( 8, (-0.5+par[0])*par[1], (7.5+par[0])*par[1]);
    hSN0 = new TSpline3(hN0,"",0.1,6.9);
    hSN0->SetNpx(1000);

    return par[2]*hSN0->Eval(xx);
  }

Double_t A_fspline(Double_t *x, Double_t *par);
Double_t A_fspline(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     // par[1] - stretch
     // par[2] - amplitude
     
    double f;
    double xx = x[0];
 
    hA0->SetBins( 8, (-0.5+par[0])*par[1], (7.5+par[0])*par[1]);
    hSA0 = new TSpline3(hA0,"",0.1,6.9);
    hSA0->SetNpx(1000);

    return par[2]*hSA0->Eval(xx);
  }

Double_t Phi_fspline(Double_t *x, Double_t *par);
Double_t Phi_fspline(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     // par[1] - stretch
     // par[2] - amplitude
     
    double f;
    double xx = x[0];
 
    hPhi0->SetBins( 8, (-0.5+par[0])*par[1], (7.5+par[0])*par[1]);
    hSPhi0 = new TSpline3(hPhi0,"",0.1,6.9);
    hSPhi0->SetNpx(1000);

    return par[2]*hSPhi0->Eval(xx);
  }


Double_t N_f(Double_t *x, Double_t *par);
Double_t N_f(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     //
     
    int is = (int)x[0];
    double f;

    if (is == 1) { // top slice = 1
      if (par[0] >= 0.0) f = (1.0-par[0])*N[0][0] + par[0]*N[0][1]; // slide upwards, gets energy from slice 2
      if (par[0] <  0.0) f = (1.0+par[0])*N[0][0]; // slides downwards, and gets nothing from region above (assumption)
      return f;
    }
    if (is == 6) {// bot slice = 6
      if (par[0] >= 0.0) f = (1.0-par[0])*N[0][5]; // sliding upwards, and gets nothing, from region below (assumption)
      if (par[0] <  0.0) f = (1.0+par[0])*N[0][5] - par[0]*N[0][4]; // sliding downwards, gets energy from slice 4
      return f;
    }

    // other slices
    if (par[0] >= 0.0) f = (1.0-par[0])*N[0][is-1] + par[0]*N[0][is]; // slides upwards, and gets energy from slice below
    if (par[0] <  0.0) f = (1.0+par[0])*N[0][is-1] - par[0]*N[0][is-2]; // slides downwards, and gets energy from slice above
    return f;
  }

Double_t Phi_f(Double_t *x, Double_t *par);
Double_t Phi_f(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     //
     
    int is = (int)x[0];
    double f;

    if (is == 1) {
      if (par[0] >= 0.0) f = (1.0-par[0])*Phi[0][0] + par[0]*Phi[0][1];
      if (par[0] <  0.0) f = (1.0+par[0])*Phi[0][0];
      return f;
    }
    if (is == 6) {
      if (par[0] >= 0.0) f = (1.0-par[0])*Phi[0][5];
      if (par[0] <  0.0) f = (1.0+par[0])*Phi[0][5] - par[0]*Phi[0][4];
      return f;
    }

    if (par[0] >= 0.0) f = (1.0-par[0])*Phi[0][is-1] + par[0]*Phi[0][is];
    if (par[0] <  0.0) f = (1.0+par[0])*Phi[0][is-1] - par[0]*Phi[0][is-2];
    return f;
  }

Double_t A_f(Double_t *x, Double_t *par);
Double_t A_f(Double_t *x, Double_t *par)
  {
     // 
     // y-distribution
     // par[0] - shift
     //
     
    int is = (int)x[0];
    double f;

    if (is == 1) {
      if (par[0] >= 0.0) f = (1.0-par[0])*A[0][0] + par[0]*A[0][1];
      if (par[0] <  0.0) f = (1.0+par[0])*A[0][0];
      return f;
    }
    if (is == 6) {
      if (par[0] >= 0.0) f = (1.0-par[0])*A[0][5];
      if (par[0] <  0.0) f = (1.0+par[0])*A[0][5] - par[0]*A[0][4];
      return f;
    }

    if (par[0] >= 0.0) f = (1.0-par[0])*A[0][is-1] + par[0]*A[0][is];
    if (par[0] <  0.0) f = (1.0+par[0])*A[0][is-1] - par[0]*A[0][is-2];
    return f;
  }

void fitslices();

int earlylatephase();
int earlylatephase(){

  int is, it;

  fix_R = 1;

  for (it = 0; it < Nt; it++){
    
    tend[it] = tstart[it] + tincr[it];
    tmean[it] = ( tstart[it] +  tend[it] )/2.;
    for (is = 0; is < 6; is++){
      
      fitCalo("horizontalslice",is,tstart[it],tend[it]);
      N[it][is]=precessf->GetParameter(0);
      dN[it][is]=precessf->GetParError(0);
      A[it][is]=precessf->GetParameter(2);
      dA[it][is]=precessf->GetParError(2);
      BR[it][is]=precessf->GetParameter(3);
      dBR[it][is]=precessf->GetParError(3);
      Phi[it][is]=precessf->GetParameter(4);
      dPhi[it][is]=precessf->GetParError(4);
    }

    fitCalo("full",0,tstart[it],tend[it]);
    Nfull[it]=precessf->GetParameter(0);
    dNfull[it]=precessf->GetParError(0);
    Afull[it]=precessf->GetParameter(2);
    dAfull[it]=precessf->GetParError(2);
    BRfull[it]=precessf->GetParameter(3);
    dBRfull[it]=precessf->GetParError(3);
    Phifull[it]=precessf->GetParameter(4);
    dPhifull[it]=precessf->GetParError(4);
    
  }

  grphPhifull = new TGraphErrors(Nt,tmean,Phifull,dxR,dPhifull);
  grphPhifull->SetTitle("PhiFull vs time");
  grphPhifull->SetMarkerStyle(8);
  grphPhifull->SetMarkerSize(1.2);
  grphPhifull->SetMarkerColor(1);
  grphPhifull->Draw("ap");

  if (doDraw2){
  sprintf(pngname,"Phipull-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  grphTs = new TGraph(Nt,xT,tstart);
  grphTe = new TGraph(Nt,xT,tend);
  grphTi = new TGraph(Nt,xT,tincr);
  grphTm = new TGraph(Nt,xT,tmean);
  grphTs->SetTitle("Ts, Te, Ti time window parameters");
  grphTs->SetMarkerStyle(8);
  grphTs->SetMarkerSize(1.2);
  grphTs->SetMarkerColor(1);
  grphTs->SetLineColor(1);
  grphTe->SetMarkerStyle(8);
  grphTe->SetMarkerSize(1.2);
  grphTe->SetMarkerColor(2);
  grphTe->SetLineColor(2);
  grphTi->SetMarkerStyle(8);
  grphTi->SetMarkerSize(1.2);
  grphTi->SetMarkerColor(3);
  grphTi->SetLineColor(3);
  grphTm->SetMarkerStyle(8);
  grphTm->SetMarkerSize(1.2);
  grphTm->SetMarkerColor(4);
  grphTm->SetLineColor(4);
  grphTs->SetMaximum(2.3e5);
  grphTs->SetMinimum(0.0);
  
  c1->Clear();
  c1->Divide(1,1);
  c1->cd(1);
  grphTs->Draw("ap");
  grphTe->Draw("psame");
  grphTi->Draw("psame");
  grphTm->Draw("psame");
  if (doDraw2){
  sprintf(pngname,"timewindows-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  double ytmp[6], dytmp[6];

  c1->Clear();
  c1->Divide(1,3);
  
  c1->cd(1);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = N[it][is];
    for (is = 0; is < 6; is++) dytmp[is] = dN[it][is];
    grphN[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphN[it]->SetTitle("N vs y-coordinate");
    grphN[it]->SetMarkerStyle(8);
    grphN[it]->SetMarkerSize(1.2);
    grphN[it]->SetMarkerColor(1+it);
    if (it==9) grphN[9]->SetMarkerColor(11);
    if (it==9) grphN[9]->SetLineColor(11);
    if (it==0){
      grphN[0]->Draw("ap");
    } else {
      grphN[it]->Draw("psame");
    }
  }

  c1->cd(2);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = N[it][is]-N[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dN[it][is];
    grphNdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphNdiff[it]->SetTitle("Ndiff vs y-coordinate");
    grphNdiff[it]->SetMarkerStyle(8);
    grphNdiff[it]->SetMarkerSize(1.2);
    grphNdiff[it]->SetMarkerColor(1+it);
    grphNdiff[it]->SetLineColor(1+it);
    grphNdiff[it]->SetLineWidth(2);
    if (it==9) grphNdiff[9]->SetMarkerColor(11);
    if (it==9) grphNdiff[9]->SetLineColor(11);
    if (it==0){
      grphNdiff[0]->SetMaximum(2.5e8);
      grphNdiff[0]->SetMinimum(-2.5e8);
      grphNdiff[0]->Draw("apl");
    } else {
      grphNdiff[it]->Draw("plsame");
    }
  }
   
  c1->cd(3);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = ( N[it][is]-N[0][is] ) / N[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dN[it][is] / N[0][is];
    grphNfdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphNfdiff[it]->SetTitle("Nfracdiff vs y-coordinate");
    grphNfdiff[it]->SetMarkerStyle(8);
    grphNfdiff[it]->SetMarkerSize(1.2);
    grphNfdiff[it]->SetMarkerColor(1+it);
    grphNfdiff[it]->SetLineColor(1+it);
    grphNfdiff[it]->SetLineWidth(2);
    if (it==9) grphNfdiff[9]->SetMarkerColor(11);
    if (it==9) grphNfdiff[9]->SetLineColor(11);
    if (it==0){
      grphNfdiff[0]->SetMaximum(0.08);
      grphNfdiff[0]->SetMinimum(-0.08);
      grphNfdiff[0]->Draw("apl");
    } else {
      grphNfdiff[it]->Draw("plsame");
    }
  }

  if (doDraw2){
  sprintf(pngname,"N-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  c1->Clear();
  c1->Divide(1,3);
  
  c1->cd(1);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = Phi[it][is];
    for (is = 0; is < 6; is++) dytmp[is] = dPhi[it][is];
    grphPhi[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphPhi[it]->SetTitle("Phi vs y-coordinate");
    grphPhi[it]->SetMarkerStyle(8);
    grphPhi[it]->SetMarkerSize(1.2);
    grphPhi[it]->SetMarkerColor(1+it);
    if (it==9) grphPhi[9]->SetMarkerColor(11);
    if (it==9) grphPhi[9]->SetLineColor(11);
    if (it==0){
      grphPhi[0]->Draw("ap");
    } else {
      grphPhi[it]->Draw("psame");
    }
  }

  c1->cd(2);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = Phi[it][is]-Phi[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dPhi[it][is];
    grphPhidiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphPhidiff[it]->SetTitle("Phidiff vs y-coordinate");
    grphPhidiff[it]->SetMarkerStyle(8);
    grphPhidiff[it]->SetMarkerSize(1.2);
    grphPhidiff[it]->SetMarkerColor(1+it);
    grphPhidiff[it]->SetLineColor(1+it);
    if (it==9) grphPhidiff[9]->SetMarkerColor(11);
    if (it==9) grphPhidiff[9]->SetLineColor(11);
    if (it==0){
      grphPhidiff[0]->SetMaximum(0.04);
      grphPhidiff[0]->SetMinimum(-0.04);
      grphPhidiff[0]->Draw("apl");
    } else {
      grphPhidiff[it]->Draw("plsame");
    }
  }

  c1->cd(3);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = ( Phi[it][is]-Phi[0][is] ) / Phi[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dPhi[it][is] / Phi[0][is];
    grphPhifdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphPhifdiff[it]->SetTitle("Phifracdiff vs y-coordinate");
    grphPhifdiff[it]->SetMarkerStyle(8);
    grphPhifdiff[it]->SetMarkerSize(1.2);
    grphPhifdiff[it]->SetMarkerColor(1+it);
    grphPhifdiff[it]->SetLineColor(1+it);
    grphPhifdiff[it]->SetLineWidth(2);
    if (it==9) grphPhifdiff[9]->SetMarkerColor(11);
    if (it==9) grphPhifdiff[9]->SetLineColor(11);
    if (it==0){
      grphPhifdiff[0]->SetMaximum(0.02);
      grphPhifdiff[0]->SetMinimum(-0.02);
      grphPhifdiff[0]->Draw("apl");
    } else {
      grphPhifdiff[it]->Draw("plsame");
    }
  }


  if (doDraw2){
  sprintf(pngname,"Phi-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  c1->Clear();
  c1->Divide(1,3);
  
  c1->cd(1);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = A[it][is];
    for (is = 0; is < 6; is++) dytmp[is] = dA[it][is];
    grphA[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphA[it]->SetTitle("A vs y-coordinate");
    grphA[it]->SetMarkerStyle(8);
    grphA[it]->SetMarkerSize(1.2);
    grphA[it]->SetMarkerColor(1+it);
    if (it==9) grphA[9]->SetMarkerColor(11);
    if (it==9) grphA[9]->SetLineColor(11);
    if (it==0){
      grphA[0]->Draw("ap");
    } else {
      grphA[it]->Draw("psame");
    }
  }

  c1->cd(2);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = A[it][is]-A[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dA[it][is];
    grphAdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphAdiff[it]->SetTitle("Adiff vs y-coordinate");
    grphAdiff[it]->SetMarkerStyle(8);
    grphAdiff[it]->SetMarkerSize(1.2);
    grphAdiff[it]->SetMarkerColor(1+it);
    grphAdiff[it]->SetLineColor(1+it);
    if (it==9) grphAdiff[9]->SetMarkerColor(11);
    if (it==9) grphAdiff[9]->SetLineColor(11);
    if (it==0){
      grphAdiff[0]->SetMaximum(0.01);
      grphAdiff[0]->SetMinimum(-0.01);
      grphAdiff[0]->Draw("apl");
    } else {
      grphAdiff[it]->Draw("plsame");
    }
  }
  grphAdiff[0]->Draw("plsame");

  c1->cd(3);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = ( A[it][is]-A[0][is] ) / A[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dA[it][is] / A[0][is];
    grphAfdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphAfdiff[it]->SetTitle("Afracdiff vs y-coordinate");
    grphAfdiff[it]->SetMarkerStyle(8);
    grphAfdiff[it]->SetMarkerSize(1.2);
    grphAfdiff[it]->SetMarkerColor(1+it);
    grphAfdiff[it]->SetLineColor(1+it);
    grphAfdiff[it]->SetLineWidth(2);
    if (it==9) grphAfdiff[9]->SetMarkerColor(11);
    if (it==9) grphAfdiff[9]->SetLineColor(11);
    if (it==0){
      grphAfdiff[0]->SetMaximum(0.04);
      grphAfdiff[0]->SetMinimum(-0.04);
      grphAfdiff[0]->Draw("apl");
    } else {
      grphAfdiff[it]->Draw("plsame");
    }
  }

  if (doDraw2){
  sprintf(pngname,"A-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  c1->Clear();
  c1->Divide(1,2);
  
  c1->cd(1);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = BR[it][is];
    for (is = 0; is < 6; is++) dytmp[is] = dBR[it][is];
    grphBR[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphBR[it]->SetTitle("N vs y-coordinate");
    grphBR[it]->SetMarkerStyle(8);
    grphBR[it]->SetMarkerSize(1.2);
    grphBR[it]->SetMarkerColor(1+it);
    if (it==9) grphNdiff[9]->SetMarkerColor(11);
    if (it==9) grphNdiff[9]->SetLineColor(11);
    if (it==0){
      grphBR[0]->Draw("ap");
    } else {
      grphBR[it]->Draw("psame");
    }
  }

  c1->cd(2);
  for (it = 0; it < Nt; it++){
    
    for (is = 0; is < 6; is++) ytmp[is] = BR[it][is]-BR[0][is];
    for (is = 0; is < 6; is++) dytmp[is] = dBR[it][is];
    grphBRdiff[it] = new TGraphErrors(6,xR,ytmp,dxR,dytmp);
    grphBRdiff[it]->SetTitle("diff vs y-coordinate");
    grphBRdiff[it]->SetMarkerStyle(8);
    grphBRdiff[it]->SetMarkerSize(1.2);
    grphBRdiff[it]->SetMarkerColor(1+it);
    if (it==9) grphNdiff[9]->SetMarkerColor(11);
    if (it==9) grphNdiff[9]->SetLineColor(11);
    if (it==0){
      grphBRdiff[0]->Draw("apl");
    } else {
      grphBRdiff[it]->Draw("plsame");
    }
  }
  grphBRdiff[0]->Draw("plsame");

  if (doDraw2){
  sprintf(pngname,"BR-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  char f_N_name[64];
  int ib;
  for (ib = 0; ib < 6; ib++) hN0->SetBinContent( ib+2, N[0][ib]);
  hN0->SetBinContent( 1, N[0][0]-(N[0][1]-N[0][0]));
  hN0->SetBinContent( 8, N[0][6]-(N[0][5]-N[0][6]));

  for (it = 0; it < Nt; it++){
    sprintf( f_N_name, "f_N%01i", it);
    printf("fit %s\n", f_N_name);

    if (Ngausfit) {

    f_N[it] = new TF1( f_N_name, "[0]*exp(-(x-[1])*(x-[1])/[2]/[2])", -5., +12.);
    f_N[it]->SetParNames("Amplitude","Mean","Sigma");
    f_N[it]->SetNpx(10000);
    f_N[it]->SetParameters(8.e8,3.5,1.5);
    grphN[it]->Fit( f_N_name, "", "NIRE", 0.5, 6.5);

    N_Ampl_fit[it] = f_N[it]->GetParameter(0);
    N_dAmpl_fit[it] = f_N[it]->GetParError(0);
    N_Cent_fit[it] = f_N[it]->GetParameter(1);
    N_dCent_fit[it]  = f_N[it]->GetParError(1);
    N_Sgma_fit[it] = f_N[it]->GetParameter(2);
    N_dSgma_fit[it]  = f_N[it]->GetParError(2);

    } else {

    f_Nspline[it] = new TF1( f_N_name, N_fspline, 0.1, 6.9, 3);
    f_Nspline[it]->SetParNames("Shft-N","Strtch-N","Norm-N");
    f_Nspline[it]->SetNpx(10000);
    f_Nspline[it]->SetParameter(0,0.0);
    if (FIXSlide) f_Phispline[it]->FixParameter(0,0.0);
    f_Nspline[it]->SetParameter(1,1.0);
    if (FIXStretch) f_Nspline[it]->FixParameter(1,1.0);
    f_Nspline[it]->SetParameter(2,1.0);
    if (FIXAmplitude) f_Nspline[it]->FixParameter(2,1.0);
    grphN[it]->Fit( f_N_name, "", "NIRE", 0.5, 6.5);

    N_Ampl_fit[it] = f_Nspline[it]->GetParameter(0);
    N_dAmpl_fit[it] = f_Nspline[it]->GetParError(0);
    N_Cent_fit[it] = f_Nspline[it]->GetParameter(1);
    N_dCent_fit[it]  = f_Nspline[it]->GetParError(1);
    N_Sgma_fit[it] = f_Nspline[it]->GetParameter(2);
    N_dSgma_fit[it]  = f_Nspline[it]->GetParError(2);

    }

  }

  c1->Clear();
  c1->Divide(Nt/2,2);
  for (it = 0; it < Nt; it++){
    c1->cd(it+1);
    grphN[it]->Draw();
  }

  if (doDraw2){
  sprintf(pngname,"N-fits-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  char f_A_name[64]; 
 for (ib = 0; ib < 6; ib++) hA0->SetBinContent( ib+2, A[0][ib]);
 hA0->SetBinContent( 1, A[0][0]-(A[0][1]-A[0][0]));
 hA0->SetBinContent( 8, A[0][6]-(A[0][5]-A[0][6]));

  for (it = 0; it < Nt; it++){
    sprintf( f_A_name, "f_A%01i", it);
    printf("fit %s\n", f_A_name);

    //f_A[it] = new TF1( f_A_name, "[0]*exp(-(x-[1])*(x-[1])/[2]/[2])", -5., +12.);
    //f_A[it]->SetParNames("Amplitude","Mean","Sigma");
    //f_A[it]->SetNpx(10000);
    //f_A[it]->SetParameters(8.e8,3.5,1.5);

    /*
    f_A[it] = new TF1( f_A_name, N_f, -5., +12.,1);
    f_A[it]->SetParNames("Shft-A");
    f_A[it]->SetNpx(10000);
    f_A[it]->SetParameter(0,0.0);
    grphA[it]->Fit( f_A_name, "", "NRE", 0.5, 6.5);
    A_Ampl_fit[it] = f_A[it]->GetParameter(0);
    A_dAmpl_fit[it] = f_A[it]->GetParError(0);
    A_Cent_fit[it] = f_A[it]->GetParameter(1);
    A_dCent_fit[it]  = f_A[it]->GetParError(1);
    A_Sgma_fit[it] = f_A[it]->GetParameter(2);
    A_dSgma_fit[it]  = f_A[it]->GetParError(2);
    */

    f_Aspline[it] = new TF1( f_A_name, A_fspline, 0.1, 6.9, 3);
    f_Aspline[it]->SetParNames("Shft-A","Strtch-A","Norm-A");
    f_Aspline[it]->SetNpx(10000);
    f_Aspline[it]->SetParameter(0,0.0);
    if (FIXSlide) f_Phispline[it]->FixParameter(0,0.0);
    f_Aspline[it]->SetParameter(1,1.0);
    if (FIXStretch) f_Aspline[it]->FixParameter(1,1.0);
    f_Aspline[it]->SetParameter(2,1.0);
    if (FIXAmplitude) f_Aspline[it]->FixParameter(2,1.0);
    grphA[it]->Fit( f_A_name, "", "NIRE", 0.5, 6.5);

    A_Ampl_fit[it] = f_Aspline[it]->GetParameter(0);
    A_dAmpl_fit[it] = f_Aspline[it]->GetParError(0);
    A_Cent_fit[it] = f_Aspline[it]->GetParameter(1);
    A_dCent_fit[it]  = f_Aspline[it]->GetParError(1);
    A_Sgma_fit[it] = f_Aspline[it]->GetParameter(2);
    A_dSgma_fit[it]  = f_Aspline[it]->GetParError(2);

  }

  c1->Clear();
  c1->Divide(Nt/2,2);
  for (it = 0; it < Nt; it++){
    c1->cd(it+1);
    grphA[it]->Draw();
  }

  sprintf(pngname,"A-fits-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);

  char f_Phi_name[64]; 
  for (ib = 0; ib < 6; ib++) hPhi0->SetBinContent( ib+2, Phi[0][ib]);
  hPhi0->SetBinContent( 1, Phi[0][0]-(Phi[0][1]-Phi[0][0]));
  hPhi0->SetBinContent( 8, Phi[0][6]-(Phi[0][5]-Phi[0][6]));

  for (it = 0; it < Nt; it++){

    sprintf( f_Phi_name, "f_Phi%01i", it);
    printf("fit %s\n", f_Phi_name);

    if (Phigausfit) {

    f_Phi[it] = new TF1( f_Phi_name, "[0]*exp(-(x-[1])*(x-[1])/[2]/[2])+[3]", -5., +12.);
    f_Phi[it]->SetParNames("Amplitude","Mean","Sigma","Offset");
    f_Phi[it]->SetNpx(10000);
    f_Phi[it]->SetParameters(-0.264,3.57,1.88,2.35);
    grphPhi[it]->Fit( f_Phi_name, "", "NIRE", 0.5, 6.5);
    grphPhi[it]->Fit( f_Phi_name, "", "NIRE", 0.5, 6.5);  // dit two times to get fit convergence

    Phi_Ampl_fit[it] = f_Phi[it]->GetParameter(0);
    Phi_dAmpl_fit[it] = f_Phi[it]->GetParError(0);
    Phi_Cent_fit[it] = f_Phi[it]->GetParameter(1);
    Phi_dCent_fit[it]  = f_Phi[it]->GetParError(1);
    Phi_Sgma_fit[it] = f_Phi[it]->GetParameter(2);
    Phi_dSgma_fit[it]  = f_Phi[it]->GetParError(2);
    Phi_Ofst_fit[it] = f_Phi[it]->GetParameter(3);
    Phi_dOfst_fit[it]  = f_Phi[it]->GetParError(3);
    } else {

    f_Phispline[it] = new TF1( f_Phi_name, Phi_fspline, 0.1, 6.9, 3);
    f_Phispline[it]->SetParNames("Shft-Phi","Strtch-Phi","Strtch-Phi");
    f_Phispline[it]->SetNpx(10000);
    f_Phispline[it]->SetParameter(0,0.0);
    if (FIXSlide) f_Phispline[it]->FixParameter(0,0.0);
    f_Phispline[it]->SetParameter(1,1.0);
    if (FIXStretch) f_Phispline[it]->FixParameter(1,1.0);
    f_Phispline[it]->SetParameter(2,1.0);
    if (FIXAmplitude) f_Phispline[it]->FixParameter(2,1.0);
    grphPhi[it]->Fit( f_Phi_name, "", "NIRE", 0.5, 6.5);

    Phi_Ampl_fit[it] = f_Phispline[it]->GetParameter(0);
    Phi_dAmpl_fit[it] = f_Phispline[it]->GetParError(0);
    Phi_Cent_fit[it] = f_Phispline[it]->GetParameter(1);
    Phi_dCent_fit[it]  = f_Phispline[it]->GetParError(1);
    Phi_Sgma_fit[it] = f_Phispline[it]->GetParameter(2);
    Phi_dSgma_fit[it]  = f_Phispline[it]->GetParError(2);

    }
  }

  c1->Clear();
  c1->Divide(Nt/2,2);
  for (it = 0; it < Nt; it++){
    c1->cd(it+1);
    grphPhi[it]->Draw();
  }

  if (doDraw2){
  sprintf(pngname,"Phi-fits-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  grphN_Ampl_fit  = new TGraphErrors(Nt,tmean,N_Ampl_fit,dxR,N_dAmpl_fit);
  grphN_Ampl_fit->SetTitle("N_Ampl_fit vs time");
  grphN_Ampl_fit->SetMarkerStyle(8);
  grphN_Ampl_fit->SetMarkerSize(1.2);
  grphN_Ampl_fit->SetMarkerColor(1);
  grphN_Cent_fit  = new TGraphErrors(Nt,tmean,N_Cent_fit,dxR,N_dCent_fit);
  grphN_Cent_fit->SetTitle("N_Cent_fit vs time");
  grphN_Cent_fit->SetMarkerStyle(8);
  grphN_Cent_fit->SetMarkerSize(1.2);
  grphN_Cent_fit->SetMarkerColor(2);
  grphN_Sgma_fit  = new TGraphErrors(Nt,tmean,N_Sgma_fit,dxR,N_dSgma_fit);
  grphN_Sgma_fit->SetTitle("N_Sgma_fit vs time");
  grphN_Sgma_fit->SetMarkerStyle(8);
  grphN_Sgma_fit->SetMarkerSize(1.2);
  grphN_Sgma_fit->SetMarkerColor(3);
  
  c1->Clear();
  c1->Divide(3,1);
  c1->cd(1);
  grphN_Ampl_fit->Draw();
  c1->cd(2);
  grphN_Cent_fit->Draw();
  c1->cd(3);
  grphN_Sgma_fit->Draw();
  if (doDraw2){
  sprintf(pngname,"N-fitpars-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  grphA_Ampl_fit  = new TGraphErrors(Nt,tmean,A_Ampl_fit,dxR,A_dAmpl_fit);
  grphA_Ampl_fit->SetTitle("A_Ampl_fit vs time");
  grphA_Ampl_fit->SetMarkerStyle(8);
  grphA_Ampl_fit->SetMarkerSize(1.2);
  grphA_Ampl_fit->SetMarkerColor(1);
  grphA_Cent_fit  = new TGraphErrors(Nt,tmean,A_Cent_fit,dxR,A_dCent_fit);
  grphA_Cent_fit->SetTitle("A_Cent_fit vs time");
  grphA_Cent_fit->SetMarkerStyle(8);
  grphA_Cent_fit->SetMarkerSize(1.2);
  grphA_Cent_fit->SetMarkerColor(2);
  grphA_Sgma_fit  = new TGraphErrors(Nt,tmean,A_Sgma_fit,dxR,A_dSgma_fit);
  grphA_Sgma_fit->SetTitle("A_Sgma_fit vs time");
  grphA_Sgma_fit->SetMarkerStyle(8);
  grphA_Sgma_fit->SetMarkerSize(1.2);
  grphA_Sgma_fit->SetMarkerColor(3);
  
  c1->Clear();
  c1->Divide(3,1);
  c1->cd(1);
  grphA_Ampl_fit->SetMaximum(0.05);
  grphA_Ampl_fit->SetMinimum(-0.05);
  grphA_Ampl_fit->Draw("ap");
  c1->cd(2);
  grphA_Cent_fit->SetMaximum(1.05);
  grphA_Cent_fit->SetMinimum(0.95);
  grphA_Cent_fit->Draw("ap");
  c1->cd(3); 
  grphA_Sgma_fit->SetMaximum(1.05);
  grphA_Sgma_fit->SetMinimum(0.95);
  grphA_Sgma_fit->Draw("ap");
  if (doDraw2){
  sprintf(pngname,"A-fitpars-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  grphPhi_Ampl_fit  = new TGraphErrors(Nt,tmean,Phi_Ampl_fit,dxR,Phi_dAmpl_fit);
  grphPhi_Ampl_fit->SetTitle("Phi_Ampl_fit vs time");
  grphPhi_Ampl_fit->SetMarkerStyle(8);
  grphPhi_Ampl_fit->SetMarkerSize(1.2);
  grphPhi_Ampl_fit->SetMarkerColor(1);
  grphPhi_Cent_fit  = new TGraphErrors(Nt,tmean,Phi_Cent_fit,dxR,Phi_dCent_fit);
  grphPhi_Cent_fit->SetTitle("Phi_Cent_fit vs time");
  grphPhi_Cent_fit->SetMarkerStyle(8);
  grphPhi_Cent_fit->SetMarkerSize(1.2);
  grphPhi_Cent_fit->SetMarkerColor(2);
  grphPhi_Sgma_fit  = new TGraphErrors(Nt,tmean,Phi_Sgma_fit,dxR,Phi_dSgma_fit);
  grphPhi_Sgma_fit->SetTitle("Phi_Sgma_fit vs time");
  grphPhi_Sgma_fit->SetMarkerStyle(8);
  grphPhi_Sgma_fit->SetMarkerSize(1.2);
  grphPhi_Sgma_fit->SetMarkerColor(3);
  /*
  grphPhi_Ofst_fit  = new TGraphErrors(Nt,tmean,Phi_Ofst_fit,dxR,Phi_dOfst_fit);
  grphPhi_Ofst_fit->SetTitle("Phi_Sgma_fit vs time");
  grphPhi_Ofst_fit->SetMarkerStyle(8);
  grphPhi_Ofst_fit->SetMarkerSize(1.2);
  grphPhi_Ofst_fit->SetMarkerColor(3);
  */

  if (doDraw2){
  c1->Clear();
  c1->Divide(3,1);
  c1->cd(1);
  grphPhi_Ampl_fit->Draw("ap");
  c1->cd(2);
  grphPhi_Cent_fit->Draw("ap");
  c1->cd(3);
  grphPhi_Sgma_fit->Draw("ap");

  sprintf(pngname,"Phi-fitpars-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  }

  fitslices();
  
  for (int it = 0; it < Nt; it++) printf("it, shift, stretch %f, %f\n", it, N_Ampl_fit[it], N_Cent_fit[it]);

  printf("start phase shift calculation\n");

  for (int it = 0; it < Nt; it++){

    double num = 0.0, den = 0.0;
    double num2 = 0.0, den2 = 0.0;
    double num3 = 0.0, den3 = 0.0;

    for (int is = 0; is < 6; is++){


      if (Ngausfit) {
	Ny[it] = f_N[it]->Eval(is);
      } else {
	Ny[it] = f_Nspline[it]->Eval(is);
      }
      //printf("is, Ny %i %f\n", is, Ny[it]);      

      if (Phigausfit) {
	Phiy[it] = f_Phi[it]->Eval(is);
      } else {
	Phiy[it] = f_Phispline[it]->Eval(is);
      }
      //printf("is, Phiy %i %f\n", is, Phiy[it]);      

      if (USEAsymmetry) {
	Ay[it] = f_Aspline[it]->Eval(is);
      } else {
 	Ay[it] = 1.0;
      }
      //printf("is, Ay %i %f\n", is, Ay[it]);      
      
      num += Ny[it]*Phiy[it]*Ay[it];
      den += Ny[it]*Ay[it];
    }

    int istop = 0.5; // o,5
    int isbot = 6.5; // 6.5
    //mean shift 0.004095 (xtal units), 0.102381 (mm units) per 100 microseconds

    double gainphinum = 0.0, lostphinum = 0.0, gainphiden = 0.0, lostphiden = 0.0;
    if (Ngausfit) {
      gainphinum = ( 0.004095e-5*tmean[it] )*f_N[0]->Eval(istop)*f_Phispline[0]->Eval(istop)*f_Aspline[0]->Eval(istop);
      lostphinum = ( 0.004095e-5*tmean[it] )*f_N[0]->Eval(isbot)*f_Phispline[0]->Eval(isbot)*f_Aspline[0]->Eval(isbot);
      gainphiden = ( 0.004095e-5*tmean[it] )*f_N[0]->Eval(istop)*f_Aspline[0]->Eval(istop); 
      lostphiden = ( 0.004095e-5*tmean[it] )*f_N[0]->Eval(isbot)*f_Aspline[0]->Eval(isbot);
    } else {
      gainphinum = ( 0.004095e-5*tmean[it] )*f_Nspline[0]->Eval(istop)*f_Phispline[0]->Eval(istop)*f_Aspline[0]->Eval(istop);
      lostphinum = ( 0.004095e-5*tmean[it] )*f_Nspline[0]->Eval(isbot)*f_Phispline[0]->Eval(isbot)*f_Aspline[0]->Eval(isbot);
      gainphiden = ( 0.004095e-5*tmean[it] )*f_Nspline[0]->Eval(istop)*f_Aspline[0]->Eval(istop); 
      lostphiden = ( 0.004095e-5*tmean[it] )*f_Nspline[0]->Eval(isbot)*f_Aspline[0]->Eval(isbot);
    }
    printf("it, gainphinum, lostphinum, gainphiden, lostphinden %i, %e, %e, %e, %e\n", 
	it, gainphinum, lostphinum, gainphiden, lostphiden );      
    
    num2 += Ny[0]*Phiy[0]*Ay[0];
    den2 += Ny[0]*Ay[0];
    num3 = num2 + gainphinum - lostphinum;
    den3 = den2 + gainphiden - lostphiden;

    meanphi[it] = num / den;
    meanphi2[it] = num2 / den2;
    meanphi3[it] = num3 / den3;

    //meanphi[it] += 2.*TMath::Pi()*75./4370.; // silly test for displacement by one 75ns bin of slice histos
    printf("time window %i, num %f, slice %i, den %f, meanphi1,2,3 %f %f %f\n", it, is, num, den, meanphi[it], meanphi2[it], meanphi3[it]);
  }

  grphPhiMean = new TGraphErrors(Nt,tmean,meanphi,0,0);
  grphPhiMean->SetTitle("Calculated phi mean vs time");
  grphPhiMean->SetMarkerStyle(8);
  grphPhiMean->SetMarkerSize(1.2);
  grphPhiMean->SetMarkerColor(3);

  grphPhiMean2 = new TGraphErrors(Nt,tmean,meanphi2,0,0);
  grphPhiMean2->SetTitle("Calculated phi mean2 vs time");
  grphPhiMean2->SetMarkerStyle(8);
  grphPhiMean2->SetMarkerSize(1.2);
  grphPhiMean2->SetMarkerColor(3);

  grphPhiMean3 = new TGraphErrors(Nt,tmean,meanphi3,0,0);
  grphPhiMean3->SetTitle("Calculated phi mean3 vs time");
  grphPhiMean3->SetMarkerStyle(8);
  grphPhiMean3->SetMarkerSize(1.2);
  grphPhiMean3->SetMarkerColor(3);

  polf = new TF1("polf",fpolname, tstart[0],tend[Nt-1]);
  grphPhiMean3->Fit( polf, "", "NRE", tstart[0], tend[Nt-1]);

  c1->Clear();
  c1->Divide(1,3);
  c1->cd(1);
  grphPhiMean2->SetMaximum(2.40);
  grphPhiMean2->SetMinimum(2.00);
  grphPhiMean2->Draw("ap");
  c1->cd(2);
  grphPhiMean2->SetMaximum(2.40);
  grphPhiMean2->SetMinimum(2.00);
  grphPhiMean2->Draw("ap");
  //grphPhifull->SetMaximum(2.276);
  //grphPhifull->SetMinimum(2.271);
  //grphPhifull->Draw("ap");
  sprintf(pngname,"PhiMeanCalc-vs-y-noA%i-fixC%i-fixS%i-fixA%i-nt%i-tstart%6.1f-%s-%s-%s.png",USEAsymmetry,FIXSlide,FIXStretch,FIXAmplitude,Nt,tstart[0],dataname,fitname,vname);
  c1->SaveAs(pngname);
  
  /*
  fPhi= new TF1("fPhi","[0]*exp(-(x-[1])*(x-[1])/[2]/[2])+[3]",0,7);
  fPhi->SetParNames("Amplitude","Mean","Sigma","offset");
  fPhi->SetNpx(10000);
  fPhi->SetParameters(-0.2,2.5,1.5,2.5);
  */

  return 0;
}

void fitslices(){

  char hslcename[64];
  int is;
  
  fitCalo("full",0,30000,215500);
  hfll=(TH1D*)htmp->Clone();
  hfll->SetName("hfull");
  hfll->SetTitle("hfull");
    if (run_2) {    
      hfll->Rebin(4);
      hfll->Scale(1./4.);
    } 

  fitCalo("horizontalslice",0,30000,215500);
  hslcesum=(TH1D*)htmp->Clone();
  hslcesum->SetName("hslicesum");
  hslcesum->SetTitle("hslicesum");
  hslcesum->Reset();
    if (run_2) {    
      hslcesum->Rebin(4);
      hslcesum->Scale(1./4.);
    } 

  for (is = 0; is < 6; is++){
    
    fitCalo("horizontalslice",is,30000,215500);
    hslce[is]=(TH1D*)htmp->Clone();
    sprintf(hslcename,"hslice%01i",is);
    hslce[is]->SetName(hslcename);
    hslce[is]->SetTitle(hslcename);
    if (run_2) {    
      hslce[is]->Rebin(4);
      hslce[is]->Scale(1./4.);
    } 
  }

  for (is = 0; is < 6; is++){

    double yscle = hfll->Integral( 241, 1478) / hslce[is]->Integral( 241, 1478); // bins for 30000,215500
    printf("is yscale %i %f\n", is, yscle);  
    hslce[is]->Scale( yscle);
    hslce[is]->Divide( hfll);
    
    fnslce[is] = new TF1("fn0","[0]*sin([1]*x+[2])+[3]+[4]*x",30000,215500);
    fnslce[is]->SetParNames("Amplitude","Omega","Phase","Offset","Slope");
    fnslce[is]->SetParameters(0.02,0.00143,0.0,1.0,0.0);
    fnslce[is]->SetNpx(10000);

    if ( !run_2 && is == 1)  fnslce[is]->FixParameter(0,0.015); // helper

    hslce[is]->Fit( fnslce[is], "", "NRE", 30000, 215500);

    omega[is] = fnslce[is]->GetParameter(1);
    domega[is] = fnslce[is]->GetParError(1);
    Offst[is] = fnslce[is]->GetParameter(3);
    dOffst[is] = fnslce[is]->GetParError(3);
    Slope[is] = fnslce[is]->GetParameter(4);
    dSlope[is] = fnslce[is]->GetParError(4);    
  } 

  double meanfshift = 0.0;
  for (is = 1; is < 6; is++){
    
    /*

      mean shift 0.004095 (xtal units), 0.102381 (mm units) per 100 microseconds

      is, N[i-1], N[i] 1, 1.27e+09, 5.35e+09, Slope +/- dSlope 5.81e-08 +/- 6.28e-10, change 1.01e+00, ratio 2.37e-01, fshift -7.62e-03
      is, N[i-1], N[i] 2, 5.35e+09, 8.06e+09, Slope +/- dSlope 1.81e-08 +/- 5.33e-10, change 1.00e+00, ratio 6.64e-01, fshift -5.40e-03
      is, N[i-1], N[i] 3, 8.06e+09, 8.34e+09, Slope +/- dSlope -1.67e-08 +/- 5.26e-10, change 9.98e-01, ratio 9.67e-01, fshift 5.01e-02
      is, N[i-1], N[i] 4, 8.34e+09, 5.84e+09, Slope +/- dSlope -5.71e-08 +/- 6.06e-10, change 9.94e-01, ratio 1.43e+00, fshift -1.34e-02
      is, N[i-1], N[i] 5, 5.84e+09, 1.44e+09, Slope +/- dSlope -9.96e-08 +/- 1.08e-09, change 9.90e-01, ratio 4.07e+00, fshift -3.25e-03

    */

    if (is >= 1) {
      double delayt = 100000.;
      double chnge = 1. + Slope[is]*delayt;
      double rtio = N[0][is-1] / N[0][is];
      double  fshift = ( chnge - 1. ) / ( rtio - 1. );
      meanfshift += fshift;
      printf(" is, N[i-1], N[i] %i, %6.2e, %6.2e, Slope +/- dSlope %6.2e +/- %6.2e, change %6.2e, ratio %6.2e, fshift %6.2e\n", 
	     is, N[0][is-1], N[0][is], Slope[is], dSlope[is], chnge, rtio, fshift);
    }
 
  }
  meanfshift /= 5.;
  printf(" mean shift %f (xtal units), %f (mm units)\n", meanfshift, 25.*meanfshift);

  c1->Clear();
  c1->Divide(2,3);
  for (is = 0; is < 6; is++){
    c1->cd(is+1);
    hslce[is]->SetMaximum(1.2);
    hslce[is]->SetMinimum(0.0);
    hslce[is]->Draw("");
  }

  return;
}
