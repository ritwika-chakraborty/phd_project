//----------------------------------------------------------------
//Copy from  "srcs/gm2midastoart/calo/ReadBankExampleQ_module.cc"
//[MODIFIED] by Ritwika Chakraborty @ May 19th 2022
//
// this version fills per calo Qmethod plots rather than per xtal, per calo
//
//Module Class: Analyzers
//Module Function: Analyze midas banks and generate Q-method histograms
//----------------------------------------------------------------


//Art components
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Art records
#include "gm2dataproducts/daq/MidasBankArtRecord.hh"
#include "gm2util/common/RootManager.hh"

// data products
#include "gm2dataproducts/calo/QHistArtRecord.hh"
#include "gm2dataproducts/daq/CCCArtRecord.hh"
// Added by Fang Han
#include "gm2dataproducts/calo/CaloCalibrationConstant.hh"
#include "gm2dataproducts/calo/QinFillGainFuncArtRecord.hh"

#include "gm2dataproducts/calo/EnergyCalibrationConstantsArtRecord.hh"
#include "gm2dataproducts/calo/OOFConstantsArtRecord.hh"
#include "gm2dataproducts/calo/IFGConstantsArtRecord.hh"
#include "gm2dataproducts/calo/ChannelStatusArtRecord.hh"



// Database service
#include "gm2util/database/Database_service.hh" //db

//local helpers
//#include "gm2midastoart/sources/MidasHelper.hh"

//C/C++ includes
#include <memory>
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include <limits>
#include <boost/range/irange.hpp>

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

namespace gm2calo
{

class QFillByFillAnalyzer : public art::EDAnalyzer
{

public:
  explicit QFillByFillAnalyzer(fhicl::ParameterSet const &pset);

  virtual ~QFillByFillAnalyzer();

  void beginRun(const art::Run &r) override; //db

  void endRun(const art::Run &r) override; //db

  void beginSubRun(art::SubRun const &sr) override;
  //void endSubRun(art::SubRun const &sr) override;

  void analyze(const art::Event &event) override;

private:
  //!root plotting members
  //std::unique_ptr<RootManager> rootManager_;
  //std::string outputRootFileName_;
  //std::unique_ptr<TFile> outputRootFile_;

  std::string processorInstanceLabel_;
  art::Selector qHistSelector_;

  // array for infill gain correction
  std::map<int, std::map<int, TF1>> gainFunctionMap_;

  //FCL params
  std::string inputModuleLabel_;
  std::string inputInstanceName_;
  bool ifOneHist_;
  int CaloNum_;
  int FlshNum_;
  int SoftRB_;
  int loNoise_;
  int hiNoise_;
  int ThresWindow_;
  int ThresSign_;
  double ThresMultiplier_;
  double ThresAbsolute_;
  double PedMaxDev_;
  int HiGap_;
  int LoGap_;

  //book histogram directory
  std::string dirName_;

  //access database and store constants
  //art::ServiceHandle<gm2util::Database> databaseHandle;								  //db
  std::vector<double> ReadCalibConstants(unsigned int runNumber, unsigned int caloNum); //db
  std::vector<std::vector<double>> xtalConstants_;
  std::vector<int> skipIndices_;
  bool useCal = true;          // switch for using calibration constants
  bool useInFillGain = true;   // switch for using infill gain correction
  bool useActiveChannelList_ = true;
  double CalSegAvg = 1.0;      // average vale of calibration constants from all xtals of all calos
  int lastsequenceindex = -99; // checking incrementing of sequence index
  double binsizefactor = 15.0; // 60.0 for run-1, changed to 15.0 for run2 
  double oldbincontent_0 = 0.0;
  double oldbincontent_[5000][54];
  int time_stamp[1000000];
  int its = 0;
  int fill_counter = 0;
  int fill_counter_subrun=0;
  int fill_counter_flush = 0;
  int fill_counter_flush_old = 0;
  int empty_data = 0 ;
  int xtal_counter = 0;
  int calo_counter = 0;
  // double ref_amp;

#define SWRBMAX 4

  //book 1D histograms
  TH1D *Hlasersniff, *Hlasersniffproj, *Hlaser_xtal_counter;
  TH2D *HpedDip;
  TH2D *Hbadfillstats;
  TH2D *HnoisyCaloXtal;
  //TH1D *HAvgNoise;
  TH1D *Hseqnumstats;
  TH1D *HFillCounter;
  TH1D *qHist_1D_[24];
  TH1D *qHist_1D_thr_[24];
  TH1D *qHist_1D_ped_all_[24];
  TH1D *qHist_1D_ped_usd_[24];
  TH1D *qHist_1D_ped_ent_[24];
  double t_lastfill_[4200*4][24*54] = {0.0};
  double t_lastfill_tmp_[4200*4];
  TH1D *qHist_1D_sig_[24];
  TH1D *qHist_1D_sig_xtal_[54];
  TH1D *qHist_1D_sig_sum;
  TH1D *qHist_1D_ydiff_sum;
  //  TH1D *qHist_1D_amp_time[200];
  TH1D *qHist_1D_puA_[24];
  TH1D *qHist_1D_puL_[24];
  TH1D *qHist_1D_puH_[24];
  TH1D *qHist_1D_ESumPerFill_[24];
  TH1D *qHist_1D_NHitPerFill_[24];
  TH1D *qHist_1D_sig_vslice_[9];
  TH1D *qHist_1D_sig_hslice_[6];
  //TH2D *qHist_2D_sig_vslice_[9];
  //TH2D *qHist_2D_sig_hslice_[6];
  TH2D *qHist_2D_puLE1E2_[24];
  TH2D *qHist_2D_puHE1E2_[24];
  TH1D *qHist_1D_sig_seqnum_[16];
  TH2D *qHist_2D_sig_lego_;
  TH2D *qHist_2D_sig_legocalo_[24];
  TH2D *qHist_2D_sig_flash_;
  TH2D *qHist_2D_sig_flashcalo_[24];
  TH1D *qHist_1D_tmp_[54][24];
  TH1D *qHist_1D_rms_[54][24];
  TH1D *qHist_1D_enrgy_p_[24];
  TH1D *qHist_1D_enrgy_e_[24];
  TH1D *qHist_1D_noise_[54][24];
  TH1D *qHist_1D_puL_hits_;
  TH1D *qHist_1D_puH_hits_;
  TH1D *qHist_1D_rej_hits_;
  TH1D *qSubHist_1D_[24][4];
  TH1D *hcalocal_;
  TH1D *hxtalcal_;
  TH1D *hthreshold_;
  TH1D *qHist_1D_time_stamp_[54];

  TH2D *qHist_2D_noise_v_signal_;
  TH2D *qHist_2D_rejectedpoints_;
  TH2D *qHist_2D_acceptedpoints_;
  TH2D *qHist_2D_pointrejection_;
  TH2D *qHist_2D_pointrejection2_;
  TH2D *qHist_2D_sigmaold_;
  TH2D *qHist_2D_sigma_;
  TH2D *qHist_2D_averageold_;
  TH2D *qHist_2D_average_;
  TH2D *qHist_2D_slpe_;
  TH2D *qHist_2D_slpeL_;
  TH2D *qHist_2D_slpeH_;
  TH2D *qHist_2D_ydifflohi_;

  // Doube_t time_stamp[10];
  //  Int_t ts=0;

  //book 2D histograms
  TH2D *qHist_2D_[24];
  TH1D *qHist_1D_ydiff_[24];
  //TH2D *qHist_2D_Coarse_[24];
  //TH2D *qSubHist_2D_[24][4];


  //calib constants from database
  std::string energyConstsInstanceLabel_;
  std::string energyConstsModuleLabel_;
  gm2calo::EnergyCalibrationConstantsArtRecord calibConstants_;

  //OOF constants from database                         
  std::string LTCorrectionInstanceLabel_;
  std::string LTCorrectionModuleLabel_;
  gm2calo::OOFConstantsArtRecord oofConstants_;

  std::string ChannelStatusInstanceLabel_;
  std::string ChannelStatusModuleLabel_;
  gm2calo::ChannelStatusSubrunArtRecord channelStatusDB_;

  //In-Fill Gain Correction from database                   
  std::string inFillConstsInstanceLabel_;
  std::string inFillConstsModuleLabel_;
  gm2calo::IFGConstantsArtRecord IFGConstants_;

  // Switches
  bool UseIFGainFuncs_;
  bool UseADC2EConsts_;
  bool UseOOFCorrects_;
  // Gain Correction
  std::vector<gm2calo::QInFillGainFuncArtRecord> IFGainFuncCol_;
  std::vector<std::vector<double>> ADC2EConsts_;
  std::vector<std::vector<double>> OOFCorrects_;
};

//! standard constructor
QFillByFillAnalyzer::QFillByFillAnalyzer(fhicl::ParameterSet const &pset)
    : art::EDAnalyzer(pset),
      processorInstanceLabel_(pset.get<std::string>("processorInstanceLabel")),
      qHistSelector_(art::ModuleLabelSelector(pset.get<std::string>("unpackerModuleLabel"))),
      inputModuleLabel_(pset.get<std::string>("inputModuleLabel", "CQUnpacker")),
      ifOneHist_(pset.get<bool>("ifOneHist", true)),
      CaloNum_(pset.get<int>("CaloNum", 23)),
      SoftRB_(pset.get<int>("SoftRB", 1)),
      loNoise_(pset.get<int>("loNoise", 3000)),
      hiNoise_(pset.get<int>("hiNoise", 3150)),
      ThresWindow_(pset.get<int>("ThresWindow", 8)),
      ThresSign_(pset.get<int>("ThresSign", +1)),
      ThresMultiplier_(pset.get<double>("ThresMultiplier", +5.0)),
      ThresAbsolute_(pset.get<double>("ThresAbsolute", +20.0)),
      PedMaxDev_(pset.get<double>("PedMaxDev", +1.0)),
      HiGap_(pset.get<int>("HiGap", +1.0)),
      LoGap_(pset.get<int>("LoGap", +1.0)),
      // Gain Correction Read In Configurations
      // Added by Fang Han
      energyConstsInstanceLabel_(pset.get<std::string>("EnergyConstsInstanceLabel", "")),
      energyConstsModuleLabel_(pset.get<std::string>("EnergyConstsModuleLabel", "EnergyCalibrationConstants")),
      LTCorrectionInstanceLabel_(pset.get<std::string>("LTCorrectionInstanceLabel", "")),
      LTCorrectionModuleLabel_(pset.get<std::string>("LTCorrectionModuleLabel", "OOFConstantsService")),
      ChannelStatusInstanceLabel_(pset.get<std::string>("ChannelStatusInstanceLabel", "")),
      ChannelStatusModuleLabel_(pset.get<std::string>("ChannelStatusModuleLabel", "ChannelStatus")),
      inFillConstsInstanceLabel_(pset.get<std::string>("inFillConstsInstanceLabel", "")),
      inFillConstsModuleLabel_(pset.get<std::string>("inFillConstsModuleLabel", "IFGnostdpConstantsService")),

      UseIFGainFuncs_(pset.get<bool>("UseIFGainFuncs", false)),
      UseADC2EConsts_(pset.get<bool>("UseADC2EConsts", false)),
      UseOOFCorrects_(pset.get<bool>("UseOOFCorrects", true))

{
  mf::LogDebug("[QFillByFillAnalyzer]") << "------------------------------------------------\n"
                                        << "CONSTRUCTING MODULE 'QFillByFillAnalyzer'";

  //Initialize all histograms as nullpointers

  Hlasersniff = nullptr;
  Hlasersniffproj = nullptr;
  Hbadfillstats = nullptr;
  Hseqnumstats = nullptr;
  HFillCounter = nullptr;
  HpedDip = nullptr;
  HnoisyCaloXtal=nullptr;
  //HAvgNoise = nullptr;
  Hlaser_xtal_counter = nullptr;
  

  hcalocal_ = nullptr;
  hxtalcal_ = nullptr;
  hthreshold_ = nullptr;

  
   qHist_1D_sig_sum = nullptr;
   qHist_1D_ydiff_sum = nullptr;
   
  qHist_2D_noise_v_signal_ = nullptr;
  qHist_2D_acceptedpoints_ = nullptr;
  qHist_2D_rejectedpoints_ = nullptr;
  qHist_2D_pointrejection_ = nullptr;
  qHist_2D_pointrejection2_ = nullptr;
  qHist_2D_sigmaold_ = nullptr;
  qHist_2D_sigma_ = nullptr;

  
    qHist_2D_averageold_ = nullptr;
    qHist_2D_average_ = nullptr;
    qHist_2D_slpe_ = nullptr;
    qHist_2D_ydifflohi_ = nullptr;
    qHist_2D_slpeL_ = nullptr;
    qHist_2D_slpeH_ = nullptr;
 

  for(int kxtal =0; kxtal < 54; kxtal++)
    {
      qHist_1D_time_stamp_[kxtal] = nullptr;
      oldbincontent_[1][kxtal] = 0.0;
    }

  for (int iSeq = 0; iSeq < 16; iSeq++)
  { // for sequence Index
    qHist_1D_sig_seqnum_[iSeq] = nullptr;
  }

  for (int iSlice = 0; iSlice < 9; iSlice++)
  { // for vertical / horizontal slices
    qHist_1D_sig_vslice_[iSlice] = nullptr;
    //qHist_2D_sig_vslice_[iSlice] = nullptr;
  }
  for (int iSlice = 0; iSlice < 6; iSlice++)
  {
    qHist_1D_sig_hslice_[iSlice] = nullptr;
    //qHist_2D_sig_hslice_[iSlice] = nullptr;
  }

  
  
    qHist_1D_puL_hits_ = nullptr;
    qHist_1D_puH_hits_ = nullptr;
    qHist_1D_rej_hits_ = nullptr;
  

  for (int iCal = 0; iCal < 24; iCal++)
  {
    qHist_2D_puLE1E2_[iCal] = nullptr;
    qHist_2D_puHE1E2_[iCal] = nullptr;
  }
  qHist_2D_sig_lego_ = nullptr;
  for (int iCal = 0; iCal < 24; iCal++)
  {
    qHist_2D_sig_legocalo_[iCal] = nullptr;
  }
  qHist_2D_sig_flash_ = nullptr;
  for (int iCal = 0; iCal < 24; iCal++)
  {
    qHist_2D_sig_flashcalo_[iCal] = nullptr;
  }

  for (int iCal = 0; iCal < 24; iCal++)
  {
      qHist_1D_[iCal] = nullptr;
      qHist_1D_ped_all_[iCal] = nullptr;
      qHist_1D_ped_usd_[iCal] = nullptr;
      qHist_1D_ped_ent_[iCal] = nullptr;
      qHist_1D_sig_[iCal] = nullptr;
      qHist_1D_puA_[iCal] = nullptr;
      qHist_1D_puL_[iCal] = nullptr;
      qHist_1D_puH_[iCal] = nullptr;
      qHist_1D_ESumPerFill_[iCal] = nullptr;
      qHist_1D_NHitPerFill_[iCal] = nullptr;
      if(iCal==0){
      for (int jxtal = 0; jxtal < 54; jxtal++) {
        qHist_1D_sig_xtal_[jxtal] = nullptr; }
     
    }

    for (int jxtal = 0; jxtal < 54; jxtal++) 
     {

      qHist_1D_tmp_[jxtal][iCal] = nullptr;
      qHist_1D_rms_[jxtal][iCal] = nullptr;
      qHist_1D_noise_[jxtal][iCal] = nullptr;
    }

    qHist_1D_enrgy_p_[iCal] = nullptr;
    qHist_1D_enrgy_e_[iCal] = nullptr;
    qHist_2D_[iCal] = nullptr;
    qHist_1D_ydiff_[iCal] = nullptr;
    //qHist_2D_Coarse_[iCal] = nullptr;

    qSubHist_1D_[iCal][0] = nullptr;
    qSubHist_1D_[iCal][1] = nullptr;
    qSubHist_1D_[iCal][2] = nullptr;
    qSubHist_1D_[iCal][3] = nullptr;


  }

  for (int caloIdx = 0; caloIdx < 24; ++caloIdx)
  {
    std::vector<double> thisUnitConsts;
    for (int chanIdx = 0; chanIdx < 54; ++chanIdx)
    {
      auto xtal_str = std::string("xtal") + std::to_string(chanIdx);
      auto f1 = new TF1("Dummy", "1", 0.0, 700.0);
      gm2calo::QInFillGainFuncArtRecord thisGainFunc;
      thisGainFunc.caloIdx = caloIdx + 1;
      thisGainFunc.chanIdx = chanIdx;
      thisGainFunc.func = std::move(*f1);
      //FlatFuncCol_.push_back(thisGainFunc);
      thisUnitConsts.push_back(1.0);
    } // End for one Calo
    //FlatConsts_.push_back(thisUnitConsts);
  } // End of setting all gain corrections to be 1
}

//! standard destructor
QFillByFillAnalyzer::~QFillByFillAnalyzer() {}

void QFillByFillAnalyzer::beginSubRun(art::SubRun const &subrun)
{
  mf::LogInfo() << "Enter QFillByFillAnalyzer::beginSubRun()";
  // If use out of fill gain correction constants
    if (UseOOFCorrects_)
  {
       art::Handle<gm2calo::OOFConstantsArtRecord> oofConstantsArtHandle;
    int returnvalue = subrun.getByLabel(LTCorrectionModuleLabel_,
					LTCorrectionInstanceLabel_,
					oofConstantsArtHandle);
    if (returnvalue == 1){
      oofConstants_ = *oofConstantsArtHandle;
    }
    else std::cout<<"gm2analyses::QHistProducer no oof corrections found!"<<std::endl;

  }   // End of UseOOFCorrects
    if(useActiveChannelList_){
      mf::LogDebug("IslandTemplateFit") << "Getting channel statuses..." << std::endl;
      //auto const tag = art::InputTag("ChannelStatus");
      // channelStatusDB_ = *event.getRun().getValidHandle<gm2calo::ChannelStatusArtRecord>(tag);       // old, before v9_65
      //channelStatusDB_ = subrun.getValidHandle<gm2calo::ChannelStatusSubrunArtRecord>(tag); // new, by subrun
      art::Handle<gm2calo::ChannelStatusSubrunArtRecord> ChannelStatusSubrunArtHandle;
      int returnvalue = subrun.getByLabel(ChannelStatusModuleLabel_,
                                          ChannelStatusInstanceLabel_,
                                          ChannelStatusSubrunArtHandle);
      if (returnvalue == 1){
	channelStatusDB_ = *ChannelStatusSubrunArtHandle;
      }
      for(unsigned int cIdx = 0; cIdx < channelStatusDB_.caloStatus.size(); cIdx++){
	if (std::find(skipIndices_.begin(), skipIndices_.end(),	 cIdx + 1) == skipIndices_.end() && !channelStatusDB_.caloStatus.at(cIdx)) {
          mf::LogDebug("IslandTemplateFit") << "Marking calo with index " << cIdx << " as bad!" << std::endl;
	  skipIndices_.push_back(cIdx + 1);
	}
      }
    }


  mf::LogInfo() << "Exit QFillByFillAnalyzer::beginSubRun()";
} // End of BeginSubRun()

void QFillByFillAnalyzer::endRun(const art::Run &r)
{
  printf("********* start endRun  ************\n");


  printf("********* end endRun ************\n");
}

void QFillByFillAnalyzer::beginRun(const art::Run &r)
{

  printf("********* start beginRun  ************\n");
  
  //Get calibration constants
  if (useCal)                                                                           
    {                                                                        
     art::Handle<gm2calo::EnergyCalibrationConstantsArtRecord> calibConstantsArtHandle;
     int returnvalue = r.getByLabel(energyConstsModuleLabel_,
                                   energyConstsInstanceLabel_,
                                   calibConstantsArtHandle);
     if (returnvalue == 1){
     calibConstants_ = *calibConstantsArtHandle;
     }
     else std::cout<<"gm2analyses::QHistProducer no energy calibration constants found!"<<std::endl;
    }

  //In-Fill Gain (IFG) correction
  if (useInFillGain)
    {
     art::Handle<gm2calo::IFGConstantsArtRecord> inFillConstantsArtHandle;
     int returnvalue = r.getByLabel(inFillConstsModuleLabel_,
				      inFillConstsInstanceLabel_,
				      inFillConstantsArtHandle);
     if (returnvalue == 1){
        IFGConstants_ = *inFillConstantsArtHandle;
     }
    }
 

  //Read in Database constants for given run
  if (hcalocal_ == nullptr)
  {
    art::ServiceHandle<art::TFileService> TFileServ;

    Hbadfillstats = TFileServ->make<TH2D>(
        Form("Hbadfillstats"),
        Form("Hbadfillstats"),
        54, 0., 54.,
        24, 0., 24.);

    Hseqnumstats = TFileServ->make<TH1D>(
        Form("Hseqnumstats"),
        Form("Hseqnumstats"),
        16, 0., 16.);
    
    HFillCounter = TFileServ->make<TH1D>(
	 Form("HFillCounter"),
	 Form("HFillCounter"),
	 1, 0., 1.);


    hcalocal_ = TFileServ->make<TH1D>(
        Form("hcalocal_"),
        Form("hcalocal_"),
        24, 0, 24);

    hxtalcal_ = TFileServ->make<TH1D>(
        Form("hxtalcal_"),
        Form("hxtalcal_"),
        54 * 24, 0, 54 * 24);

    hthreshold_ = TFileServ->make<TH1D>(
        Form("hthreshold_"),
        Form("hthreshold_"),
        500, 0., 23.63*100.);

    // only fill array and histos on first call

    printf(" fill, entries, integral bin 100, %i, not valid first sub-run\n", FlshNum_);
  }
  else
  {
    printf(" fill, entries, integral bin 100, %i\n", FlshNum_);
  }
                               

  printf("********* end beginRun ************\n");
}

//! analyzes data
void QFillByFillAnalyzer::analyze(const art::Event &event)
{

  printf("start Q-method analyzer Production module\n"); // create timestamps
     
  struct timeval t_start, t_end;
  gettimeofday(&t_start, NULL);

  // gaussian function for noise fitting
  //TF1 *fg;
  //fg = new TF1( "fg", "[0]*TMath::Gaus(x,[1],[2])", -10., 10.);

  mf::LogDebug("[QFillByFillAnalyzer]") << "START Tim's QFillByFillAnalyzer::analyze(" << event.id() << ")";
  art::ServiceHandle<art::TFileService> TFileServ;

  std::vector<art::Handle<gm2calo::CaloQHistArtRecord>> qHistHandleVec;
  event.getMany(qHistSelector_, qHistHandleVec);

  // create a vector of ptrs to the NON EMPTY QHistArtRecords
  std::vector<const CaloQHistArtRecord *> nonemptyQData;

  // get sequence index from fc7
  const auto &encodercol =
      //*event.getValidHandle<gm2ccc::EncoderFC7ArtRecordCollection>({fc7ModuleLabel_, fc7InstanceLabel_});
      *event.getValidHandle<gm2ccc::EncoderFC7ArtRecordCollection>({"cccUnpacker", "unpacker"});
  int sequenceIndex = encodercol.size() ? encodercol.front().sequenceIndex : -1;

  if (qHistHandleVec.size() > 0)
    std::cout << "QAnalyzer: CQ event.size() = " << event.size() << " : handle.size() = " << qHistHandleVec.size() << std::endl;




  

  // test finding calo data
  nonemptyQData.reserve(qHistHandleVec.size());


  // sniff laser calibration fills
  bool laserON = false;
  bool pedDip= false;

  TF1 *IFGfunc = new TF1("IFGfunc",IFGConstants_.fcnstring.c_str(),0.0,500.0);//define IFG function 
  //int jj = 0;
  //int nq = 0;

  for (const auto &handle : qHistHandleVec)
  {
    const auto &qhand = *handle;
    auto iCalo = qhand.caloNum;
    

     

    //       printf("__________________ ######   1st icalo definition  ######  ____ fill_counter, iCalo %i, %i \n ",fill_counter, iCalo);

    if (iCalo >= 1 && iCalo <= 24)
    { // search all calos for laser pulses (just in case a callo is turned off)

      if (!qhand.qHists.empty())
      {
        nonemptyQData.push_back(&qhand);

        if (nonemptyQData.size() > 0)
        {
          int nBins = nonemptyQData[0]->length;
          
            if(nBins<16800)
	    {
              printf("Incorrect sample length");
	      return;
	    }

          if (iCalo == 1)
          { // put here to print sequence index once for each muon type fill
            printf("sequence index %i\n", sequenceIndex);
            if (lastsequenceindex < 15 && (lastsequenceindex + 1) != sequenceIndex)
              printf("last, current seq index %i, %i\n", lastsequenceindex, sequenceIndex); // sequence increment not sensible
            if (lastsequenceindex == 15 && sequenceIndex != 0)
              printf("last, current seq index %i, %i\n", lastsequenceindex, sequenceIndex); // sequence increment not sensible
            lastsequenceindex = sequenceIndex;
          }

          if (Hlasersniff == nullptr)
          {
            Hlasersniff = TFileServ->make<TH1D>(
                Form("Hlasersniff"),
                Form("Hlasersniff"),
                nBins, 0, nBins);
           
            Hlasersniffproj = TFileServ->make<TH1D>(
                Form("Hlasersniffproj"),
                Form("Hlasersniffproj"),
                10000, 0, 10000);
          }
          Hlasersniff->Reset();

	  if(HpedDip == nullptr)
	    {
	      HpedDip = TFileServ->make<TH2D>(
						  Form("HpedDip"),
						  Form("HpedDip"),
						  500, 0, 500, 500, 0, 500);
            }

	  if(HnoisyCaloXtal == nullptr)
            {
              HnoisyCaloXtal = TFileServ->make<TH2D>(
					      Form("HnoisyCaloXtal"),
					      Form("HnoisyCaloXtal"),
					      54, 0, 54, 24, 0, 24);
            }



                  
	  if(Hlaser_xtal_counter == nullptr)
            {
              Hlaser_xtal_counter = TFileServ->make<TH1D>(
							  Form("Hlaser_xtal_counter"),
							  Form("Hlaser_xtal_counter"),
							  nBins, 0, nBins);
            }
          Hlaser_xtal_counter->Reset();



          for (auto art_qHist : nonemptyQData[0]->qHists)
          { // start xtal loop

            int iADC = 0;
            for (double ADC_Value : art_qHist.trace)
            { // start sample loop
              ADC_Value = ADC_Value / binsizefactor;
              iADC++;
              Hlasersniff->Fill(iADC, ADC_Value);
	      if( iADC > 400 && ADC_Value*binsizefactor > 1000 )
                {
		  Hlaser_xtal_counter->Fill(iADC);
                }
            } // end sample loop
          }   // end segment loop
        }     // end empty qHist check
        

	for(int iADC = 1; iADC <= Hlaser_xtal_counter->GetNbinsX(); iADC++)
          {
            if(Hlaser_xtal_counter->GetBinContent(iADC)>=50)
              {
                xtal_counter=1;
		std::cout<<"pulse found in "<<Hlaser_xtal_counter->GetBinContent(iADC)<<" xtals;\
 LASER DETECTED!! at calo "<<iCalo<<endl;
              }
          }


	// int minib = 100; // cut region of beam injection to avoid false lasers
        int minib = 400; // cut the region of beam injection for run 2
        for (int ib = minib; ib <= Hlasersniff->GetNbinsX() - 1; ib++)
        {
          Hlasersniffproj->Fill(Hlasersniff->GetBinContent(ib));
	  // if (Hlasersniff->GetBinContent(ib) >= 2500.)
	  if (Hlasersniff->GetBinContent(ib)*binsizefactor >= 1000. && xtal_counter == 1)// factor of 4 comes from run-2 time binning.
          {
            laserON = true;
	    // printf("Laser On in iCalo, fill_counter, %i, %i \n", fill_counter, iCalo);
          } // set laser on
        }   // loop over bins
        xtal_counter=0;
	int minic = 1920;
        for (int ic = minic; ic <= Hlasersniff->GetNbinsX() - 1; ic++)
	  {
	    
	    if (Hlasersniff->GetBinContent(ic) <= -3000. || Hlasersniff->GetBinContent(ic) >= 10000)
              {
		pedDip = true;
                HpedDip->Fill(event.event(),event.subRun());
	      }
           }   // loop over bins                                                                                                 

        nonemptyQData.clear();

      } // calo data empty
    }   // calo 1-24
  }     // calo data handle loop
 
  fill_counter++;
  fill_counter_subrun++;
  printf("Fill number is %i, %i \n", fill_counter, fill_counter_subrun);  

  printf("laser on = %i\n", laserON);
  
if (laserON)
  {
    std::cout<<"Laser is on!!"<<endl;
    return;
  }

  if(pedDip)
    {
      printf("Negative Pulse found!"); //these sould have been rejected by QSpikeCheckFilter_module.cc
    //return;
    }

  std::cout<<" Run number is "<<  event.run() <<" Subrun number is "<< event.subRun() << " Event number is "<< event.event() << std::endl;
 


  printf(" histogram sequence number %i\n", sequenceIndex);
  Hseqnumstats->Fill(sequenceIndex); // store sequence index if laser off

  // Int_t iLaser;
  calo_counter=0;

  for (const auto &handle : qHistHandleVec)
  {
     
    const auto &qhand = *handle;
    auto iCalo = qhand.caloNum;
    if( iCalo<1 || iCalo>24)
      {
	empty_data = 1;
      }
    else
      {
	empty_data = 0;
      }

    //  printf(" ######    Before !qhand.qHists.empty() check: fill_counter, iCalo %i, %i \n ",fill_counter, iCalo);
   
   if (!qhand.qHists.empty())
    {
 
      //   printf(" ######    After !qhand.qHists.empty() check: fill_counter, iCalo %i, %i \n ",fill_counter, iCalo);

      nonemptyQData.push_back(&qhand);
      
      if (qhand.caloNum >= 1 && qhand.caloNum <= 24)
      {
	
        
        if (nonemptyQData.size() > 0)
        {
          
          fill_counter_flush_old=fill_counter_flush;
          fill_counter_flush=fill_counter;
          printf("_________________________check here_______________________calo counter, icalo, fill_counter, fill_counter_flush, fill_counter_flush_old %i, %i, %i, %i, %i \n", iCalo, calo_counter, fill_counter, fill_counter_flush_old, fill_counter_flush);

          if( fill_counter_flush-fill_counter_flush_old < 4 && calo_counter==0 )
	    {
	       if(fill_counter_subrun==1 || fill_counter_subrun==2 || fill_counter_subrun==3)
		{
		  fill_counter_flush=fill_counter;
		}
	       else
		{
	          fill_counter_flush=fill_counter_flush_old;
		  std::cout<<"Spurious fill encountered"<<endl; //These should have been rejected by QFlushQualityFilter_module.cc
	          //return;
		}
	    }
	    calo_counter++;

          int nBins = nonemptyQData[0]->length;
          int iCal = nonemptyQData[0]->caloNum;
          int nSegs = nonemptyQData[0]->qHists.size();
          int nIntv = nonemptyQData[0]->intervalNum;
          int multplr = nonemptyQData[0]->multiplier;
          int rebin_0 = nonemptyQData[0]->initRebinConst;
          int first_bin_ = nonemptyQData[0]->firstSampleNum;
          int last_bin_ = nonemptyQData[0]->lastSampleNum;
          
	 
          
          iCal += 0; // avoid compiler complaints if not printing below
          nSegs += 0;
          nIntv += 0;
          multplr += 0;

          // fill calibration constant array for given calo
          int countSeg = 0;  // xtal counter gor given calo
          double CalSeg[54]; // calibration constant array
          if (useCal)
          {
 
	    //    for (int caloIdx = 0; caloIdx < 24; ++caloIdx)
	    //  {
		for (int chanIdx = 0; chanIdx < 54; ++chanIdx)
		  {
                    if((oofConstants_.oofcorrection.at((iCalo-1)*54 + chanIdx)) != 0 )
		      {
			//	std::cout<<calibConstants_.energyCalibrationConstant.at((iCalo-1)*54 + chanIdx)<<endl;
		       CalSeg[countSeg] = (rebin_0) *  (calibConstants_.energyCalibrationConstant.at((iCalo-1)*54 + chanIdx)) / (oofConstants_.oofcorrection.at((iCalo-1)*54 + chanIdx));
		      }//CalSeg[countSeg] = oofConstants_.oofcorrection.at((iCalo-1)*54 + chanIdx);
		    else
		      {
			CalSeg[countSeg] = (rebin_0) *  (calibConstants_.energyCalibrationConstant.at((iCalo-1)*54 + chanIdx));
		      }
		    //std::cout<<countSeg<<" "<<CalSeg[countSeg]<<std::endl;
                    countSeg++;
		  }
		//  }

          }
     

          // fix histogram binning so exaxtcy devisible by max rebinning factor
          int nBinsNew = nBins, iBinSize = 1;
       	  // nBinsNew = 3200; // quick fix different odb pars in run 2 FIX ME!!!!!!!!!!!!
          nBinsNew = 4200*4;
          iBinSize = (last_bin_ + 1 - first_bin_) / nBins;
          last_bin_ += (nBinsNew - nBins) * iBinSize;
          //printf("after: nBins, nBinsNew, iBinSize, first_bin_, last_bin_ %i, %i, %i, %i %i\n", nBins, nBinsNew, iBinSize, last_bin_, last_bin_);

          // start xtal loop
          for (auto art_qHist : nonemptyQData[0]->qHists)
          {
            auto iSeg = art_qHist.xtalNum;
	    
            unsigned int entriesN;
            unsigned int entriesNped;
            unsigned int entriesNrolling;

           
              entriesN = 0;
              entriesNped = 0;
              entriesNrolling = 0;
            
            

            unsigned int entriesNrollingvslice[9], entriesNrollinghslice[6];
            for (int islice = 0; islice < 9; islice++)
            {
              entriesNrollingvslice[islice] = 0;
            }
            for (int islice = 0; islice < 6; islice++)
            {
              entriesNrollinghslice[islice] = 0;
            }

            int iADC; // ADC sample index

            switch (ifOneHist_) // switch for single / staggered time decimation
            {
	       
            case true: // single time decimation
            {
	      
              // define 1D histograms if undefined
              if (qHist_1D_[iCalo - 1] == nullptr)
              {
                FlshNum_ = 0; // fill counter reset
                              //qHist_1D_thr_[iCalo-1] = TFileServ->make<TH1D>(
                              //    Form("qHist1D_thr_%i_%i", iCalo, iCalo),
                              //    Form("qHist1D_thr_%i_%i", iCalo, iCalo),
                              //    nBins,
                              //    edges);

                if (qHist_1D_sig_seqnum_[0] == nullptr)
                {
	      
                  for (int iSeq = 0; iSeq < 8; iSeq++)
                  {
                    //printf("define qHist1D_sig_seqnum_%i", iSeq);
                    qHist_1D_sig_seqnum_[iSeq] = TFileServ->make<TH1D>(
                        Form("qHist1D_sig_seqnum_%i", iSeq),
                        Form("qHist1D_sig_seqnum_%i", iSeq),
                        nBinsNew, first_bin_, last_bin_ + 1);
                  }

                  qHist_2D_noise_v_signal_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_noise_v_signal_"),
                      Form("qHist_2D_noise_v_signal_"),
                      150, -50., 100.,
                      150, -50., 100.);

                  qHist_2D_acceptedpoints_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_acceptedpoints_"),
                      Form("qHist_2D_acceptedpoints_"),
                      3168 / 32, 0.0, 3168.,
                      55, -100., 1100.);

                  qHist_2D_rejectedpoints_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_rejectedpoints_"),
                      Form("qHist_2D_rejectedpoints_"),
                      3168 / 32, 0.0, 3168.,
                      55, -100., 1100.);

                  qHist_2D_pointrejection_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_pointrejection_"),
                      Form("qHist_2D_pointrejection_"),
                      3 * 32, 0., 3 * 32.,
                      550, -100., 1100.);

                  qHist_2D_pointrejection2_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_pointrejection2_"),
                      Form("qHist_2D_pointrejection2_"),
                      3 * 32, 0., 3 * 32.,
                      275, -100., 1100.);

                  qHist_2D_sigmaold_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_deltasigmaold_"),
                      Form("qHist_2D_deltasigmaold_"),
                      //							    nBinsNew, first_bin_, last_bin_+1,
                      3168 / 32, 0.0, 3168.,
                      40, -0.2, 0.2);
                  qHist_2D_sigma_ = TFileServ->make<TH2D>(
                      Form("qHist_2D_sigma_"),
                      Form("qHist_2D_sigma_"),
                      //							    nBinsNew, first_bin_, last_bin_+1,
                      3168 / 32, 0.0, 3168.,
                      25, 0., 5.);

                 
                  
                    qHist_2D_averageold_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_deltaaverageold_"),
                        Form("qHist_2D_deltaaverageold_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        500, -25., 25.);
                    qHist_2D_average_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_average_"),
                        Form("qHist_2D_average_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        500, -25., 25.);
                    qHist_2D_slpe_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_slpe_"),
                        Form("qHist_2D_slpe_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        500, -5., 5.);
                    qHist_2D_ydifflohi_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_ydifflohi_"),
                        Form("qHist_2D_ydifflohi_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        400, -5., 5.);
                    qHist_2D_slpeL_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_slpeL_"),
                        Form("qHist_2D_slpeL_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        500, -5., 5.);
                    qHist_2D_slpeH_ = TFileServ->make<TH2D>(
                        Form("qHist_2D_slpeH_"),
                        Form("qHist_2D_slpeH_"),
                        //							    nBinsNew, first_bin_, last_bin_+1,
                        3168 / 32, 0.0, 3168.,
                        500, -5., 5.);
                  
                }

                if (qHist_1D_puL_hits_ == nullptr)
                {
                    qHist_1D_puL_hits_ = TFileServ->make<TH1D>(
                        Form("qHist_1D_puL_hits_"),
                        Form("qHist_1D_puL_hits_"),
                        3168, 0.0, 3168.);

                    qHist_1D_puH_hits_ = TFileServ->make<TH1D>(
                        Form("qHist_1D_puH_hits_"),
                        Form("qHist_1D_puH_hits_"),
                        3168, 0.0, 3168.);

                    qHist_1D_rej_hits_ = TFileServ->make<TH1D>(
                        Form("qHist_1D_rej_hits_"),
                        Form("qHist_1D_rej_hits_"),
                        3168, 0.0, 3168.); 
                }

                if (qHist_1D_sig_vslice_[0] == nullptr)
                {
                  for (int iSlice = 0; iSlice < 9; iSlice++)
                  {
                    qHist_1D_sig_vslice_[iSlice] = TFileServ->make<TH1D>(
                        Form("qHist1D_sig_vslice_%i", iSlice),
                        Form("qHist1D_sig_vslice_%i", iSlice),
                        nBinsNew, first_bin_, last_bin_ + 1);

		    /* qHist_2D_sig_vslice_[iSlice] = TFileServ->make<TH2D>(
                        Form("qHist2D_sig_vslice_%i", iSlice),
                        Form("qHist2D_sig_vslice_%i", iSlice),
                        nBinsNew, first_bin_, last_bin_ + 1,
                        24, 0, 24 + 1);*/
                  }
                }

                if (qHist_1D_sig_hslice_[0] == nullptr)
                {
                  for (int iSlice = 0; iSlice < 6; iSlice++)
                  {
                    qHist_1D_sig_hslice_[iSlice] = TFileServ->make<TH1D>(
                        Form("qHist1D_sig_hslice_%i", iSlice),
                        Form("qHist1D_sig_hslice_%i", iSlice),
                        nBinsNew, first_bin_, last_bin_ + 1);

                    /*qHist_2D_sig_hslice_[iSlice] = TFileServ->make<TH2D>(
                        Form("qHist2D_sig_hslice_%i", iSlice),
                        Form("qHist2D_sig_hslice_%i", iSlice),
                        nBinsNew, first_bin_, last_bin_ + 1,
                        24, 0, 24 + 1);*/
                  }
                }

                if (qHist_2D_sig_lego_ == nullptr)
                {
                  qHist_2D_sig_lego_ = TFileServ->make<TH2D>(
                      Form("qHist2D_sig_lego"),
                      Form("qHist2D_sig_lego"),
                      9, 0., 9.,
                      6, 0., 6.);
                }

                if (qHist_2D_sig_flash_ == nullptr)
                {
                  qHist_2D_sig_flash_ = TFileServ->make<TH2D>(
                      Form("qHist2D_sig_flash"),
                      Form("qHist2D_sig_flash"),
                      9, 0., 9.,
                      6, 0., 6.);
                }



		/*	if (time_stamp == nullptr)
		  {
		    time_stamp = TFileServ->make<TH1D>(
					  Form("time_stamp"),
					  Form("time_stamp"),
					  2000, 0, 2000);

		  }
 
		*/
                if (qHist_2D_puHE1E2_[0] == nullptr)
                {
                  for (int iCal = 0; iCal < 24; iCal++)
                  {
                    qHist_2D_puLE1E2_[iCal] = TFileServ->make<TH2D>(
                        Form("qHist2D_puLE1E2_%i", iCal + 1),
                        Form("qHist2D_puLE1E2_%i", iCal + 1),
                        100, -23.63*10., 23.63*110.,
                        100, -23.63*10., 23.63*110.);
                    qHist_2D_puHE1E2_[iCal] = TFileServ->make<TH2D>(
                        Form("qHist2D_puHE1E2_%i", iCal + 1),
                        Form("qHist2D_puHE1E2_%i", iCal + 1),
                        100, -23.63*10., 23.63*110.,
                        100, -23.63*10., 23.63*110.);
                  }
                }

                if (qHist_2D_sig_legocalo_[0] == nullptr)
                {
                  for (int iCal = 0; iCal < 24; iCal++)
                  {
                    qHist_2D_sig_legocalo_[iCal] = TFileServ->make<TH2D>(
                        Form("qHist2D_sig_legocalo_%i", iCal + 1),
                        Form("qHist2D_sig_legocalo_%i", iCal + 1),
                        9, 0., 9.,
                        6, 0., 6.);
                  }
                }

                if (qHist_2D_sig_flashcalo_[0] == nullptr)
                {
                  for (int iCal = 0; iCal < 24; iCal++)
                  {
                    qHist_2D_sig_flashcalo_[iCal] = TFileServ->make<TH2D>(
                        Form("qHist2D_sig_flashcalo_%i", iCal + 1),
                        Form("qHist2D_sig_flashcalo_%i", iCal + 1),
                        9, 0., 9.,
                        6, 0., 6.);
                  }
                }

                
                  qHist_1D_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_%i", iCalo),
                      Form("qHist1D_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);

                  qHist_1D_ped_all_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_ped_all_%i", iCalo),
                      Form("qHist1D_ped_all_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_ped_usd_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_ped_usd_%i", iCalo),
                      Form("qHist1D_ped_usd_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_ped_ent_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_ped_ent_%i", iCalo),
                      Form("qHist1D_ped_ent_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_sig_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_sig_%i", iCalo),
                      Form("qHist1D_sig_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);

		    if(iCalo-1==0){
		  for (int jxtal = 0; jxtal < 54; jxtal++){ 
		    qHist_1D_sig_xtal_[jxtal] = TFileServ->make<TH1D>(
			     Form("qHist1D_sig_xtal_%i", jxtal),
			     Form("qHist1D_sig_xtal_%i", jxtal),
			     nBinsNew, first_bin_, last_bin_ + 1);
		  //		  printf("chrystal number is %i\n",iSeg);

		              }
                         }
               

                  qHist_1D_puA_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_puA_%i", iCalo),
                      Form("qHist1D_puA_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_puL_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_puL_%i", iCalo),
                      Form("qHist1D_puL_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_puH_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_puH_%i", iCalo),
                      Form("qHist1D_puH_%i", iCalo),
                      //nBins,
                      //edges);
                      nBinsNew, first_bin_, last_bin_ + 1);

                  qHist_1D_ESumPerFill_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_ESumPerFill_%i", iCalo),
                      Form("qHist1D_ESumPerFill_%i", iCalo),
                      nBinsNew, first_bin_, last_bin_ + 1);
                  qHist_1D_NHitPerFill_[iCalo - 1] = TFileServ->make<TH1D>(
                      Form("qHist1D_NHitPerFill_%i", iCalo),
                      Form("qHist1D_NHitPerFill_%i", iCalo),
                      nBinsNew, first_bin_, last_bin_ + 1);
                


	
	       for (int jxtal = 0; jxtal < 54; jxtal++)
               {
                qHist_1D_rms_[jxtal][iCalo - 1] = TFileServ->make<TH1D>(
	        	Form("qHist1D_rms_%i_%i", jxtal, iCalo),
			Form("qHist1D_rms_%i_%i", jxtal, iCalo),
                        400, -50., +50.);
                qHist_1D_tmp_[jxtal][iCalo - 1] = TFileServ->make<TH1D>(
			Form("qHist1D_tmp_%i_%i", jxtal, iCalo),
			Form("qHist1D_tmp_%i_%i", jxtal, iCalo),
                        500, -50., +50.);
                qHist_1D_noise_[jxtal][iCalo - 1] = TFileServ->make<TH1D>(
			Form("qHist1D_noise_%i_%i", jxtal, iCalo),
			Form("qHist1D_noise_%i_%i", jxtal, iCalo),
                        160, 0.0, 40.0);
	       }

                qHist_1D_enrgy_p_[iCalo - 1] = TFileServ->make<TH1D>(
			Form("qHist1D_enrgy_p_%i",  iCalo),
			Form("qHist1D_enrgy_p_%i",  iCalo),
                        110, -23.63*40.0, 23.63*400.0);
                qHist_1D_enrgy_e_[iCalo - 1] = TFileServ->make<TH1D>(
			Form("qHist1D_enrgy_e_%i",  iCalo),
			Form("qHist1D_enrgy_e_%i",  iCalo),
                        110, -23.63*40.0, 23.63*400.0);
	       
                // time decimation rebinning for pileup studies
                int irb = 1;
               
                  qHist_1D_[iCalo - 1]->Rebin(irb);
                  qHist_1D_ped_all_[iCalo - 1]->Rebin(irb);
                  qHist_1D_ped_usd_[iCalo - 1]->Rebin(irb);
                  qHist_1D_ped_ent_[iCalo - 1]->Rebin(irb);
                  qHist_1D_sig_[iCalo - 1]->Rebin(irb);
                  qHist_1D_puA_[iCalo - 1]->Rebin(irb);
                  qHist_1D_puL_[iCalo - 1]->Rebin(irb);
                  qHist_1D_puH_[iCalo - 1]->Rebin(irb);
                  qHist_1D_ESumPerFill_[iCalo - 1]->Rebin(irb);
                  qHist_1D_NHitPerFill_[iCalo - 1]->Rebin(irb);
		  if(iCalo-1==16){
		  for (int jxtal = 0; jxtal < 54; jxtal++) {                                                                                                                                                    
		    qHist_1D_sig_xtal_[jxtal]->Rebin(irb);                                                                                                                                      
		   }    
		  }
		  
                 
                
                       
                
                // define 2D histograms if undefined
                if (qHist_2D_[iCalo - 1] == nullptr)
                {
                  		    
		    double ymin = -100*23.64;
              double ymax = +250.*23.64;
              int nBinsY = 5*(ymax - ymin);
              
	   			   	      

	      printf("_________________________check here_______________________ %i, %f, %f, %i, %i, %i, %i, %i \n", iCalo, ymin, ymax, nBinsY, nBins, first_bin_, last_bin_+1, nBinsNew); 
		  qHist_1D_ydiff_[iCalo-1] = TFileServ->make<TH1D>(Form("qHist1D_ydiff_%i", iCalo), Form("qHist1D_ydiff_%i", iCalo), 8077, -100*23.64, 4000*23.64);
		 
		  // qHist_1D_ydiff_[iCalo-1]->SetStats(kFALSE);      
	  
                } //End of one-time Initializing 2D

                //printf("finish defining histograms\n");

              } //End of one-time Initializing Histograms.

              if (iCalo == 1 && iSeg == 0)
		{
                 FlshNum_++; // increment flush counter on first calo, first xtal
	         HFillCounter->Fill(0.5); // store sequence index if laser off 
           	}
              // for each xtal calculate the noise RMS in region before flash in order to set threshold
              //int iNoise = 0, lobin = 1500, hibin = 3000; // lobin, hibin are range for noise estimate, June 17 commissioning
              //int iNoise = 0, lobin = 3000, hibin = 3150; // lobin, hibin are range for noise estimate, fall 17 commissioning
              int iNoise = 0, lobinNoise = loNoise_, hibinNoise = hiNoise_; // lobinNoise, hibin are range for noise estimate
              double sumNoise = 0.0, averageNoise = 0.0, noise_ADC_range = 10.0;

              // calculate average value in noise window to help with window in computing the rms
              for (double ADC_Value : art_qHist.trace)
              {
                ADC_Value = ADC_Value / binsizefactor;
                iNoise++;
                if (iNoise < lobinNoise || iNoise > hibinNoise)
                  continue; // noise time window
                sumNoise += ADC_Value;
              }
	      //std::cout<<"inoise is "<<iNoise<<endl;
              averageNoise = sumNoise / iNoise;
              iNoise = 0;

              for (double ADC_Value : art_qHist.trace)
              {
                ADC_Value = ADC_Value / binsizefactor;
                iNoise++;
                if (iNoise < lobinNoise || iNoise > hibinNoise)
                  continue; // noise time window
                if (ADC_Value - averageNoise < -noise_ADC_range || ADC_Value - averageNoise > +noise_ADC_range)
                  continue;                                               // noise ADC window
                qHist_1D_tmp_[iSeg][iCalo - 1]->Fill(ADC_Value - averageNoise); // stores noise relative to averageNoise value in range lobinNoise, hibinNoise
                qHist_1D_rms_[iSeg][iCalo - 1]->Fill(ADC_Value - averageNoise); // stores sample values in range lobinNoise, hibinNoise
		//	std::cout<<"adc-averagenoise,iseg,icalo-1"<<ADC_Value-averageNoise<<", "<<iSeg<<", "<<iCalo-1<<endl;
              }
             
	      // double noiseRMS = qHist_1D_tmp_[iCalo - 1]->GetRMS();
              double noiseEntries = qHist_1D_tmp_[iSeg][iCalo - 1]->GetEntries();
              double noiseIntegral = qHist_1D_tmp_[iSeg][iCalo - 1]->Integral();
              double noiseRMS[54];

              noiseRMS[iSeg] = qHist_1D_tmp_[ iSeg][iCalo - 1]->GetStdDev();

	      //std::cout<<"uncalibrated noiserms in xtal, calo"<<noiseRMS[iSeg]<<", "<<iSeg<<", "<<iCalo-1<<endl;             

              noiseRMS[iSeg] *= CalSeg[iSeg];
               
	      //std::cout<<"calibrated noiserms in xtal, calo"<<noiseRMS[iSeg]<<", "<<iSeg<<", "<<iCalo-1<<endl;
              //HAvgNoise->Fill( (iCalo-1)*54 + iSeg, noiseRMS[iSeg] );

              if (noiseRMS[iSeg] == 0.0)
              {
                noiseRMS[iSeg] = -999.; //store in underflow
                printf("NO VALID DATA, Flush %i, ISeg, ICalo %i, %i, averageNoise, noiseRMS, noiseEntries, noiseIntegral %f, %f, %f, %f. lobin, hibin ADC %f, %f\n",
                       FlshNum_, iSeg, iCalo, averageNoise, noiseRMS[iSeg], noiseEntries, noiseIntegral, art_qHist.trace.at(lobinNoise) / binsizefactor, art_qHist.trace.at(hibinNoise) / binsizefactor);
                Hbadfillstats->Fill(iSeg - 1, iCalo - 1);
                qHist_1D_noise_[iSeg][iCalo - 1]->Fill(noiseRMS[iSeg]); // stores noise RMS from lobinNoise, hibinNoise
                continue;                                   // no data in CQ bank? happens with last bank?
              }
              qHist_1D_noise_[iSeg][iCalo - 1]->Fill(noiseRMS[iSeg]); // stores noise RMS from lobinNoise, hibinNoise
              qHist_1D_tmp_[iSeg][iCalo - 1]->Reset();          // temporary store of calo noise for rms calculation

              if (iSeg == 0)
              {
                  qHist_1D_ESumPerFill_[iCalo - 1]->Reset(); // fill-by-fill xtal-sum, energy distrbution
                  qHist_1D_NHitPerFill_[iCalo - 1]->Reset(); // fill-by-fill xtal-sum, hit distribution
              } // if iSeg = 0
              
              //get the parameters for the IFG function              
	      for(int par = 0; par < IFGConstants_.nparams; par++){
		IFGfunc->SetParameter(par,IFGConstants_.params[par].at((iCalo-1)*54 + iSeg));
	      }
              int flashbin = std::max_element(art_qHist.trace.begin(),art_qHist.trace.end()) - art_qHist.trace.begin();


              double rollingthreshold;                        // for rolling threshold Q-method histogram
              rollingthreshold = ThresMultiplier_ * noiseRMS[iSeg]; // set thresholds from noise RMSdetermination
              if (rollingthreshold < ThresAbsolute_)
                {
                 rollingthreshold = ThresAbsolute_;              // override noise threshold with absolute threshold unless noise too large
		}
             else
		{
		  HnoisyCaloXtal->Fill(iSeg,iCalo-1);
		  rollingthreshold = ThresSign_ * rollingthreshold; // set sign for negative, positive going pulses (we changed the polarity in 2017)
		}
              hthreshold_->Fill(rollingthreshold);              // build rolling threshold histogram

              double newbinerror;      // for energy-weighted bin errors
              int wndw = ThresWindow_; // pedestal sample window for rolling threshold is +/-wndw

              // loop over software rebins SoftRB_
              int irb = 1; // software rebinning factor

                if (wndw < irb)
                  printf("error - pedestal window %i less than software rebinning %i\n", wndw, irb); // note  doesnt make sense if wndw is smaller than SOFTRB_

                // loop over samples for given xtal, calo (iSeg, iCalo)
                iADC = 0; // zero sample index

       
 
                for (int ilastfill = 0; ilastfill < nBinsNew; ilastfill++)
                {                                                                                      // memory problems with histograms? try array
                  t_lastfill_tmp_[ilastfill] = t_lastfill_[ilastfill][iSeg + 54 * (iCalo - 1)]; // make temp copy of time distribution for iSeg, iCalo from last fill array
                  t_lastfill_[ilastfill][iSeg + 54 * (iCalo - 1)] = 0.0;                        // zero time distribution for iSeg, iCalo in last fill array
                }
                
                for (double ADC_Value : art_qHist.trace)
                {
                  iADC++;           // increment sample counter
                  ADC_Value += 0.0; // use ADC_value to avoid compiler error
                  ADC_Value = ADC_Value / binsizefactor;

                  if (iADC % irb)
                    continue; // apply software rebinning

                  double ADC_Value_rebinned = 0.0; // calculate the software rebinned ADC, ADC_Value_rebinned
                  for (int jADC = iADC - irb; jADC <= iADC - 1; jADC++)
                    ADC_Value_rebinned += art_qHist.trace.at(jADC) / binsizefactor;

		  // int tskipflash = 200; // used to skip the region of flash in histogramming rejected points in pedestal calculation
                  int tskipflash = 4*200; // factor of 4 due to time bin size of run2
                  int windoffsetlo = LoGap_, windoffsethi = HiGap_;
                  int shift_puL = 2 * wndw + 1 + windoffsetlo, shift_puH = 2 * wndw + 1 + windoffsethi; // calculate the relative positions of lo, hi pileup windows
                  if ((iADC >= wndw + shift_puL + irb + windoffsetlo) && (iADC <= nBins - wndw - shift_puH - irb - windoffsethi))
                  { // bookend samples to out-of bounds lo, hi pileup windows

                    double yaveragelesskthsample = 0.0, yaverage = 0.0, ysum = 0.0, ysum2 = 0.0, ydiff = 0.0, ysigma = 0.0; // rolling threshold variables - diff sum, sigma, sigma

                    //if (iADC == 500 && irb == 1) printf("iADC-1, ADC_Value_rebinned %i %f\n", iADC-1, ADC_Value_rebinned);
                    // calculate the rolling average (i.e. pedestal) of neighboring (+/- wndw) )ADC samples
                    for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      if (kADC != iADC / irb)
                      {
                        for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                        {
                          if (jADC > iADC)
                          {
                            ysum += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                            //if (iADC == 500 && irb == 1) printf(" art_qHist.trace.at(jADC-1+windoffset)/binsizefactor, jADC-1+windoffset, ysum %f %i %f\n", art_qHist.trace.at(jADC-1+windoffset)/binsizefactor, jADC-1+windoffset, ysum);
                          }
                          if (jADC < iADC)
                          {
                            ysum += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                            //if (iADC == 500 && irb == 1) printf(" art_qHist.trace.at(jADC-1+windoffset), jADC-1+windoffset, ysum %f %i %f\n", art_qHist.trace.at(jADC-1-windoffset)/binsizefactor, jADC-1-windoffset, ysum);
                          }
                        }                                    // loop over raw samples of rebinned samples, 1 sample if irb=1, 2 samples if irb=2, etc
                      }                                      // avoid the trigger sample
                    }                                        // loop over rebinned samples in pedestal window
                    yaverage = ysum / (2. * ((double)wndw)); // rolling average pedestal around trigger sample w/o ADC calibration
                    //if (iADC == 500 && irb == 1) printf("ysum, yaverage %f %f\n", ysum, yaverage);
                    qHist_2D_averageold_->Fill(iADC, yaverage);

                    // calculate the sigma of rolling average of ADC samples
                    for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      if (kADC != iADC / irb)
                      {
                        for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                        {
                          if (jADC > iADC)
                            ysum2 += (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaverage) * (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaverage);
                          if (jADC < iADC)
                            ysum2 += (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaverage) * (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaverage);
                        }                                         // loop over raw samples of rebinned samples, 1 sample if irb=1, 2 samples if irb=2, etc
                      }                                           // avoid the trigger sample
                    }                                             // loop over rebinned samples in pedestal window
                    ysigma = sqrt(ysum2 / (2. * ((double)wndw))); // rolling average pedestal sigma around trigger sample w/o calibration

                    // use infill gain correction
		       double infillgaincorrection = 1.0;
		                         double IFGcorrection = 1.0;
		       if(useInFillGain)
			 {
			   
			 			   

			   if(iADC - flashbin >= 0 && irb == 1)
			       {
				 IFGcorrection = IFGfunc->Eval((iADC - flashbin)/1000.0 * binsizefactor);
                                 //printf("iADC is %i, %i, %f\n", iADC, iADC - flashbin, IFGfunc->Eval(iADC - flashbin));
			       }
			   			     //v2.at(timebin) = v1.at(timebin)/IFGcorrection;
						      // }
			      //printf("IFG correction is %f, %i, %i, %f\n",IFGfunc->Eval((iADC - flashbin)/1000.0), iSeg, flashbin, binsizefactor);
                           //if(irb==1 && iADC>=320 && iADC<=400) std::cout<<"ifg correctiom"<<iADC<<" "<<IFGcorrection<<std::endl;

			 }//end of if(useInFillGain)
		       
		        
                    // pileup rejection based on dropping above-threshold samples. PedMaxDev_ is multiplicative factor for dropping pedestal samples and switching dropping on, off
                    int ndrop = 0;                 // ndrop is number of  pedestal samples dropped from pedestal average calculation
                    bool dropsample[33] = {false}; // true/false array of dropped pedestal samples (dimension of 33 for max of 16 samples in lo, hi windows)
                    double pedhit = 0.0;           // pedestal-subtracted, sign-corrected, energy-calibrated, pedestal value for testing against threshold
                    for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      if (kADC != iADC / irb)
                        for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                        {
                          if (jADC > iADC)
                            yaveragelesskthsample = 2. * ((double)wndw) / (2. * ((double)wndw) - 1.) * (yaverage - (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor / (2. * ((double)wndw)))); // average pedestal without contribution of test sample
                          if (jADC < iADC)
                            yaveragelesskthsample = 2. * ((double)wndw) / (2. * ((double)wndw) - 1.) * (yaverage - (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor / (2. * ((double)wndw)))); // average pedestal without contribution of test sample
                          if (jADC > iADC)
                            pedhit = ThresSign_ * (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaveragelesskthsample); // calc signal on pedestal sample
                          if (jADC < iADC)
                            pedhit = ThresSign_ * (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaveragelesskthsample); // calc signal on pedestal sample
                          if (useCal) pedhit = CalSeg[iSeg] * pedhit; // use calibration if requested
			  // printf("chrystal number is %i\n",iSeg);

			  //  if (useInFillGain) pedhit = pedhit / infillgaincorrection;  //use infill gain correction
                          if (useInFillGain) pedhit = pedhit / IFGcorrection; //IFG from database
                          
                          if (pedhit > rollingthreshold * PedMaxDev_)
                          { // reject pedestal sample if exceeding threshold setting

                            dropsample[jADC - iADC + wndw + irb - 1] = true; // record rejected sample
                            ndrop++;                                         // count rejected samples
                            // fill diagnostic histograms for rejected pedestal samples
                            if (iADC > tskipflash)
                            {
                              if (jADC > iADC)
                                qHist_2D_pointrejection2_->Fill(32 + jADC - iADC + wndw + irb - 1, art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaveragelesskthsample); // store pedestal subtracted value of dropped sample
                              if (jADC < iADC)
                                qHist_2D_pointrejection2_->Fill(32 + jADC - iADC + wndw + irb - 1, art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaveragelesskthsample); // store pedestal subtracted value of dropped sample
                            }
                            if (irb == 1 && iADC > tskipflash)
                            { // avoid flash
                              if (jADC > iADC)
                                qHist_2D_noise_v_signal_->Fill(ADC_Value_rebinned, art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor); // store diagnostics on energy correlations
                              if (jADC < iADC)
                                qHist_2D_noise_v_signal_->Fill(ADC_Value_rebinned, art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor); // store diagnostics on energy correlations
                              if (jADC > iADC)
                                qHist_2D_rejectedpoints_->Fill(jADC, art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaveragelesskthsample); // store diagnostics data on rejected sample (compare to low-side, high-side pu histograms
                              if (jADC < iADC)
                                qHist_2D_rejectedpoints_->Fill(jADC, art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaveragelesskthsample); // store diagnostics data on rejected sample (compare to low-side, high-side pu histograms
                            }
                          }
                          else
                          { // fill diagnostic histograms for accepted pedestal samples
                            if (irb == 1 && iADC > tskipflash)
                            {
                              if (jADC > iADC)
                                qHist_2D_acceptedpoints_->Fill(jADC, art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor); // store diagnostics on accepted sample
                              if (jADC < iADC)
                                qHist_2D_acceptedpoints_->Fill(jADC, art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor); // store diagnostics on accepted sample
                            }
                          }
                        } // if to avoid signal sample
                    }     // loop over all pedestal samples
		    
                    // calculate the new rolling average of ADC samples with data point rejection
                    int cntr = 0;
                    int cntrlo = 0, cntrhi = 0;
                    double xysum = 0.0, x2sum = 0.0;
                    double yaveragelo = 0.0, yaveragehi = 0.0, ydifflohi = 0.0;
                    ysum = 0.0;
                    for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      if (kADC != iADC / irb)
                      {
                        for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                        {
                          if (!dropsample[jADC - iADC + wndw + irb - 1])
                          {
                            if (jADC > iADC)
                              ysum += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                            if (jADC < iADC)
                              ysum += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                            if (jADC > iADC)
                              xysum += (jADC - iADC) * art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                            if (jADC < iADC)
                              xysum += (jADC - iADC) * art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                            x2sum += (jADC - iADC) * (jADC - iADC);
                            cntr++;
                            if (kADC < iADC / irb)
                            {
                              if (jADC > iADC)
                                yaveragelo += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                              if (jADC < iADC)
                                yaveragelo += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                              cntrlo++;
                            }
                            if (kADC > iADC / irb)
                            {
                              if (jADC > iADC)
                                yaveragehi += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                              if (jADC < iADC)
                                yaveragehi += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                              cntrhi++;
                            }
                          } // if pedestal sample isn't dropped sample
                        }   // loop over raw samples of rebinned samples, 1 sample if irb=1, 2 samples if irb=2, etc
                      }     // avoid the trigger sample
                    }       // loop over rebinned samples in pedestal window

                    double slpe = xysum / x2sum;   // straight-line slope of pedestal samples
                    double yaverageold = yaverage; // store old average pedestal before pedestal sample rejection
                    if (cntrlo > 0)
                      yaveragelo /= cntrlo;
                    if (cntrlo > 0)
                      yaveragehi /= cntrhi;
                    ydifflohi = yaveragehi - yaveragelo;
                    yaverage = ysum / (2. * (double)wndw - (double)ndrop); // w/o calibration

                    qHist_2D_average_->Fill(iADC, yaverage);
                    qHist_2D_slpe_->Fill(iADC, slpe);
                    qHist_2D_ydifflohi_->Fill(iADC, ydifflohi);

                    // calculate the slope for low-side window of trigger sample
                    ysum = 0.0;
                    xysum = 0.0;
                    x2sum = 0.0;
                    double nsum = 0.0, xsum = 0.0;
                    for (int kADC = iADC / irb - wndw / irb; kADC < iADC / irb; kADC++)
                    {
                      for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                      {
                        if (!dropsample[jADC - iADC + wndw + irb - 1])
                        {
                          nsum += 1.0;
                          xsum += jADC;
                          if (jADC > iADC)
                            ysum += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                          if (jADC < iADC)
                            ysum += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                          if (jADC > iADC)
                            xysum += jADC * art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                          if (jADC < iADC)
                            xysum += jADC * art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                          x2sum += jADC * jADC;
                        }
                      }
                    }
                    double slpeL = (nsum * xysum - xsum * ysum) / (nsum * x2sum - xsum * xsum);
                    qHist_2D_slpeL_->Fill(iADC, slpeL); // diagnostic histogram of slope of pedestal

                    // calculate the slope for high-side window of trigger sample
                    ysum = 0.0;
                    xysum = 0.0;
                    x2sum = 0.0;
                    nsum = 0.0;
                    xsum = 0.0;
                    for (int kADC = iADC / irb + 1; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                      {
                        if (!dropsample[jADC - iADC + wndw + irb - 1])
                        {
                          nsum += 1.0;
                          xsum += jADC;
                          if (jADC > iADC)
                            ysum += art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                          if (jADC < iADC)
                            ysum += art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                          if (jADC > iADC)
                            xysum += jADC * art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor;
                          if (jADC < iADC)
                            xysum += jADC * art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor;
                          x2sum += jADC * jADC;
                        }
                      }
                    }
                    double slpeH = (nsum * xysum - xsum * ysum) / (nsum * x2sum - xsum * xsum);
                    qHist_2D_slpeH_->Fill(iADC, slpeH); // diagnostic histogram of slope of pedestal

                    // calculate the new sigma of rolling average of ADC samples with point rejection
                    ysum2 = 0.0;
                    for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                    {
                      if (kADC != iADC / irb)
                      {
                        for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                        {
                          if (!dropsample[jADC - iADC + wndw + irb - 1])
                          {
                            if (jADC > iADC)
                              ysum2 += (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaverage) * (art_qHist.trace.at(jADC - 1 + windoffsethi) / binsizefactor - yaverage); // if pedestal sample isn't dropped sample
                            if (jADC < iADC)
                              ysum2 += (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaverage) * (art_qHist.trace.at(jADC - 1 - windoffsetlo) / binsizefactor - yaverage); // if pedestal sample isn't dropped sample
                          }
                        }                                                         // loop over raw samples of rebinned samples, 1 sample if irb=1, 2 samples if irb=2, etc
                      }                                                           // avoid the trigger sample
                    }                                                             // loop over rebinned samples in pedestal window
                    

                    double ysigmaold = ysigma;                                    // store the old sigma of pedestal samples without sample rejection
                    ysigma = sqrt(ysum2 / (2. * (double)(wndw) - (double)ndrop)); // w/o calibration

                    if (irb == 1 && iADC > tskipflash)
                    { // diagnostic histograms of sigma's of pedestal samples
                      qHist_2D_sigma_->Fill(iADC, ysigma);
                      qHist_2D_sigmaold_->Fill(iADC, ysigma - ysigmaold);
                    }
                    if (iADC > tskipflash)
                      qHist_2D_pointrejection_->Fill(32 + ndrop, yaverageold - yaverage); // store diagnostics data on rejected samples

                    if (useCal)
                    {                                                         //if using calibration we multiplying everything by cal. constant (so be careful on the meaning of the pedestal)
		      //std::cout<<"before calibration, ADC_value_rebinned "<<ADC_Value_rebinned<<endl;
                      yaverage = CalSeg[iSeg] * yaverage;                     // w/ calibration
                      ADC_Value_rebinned = CalSeg[iSeg] * ADC_Value_rebinned; // w/ calibration
		      //std::cout<<"after calibration, ADC_value_rebinned "<<ADC_Value_rebinned<<endl;
                    }

                    if (yaverage < -1.e6 || yaverage > 1.e6)
                      printf("warning yaverage out-of-range, iADC %i, irb %i, yaverage %f, ndrop %i\n", iADC, irb, yaverage, ndrop); // warn for inifinities in calculation of pedestal

                    ydiff = ADC_Value_rebinned - irb * yaverage;       // calculate pedestal subtracted signal for trigger sample
                    //ydiff = ThresSign_ * ydiff / infillgaincorrection; // handle signal polarities - legacy issue
                    ydiff = ThresSign_ * ydiff / IFGcorrection;
            
                    // fill diagnostic histograms
		    
		      qHist_1D_[iCalo - 1]->AddBinContent(iADC / irb, ADC_Value_rebinned); // store all raw samples
		      
                    if (iADC / irb >= 120 && iADC / irb <= 125)
                      qHist_2D_sig_flash_->Fill(iSeg % 9, iSeg / 9, -ADC_Value_rebinned); // diagnostics - note sign change to make pedestal from flash go positive
                    if (iADC / irb >= 120 && iADC / irb <= 125)
                      qHist_2D_sig_flashcalo_[iCalo - 1]->Fill(iSeg % 9, iSeg / 9, -ADC_Value_rebinned); // diagnostics - note sign change to make pedestal from flash go positive
                    qHist_1D_ped_all_[iCalo - 1]->AddBinContent(iADC / irb, yaverage);            // store calculated rolling-average pedestal
                    entriesNped++;                                                                // counter for histogram entries

                    if (ydiff >= rollingthreshold)
                    { // pedestal-subtracted value of trigger sample exceeds threshold

                      
                      t_lastfill_[iADC / irb][iSeg + 54 * (iCalo - 1)] = ydiff; // fill above-threshold signals into last fill array for iSeg, iCalo

                      // calculate the rolling average of ADC samples for lower pileup window
                      double ADC_Value_rebinned_puL = 0.0, yaverage_puL = 0.0, ysum_puL = 0.0, ydiff_puL = 0.0; // low-side pileup sample sum, average, difference
                      for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                      {
                        if (kADC != iADC / irb)
                        { // avoid trigger sample
                          for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                          {
                            ysum_puL += art_qHist.trace.at(jADC - shift_puL - 1 - windoffsetlo) / binsizefactor; // lo-side pedestal sum
                                                                                                                 //if (irb == 1) printf("jADC, jADC-shift_puL-1, ysum_puL %i, %i, %f\n", jADC, jADC-shift_puL-1, ysum_puL); //debug
                          }
                        }
                        else
                        { // get trigger sample
                          for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                          {
                            ADC_Value_rebinned_puL += art_qHist.trace.at(jADC - shift_puL - 1 - windoffsetlo) / binsizefactor; // lo-side pileup sample (software-rebinned)
                                                                                                                               //if (irb == 1) printf("jADC, jADC-shift_puL-1, ADC_Value_rebinned_puL %i, %i, %f\n", jADC, jADC-shift_puL-1, ADC_Value_rebinned_puL); //debug
                          }                                                                                                    // loop over software rebinning
                        }
                      }
                      yaverage_puL = ysum_puL / (2. * ((double)wndw)); // average pedestal for lo-side pu sample w/o calibration

                      // calculate the rolling average of ADC samples for higher pileup window
                      double ADC_Value_rebinned_puH = 0.0, yaverage_puH = 0.0, ysum_puH = 0.0, ydiff_puH = 0.0; // low-side pileup sample sum, average, difference
                      for (int kADC = iADC / irb - wndw / irb; kADC <= iADC / irb + wndw / irb; kADC++)
                      {
                        if (kADC != iADC / irb)
                        { // avoid trigger sample
                          for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                          {
                            ysum_puH += art_qHist.trace.at(jADC + shift_puH - 1 + windoffsethi) / binsizefactor; // hi-side pedestal sum
                          }
                        }
                        else
                        { // get trigger sample
                          for (int jADC = kADC * irb - irb + 1; jADC <= kADC * irb; jADC++)
                          {
                            ADC_Value_rebinned_puH += art_qHist.trace.at(jADC + shift_puH - 1 + windoffsethi) / binsizefactor; // hi-side pileup sample (software-rebinned)
                          }                                                                                                    // loop over software rebinning
                        }
                      }
                      yaverage_puH = ysum_puH / (2. * ((double)wndw)); // average pedestal for lo-side pu sample w/o calibration

                      if (useCal)
                      {                                                                 //if using calibration we multiplying everything by cal. constant (so be careful on the meaning of the pedestal)
                        yaverage_puL = CalSeg[iSeg] * yaverage_puL;                     // w/ calibration
                        yaverage_puH = CalSeg[iSeg] * yaverage_puH;                     // w/ calibration
                        ADC_Value_rebinned_puL = CalSeg[iSeg] * ADC_Value_rebinned_puL; // w/ calibration
                        ADC_Value_rebinned_puH = CalSeg[iSeg] * ADC_Value_rebinned_puH; // w/ calibration
                      }

                      ydiff_puL = ADC_Value_rebinned_puL - irb * yaverage_puL;   // calculate pedestal subtracted signal for lo-side pu sample
                      ydiff_puH = ADC_Value_rebinned_puH - irb * yaverage_puH;   // calculate pedestal subtracted signal for hi-side pu sample
                      ydiff_puL = ThresSign_ * ydiff_puL / infillgaincorrection; // handle signal polarities - legacy issue
                      ydiff_puH = ThresSign_ * ydiff_puH / infillgaincorrection; // handle signal polarities - legacy issue

                      //double ydiff_puA = qHist_1D_lastfill_tmp_->GetBinContent( iADC/irb); // pileup signal from prior fill (error handling done later)
                      double ydiff_puA = t_lastfill_tmp_[iADC / irb] / infillgaincorrection; // grab signal in last fill shadow window at time iADC/irb from temp array with iSeg, iCalo data

                      qHist_1D_rej_hits_->Fill(iADC, ndrop); // record the trigger time for rejected pedestal samples
                      if (ydiff_puL >= rollingthreshold)
                        qHist_1D_puL_hits_->Fill(iADC); // record the trigger time for low-side shadow window pu hits
                      if (ydiff_puH >= rollingthreshold)
                        qHist_1D_puH_hits_->Fill(iADC); // record the trigger time for low-side shadow window pu hits

                      // pileup histograms for adjacent fill pileup window
                      if (ydiff_puA >= rollingthreshold)
                      {
                        
                        if ((ydiff - ydiff_puA / (2. * ((double)wndw))) < rollingthreshold)
                        {                                                                                                 // case where pileup causes complete loss of trigger sample signal
                          qHist_1D_puA_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff);                  // lower pile-up window fill time distribution
                          qHist_1D_puA_[iCalo - 1]->SetEntries(1 + qHist_1D_puA_[iCalo - 1]->GetEntries()); // because I'm using AddBinContent()
                        }
                        if ((ydiff - ydiff_puA / (2. * ((double)wndw))) >= rollingthreshold)
                        {                                                                                                  // case where pileup causes reduced value for trigger sample signal
                          qHist_1D_puA_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff_puA / (2. * wndw)); // lower pile-up window fill time distribution
                          qHist_1D_puA_[iCalo - 1]->SetEntries(1 + qHist_1D_puA_[iCalo - 1]->GetEntries());  // because I'm using AddBinContent()
                        }
                      } // end adjacent fill pu window

                      // pileup histograms for lo-side pileup window
                      if (ydiff_puL >= rollingthreshold)
                      {
                    
                        if ((ydiff - ydiff_puL / (2. * ((double)wndw))) < rollingthreshold)
                        {                                                                                                 // case where pileup causes complete loss of trigger sample signal
                          qHist_1D_puL_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff);                  // lower pile-up window fill time distribution
                          qHist_1D_puL_[iCalo - 1]->SetEntries(1 + qHist_1D_puL_[iCalo - 1]->GetEntries()); // because I'm using AddBinContent()
                                                                                                                          //printf("puL killed ydiff, iCalo, iSeg, %i, %i iADC %i, ydf %f, ydfpuL %f, yav %f yavpuL %f\n", iCalo, iSeg, iADC, ydiff, ydiff_puL, yaverage, yaverage_puL);
                        }
                        if ((ydiff - ydiff_puL / (2. * ((double)wndw))) >= rollingthreshold)
                        {                                                                                                  // case where pileup causes reduced value for trigger sample signal
                          qHist_1D_puL_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff_puL / (2. * wndw)); // lower pile-up window fill time distribution
                          qHist_1D_puL_[iCalo - 1]->SetEntries(1 + qHist_1D_puL_[iCalo - 1]->GetEntries());  // because I'm using AddBinContent()
                                                                                                                           //printf("puL wounded ydiff, iCalo, iSeg, %i, %i iADC %i, ydf %f, ydfpuL %f, yav %f yavpuL %f\n", iCalo, iSeg, iADC, ydiff, ydiff_puL, yaverage, yaverage_puL);
                        }
                      } // end lo-side pu window

                      // pileup histograms for hi-side pileup window
                      if (ydiff_puH >= rollingthreshold)
                      {
                        if ((ydiff - ydiff_puH / (2. * ((double)wndw))) < rollingthreshold)
                        {                                                                                                 // case where pileup causes complete loss of trigger sample signal
                          qHist_1D_puH_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff);                  // higher pile-up window fill time distribution
                          qHist_1D_puH_[iCalo - 1]->SetEntries(1 + qHist_1D_puH_[iCalo - 1]->GetEntries()); // because I'm using AddBinContent()
                                                                                                                          //printf("puH killed ydiff, iCalo, iSeg, %i, %i iADC %i, ydf %f, ydfpuH %f, yav %f yavpuH %f\n", iCalo, iSeg, iADC, ydiff, ydiff_puH, yaverage, yaverage_puH);
                        }
                        if ((ydiff - ydiff_puH / (2. * ((double)wndw))) >= rollingthreshold)
                        {                                                                                                  // case where pileup causes reduced value for trigger sample signal
                          qHist_1D_puH_[iCalo - 1]->AddBinContent(iADC / irb, 2. * wndw * ydiff_puH / (2. * wndw)); // higher pile-up window fill time distribution
                          qHist_1D_puH_[iCalo - 1]->SetEntries(1 + qHist_1D_puH_[iCalo - 1]->GetEntries());  // because I'm using AddBinContent()
                                                                                                                           //printf("puH wounded ydiff, iCalo, iSeg, %i, %i iADC %i, ydf %f, ydfpuH %f, yav %f yavpuH %f\n", iCalo, iSeg, iADC, ydiff, ydiff_puH, yaverage, yaverage_puH);
                        }
                      } // end hi-side pu window
                      
		      if(useActiveChannelList_ && channelStatusDB_.xtalStatus.at(iCalo - 1).at(iSeg))
		      {
                 	 if (iADC > tskipflash) qHist_1D_ydiff_[iCalo-1]->Fill( ydiff);
		      }    
                    
                      int lobin_p = 110, hibin_p = 2950; // lobin, hibin for proton launch, run 15960
                      int lobin_e = 660, hibin_e = 2950; // lobin, hibin for electron signal, run 15960
                      if (iADC >= lobin_p && iADC < hibin_p && irb == 1)
                        qHist_1D_enrgy_p_[iCalo - 1]->Fill(ydiff); // fill energy distribution
                      if (iADC >= lobin_e && iADC < hibin_e && irb == 1)
                        qHist_1D_enrgy_e_[iCalo - 1]->Fill(ydiff); // fill energy distribution

                      qHist_1D_ped_usd_[iCalo - 1]->AddBinContent(iADC / irb, yaverage); // records pedestal for above threshold signals
                      qHist_1D_ped_ent_[iCalo - 1]->AddBinContent(iADC / irb, 1.0);      // counts entries of above threshold signals
		  
        
		      
			
		      if(iCalo-1==16){
			qHist_1D_sig_xtal_[iSeg]->AddBinContent(iADC / irb, ydiff);
			//  printf("chrystal number is %i\n",iSeg);   
		      }

		      if(useActiveChannelList_ && channelStatusDB_.xtalStatus.at(iCalo - 1).at(iSeg))
                      {
			// mf::LogDebug("IslandTemplateFit") << "Marking Fit in Calo " << cn << " xtal "
			//                                   << xn << " as bad!" << std::endl;
			qHist_1D_sig_[iCalo - 1]->AddBinContent(iADC / irb, ydiff);
		      }//hijacking this loop to mark the fits with a negative status if in the ignore list
                      /*else
		      {
		        qHist_1D_sig_[iCalo - 1]->AddBinContent(iADC / irb, ydiff);
		      }
		      */
                             
		    
		      // printf("chrystal number is %i\n",iSeg);
                      // diagnostic histograms for time distribution by fill sequence index
                      if (irb == 1 && sequenceIndex >= 0 && sequenceIndex < 8)
                      {
                        qHist_1D_sig_seqnum_[sequenceIndex]->AddBinContent(iADC / irb, ydiff);                                                                                           // fill time distribution versus seq number
                        newbinerror = sqrt(ydiff * ydiff + qHist_1D_sig_seqnum_[sequenceIndex]->GetBinError(iADC / irb) * qHist_1D_sig_seqnum_[sequenceIndex]->GetBinError(iADC / irb)); // this is improper error handling
                        qHist_1D_sig_seqnum_[sequenceIndex]->SetBinError(iADC / irb, newbinerror);                                                                                       // this is improper error handling
                        qHist_1D_sig_seqnum_[sequenceIndex]->SetEntries(1 + qHist_1D_sig_seqnum_[sequenceIndex]->GetEntries());                                                          // set entries as number of fills for particular sequence cycle
                      }

                      if (irb == 1)
                      { // diagnostic histograms of energy correlations between trigger sample and pu samples
                        qHist_2D_puLE1E2_[iCalo - 1]->Fill(ydiff, ydiff_puL);
                        qHist_2D_puHE1E2_[iCalo - 1]->Fill(ydiff, ydiff_puH);
                      }

                      double yerrwithpedsigma = sqrt( ydiff*ydiff + ysigma*ysigma / (2. * (double)(wndw) - (double)ndrop)); // w/o calibration
                    
		      qHist_1D_ESumPerFill_[iCalo - 1]->AddBinContent( iADC / irb, yerrwithpedsigma); // diagnostic fill-by-fill xtal-sum, energy distrbution
		
                      qHist_1D_NHitPerFill_[iCalo - 1]->AddBinContent(iADC / irb, 1);     // diagnostic fill-by-fill xtal-sum, hit distribution
                      entriesNrolling++;

                      if (irb == 1)
                      { // vertical slice, horizontal slice histograms (only filled for software rebinning = 1)

		                                                                                          
                          {
                            //std::cout<<"skipping noisy crystal!! iseg,icalo "<<iSeg<<" , "<<iCalo<<endl;                                                                                                
                           qHist_1D_sig_vslice_[iSeg % 9]->AddBinContent(iADC / irb, 0);      // fill horizontal slice time distribution avoiding noisy xtal                                               
			  }        // time distribution of above threshold signals (error handling done later)      
			*/
			//else
			//{
			if(useActiveChannelList_ && channelStatusDB_.xtalStatus.at(iCalo - 1).at(iSeg))
                        {
			  // mf::LogDebug("IslandTemplateFit") << "Marking Fit in Calo " << cn << " xtal "
			  //     << xn << " as bad!" << std::endl;
			  qHist_1D_sig_vslice_[iSeg % 9]->AddBinContent(iADC / irb, ydiff);
			}//hijacking this loop to mark the fits with a negative status if in the ignore list                                                                                               
		                           newbinerror = sqrt(ydiff * ydiff + qHist_1D_sig_vslice_[iSeg % 9]->GetBinError(iADC / irb) * qHist_1D_sig_vslice_[iSeg % 9]->GetBinError(iADC / irb));
                        qHist_1D_sig_vslice_[iSeg % 9]->SetBinError(iADC / irb, newbinerror); // incorrect error handling
                        entriesNrollingvslice[iSeg % 9]++;

        	                        //else
			//{
			if(useActiveChannelList_ && channelStatusDB_.xtalStatus.at(iCalo - 1).at(iSeg))
                        {
                          // mf::LogDebug("IslandTemplateFit") << "Marking Fit in Calo " << cn << " xtal "                                                                                                 
                          //     << xn << " as bad!" << std::endl;                                                                                                                                         
                          qHist_1D_sig_hslice_[iSeg / 9]->AddBinContent(iADC / irb, ydiff);
			}//hijacking this loop to mark the fits with a negative status if in the ignore list                                                                                               
                                         newbinerror = sqrt(ydiff * ydiff + qHist_1D_sig_hslice_[iSeg / 9]->GetBinError(iADC / irb) * qHist_1D_sig_hslice_[iSeg / 9]->GetBinError(iADC / irb));
                        qHist_1D_sig_hslice_[iSeg / 9]->SetBinError(iADC / irb, newbinerror); // incorrect error handling
                        entriesNrollinghslice[iSeg / 9]++;

                        if (iADC / irb >= 180)
                          qHist_2D_sig_lego_->Fill(iSeg % 9, iSeg / 9, ydiff); // xy distributions avoiding the flash region
                        if (iADC / irb >= 180)
                          qHist_2D_sig_legocalo_[iCalo - 1]->Fill(iSeg % 9, iSeg / 9, ydiff); // xy distributions avoiding the flash region
                      }
                    } // end exceeding rolling threshold
                    
		    
                    // does handle the sharing of positron energy between xtals, but suffers from pileup
                    if (iSeg == 53)
                    { // only do for single xtal
                      newbinerror = sqrt(qHist_1D_ESumPerFill_[iCalo - 1]->GetBinContent(iADC / irb) * qHist_1D_ESumPerFill_[iCalo - 1]->GetBinContent(iADC / irb) + qHist_1D_sig_[iCalo - 1]->GetBinError(iADC / irb) * qHist_1D_sig_[iCalo - 1]->GetBinError(iADC / irb));
                       qHist_1D_sig_[iCalo - 1]->SetBinError(iADC / irb, newbinerror); // energy-weighted error bars
		       //std::cout<<"error bar iseg,icalo "<<iSeg<<", "<<iCalo<<", "<<newbinerror<<endl;
                    }
		    
                  } // end book-ending
                  entriesN++;
		  //            printf(" ______######__________  fill_counter, iCalo, iSeg ,iADC , %i, %i, %i, %i \n ",fill_counter, iCalo, iSeg, iADC);

		} // end loop over trace samples
		//               printf(" ______######__________  fill_counter, iCalo, iSeg %i, %i, %i \n ",fill_counter, iCalo, iSeg);

                
                qHist_1D_[iCalo - 1]->SetEntries(entriesN + qHist_1D_[iCalo - 1]->GetEntries());
                qHist_1D_sig_[iCalo - 1]->SetEntries(entriesNrolling + qHist_1D_sig_[iCalo - 1]->GetEntries());
                if(iCalo-1==16){qHist_1D_sig_xtal_[iSeg]->SetEntries(entriesNrolling + qHist_1D_sig_xtal_[iSeg]->GetEntries());}
                qHist_1D_ped_usd_[iCalo - 1]->SetEntries(entriesNrolling + qHist_1D_ped_usd_[iCalo - 1]->GetEntries());
                qHist_1D_ped_ent_[iCalo - 1]->SetEntries(entriesNrolling + qHist_1D_ped_ent_[iCalo - 1]->GetEntries());
                qHist_1D_ped_all_[iCalo - 1]->SetEntries(entriesNped + qHist_1D_ped_all_[iCalo - 1]->GetEntries());
                //qHist_1D_thr_[iCalo-1]->SetEntries( entriesNfixed + qHist_1D_thr_[iCalo-1]->GetEntries() );
                if (irb == 1)
                {
                  for (int islice = 0; islice < 9; islice++)
                  {
                    qHist_1D_sig_vslice_[islice]->SetEntries(entriesNrollingvslice[islice] + qHist_1D_sig_vslice_[islice]->GetEntries());
                  }
                  for (int islice = 0; islice < 6; islice++)
                  {
                    qHist_1D_sig_hslice_[islice]->SetEntries(entriesNrollinghslice[islice] + qHist_1D_sig_hslice_[islice]->GetEntries());
                  }
                }

                //qHist_2D_[iCalo-1]->SetEntries( entriesN + qHist_1D_[iCalo-1]->GetEntries() ); // large 2D histograms of E vs T
                //qHist_2D_Coarse_[iCalo-1]->SetEntries(entriesN); // large 2D histograms of E vs T
                

   
              

              break;
            } //end of case: true.
	    
            
            } //End of switch
	             
          } //End loop over segments/channels/crystals

        } //End loop over calorimeters
         
    
         
      } // empty Qdata check
    }   // calo 1-24 check
   

 
    nonemptyQData.clear();
  } // empty Qhist check
  
  //} //End of If return_success
  //  printf("check value %i \n", qHist_1D_sig_[0][0]->GetNbinsX());
  // fill_counter++;
  // printf("Fill number is %i \n", fill_counter);
      if (empty_data == 0)
    {
      if(qHist_1D_ydiff_sum == nullptr)
	{
	  qHist_1D_ydiff_sum = TFileServ->make<TH1D>( Form("qHist_1D_ydiff_sum"),
						      Form("qHist_1D_ydiff_sum"),
						      8077, -100*23.64, 4000*23.64);
	}
  
      if(qHist_1D_sig_sum == nullptr)
        {
          qHist_1D_sig_sum = TFileServ->make<TH1D>( Form("qHist_1D_sig_sum"), 
                                                      Form("qHist_1D_sig_sum"), 
						    16800, 100001, 352001);
        }
       }
 
  gettimeofday(&t_end, NULL);
  printf("QFillByFillAnalyzer duration: %f s \n", (t_end.tv_sec - t_start.tv_sec) + ((t_end.tv_usec - t_start.tv_usec) / 1000000.0));

  mf::LogDebug("[QFillByFillAnalyzer]") << "Exit QFillByFillAnalyzer::analyze(" << event.id() << ")";
  //  printf("calo num and xtal num are: %i \n",iCalo,iSeg);

  printf("End of analyzer \n");
} //analyze

// These are some necessary boilerplate for the ROOT persistency system
using gm2calo::QFillByFillAnalyzer;
DEFINE_ART_MODULE(QFillByFillAnalyzer)

} // End of namespace gm2calo
