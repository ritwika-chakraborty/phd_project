/*****************************************************************************\
Name:   CaloSC_frontend.cxx
Author: Aaron Fienberg
Email:  fienberg@uw.edu
About:  Calo slow control midas frontend
\*****************************************************************************/

//--- std includes ----------------------------------------------------------//
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <random>
#include <sys/time.h>
#include <cassert>
#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>


//--- other includes --------------------------------------------------------//
#include "midas.h"

//--- project includes ------------------------------------------------------//
#include "CaloSCODB.hh"
#include "CaloSCHelpers.hh"

//--- globals ---------------------------------------------------------------//

extern "C" {

// The frontend name (client name) as seen by other MIDAS clients
const char *frontend_name = "CaloSC";

// The frontend file name, don't change it.
char *frontend_file_name = (char *)__FILE__;

// frontend_loop is called periodically if this variable is TRUE
BOOL frontend_call_loop = TRUE;

// A frontend status page is displayed with this frequency in ms.
INT display_period = 0;

// maximum event size produced by this frontend
INT max_event_size = 0x8000;  // 32 kB @EXAMPLE - adjust if neeeded

// maximum event size for fragmented events (EQ_FRAGMENTED)
INT max_event_size_frag = 0x800000;  // DEPRECATED

// buffer size to hold events
INT event_buffer_size = DEFAULT_MAX_EVENT_SIZE;

  //extern INT frontend_index;  // frontend index from command line

extern INT run_state;

// Function declarations
INT frontend_init();
INT frontend_exit();
INT begin_of_run(INT run_number, char *error);
INT end_of_run(INT run_number, char *error);
INT pause_run(INT run_number, char *error);
INT resume_run(INT run_number, char *error);

INT frontend_loop();
INT read_event(char *pevent, INT off);
INT poll_event(INT source, INT count, BOOL test);
INT interrupt_configure(INT cmd, INT source, PTYPE adr);

// Equipment list

EQUIPMENT equipment[] = {{
    //                          "CaloSC%02d", /* equipment name */
           "CaloSlowControls", /*equipment name */
                          {
                           // EQUIPMENT_INFO
                           2,           /* event ID */
                           0,           /* trigger mask */
                           "SYSTEM",    /* event buffer */
			   //                           EQ_PERIODIC, /* equipment type */
                           EQ_POLLED, /* equipment type */

                           0,           /* event source */
                           "MIDAS",     /* format */
                           TRUE,        /* enabled */
			   //                     RO_ALWAYS,   /* read all the time */
                           RO_RUNNING,   /* read during a run */  
			   //     1000,       /* time of periodic reads (ms) */
			   500,         /* poll for 500 (ms) */
                           0,           /* stop run after this event limit */
                           0,           /* number of sub events */
                           0,           /* log history */
                           "",          /* frontend host */
                           "",          /* frontend name */
                           "",          /* Source file  */
                           "",          /* Current status of equipment */
                           "", /* Color to be used by mhttpd for status */
                          },
	                  // read_periodic_event, /* pointer to readout routine */
                	  read_event, /* pointer to readout routine */

                         },

                         {""}};

}  // extern C

// anonymous namespace for non midas stuff used here
namespace {
// odb structure
std::string odbDir;
  //std::vector<bk_folder> bk_folders(4);
  //std::vector<sipm_folder> sipm_folders(54);
BOOL temp_alarm_enabled[54];
calosc_settings settings;
alarm_monitors monitors = {0, 0};
int n_consecutive_timeouts = 0;

// copy of readings at begin of run
std::vector<bk_folder> bor_bk_folders;
std::vector<sipm_folder> bor_sipm_folders;

// helper classes
//std::unique_ptr<CaloSC::BeagleCommunicator> beagle_comm;
//std::unique_ptr<CaloSC::PSQLCommunicator> psql_comm;

  // Thread buffers                                                                     
  std::vector<std::vector<std::vector<sipm_folder>>> sipm_buffer(54); //i: crate, j\  :event, k: slot

// Thread control                                                                        
bool FrontendActive=false;
bool RunActive=false;

// Mutex locks                                                                           
std::mutex mlock;
std::mutex data_lock;
std::mutex database_lock;


// reload all info into the sipm and bk folders from db and beagles
bool reload_all_info() {
  if (!psql_comm->fill_sipm_maps(sipm_folders)) {
    cm_msg(MERROR, __FILE__, "Could not load SiPM mapping from database!");
    monitors.db_failure = 1;
    return false;
  }
  // try twice to read beagle info before raising alarm
  if (!beagle_comm->read_beagle_info(bk_folders, sipm_folders) && !beagle_comm->read_beagle_info(bk_folders, sipm_folders)) {
    cm_msg(MERROR, __FILE__, "Could not read all values from BeagleBone!");
    monitors.beagle_failure = 1;
    return false;
  }
  return true;
}

}  // end namespace

// check for alarm conditions (checked after each read)
void check_and_raise_alarms();

//--- Frontend Init ---------------------------------------------------------//
INT frontend_init() {
  INT frontend_index = get_frontend_index();  // frontend index from command line 
  // Check that an index has been provided
  if (frontend_index < 0) {
    cm_msg(MERROR, __FUNCTION__,
           "No frontend index supplied (must use -i <index> when starting "
           "frontend)");
    return FE_ERR_HW;
  }

  // init helpers
  //std::string connection_addr =
  //   "tcp://192.168." + std::to_string(frontend_index) + ".21:6669";
  // d0 address
  //"tcp://192.168.1.21:6669";
  // localhost for testing on my laptop
  // std::string connection_addr = "tcp://localhost:6669";
  //beagle_comm = std::make_unique<CaloSC::BeagleCommunicator>(connection_addr);

  // try to turn on JMU supply
  //  if (!beagle_comm->turn_on_JMU()) {
  //  cm_msg(MERROR, __FILE__, "could not connect to BeagleBone!");
  // return FE_ERR_HW;
  //}

  //std::string config_file = "config/dbconnection.json";
  //psql_comm =
  //    std::make_unique<CaloSC::PSQLCommunicator>(config_file, frontend_index);

  //
  // INIT ODB
  //
  //  char indexedName[16];
  // sprintf(indexedName, "CaloSC%02i", frontend_index);

  // get db handle
  HNDLE hDB, hKey;
  cm_get_experiment_database(&hDB, NULL);
  // odbDir = "/Equipment/" + std::string(indexedName);
  odbDir = "/Equipment/CaloSlowControls";

  // settings init
  auto dirStr = odbDir + "/Settings/Globals";
  //  db_create_record(hDB, 0, dirStr.c_str(), settings_str);
  db_check_record(hDB, 0, dirStr.c_str(), settings_str,TRUE);
  db_find_key(hDB, 0, dirStr.c_str(), &hKey);
  /*auto ret_code = db_open_record(hDB, hKey, &settings, sizeof(settings),
                                 MODE_READ, NULL, NULL);
  if (ret_code != SUCCESS) {
    std::string errormsg =
        "failed to open record with error code " + std::to_string(ret_code);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
    }*/

  int size_settings = sizeof(settings);
  auto status  = db_get_record(hDB, hKey, &settings, &size_settings, 0);
  if (status != SUCCESS) {
    std::string errormsg =
      "failed to get record with error code " + std::to_string(status);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
  }
  NumberOfCrates = settings.number_of_crates;
  ReadPeriod = settings.read_period;


  std::string tempal_enabled_str = "[.]";
  for (unsigned int i = 0; i < 54; ++i) {
    tempal_enabled_str += "\nSiPM " + std::to_string(i) + " = BOOL : y";
  }
  dirStr = odbDir + "/Settings/Enabled temperature checks";
  //db_create_record(hDB, 0, dirStr.c_str(), tempal_enabled_str.c_str());
  db_check_record(hDB, 0, dirStr.c_str(), tempal_enabled_str.c_str(), TRUE);
  db_find_key(hDB, 0, dirStr.c_str(), &hKey);
  /*ret_code = db_open_record(hDB, hKey, temp_alarm_enabled, sizeof(temp_alarm_enabled),
                                 MODE_READ, NULL, NULL);
  if (ret_code != SUCCESS) {
    std::string errormsg =
        "failed to open record with error code " + std::to_string(ret_code);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
    }*/
  
  int size_temp_alarm_enabled = sizeof(temp_alarm_enabled);
  auto status  = db_get_record(hDB, hKey, &temp_alarm_enabled, &size_temp_alarm_enabled, 0);
  if (status != SUCCESS) {
    std::string errormsg =
      "failed to get record with error code " + std::to_string(status);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
  }

  

  // monitors init
  dirStr = odbDir + "/Monitors";
  //db_create_record(hDB, 0, dirStr.c_str(), monitors_str);
  db_check_record(hDB, 0, dirStr.c_str(), monitors_str, TRUE);  
  db_find_key(hDB, 0, dirStr.c_str(), &hKey);
  /* ret_code = db_open_record(hDB, hKey, &monitors, sizeof(monitors), MODE_WRITE,
                            NULL, NULL);
  if (ret_code != SUCCESS) {
    std::string errormsg =
        "failed to open record with error code " + std::to_string(ret_code);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
    }*/

  int size_monitors= sizeof(monitors);
  auto status  = db_get_record(hDB, hKey, &monitors, &size_monitors, 0);
  if (status != SUCCESS) {
    std::string errormsg =
      "failed to get record with error code " + std::to_string(status);
    cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
    return FE_ERR_ODB;
  }


  odbDir += "/Hardware";
  
  for (int iCrate = 0 ; iCrate<NumberOfCrates ; iCrate++)
    {
      char crate_name[100];
      sprintf(crate_name,"/Crate%02d",iCrate+1);
      std::string crateDir = odbDir + std::string(crate_name);

  
      // bk folder init
      for (unsigned int i = 0; i < 4; ++i) {
         dirStr = odbDir + "/BK" + std::to_string(i + 1);
         //db_create_record(hDB, 0, dirStr.c_str(), bk_folder_str);
         db_check_record(hDB, 0, dirStr.c_str(), bk_folder_str,TRUE);   
      /*  db_find_key(hDB, 0, dirStr.c_str(), &hKey);

         // link to odb
         bk_folder &this_folder = bk_folders[i];

         ret_code = db_open_record(hDB, hKey, &this_folder, sizeof(this_folder),
                              MODE_WRITE, NULL, NULL);
         if (ret_code != SUCCESS) {
           std::string errormsg =
          "failed to open record with error code " + std::to_string(ret_code);
           cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
           return FE_ERR_ODB;
          }*/
     }

    // sipm folder init
    for (unsigned int i = 0; i < 54; ++i) {
      dirStr = odbDir + "/SiPM" + std::to_string(i);
      //db_create_record(hDB, 0, dirStr.c_str(), sipm_folder_str);
      db_check_record(hDB, 0, dirStr.c_str(), sipm_folder_str,TRUE);

     /* db_find_key(hDB, 0, dirStr.c_str(), &hKey);

       // link to odb
       sipm_folder &this_folder = sipm_folders[i];
       ret_code = db_open_record(hDB, hKey, &this_folder, sizeof(this_folder),
                              MODE_WRITE, NULL, NULL);
       if (ret_code != SUCCESS) {
       std::string errormsg =
          "failed to open record with error code " + std::to_string(ret_code);
       cm_msg(MERROR, __FUNCTION__, errormsg.c_str());
       return FE_ERR_ODB;
       }*/
  }
   
  //create odb keys for thread monitor                                                        
  std::string keyStr = "/Equipment/CaloSlowControls/ThreadMonitors/Crate" + std::to_string(iCrate+1);
  int thread_status = 0;
  db_set_value(hDB,0,keyStr.c_str(),&(thread_status),sizeof(thread_status), 1 ,TID_INT);

 }

  //Start Threads                                                                               
  mlock.lock();
  FrontendActive=true;
  RunActive=false;
  mlock.unlock();
  ReadingThreads.clear();
  for (int iCrate = 0 ; iCrate<NumberOfCrates ; iCrate++)
    {
      ReadingThreads.push_back(std::thread(ReadFromDevice,iCrate));
    }



  // read initial info
  if (!reload_all_info()) {
    return FE_ERR_HW;
  }

  bor_bk_folders = bk_folders;
  bor_sipm_folders = sipm_folders;

  return SUCCESS;
}

//--- Frontend Exit
//---------------------------------------------------------//
INT frontend_exit() { 
  mlock.lock();
  RunActive=false;
  FrontendActive=false;
  mlock.unlock();
  for (int iCrate = 0 ; iCrate<NumberOfCrates ; iCrate++)
    {
      ReadingThreads[iCrate].join();
    }
  cm_msg(MINFO,"exit","All threads joined.");

  return SUCCESS; 
}

//--- Begin of Run
//----------------------------------------------------------//
INT begin_of_run(INT run_number, char *error) {
  // turn on JMU
  //  if (!beagle_comm->turn_on_JMU()) {
  //  cm_msg(MERROR, __FILE__, "could not connect to BeagleBone!");
  //  return FE_ERR_HW;
  // }

  // reload all info at begin of run
<<<<<<< Updated upstream
  // if (!reload_all_info()) {
  //   return FE_ERR_HW;
  // }
  // will use last values read back
  // if major issue (changed variable, bad communication); will alarm at next
  // read attempt
=======
//  if (!reload_all_info()) {
//    return FE_ERR_HW;
//  }
// trust last read that's stored in ODB -- if something isn't communicating, it'll throw an error at the next read back 
>>>>>>> Stashed changes

  // copy begin of run values
  bor_sipm_folders = sipm_folders;
  bor_bk_folders = bk_folders;
  
  mlock.lock();
  RunActive = true;
  mlock.unlock(); 
 
  return SUCCESS;
}

//--- End of Run
//------------------------------------------------------------//
INT end_of_run(INT run_number, char *error) {
  // turn on JMU
  //  if (!beagle_comm->turn_on_JMU()) {
  //  cm_msg(MERROR, __FILE__, "could not connect to BeagleBone!");
  //  return FE_ERR_HW;
  //}

  // reload all info at end of run, too
<<<<<<< Updated upstream
  //if (!reload_all_info()) {
  //  return FE_ERR_HW;
  //}
  // keep the last read back value as the one stored in the ODB -- don't need
  // a new set 
=======
//  if (!reload_all_info()) {
//    return FE_ERR_HW;
//  }
// does not need to reload settings at end of run -- will continue to store last read
>>>>>>> Stashed changes

  mlock.lock();
  RunActive = false;
  mlock.unlock();


  return SUCCESS;
}

//--- Pause Run
//-------------------------------------------------------------//
INT pause_run(INT run_number, char *error) { return SUCCESS; }

//--- Resume Run
//------------------------------------------------------------//
INT resume_run(INT run_number, char *error) { return SUCCESS; }

//--- Frontend Loop -----------------------------------------------------//
INT frontend_loop() { 
  // try to limit CPU usage... what is midas doing?
  cm_yield(50);
  return SUCCESS; 
}

//---------------------------------------------------------------------------//

/****************************************************************************\

  Readout routines for different events
\*****************************************************************************/

//--- Trigger event routines
//------------------------------------------------//

INT poll_event(INT source, INT count, BOOL test)
 { 
   // fake calibration                                                                               
   if (test) {
     for (int i = 0; i < count; i++) {
       usleep(10);
     }
     return 0;
   }

   INT retval = 1;

   data_lock.lock();
   for (int i=0;i<NumberOfCrates;i++)
     {
       if (sipm_buffer[i].size()==0)
	 {
	   retval = 0;
	   break;
	 }
     }
   data_lock.unlock();
   return retval;

   //return SUCCESS;  
 }

//--- Interrupt configuration ---------------------------------------*/

INT interrupt_configure(INT cmd, INT source, PTYPE adr) {
  switch (cmd) {
    case CMD_INTERRUPT_ENABLE:
      break;
    case CMD_INTERRUPT_DISABLE:
      break;
    case CMD_INTERRUPT_ATTACH:
      break;
    case CMD_INTERRUPT_DETACH:
      break;
  }
  return SUCCESS;
}

//--- Event readout -------------------------------------------------*/

INT read_event(char *pevent, INT off) {
  // cm_msg(MDEBUG, __FILE__, "in read periodic event");

  // read values from beagle board
  /*  auto start = std::chrono::high_resolution_clock::now();
  //try twice to read beagle info before raising alarm
  if (!beagle_comm->read_beagle_info(bk_folders, sipm_folders) && !beagle_comm->read_beagle_info(bk_folders, sipm_folders)) {
    monitors.beagle_failure = 1;
    std::cout << ++n_consecutive_timeouts << " consecutive timeouts." << std::endl;
    if (n_consecutive_timeouts >= 10) {
      char alarm_message[80];
      sprintf(alarm_message, "Cannot communicate with Calo %d BeagleBone!",
	      frontend_index);

      int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
				      "communication failure", AT_INTERNAL);
      if (ret_code != AL_SUCCESS) {
	cm_msg(MERROR, __FILE__,
	       "Failure raising alarm. Error num %d, alarm message: %s", ret_code,
	       alarm_message);
      }
    }
    return 0;
  } else {
    check_and_raise_alarms();
    monitors.beagle_failure = 0;
    n_consecutive_timeouts = 0;
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::string timestr =
      "reading all values took " +
      std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start).count()) +
      " ms";
  std::cout << timestr << std::endl;

  // write to db
  start = std::chrono::high_resolution_clock::now();
  if (settings.write_temps) {
    if (!psql_comm->write_sipm_temps(sipm_folders)) {
      monitors.db_failure = 1;
      char alarm_message[80];
      sprintf(alarm_message, "Cannot write Calo %d values to database!",
              frontend_index);

      int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
                                      "database failure", AT_INTERNAL);
      if (ret_code != AL_SUCCESS) {
        cm_msg(MERROR, __FILE__,
               "Failure raising alarm. Error num %d, alarm message: %s",
               ret_code, alarm_message);
      }
      return 0;
    } else {
      monitors.db_failure = 0;
    }
  }
  if (settings.write_gains) {
    if (!psql_comm->write_sipm_gains(sipm_folders)) {
      monitors.db_failure = 1;
      return 0;
    } else {
      monitors.db_failure = 0;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  timestr =
      "queries took " +
      std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start).count()) +
      " ms";
  std::cout << timestr << std::endl;
  */
  if (settings.write_midas_bank) {
    // cm_msg(MDEBUG, __FILE__, "writing midas temperature bank");
    printf(" Writing midas temperature bank\n");
    data_lock.lock();
    for (int i=0;i<NumberOfCrates;i++)
      {
	INT frontend_index = i + 1;  // this corresponds to the old frontend index set by commandline
        char name[16];
        sprintf(name, "SI%02d", frontend_index);
        float *pdata;
        bk_init(pevent);
        bk_create(pevent, name, TID_FLOAT, (void **)&pdata);
        // write temps to the midas file
        //for (const auto &sipm : sipm_folders) {
        for (const auto &sipm : sipm_buffer[i][0]) { 
          *pdata++ = sipm.temp;
         }
        bk_close(pevent, pdata);
      }

    for (int i=0;i<NumberOfCrates;i++)
      {
	sipm_buffer[i].erase(sipm_buffer[i].begin());
      }
    data_lock.unlock();
    printf(" Done writing to midas temperature bank\n");

    return bk_size(pevent);
  }

  return 0;
}

//void check_and_raise_alarms() {
void check_and_raise_alarms(std::vector<sipm_folder> & sipm_folders, std::vector<bk_folder> & bk_folders){
  std::vector<int> sipm_nums;
  for (unsigned int i = 0; i < sipm_folders.size(); ++i) {
    // don't raise alarms for temperatures > 120 degrees
    // those are definitely communication failures
    if (temp_alarm_enabled[i] && sipm_folders[i].temp > settings.temp_al_thresh&& sipm_folders[i].temp < 120) {
      sipm_nums.push_back(i);
    }
  }
  if (!sipm_nums.empty()) {
    char alarm_message[80];
    if (sipm_nums.size() == 1) {
      sprintf(alarm_message, "SiPM %d in Calo %d has high temperature!",
              sipm_nums.front(), frontend_index);
    } else {
      sprintf(alarm_message, "Multiple high SiPM temperatures in Calo %d",
              frontend_index);
    }

    int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
                                    "high temperatures", AT_INTERNAL);
    if (ret_code != AL_SUCCESS) {
      cm_msg(MERROR, __FILE__,
             "Failure raising alarm. Error num %d, alarm message: %s", ret_code,
             alarm_message);
    }
  }

  std::vector<int> bk_nums;
  for (unsigned int i = 0; i < bk_folders.size(); ++i) {
    if (bk_folders[i].curr_out > settings.curr_al_thresh) {
      bk_nums.push_back(i);
    }
  }
  if (!bk_nums.empty()) {
    char alarm_message[80];
    if (bk_nums.size() == 1) {
      sprintf(alarm_message, "BK %d in Calo %d has high current!",
              bk_nums.front() + 1, frontend_index);
    } else {
      sprintf(alarm_message, "Multiple BKs in Calo %d have high currents!",
              frontend_index);
    }

    int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
                                    "high current", AT_INTERNAL);
    if (ret_code != AL_SUCCESS) {
      cm_msg(MERROR, __FILE__,
             "Failure raising alarm. Error num %d, alarm message: %s", ret_code,
             alarm_message);
    }
  }

  // rest of checks only occur if a run is in progress
  if (run_state == STATE_RUNNING) {
    assert(sipm_folders.size() == bor_sipm_folders.size());
    sipm_nums.clear();
    if (settings.pga_alarm) {
      // check if any SiPMs had their gain setting change
      for (unsigned int i = 0; i < sipm_folders.size(); ++i) {
        if (sipm_folders[i].gain != bor_sipm_folders[i].gain) {
          sipm_nums.push_back(i);
        }
      }
      if (!sipm_nums.empty()) {
        char alarm_message[80];
        if (sipm_nums.size() == 1) {
          sprintf(alarm_message,
                  "PGA setting for SiPM %d in Calo %d changed during the run!",
                  sipm_nums.front(), frontend_index);
        } else {
          sprintf(
              alarm_message,
              "PGA settings for multiple SiPMs in Calo %d changed during the "
              "run!",
              frontend_index);
        }

        int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
                                        "PGA changed", AT_INTERNAL);
        if (ret_code != AL_SUCCESS) {
          cm_msg(MERROR, __FILE__,
                 "Failure raising alarm. Error num %d, alarm message: %s",
                 ret_code, alarm_message);
        }
      }
    }

    assert(bk_folders.size() == bor_bk_folders.size());
    bk_nums.clear();
    if (settings.bk_off_alarm) {
      // check if any BKs turned off during the run
      std::vector<int> bk_nums;
      for (unsigned int i = 0; i < bk_folders.size(); ++i) {
        if (std::string(bor_bk_folders[i].output_stat) == "ON" &&
            std::string(bk_folders[i].output_stat) == "OFF") {
          bk_nums.push_back(i);
        }
      }
      if (!bk_nums.empty()) {
        char alarm_message[80];
        if (bk_nums.size() == 1) {
          sprintf(alarm_message, "BK %d in Calo %d turned off during the run!",
                  bk_nums.front() + 1, frontend_index);
        } else {
          sprintf(alarm_message,
                  "Multiple BKs in Calo %d turned off during the run!",
                  frontend_index);
        }

        int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm",
                                        "BK turned off", AT_INTERNAL);
        if (ret_code != AL_SUCCESS) {
          cm_msg(MERROR, __FILE__,
                 "Failure raising alarm. Error num %d, alarm message: %s",
                 ret_code, alarm_message);
        }
      }
    }
  }
}

//The actual read function                                                                       

void ReadFromDevice(int CrateIndex)
{
  std::vector<bk_folder> bk_folders(4);                                                            
  std::vector<sipm_folder> sipm_folders(54); 

  INT frontend_index = CrateIndex + 1;  // this corresponds to the old frontend index set by command line                                                                                        

  //Set thread monitor to 1                                                                      
  std::string keyStr = "/Equipment/uTCASlowControls/ThreadMonitors/Crate" + std::to_string(frontend_index);
 int thread_status = 1;
 mlock.lock();
 db_set_value(hDB,0,keyStr.c_str(),&(thread_status),sizeof(thread_status), 1 ,TID_INT);
 mlock.unlock();

 //Get hDB;                                                                                     
 HNDLE hDB, hKey;
 cm_get_experiment_database(&hDB, NULL);

 // init helpers                                                                                  
 std::string connection_addr =
   "tcp://192.168." + std::to_string(frontend_index) + ".21:6669";
 // d0 address                                                                                    
 //"tcp://192.168.1.21:6669";                                                                     
 // localhost for testing on my laptop                                                            
 // std::string connection_addr = "tcp://localhost:6669";                                         
 std::unique_ptr<CaloSC::BeagleCommunicator> beagle_comm = std::make_unique<CaloSC::BeagleCommunicator>(connection_addr);

 // try to turn on JMU supply                                                                     
 //  if (!beagle_comm->turn_on_JMU()) {                                                           
 //  cm_msg(MERROR, __FILE__, "could not connect to BeagleBone!");                                
 // return FE_ERR_HW;                                                                             
 //}                                                                                              

 std::string config_file = "config/dbconnection.json";
 std::unique_ptr<CaloSC::PSQLCommunicator> psql_comm = std::make_unique<CaloSC::PSQLCommunicator>(config_file, frontend_index);



 //psql comm handle                                                                             
 // std::unique_ptr<utcaSC::PSQLCommunicator> psql_comm = std::make_unique<utcaSC::PSQLCommunicator>(db_config_file, frontend_index);

 int ReadPeriodSecond = static_cast<int>(ReadPeriod);
 int ReadCycle = 1;// in s                                                                      

 int count = 0;
 while (1){
   bool localFrontendActive;
   bool localRunActive;
   mlock.lock();
   localFrontendActive = FrontendActive;
   localRunActive = RunActive;
   mlock.unlock();
   if (!localFrontendActive)break;

   if (count % (ReadPeriodSecond/ReadCycle) == 0)
     {
         auto start = std::chrono::high_resolution_clock::now();                                             //try twice to read beagle info before raising alarm                                              
         if (!beagle_comm->read_beagle_info(bk_folders, sipm_folders) && !beagle_comm->read_beagle_info(bk_folders, sipm_folders)) {                                                                            monitors.beagle_failure = 1;                                                                        std::cout << ++n_consecutive_timeouts << " consecutive timeouts." << std::endl;                                                             
        if (n_consecutive_timeouts >= 10) {                                                                     char alarm_message[80];                                                                             sprintf(alarm_message, "Cannot communicate with Calo %d BeagleBone!", frontend_index);                                                                                                                  int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm", "communication failure", AT_INTERNAL);                                                                                                                                
      if (ret_code != AL_SUCCESS) {                                                                                cm_msg(MERROR, __FILE__,"Failure raising alarm. Error num %d, alarm message: %s", ret_code, alarm_message);                                                                                    }                                                                                                 }                                                                                                   return 0;                                                                                         }
	 else {                                                                                                     check_and_raise_alarms(sipm_folders, bk_folders);                                                   monitors.beagle_failure = 0;                                                                        n_consecutive_timeouts = 0;                                                                        }                                                                                                  auto end = std::chrono::high_resolution_clock::now();                                               std::string timestr = "reading all values took " +   std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) + " ms";                                                std::cout << timestr << std::endl;

  // write to db                                                                                       
	 start = std::chrono::high_resolution_clock::now();                                                     if (settings.write_temps) {                                                                             if (!psql_comm->write_sipm_temps(sipm_folders)) {                                                    monitors.db_failure = 1;                                                                            char alarm_message[80];                                                                             sprintf(alarm_message, "Cannot write Calo %d values to database!", frontend_index);                                                                       
                 int ret_code = al_trigger_alarm("CaloSC", alarm_message, "CaloSC Alarm", "database failure", AT_INTERNAL);                              
	     if (ret_code != AL_SUCCESS) { cm_msg(MERROR, __FILE__, "Failure raising alarm. Error num %d, alarm message: %s", ret_code, alarm_message);  }                                                                    
                return 0;                                                                                          } 
          
          else {      
                 monitors.db_failure = 0;        
               }                                                                           
	 }                                                                                                  
	 if (settings.write_gains) {                                                                            if (!psql_comm->write_sipm_gains(sipm_folders)) {                                                        monitors.db_failure = 1;                                                                            return 0;                                                                                               } 
                    else{                                                                                                    monitors.db_failure = 0;                                                                             }                                                                                                 }                                                                                        end = std::chrono::high_resolution_clock::now();                                                   
	 timestr = "queries took " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>( end - start).count()) + " ms";                                                                        
         std::cout << timestr << std::endl;                                                              }
   sleep(ReadCycle);
   count++;
  }

 thread_status = 0;
 mlock.lock();
 db_set_value(hDB,0,keyStr.c_str(),&(thread_status),sizeof(thread_status), 1 ,TID_INT);
 mlock.unlock();

 return;

}
