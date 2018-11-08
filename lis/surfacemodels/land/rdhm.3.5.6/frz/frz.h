#ifndef FRZ_H
#define FRZ_H

#include <iostream>
#include "ModelBase.hpp"
#include "PixelGraph.hpp"
#include "default_model_func.h"

using namespace Ohd_Hydro;

int frzFuncBefore( PixelGraph& g );
int frzInCalLoop( PixelGraph& g );

DEFINE_MODEL( frz, 
	      frzFuncBefore, 
	      default_model_func, 
	      default_model_func,
              default_model_func /*frzInCalLoop*/ );

const ModelDataDescript frz::inputBeforeLoop[] = {
	//
	//
	// Frozen ground parameters
	//
	//
	// 1st Frozen ground parameter
	//   
	//   soil texture
	//
     {
	  "frz_STXT", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DO_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// second Frozen Ground parameter
	//
     {
	  "frz_TBOT", //name
	  "K",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  0,           //max number of missing
	  DO_FILLMISSING,    //if fill missing
          60.f,                 //missing value
          60.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 3rd Frozen Ground parameter
	//
     {
	  "frz_RSMAX", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 4th Frozen Ground parameter
	//
     {
	  "frz_CKSL", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 5th Frozen Ground parameter
	//
     {
	  "frz_ZBOT", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 6th Frozen Ground parameter
	 // calculated
	//
     {
	  "frz_RTUP", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 7th Frozen Ground parameter
	//
	 //calculated
     {
	  "frz_RTLW", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 8th Frozen Ground parameter
	//
     {
	  "frz_PSISAT", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 9th Frozen Ground parameter
	//
     {
	  "frz_SWLT", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 10th Frozen Ground parameter
	//
     {
	  "frz_Z0", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 11th Frozen Ground parameter
	//
     {
	  "frz_Z1", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 12th Frozen Ground  parameter
	//
     {
	  "frz_Z2", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 13th Frozen Ground  parameter
	//
     {
	  "frz_Z3", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 14th Frozen Ground  parameter
	//
     {
	  "frz_Z4", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 15th Frozen Ground  parameter
	//
     {
	  "frz_nsoil", //name, number of soil layers
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 16th Frozen Ground  parameter
	//
     {
	  "frz_nupl", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
	//
	// 17th Frozen Ground  parameter
	//
     {
	  "frz_nsac", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -99.f,                 //missing value
          -99.f,                 //missing value for the first time step
         -99.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
     },
         //
	 // Frozen Ground states
	 //
	//
	// 1st Frozen Ground state
	//
     {
	  "ts0", //name
	  "C",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 2nd Frozen Ground state
	//
     {
	  "ts1", //name
	  "C",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 3rd Frozen Ground state
	//
     {
	  "ts2", //name
	  "C",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 4th Frozen Ground state
	//
     {
	  "ts3", //name
	  "C",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 5th Frozen Ground state
	//
     {
	  "ts4", //name
	  "C",        //unit
	  "T",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 6th Frozen Ground state
	//
     {
	  "uztwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 7th Frozen Ground state
	//
     {
	  "uzfwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 8th Frozen Ground state
	//
     {
	  "lztwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 9th Frozen Ground state
	//
     {
	  "lzfsh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 10th Frozen Ground state
	//
     {
	  "lzfph", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 6th Frozen Ground state real value
	//
     {
	  "real_uztwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 7th Frozen Ground state real value
	//
     {
	  "real_uzfwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 8th Frozen Ground state real value
	//
     {
	  "real_lztwh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 9th Frozen Ground state real value
	//
     {
	  "real_lzfsh", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 10th Frozen Ground state real value
	//
     {
	  "real_lzfph", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,                 //missing value
          DEFAULT_MISSING_VALUE,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 1st SAC_SMA state previous
	//
     {
	  "uztwc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 2nd SAC_SMA state previous
	//
     {
	  "uzfwc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 3rd SAC_SMA state previous
	//
     {
	  "lztwc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 4th SAC_SMA state previous
	//
     {
	  "lzfsc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 5th SAC_SMA state previous
	//
     {
	  "lzfpc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },
	//
	// 6th SAC_SMA state previous
	//
     {
	  "adimpc_prv", //name
	  "mm",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
     },

     {
	  "smc0", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },

     {
	  "smc1", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "smc2", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "smc3", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },

     {
	  "smc4", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "smc5", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "sh2o0", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },

     {
	  "sh2o1", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "sh2o2", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },

     {
	  "sh2o3", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "sh2o4", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
     {
	  "sh2o5", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
       {
	  "surf_water", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "MAP",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
       {
	  "surf_water_ratio", //name
	  "MM",        //unit
	  "L",         //dimension
	  "MAP",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       }
	       };

const int frz::numOfInputBeforeLoop = 52;
const ModelDataDescript frz::inputInsideLoop[] = {
                                 	{"tair", "DEGF", "T", "TEMP",
	                                 AVERAGE,
						       12,
					 DO_FILLMISSING,
					 USE_PREVIOUS,
                                                       40.f,
						       -99.f,
						       "1:00:00"
//                                	"xmrg", "mm", "L", "MAPX",
//	                                 SUM,
//					 UNLIMITED_MISSING,
//					 DONT_FILLMISSING,
//						       0.f,
//						       0.f,
//						       -1.f,
//						       "1:00:00"
                                   },
       {
	  "surf_water", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "MAP",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          USE_INIT_STATE,        //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       }
						        };

const int frz::numOfInputInsideLoop = 2;

const ModelDataDescript frz::inputAfterLoop[] = {};
const int frz::numOfInputAfterLoop = 0;
#endif//#ifndef frz
