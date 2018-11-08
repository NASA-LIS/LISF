#ifndef SAC_H
#define SAC_H

#include <iostream>
#include "ModelBase.hpp"
#include "PixelGraph.hpp"
#include "default_model_func.h"

using namespace Ohd_Hydro;

int sacFuncBefore( PixelGraph& g );

int sacFuncInside( PixelGraph& g );
int sacInCalLoop( PixelGraph& g );

DEFINE_MODEL( sac, 
	      sacFuncBefore, 
	      sacFuncInside, 
	      default_model_func,
              default_model_func /*sacInCalLoop*/ );

const ModelDataDescript sac::inputBeforeLoop[] = {
	//
	//
	// SAC_SMA parameters
	//
	//
	// first SAC_SMA parameter
	//
       {
	  "sac_UZTWM", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// second SAC_SMA parameter
	//
       {
	  "sac_UZFWM", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 3rd SAC_SMA parameter
	//
       {
	  "sac_UZK", //name
	  "1/deg",        //unit
	  "L",         //dimension
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
	// 4th SAC_SMA parameter
	//
       {
	  "sac_PCTIM", //name
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
	// 5th SAC_SMA parameter
	//
       {
	  "sac_ADIMP", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 6th SAC_SMA parameter
	//
       {
	  "sac_RIVA", //name
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
	// 7th SAC_SMA parameter
	//
       {
	  "sac_ZPERC", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 8th SAC_SMA parameter
	//
       {
	  "sac_REXP", //name
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
	// 9th SAC_SMA parameter
	//
       {
	  "sac_LZTWM", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 10th SAC_SMA parameter
	//
       {
	  "sac_LZFSM", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 11th SAC_SMA parameter
	//
       {
	  "sac_LZFPM", //name
	  "MM",        //unit
	  "L",         //dimension
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
	// 12th SAC_SMA parameter
	//
       {
	  "sac_LZSK", //name
	  "1/deg",        //unit
	  "L",         //dimension
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
	// 13th SAC_SMA parameter
	//
       {
	  "sac_LZPK", //name
	  "1/deg",        //unit
	  "L",         //dimension
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
	// 14th SAC_SMA parameter
	//
       {
	  "sac_PFREE", //name
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
	// 15th SAC_SMA parameter
	//
       {
	  "sac_SIDE", //name
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
	// 16th SAC_SMA parameter
	//
       {
	  "sac_RSERV", //name
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
	// 17th SAC_SMA parameter --- fraction of forest cover
	//
       {
	  "sac_EFC", //name
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
	 //PE
	 //
	//
	// Jan PE
	//
       {
	  "pe_JAN", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// Feb PE
	//
       {
	  "pe_FEB", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// MAR PE
	//
       {
	  "pe_MAR", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// APR PE
	//
       {
	  "pe_APR", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// MAY PE
	//
       {
	  "pe_MAY", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// JUN PE
	//
       {
	  "pe_JUN", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// JUL PE
	//
       {
	  "pe_JUL", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// AUG PE
	//
       {
	  "pe_AUG", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// SEP PE
	//
       {
	  "pe_SEP", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// OCT PE
	//
       {
	  "pe_OCT", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// NOV PE
	//
       {
	  "pe_NOV", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	// DEC PE
	//
       {
	  "pe_DEC", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
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
	 //PE Adj
	 //
	//
	// Jan PE adj
	//
       {
	  "peadj_JAN", //name
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
	// Feb PE adj
	//
       {
	  "peadj_FEB", //name
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
	// MAR PE adj
	//
       {
	  "peadj_MAR", //name
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
	// APR PE adj
	//
       {
	  "peadj_APR", //name
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
	// MAY PE adj
	//
       {
	  "peadj_MAY", //name
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
	// JUN PE adj
	//
       {
	  "peadj_JUN", //name
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
	// JUL PE adj
	//
       {
	  "peadj_JUL", //name
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
	// AUG PE adj
	//
       {
	  "peadj_AUG", //name
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
	// SEP PE adj
	//
       {
	  "peadj_SEP", //name
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
	// OCT PE adj
	//
       {
	  "peadj_OCT", //name
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
	// NOV PE adj
	//
       {
	  "peadj_NOV", //name
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
	// DEC PE adj
	//
       {
	  "peadj_DEC", //name
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
	 // SAC_SMA states in precentage
	 //
	//
	// 1st SAC_SMA state
	//
       {
	  "uztwc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
	//
	// 2nd SAC_SMA state
	//
       {
	  "uzfwc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
	//
	// 3rd SAC_SMA state
	//
       {
	  "lztwc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
	//
	// 4th SAC_SMA state
	//
       {
	  "lzfsc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
	//
	// 5th SAC_SMA state
	//
       {
	  "lzfpc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },
	//
	// 6th SAC_SMA state
	//
       {
	  "adimpc", //name
	  "PCTD",        //unit
	  "DLESS",         //dimension
	  "SASC",         //Data Type Code
          AVERAGE,     //accumulation policy
	  0,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          -1.f,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       },

         //
	 // SAC_SMA states in real value not precentage
	 //
	//
	// 1st SAC_SMA state
	//
       {
	  "real_uztwc", //name
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
       },
	//
	// 2nd SAC_SMA state
	//
       {
	  "real_uzfwc", //name
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
       },
	//
	// 3rd SAC_SMA state
	//
       {
	  "real_lztwc", //name
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
       },
	//
	// 4th SAC_SMA state
	//
       {
	  "real_lzfsc", //name
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
       },
	//
	// 5th SAC_SMA state
	//
       {
	  "real_lzfpc", //name
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
       },
	//
	// 6th SAC_SMA state
	//
       {
	  "real_adimpc", //name
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

const int sac::numOfInputBeforeLoop = 53;
const ModelDataDescript sac::inputInsideLoop[] = {
//                                 	"tair", "C", "T",
//	                                 AVERAGE,
//						       12,
//					 DO_FILLMISSING,
//						       40.f,
//					 USE_PREVIOUS,
//						       -99.f,
//						       "1:00:00"
                                	{"xmrg", "MM", "L", "MAPX",
	                                 SUM,
					 UNLIMITED_MISSING,
					 DONT_FILLMISSING,
						       0.f,
						       0.f,
						       -1.f,
						       "1:00:00" },
/*       {
	  "surf_water", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "MAP",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          USE_PREVIOUS,                 //missing value
          -1.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
       }
*/
						        };

const int sac::numOfInputInsideLoop = 1;

const ModelDataDescript sac::inputAfterLoop[] = {};
const int sac::numOfInputAfterLoop = 0;
#endif//#ifndef SAC_H
