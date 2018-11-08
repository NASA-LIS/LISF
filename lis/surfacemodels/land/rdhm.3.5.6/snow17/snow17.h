#ifndef SNOW17_H
#define SNOW17_H

#include <iostream>
#include "ModelBase.hpp"
#include "PixelGraph.hpp"
#include "default_model_func.h"

using namespace Ohd_Hydro;

int snow17FuncBeforeWrap( PixelGraph& g );

int snow17FuncInside( PixelGraph& g );

DEFINE_MODEL( snow17, 
	      snow17FuncBeforeWrap, 
	      snow17FuncInside, 
	      default_model_func,
              default_model_func );
#if 0
const ModelDataDescript snow17::inputBeforeLoop[] = {
	//
	//
	// SNOW 17parameters
	//
	//
	// first SNOW 17  parameter
	//
      {
	  "snow_ALAT", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// second SNOW 17 parameter
	//
      {
	  "snow_SCF", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 3rd SNOW 17 parameter, input as mm/hr but 'do_snow' converts to 
	 //                       mm/dthr
	//
      {
	  "snow_MFMAX", //name
	  "mm/DEGC/6hr",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 4th SNOW 17 parameter, input as mm/hr but 'do_snow' converts to 
	 //                       mm/dthr
	//
      {
	  "snow_MFMIN", //name
	  "mm/DEGC/6hr",        //unit
	  "L/D/T",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 5th SNOW 17 parameter
	//
      {
	  "snow_NMF", //name
	  "mm/DEGC/6hr",        //unit
	  "L/D/T",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 6th SNOW 17 parameter
	//
      {
	  "snow_UADJ", //name
	  "mm/mb",        //unit
	  "L/P",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 7th SNOW 17 parameter
	//
      {
	  "snow_SI", //name
	  "MM",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 8th SNOW 17 parameter
	//
      {
	  "snow_MBASE", //name
	  "DEGC",        //unit
	  "D",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 9th SNOW 17 parameter
	//
      {
	  "snow_PXTMP", //name
	  "DEGC",        //unit
	  "D",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 10th SNOW 17 parameter
	//
      {
	  "snow_PLWHC", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 11th SNOW 17 parameter
	//
      {
	  "snow_TIPM", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 12th SNOW 17 parameter
	//
      {
	  "snow_PGM", //name
	  "mm/day",        //unit
	  "L/T",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 13th SNOW 17 parameter
	//
      {
	  "snow_ELEV", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 14th SNOW 17 parameter
	//
      {
	  "snow_LAEC", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 15th SNOW 17 parameter
	//
      {
	  "snow_ADC1", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 16th SNOW 17 parameter
	//
      {
	  "snow_ADC2", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 17th SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC3", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 18th SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC4", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 19th SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC5", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 20th SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC6", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 21st SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC7", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 22nd SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC8", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 23nd SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC9", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 24nd SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC10", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//
	// 25nd SNOW 17 parameter --- fraction of forest cover
	//
      {
	  "snow_ADC11", //name
	  "N/A",        //unit
	  "DLESS",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
         //
	 // SNOW 17 states 
	 //
	//
	// 1st SNOW 17 state
	//
      {
	  "we", //name
	  "MM",        //unit
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
	// 2nd SNOW 17 state
	//
      {
	  "neghs", //name
	  "MM",        //unit
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
	// 3rd SNOW 17 state
	//
      {
	  "liqw", //name
	  "MM",        //unit
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
	// 4th SNOW 17 state
	//
      {
	  "tindex", //name
	  "DEGC",        //unit
	  "D",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
      },
	//
	// 5th SNOW 17 state
	//
      {
	  "accmax", //name
	  "MM",        //unit
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
	// 6th SNOW 17 state
	//
      {
	  "sndpt", //name
	  "cm",        //unit
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
	// 7th SNOW 17 state
	//
      {
	  "sntmp", //name
	  "DEGC",        //unit
	  "D",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -99.f,                 // no data value
         "1:00:00"     //time interval
      },
	//
	// 8th SNOW 17 state
	//
      {
	  "sb", //name
	  "MM",        //unit
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
	// 9th SNOW 17 state
	//
      {
	  "sbaesc", //name
	  "fraction",        //unit
	  "N/A",         //dimension
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
	// 10th SNOW 17 state
	//
      {
	  "sbws", //name
	  "MM",        //unit
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
	// 11th SNOW 17 state
	//
      {
	  "storge", //name
	  "MM",        //unit
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
	// 12th SNOW 17 state
	//
      {
	  "aeadj", //name
	  "MM",        //unit
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
	// 13th SNOW 17 state
	//
      {
	  "sxlag1", //name
	  "MM",        //unit
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
	// 14th SNOW 17 state
	//
      {
	  "sxlag2", //name
	  "MM",        //unit
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
	// 15th SNOW 17 state
	//
      {
	  "sxlag3", //name
	  "MM",        //unit
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
	// 16th SNOW 17 state
	//
      {
	  "sxlag4", //name
	  "MM",        //unit
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
	// 17th SNOW 17 state
	//
      {
	  "sxlag5", //name
	  "MM",        //unit
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
	// 18th SNOW 17 state
	//
      {
	  "sxlag6", //name
	  "MM",        //unit
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
	// 19th SNOW 17 state
	//
      {
	  "sxlag7", //name
	  "MM",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          AVERAGE,     //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          0.f,                 //missing value
          0.f,                 //missing value for the first time step
         -1.f,                 // no data value
         "1:00:00"     //time interval
      }
	       };

const int snow17::numOfInputBeforeLoop = 44;
const ModelDataDescript snow17::inputInsideLoop[] = {
                                     {
                                	"xmrg", "MM", "L", "MAPX",
	                                 SUM,
					 UNLIMITED_MISSING,
					 DONT_FILLMISSING,
						       0.f,
						       0.f,
						       -1.f,
						       "1:00:00"
                                      },
                                      {
                                 	"tair", "DEGF", "T", "TEMP",
	                                 AVERAGE,
						       12,
					 DO_FILLMISSING,
						       40.f,
					 USE_PREVIOUS,
						       -99.f,
						       "1:00:00"
                                     },
                                     {
                                	"psfrac", "N/A", "DLESS",  "TEMP",
	                                 AVERAGE,
					 UNLIMITED_MISSING,
					 DONT_FILLMISSING,
						       -1.f,
						       -1.f,
						       -1.f,
						       "1:00:00"
                                     }
						        };

const int snow17::numOfInputInsideLoop = 3;

const ModelDataDescript snow17::inputAfterLoop[] = {};

const int snow17::numOfInputAfterLoop = 0;
#endif //#if 0

#endif//#ifndef SNOW17_H
