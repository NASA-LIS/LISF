#ifndef SACHTET_H
#define SACHTET_H

#include <iostream>
#include "ModelBase.hpp"
#include "PixelGraph.hpp"
#include "default_model_func.h"

using namespace Ohd_Hydro;

int sachtetFuncBefore( PixelGraph& g );

int sachtetFuncInside( PixelGraph& g );
//int sachtetInCalLoop( PixelGraph& g );

DEFINE_MODEL( sachtet, 
	      sachtetFuncBefore, 
	      sachtetFuncInside, 
	      default_model_func,
              default_model_func /*sacInCalLoop*/ );

const ModelDataDescript sachtet::inputBeforeLoop[] = {
	//
	//
	// SAC_HZTET parameters
	//
	//
	// first SAC_HTET parameter
	//

/////////////////////////////////////////////////////
//Vegetation dependent parameter grids:
///////////////////////////////////////////////////
//RCMIN minimal stomatal resistance, s/m
//
{
      "sachtet_RCMIN",			//name
      "S/M",	  	                //unit
      "T/L",				//dimension
      "SQIN",				//DataTypeCode
      UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
//      5000.0f,				//missing value
//      5000.0f,				//missing value for the first time step
      DEFAULT_MISSING_VALUE,	       //missing value
      DEFAULT_MISSING_VALUE,           //missing value for the first time step
      DEFAULT_MISSING_VALUE,     	//no data value
      NO_TIME_INTERVAL			//time interval
},
	// 2nd SAC_HTET parameter
//RGL,solar radiation threshold for which resistance factor Rsr
//is about to double its minimum value
{
      "sachtet_RGL",			//name
      "W/M2",				//unit
      "E/L2",				//dimension
      "SQIN",				//DataTypeCode
      UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
      NO_TIME_INTERVAL			//time interval
},
	// 3rd SAC_HTET parameter
//HS,parameter in the vapor pressure resistance factor
{
      "sachtet_HS",			//name
      "PCTD",				//unit
      "DLES",				//dimension
      "PTPS",				//DataTypeCode
      UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
      NO_TIME_INTERVAL			//time interval
},
	// 4th SAC_HTET parameter
//                                                
//ZR0(Z0),the roughness length, m
{
      "sachtet_ZR0",			//name
      "M",				//unit
      "L",				//dimension
      "SQIN",				//DataTypeCode
      UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
      NO_TIME_INTERVAL			//time interval
},
	// 5th SAC_HTET parameter
//LAI, the leaf area index presently set to universal value of 5.0
{
      "sachtet_LAI",		        //name
      "PCTD",				//unit
      "DLES",				//dimension
      "PTPS",				//DataTypeCode
      UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
      NO_TIME_INTERVAL			//time interval
},
	// 6th SAC_HTET parameter
//                                                                              //D50, the depth(cm) at which 50% roots are allocated
{
      "sachtet_D50",		        //name
      "CM",				//unit
      "L",				//dimension
      "SQIN",				//DataTypeCode
      UNKNOWN,				//accumulation policy
     UNLIMITED_MISSING,   //max number of missing
     DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,	       //missing value
      DEFAULT_MISSING_VALUE,           //missing value for the first time step
      DEFAULT_MISSING_VALUE,     	//no data value
      NO_TIME_INTERVAL			//time interval
},
	// 7th SAC_HTET parameter
//RPOW(rn),a dimensionless shape-parameter of root distribution &#x2198;
//function
{
     "sachtet_CROOT",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     UNLIMITED_MISSING,   //max number of missing
     DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,	       //missing value
      DEFAULT_MISSING_VALUE,           //missing value for the first time step
      DEFAULT_MISSING_VALUE,     	//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 8th SAC_HTET parameter
/////////////////////////////////////////////////////
//Soil texture dependent parameter grids
/////////////////////////////////////////////////////
//
//dksat, soil moisture conductivity
//
{
     "sachtet_DKSAT",			//name
     "M/S",				//unit
     "L/T",				//dimension
     "SQIN",				//DataTypeCode
     UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
     NO_TIME_INTERVAL			//time interval
},
	// 9th SAC_HTET parameter
//
//dwsat, soil moisture diffusivity
//
{
     "sachtet_DWSAT",			//name
     "M2/S",				//unit
     "L2/T",				//dimension
     "SQIN",				//DataTypeCode
     UNKNOWN,				//accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
     NO_TIME_INTERVAL			//time interval
},
//////////////////////////////////////////////////////////
//Not vegetation or soil texture dependent parameters:
//////////////////////////////////////////////////////////
//
	// 10th SAC_HTET parameter
//soil Albedo, soilsurface albedo
//missing=0.15
//
{
     "sachtet_soilAlbedo",		//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     UNLIMITED_MISSING,	        //max number of missing
     DONT_FILLMISSING,			//if fill missing
     0.15f,				//missing value
     0.15f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 11th SAC_HTET parameter
//
//snowAlbedo, soil surface albedo
//missing=0.7
//
{
     "sachtet_snowAlbedo",		//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     UNLIMITED_MISSING,          //max number of missing
     DONT_FILLMISSING,			//if fill missing
     0.7f,				//missing value
     0.7f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 12th SAC_HTET parameter
//12 Grinness fraction Jan to Dec
//
//
//Monthly grinness fraction grids: Jan
//
{
     "sachtet_GRN_Jan",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 13th SAC_HTET parameter
//
//Monthly grinness fraction grids: Feb
//
{
     "sachtet_GRN_Feb",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 14th SAC_HTET parameter
//
//Monthly grinness fraction grids: Mar 
//
{
     "sachtet_GRN_Mar",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 15th SAC_HTET parameter
//
//Monthly grinness fraction grids: Apr 
//
{
     "sachtet_GRN_Apr",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 16th SAC_HTET parameter
//
//Monthly grinness fraction grids: May 
//
{
     "sachtet_GRN_May",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 17th SAC_HTET parameter
//
//Monthly grinness fraction grids: Jun 
//
{
     "sachtet_GRN_Jun",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 18th SAC_HTET parameter
//
//Monthly grinness fraction grids: Jul 
//
{
     "sachtet_GRN_Jul",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 19th SAC_HTET parameter
//
//Monthly grinness fraction grids: Aug
//
{
     "sachtet_GRN_Aug",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 20st SAC_HTET parameter
//
//Monthly grinness fraction grids: Sep 
//
{
     "sachtet_GRN_Sep",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 21nd SAC_HTET parameter
//
//Monthly grinness fraction grids: Oct
//
{
     "sachtet_GRN_Oct",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 22nd SAC_HTET parameter
//
//Monthly grinness fraction grids: Nov 
//
{
     "sachtet_GRN_Nov",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	// 23rd SAC_HTET parameter
//
//Monthly grinness fraction grids: Dec
//
{
     "sachtet_GRN_Dec",			//name
     "PCTD",				//unit
     "DLES",				//dimension
     "PTPS",				//DataTypeCode
     UNKNOWN,				//accumulation policy
     0,					//max number of missing
     DO_FILLMISSING,			//if fill missing
     -1.0f,				//missing value
     -1.0f,				//missing value for the first time step
     -1.0f,				//no data value
     NO_TIME_INTERVAL			//time interval
},
	//24th SAC_HTET parameter
////////////////////////////////////////////////////////
//vegetation class. option: -1 variable, other input constant
//veg_type available grid;
////////////////////////////////////////////////////////
//
//
{
    "sachtet_veg_type",		//name
    "N/A",				//unit
    "DLES",				//dimension
    "PTPS",				//DataType Code
    UNKNOWN,				//accumulation policy
    0,					//max number of missing
    DONT_FILLMISSING,			//if fill missing
      DEFAULT_MISSING_VALUE,		//missing value
      DEFAULT_MISSING_VALUE,		//missing value for the first time step
      DEFAULT_MISSING_VALUE,            //no data value
    NO_TIME_INTERVAL			//time interval
},
////////////////////////////////////////////////
// Globally constant
	//25th SAC_HTET parameter
/////////////////////////////////////////////////////////
//pc: plant coef. default pc = -1, 0.6 - 0.8
////////////////////////////////////////////////////////
//
//{
//    "sachtet_pc",		//name
//    "N/A",				//unit
//    "DLES",				//dimension
//    "PTPS",				//DataType Code
//    UNKNOWN,				//accumulation policy
//    0,					//max number of missing
//    DO_FILLMISSING,			//if fill missing
//    -1.0f,				//missing value
//    -1.0f,				//missing value for the first time step
//    -1.0f,				//no data value
//    NO_TIME_INTERVAL			//time interval
//},
	//26th SAC_HTET parameter
////////////////////////////////////////////////////////
//Internal estimated variable which may vary in space:
////////////////////////////////////////////////////////
//offsetTime, localtime offset from Z-time
//CALCULATED
//
{
    "sachtet_offsetTime",		//name
    "HR",				//unit
    "T",				//dimension
    "TIME",				//DataType Code
    UNKNOWN,				//accumulation policy
    0,					//max number of missing
    DONT_FILLMISSING,			//if fill missing
    -1.0f,				//missing value
    -1.0f,				//missing value for the first time step
    -1.0f,				//no data value
    NO_TIME_INTERVAL			//time interval
},
	//27th SAC_HTET parameter
      {
	  "snow_ELEV", //name
	  "m",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          0,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          DEFAULT_MISSING_VALUE,	//missing value
          DEFAULT_MISSING_VALUE,        //missing value for the first time step
          DEFAULT_MISSING_VALUE,        //no data value
         NO_TIME_INTERVAL     //time interval
      },

	//28th SAC_HTET parameter
//
// refernece wind speed for PET adj (parameter), replacement for reference
// wind speed grid grid.
//
      {
	  "sachtet_sfcref", //name
	  "M/S",        //unit
	  "L/T",         //dimension
	  "UAVG",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,	//max number of missing
          DO_FILLMISSING,		//if fill missing
          4,	                        //missing value
          4,	                        //missing value for the first time step
          -1.f,        //no data value
          NO_TIME_INTERVAL     //time interval
      },

	//29th SAC_HTET parameter
//
//air temperature measurement hight, m; default is 2 m
//variable at each gauge; assume a country-wide constant?
//
      {
	  "sachtet_ztmp", //name
	  "M",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
      UNLIMITED_MISSING,		//max number of missing
      DONT_FILLMISSING,			//if fill missing
          2.f,                 //missing value
          2.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//30th SAC_HTET parameter
//
//wind speed measurement hight, m; default is 10 m
//variable at each gauge; assume a country-wide constant?
//
      {
	  "sachtet_zspd", //name
	  "M",        //unit
	  "L",         //dimension
	  "SQIN",         //Data Type Code
          UNKNOWN,     //accumulation policy
          UNLIMITED_MISSING,           //max number of missing
	  DONT_FILLMISSING,    //if fill missing
          10.f,                 //missing value
          10.f,                 //missing value for the first time step
         -1.f,                 // no data value
         NO_TIME_INTERVAL     //time interval
      },
	//31st SAC_HTET parameter
// //
 // noonTime, noontime, constant in US
 // 
 //
      { 
         "sachtet_noonTime", //name
	  "HR",               //unit
	  "T",           //dimension
	  "TIME",     //dataTypeCode
          UNKNOWN,       //accumulation policy
	  UNLIMITED_MISSING,   //max number of missing
	 DONT_FILLMISSING,    //if fill missing
         12.f,         //missing value
         12.f,    //missing value for the first time step
         -1.f,           // no data value
         NO_TIME_INTERVAL  //time interval
      },
	//32nd SAC_HTET parameter
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
	// 1st SAC_HTET state
//
//Spatiallyvariable,state
//chisdefinedinfland1.itmayberecalculated
//infland1ateachtimestep
{
    "sachtet_ch",			//name
    "PCTD",				//unit
    "DLES",				//dimension
    "PTPS",				//DataType Code
    UNKNOWN,				//accumulation policy
    UNLIMITED_MISSING,           //max number of missing
    DONT_FILLMISSING,			//if fill missing
    1.e-4,				//missing value
    1.e-4,				//missing value for the first time step
    -1.f,				//no data value
    "1:00:00"    		        //time interval
},
	//2nd SAC_HTET state
{
    "sachtet_cm",			//name
    "PCTD",				//unit
    "DLES",				//dimension
    "PTPS",				//DataType Code
    UNKNOWN,				//accumulation policy
    UNLIMITED_MISSING,           //max number of missing
    DONT_FILLMISSING,			//if fill missing
    1.e-4,				//missing value
    1.e-4,				//missing value for the first time step
    -1.f,				//no data value
    "1:00:00"    		        //time interval
},
	//3rd SAC_HTET state
//
//windadjp(in fland1) is a state global variable grid with default
//initialvalues=1.0: it may be recalculated in time loop in '&#x2198;
// fland1'
{
    "sachtet_windadjpc",			//name
    "PCTD",				//unit
    "DLES",				//dimension
    "PTPS",				//DataTypeCode
    UNKNOWN,				//accumulation policy
    UNLIMITED_MISSING,           //max number of missing
    DONT_FILLMISSING,			//if fill missing
    1.0,				//missing value
    1.0,				//missing value for the first time step
    -1.f,				//no data value
    "1:00:00"			        //time interval
},
//////////////////////////////////////
//     Globally constant
//
	//3rd SAC_HTET state
//////////////////////////////////////	
//
//bareadj – Ek-Chen evaporation threshold switch
//changed bare soil evaporation options depending on greenness:
//-bareadj is switch direct ET from Ek to Chen depending
//on greenness: if grn < bareadj, Ek option;
//if grn > bareadj, Chen option
//{
//    "sachtet_BAREADJ",			//name
//    "PCTD",				//unit
//    "DLES",				//dimension
//    "PTPS",				//DataType Code
//    UNKNOWN,				//accumulation policy
//    UNLIMITED_MISSING,           //max number of missing
//    DO_FILLMISSING,			//if fill missing
//    0.23f,				//missing value
//    0.23f,				//missing value for the first time step
//    -1.f,				//no data value
//    "1:00:00"	    		        //time interval
//},
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
       },
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
         -999.f,                 // no data value
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
         -999.f,                 // no data value
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
         -999.f,                 // no data value
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
         -999.f,                 // no data value
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
         -999.f,                 // no data value
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

const int sachtet::numOfInputBeforeLoop = 139; // = 34 + 53 + 52
const ModelDataDescript sachtet::inputInsideLoop[] = {
 //////////////////////////////////////////////////
 //SAC-HTET Input data
 //////////////////////////////////////////////////
 //daily maximum air temperature
 {
    "tairmax",				//name
    "DEGF",				//unit
    "Temp",				//dimension
    "SQIN",				//DataType Code
    UNKNOWN,				//accumulation policy
    0,					//max number of missing
    DO_FILLMISSING,			//if fill missing
    USE_PREVIOUS,				//missing value
    -999.0f,				//missing value for the first time step
    -999.0f,				//no data value
    "24:00:00"				//time interval
 },
 //daily minimum air temperature
 {
    "tairmin",				//name
    "DEGF",				//unit
    "TEMP",				//dimension
    "TAMX",				//DataType Code
    UNKNOWN,				//accumulation policy
    0,					//max number of missing
    DO_FILLMISSING,			//if fill missing
    USE_PREVIOUS,				//missing value
    -999.0f,				//missing value for the first time step
    -999.0f,				//no data value
    "24:00:00"				//time interval
 },
 //wind speed(at simulation time step) if available
 //if wind speed not available ,reference wind speed will be used
 {
    "windspeed",			//name
    "M/S",				//unit
    "L/T",				//dimension
    "UAVG",				//DataType Code
    UNKNOWN,				//accumulation policy
    UNLIMITED_MISSING,			//max number of missing
    DO_FILLMISSING,			//if fill missing
    -1.0f,				//missing value
    -1.0f,				//missing value for the first time step
    -1.0f,				//no data value
    "1:00:00"				//time interval
 },
 {
    "xmrg", 
    "MM", 
    "L", 
    "MAPX",
    SUM,
    UNLIMITED_MISSING,
    DONT_FILLMISSING,
    0.f,
    0.f,
    -1.f,
    "1:00:00" 
},
                                      {
                                 	"tair", "DEGF", "TEMP", "MAT",
	                                 AVERAGE,
						       12,
					 DO_FILLMISSING,
						       USE_PREVIOUS,
					 40.f,
						       -99.f,
						       "1:00:00"
                                     }
 
};

const int sachtet::numOfInputInsideLoop = 5;

const ModelDataDescript sachtet::inputAfterLoop[] = {};
const int sachtet::numOfInputAfterLoop = 0;
#endif//#ifndef SACHTET_H
