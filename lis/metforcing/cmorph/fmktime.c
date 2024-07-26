//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

/* the reverse of gmtime(), output the same as 'date -u -d MM/DD/YYY +%s' */
/* input yr, mon, day, hr, min, sec, output GMT epoch seconds */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ftn.h"


int FTN(fmktime)(ta)
  int ta[]; 
{

  struct tm time_str;

    time_str.tm_year = ta[5]; 
    time_str.tm_mon = ta[4]; 
    time_str.tm_mday = ta[3]; 
    time_str.tm_hour = ta[2]; 
    time_str.tm_min = ta[1]; 
    time_str.tm_sec = ta[0]; 
    time_str.tm_isdst = 0;
#if ( !defined AIX )
    time_str.tm_zone = "GMT"; 
#endif

    putenv("TZ=GMT");
    return mktime(&time_str); 

}



