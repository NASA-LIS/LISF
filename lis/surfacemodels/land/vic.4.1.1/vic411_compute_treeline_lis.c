//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_compute_treeline.c,v 1.3.2.1 2009/01/24 01:15:42 vicadmin Exp $";

void vic411_compute_treeline_lis(double                   avgJulyAirTemp,
                                 double                   *Tfactor,
                                 char                     *AboveTreeLine)
/**********************************************************************
  vic411_compute_treeline.c        Keith Cherkauer          March 12, 2003

  This routine computes the annual average July temperature for the 
  current gridcell.  The temperature is than lapsed to determine the 
  elevation at which the annual average temperature is equal to 10C.
  Snow elevation bands above this elevation are considered to be above
  the treeline.  When gridcell data is output at the end of each time 
  step, vegetation types with overstory will be excluded from the 
  variable averages of snow bands higher than the treeline (e.g. decidous
  trees will be removed from high elevation snow bands, while grass
  and shrubs will remain).  This is to serve as a preliminary fix for
  high elevation "glaciers", a more permanent vic411_version would actually 
  allow for vegetation types to be excluded from various snow bands.

  Modifications:
  2009-Jan-12 Added extra input parameter (avgJulyAirTemp) and logic to
	      handle it in the event JULY_TAVG_SUPPLIED is TRUE.		TJB
  2012-Jun-06 Modified by Shugong Wang for use in LIS-VIC implementation.
          Computation of average July air temperature has been removed. 
************************************************************************/
{

  extern vic411_option_struct       vic411_options;
  extern vic411_global_param_struct vic411_global_param;
  extern int                 vic411_NR, vic411_NF;

  double MonthSum;
  double AnnualSum;
  int    MonthCnt;
  int    AnnualCnt;
  int    rec;
  int    band;
  int    i;

  if (vic411_options.JULY_TAVG_SUPPLIED) 
  {
    // use supplied average annual July air temperature
    AnnualSum = avgJulyAirTemp;
  }
  else 
  {
    printf("Fatal error: \n");
    printf("Average July air temperature should be provided in VIC soil file!\n ");
    printf("LIS-VIC (4.1.1) doesn't support compuation of average July air temparature.\n");
    exit(1);
  }

  // Lapse average annual July air temperature to 10C and determine elevation
  for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) 
  {
    if ( AnnualSum + Tfactor[band] <= 10. )
      // Band is above treeline
      AboveTreeLine[band] = TRUE;
    else
      AboveTreeLine[band] = FALSE;
  }
}

