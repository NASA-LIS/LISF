#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_compute_treeline.c,v 1.3.2.1 2009/01/24 01:15:42 vicadmin Exp $";

void vic411_compute_treeline(vic411_atmos_data_struct        *atmos,
                      vic411_dmy_struct               *dmy,
		      double                   avgJulyAirTemp,
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

  if (vic411_options.JULY_TAVG_SUPPLIED) {

    // use supplied average annual July air temperature
    AnnualSum = avgJulyAirTemp;

  }
  else {

    // compute average annual July air temperature from forcing
    AnnualSum = 0;
    AnnualCnt = 0;
    rec = 0;
    while ( rec < vic411_global_param.nrecs ) {
      if ( dmy[rec].month == 7 ) {
        MonthSum = 0;
        MonthCnt = 0;
        while ( dmy[rec].month == 7 ) {
          for (i = 0; i < vic411_NF; i++) {
            MonthSum += atmos[rec].air_temp[i];
            MonthCnt++;
          }
          rec++;
        }
        if ( MonthCnt > 0 ) {
          // Sum monthly average July temperature
          AnnualSum += MonthSum / (double)MonthCnt;
          AnnualCnt++;
        }
      }
      rec++;
    }

    // Compute average annual July air temperature
    if ( AnnualCnt > 0 )
      AnnualSum /= (double)AnnualCnt;

  }

  // Lapse average annual July air temperature to 10C and determine elevation
  for ( band = 0; band < vic411_options.SNOW_BAND; band++ ) {
    if ( AnnualSum + Tfactor[band] <= 10. )
      // Band is above treeline
      AboveTreeLine[band] = TRUE;
    else
      AboveTreeLine[band] = FALSE;
  }

}

