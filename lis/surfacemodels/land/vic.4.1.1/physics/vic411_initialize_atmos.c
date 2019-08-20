#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_initialize_atmos.c,v 5.9.2.7 2009/10/13 21:35:10 vicadmin Exp $";

void vic411_initialize_atmos(vic411_atmos_data_struct        *atmos,
                      vic411_dmy_struct               *dmy,
		      FILE                    **infile,
                      double                    theta_l,
                      double                    theta_s,
                      double                    phi,
		      double                    elevation,
		      double                    annual_prec,
		      double                    wind_h,
		      double                    roughness,
		      double                    avgJulyAirTemp,
		      double                   *Tfactor,
#if OUTPUT_FORCE
                      char                     *AboveTreeLine,
                      vic411_out_data_file_struct     *out_data_files,
                      vic411_out_data_struct          *out_data)
#else /* OUTPUT_FORCE */
                      char                     *AboveTreeLine)
#endif /* OUTPUT_FORCE */
/**********************************************************************
  vic411_initialize_atmos	Keith Cherkauer		February 3, 1997

  This routine initializes atmospheric variables for both the model
  time step, and the time step used by the snow algorithm (if different).
  Air temperature is estimated using MTCLIM (see routine for reference),
  atmospheric moisture is estimated using Kimball's algorithm (see 
  routine for reference), and radiation is estimated using Bras's algorithms
  (see routines for reference).

  WARNING: This subroutine is site specific.  Location parameters
    must be changed before compilation.

  UNITS: mks
	energy - W/m^2

  Modifications:
  11-18-98  Removed variable array yearly_epot, since yearly potential
            evaporation is no longer used for estimating the dew
            point temperature from daily minimum temperature.   KAC
  11-25-98  Added second check to make sure that the difference 
            between tmax and tmin is positive, after being reset
            when it was equal to 0.                        DAG, EFW
  12-1-98   Changed relative humidity computations so that they 
            use air temperature for the time step, instead of average
            daily temperature.  This allows relative humidity to
            change during the day, when the time step is less than
            daily.                                              KAC
  8-19-99   MIN_TDEW was added to prevent the dew point temperature
            estimated by Kimball's equations from becoming so low
            that vic411_svp() fails.							Bart
  9-4-99    Code was largely rewritten to change make use of the MTCLIM
            meteorological preprocessor which estimates sub-daily 
	    met forcings for all time steps.  The vic411_atmos_data_struct was
	    also reconfigured so that it has a new record for each
	    model time step, but stores sub-time step forcing data
	    (that might be needed for the snow model) within each
	    record, eliminating the on the fly estimations used in
	    previous versions of the model.					Bart and Greg
  01-17-01  Pressure and vapor pressure read from a forcing file are
            converted from kPa to Pa.  This preserves the original
            format of the forcing files (where pressure was supposed 
            to be in kPa, but allows VIC to use Pa internally, eliminating
            the need to convert to Pa every time it is used.			KAC
  03-12-03 Modifed to add AboveTreeLine to vic411_soil_con_struct so that
           the model can make use of the computed treeline.			KAC
  04-Oct-04 Changed logic to allow VP to be supplied without
	    SHORTWAVE.								TJB
  2005-Mar-24 Modified to handle ALMA forcing variables.			TJB
  2005-Apr-30 Fixed typo in QAIR calculation.					TJB
  2005-May-01 Added logic for CSNOWF and LSSNOWF.				TJB
  2005-May-02 Added logic for WIND_E and WIND_N.				TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
	      out_data and out_data_files structures.				TJB
  2006-Dec-20 Replaced 1000.0 with kPa2Pa in pressure conversion.		TJB
  2006-Dec-29 Added REL_HUMID to the list of supported met input variables.	TJB
  2007-Jan-02 Added ALMA_INPUT option; removed TAIR and PSURF from list of
	      supported met input variables.					TJB
  2008-Jan-25 Fixed conditions under which net longwave replaces incoming
	      longwave in atmos[rec].longwave[vic411_NR].  Previously, net longwave
	      was stored if SNOW_STEP != global.dt.  Now, net longwave is
	      stored if vic411_options.FULL_ENERGY and vic411_options.FROZEN_SOIL are both
	      FALSE, i.e. for a water balance mode run.				TJB
  2009-Jan-12 Modified to pass avgJulyAirTemp argument to
	      vic411_compute_treeline(). 						TJB
  2009-May-18 Added vic411_options.PLAPSE, which when TRUE changes pressure
	      calculation to be a function of elevation and air temperature
	      (as opposed to a constant 95.5 kPa, as it was previously).
	      Made similar change to density calculation.			TJB
  2009-Jun-10 Fixed incorrect handling of cases when incoming longwave and
	      shortwave radiation are supplied.					TJB
  2009-Jul-26 Removed the special logic for the water balance mode, in
	      which net longwave is stored in the "longwave" variable.		TJB
  2009-Oct-13 Removed condition if(vic411_options.SNOW_BAND) for call to
	      vic411_compute_treeline(), since vic411_options.SNOW_BAND is always > 0.	TJB

**********************************************************************/
{
  extern vic411_option_struct       vic411_options;
  extern vic411_param_set_struct    vic411_param_set;
  extern vic411_global_param_struct vic411_global_param;
  extern int                 vic411_NR, vic411_NF;

  int     i;
  int     j;
  int     band;
  int     day;
  int     hour;
  int     rec;
  int     step;
  int     idx;
  int    *tmaxhour;
  int    *tminhour;
  double  deltat;
  double  min_Tfactor;
  double  shortwave;
  double  svp_tair;
  double *hourlyrad;
  double *prec;
  double *tmax;
  double *tmin;
  double *tair;
  double *tskc;
  double *vp;
  double  min, max;
  double  rainonly;
  int     Ndays;
  int     stepspday;
  double  sum;
  double **forcing_data;
  int     type;
  double  air_temp;
  double  factor;

//<devel -- bug.patch>
  // step was uninitialized
  step = 0;
//</devel -- bug.patch>
  /* compute number of simulation days */
  Ndays = ( vic411_global_param.nrecs * vic411_global_param.dt) / 24;

  /* compute number of full model time steps per day */
  stepspday = 24/vic411_global_param.dt;
  
  if ( !vic411_param_set.TYPE[PREC].SUPPLIED
    && ( ( !vic411_param_set.TYPE[RAINF].SUPPLIED && ( !vic411_param_set.TYPE[LSRAINF].SUPPLIED || !vic411_param_set.TYPE[CRAINF].SUPPLIED ) )
      || ( ( !vic411_param_set.TYPE[SNOWF].SUPPLIED && ( !vic411_param_set.TYPE[LSSNOWF].SUPPLIED || !vic411_param_set.TYPE[CSNOWF].SUPPLIED ) ) ) ) )
    vic411_nrerror("Precipitation (PREC, or { {RAINF or {LSRAINF and CRAINF}} and {SNOWF or {LSSNOWF and CSNOWF}} }) must be given to the model, check input files\n");

  if ((!vic411_param_set.TYPE[TMAX].SUPPLIED || !vic411_param_set.TYPE[TMIN].SUPPLIED) && !vic411_param_set.TYPE[AIR_TEMP].SUPPLIED )
    vic411_nrerror("Daily maximum and minimum air temperature or sub-daily air temperature must be given to the model, check input files\n");
  
  if ( vic411_param_set.TYPE[AIR_TEMP].SUPPLIED && vic411_param_set.FORCE_DT[vic411_param_set.TYPE[AIR_TEMP].SUPPLIED-1] == 24 )
    vic411_nrerror("Model cannot use daily average temperature, must provide daily maximum and minimum or sub-daily temperatures.");
  
  if ( vic411_param_set.TYPE[SHORTWAVE].SUPPLIED
        && !(vic411_param_set.TYPE[VP].SUPPLIED || vic411_param_set.TYPE[QAIR].SUPPLIED || vic411_param_set.TYPE[REL_HUMID].SUPPLIED) )
    vic411_nrerror("Sub-daily shortwave and vapor pressure forcing data must be supplied together.");

  /* mtclim routine memory allocations */

  hourlyrad  = (double *) calloc(Ndays*24, sizeof(double));
  prec       = (double *) calloc(Ndays*24, sizeof(double));
  tair       = (double *) calloc(Ndays*24, sizeof(double));
  tmax       = (double *) calloc(Ndays, sizeof(double));
  tmaxhour   = (int *)    calloc(Ndays, sizeof(double));
  tmin       = (double *) calloc(Ndays, sizeof(double));
  tminhour   = (int *)    calloc(Ndays, sizeof(double));
  tskc       = (double *) calloc(Ndays*24, sizeof(double));
  vp         = (double *) calloc(Ndays*24, sizeof(double));
  
  if (hourlyrad == NULL || prec == NULL || tair == NULL || tmax == NULL ||
      tmaxhour == NULL ||  tmin == NULL || tminhour == NULL || tskc == NULL ||
      vp == NULL)
    vic411_nrerror("Memory allocation failure in vic411_initialize_atmos()");
  
  /*******************************
    read in meteorological data 
  *******************************/

  forcing_data = vic411_read_forcing_data(infile, vic411_global_param);
  
  fprintf(stderr,"\nRead meteorological forcing file\n");

  /*************************************************
    Pre-processing
  *************************************************/

  /*************************************************
    Convert units from ALMA to VIC standard, if necessary
  *************************************************/
  if (vic411_options.ALMA_INPUT) {
    for (type=0; type<N_FORCING_TYPES; type++) {
      if (vic411_param_set.TYPE[type].SUPPLIED) {
        /* Convert moisture flux rates to accumulated moisture flux per time step */
        if (   type == PREC
            || type == RAINF
            || type == CRAINF
            || type == LSRAINF
            || type == SNOWF
            || type == CSNOWF
            || type == LSSNOWF
           ) {
          for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
            forcing_data[type][idx] *= vic411_global_param.dt * 3600;
          }
        }
        /* Convert temperatures from K to C */
        else if (   type == AIR_TEMP
                 || type == TMIN
                 || type == TMAX
                ) {
          for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
            forcing_data[type][idx] -= KELVIN;
          }
        }
      }
    }
  }
  else {
    for (type=0; type<N_FORCING_TYPES; type++) {
      if (vic411_param_set.TYPE[type].SUPPLIED) {
        /* Convert pressures from kPa to Pa */
        if (   type == PRESSURE
            || type == VP
           ) {
          for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
            forcing_data[type][idx] *= kPa2Pa;
          }
        }
      }
    }
  }

  /*************************************************
    Precipitation
  *************************************************/

  /*************************************************
    If provided, translate rainfall and snowfall
    into total precipitation
    NOTE: this overwrites any PREC data that was supplied
  *************************************************/

  if(vic411_param_set.TYPE[RAINF].SUPPLIED && vic411_param_set.TYPE[SNOWF].SUPPLIED) {
    /* rainfall and snowfall supplied */
    if (forcing_data[PREC] == NULL) {
      forcing_data[PREC] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),sizeof(double));
    }
    for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
      forcing_data[PREC][idx] = forcing_data[RAINF][idx] + forcing_data[SNOWF][idx];
    }
    vic411_param_set.TYPE[PREC].SUPPLIED = vic411_param_set.TYPE[RAINF].SUPPLIED;
  }
  else if(vic411_param_set.TYPE[CRAINF].SUPPLIED && vic411_param_set.TYPE[LSRAINF].SUPPLIED
    && vic411_param_set.TYPE[CSNOWF].SUPPLIED && vic411_param_set.TYPE[LSSNOWF].SUPPLIED) {
    /* convective and large-scale rainfall and snowfall supplied */
    if (forcing_data[PREC] == NULL) {
      forcing_data[PREC] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),sizeof(double));
    }
    for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
      forcing_data[PREC][idx] = forcing_data[CRAINF][idx] + forcing_data[LSRAINF][idx]
                               + forcing_data[CSNOWF][idx] + forcing_data[LSSNOWF][idx];
    }
    vic411_param_set.TYPE[PREC].SUPPLIED = vic411_param_set.TYPE[LSRAINF].SUPPLIED;
  }

  /*************************************************
    Create sub-daily precipitation if not provided
  *************************************************/

  if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[PREC].SUPPLIED-1] == 24) {
    /* daily prec provided */
    rec = 0;
    for (day = 0; day < Ndays; day++) {
      for (i = 0; i < stepspday; i++) {
	sum = 0;
	for (j = 0; j < vic411_NF; j++) {
	  atmos[rec].prec[j] = forcing_data[PREC][day] 
	    / (float)(vic411_NF * stepspday);
	  sum += atmos[rec].prec[j];
	}
	if(vic411_NF>1) atmos[rec].prec[vic411_NR] = sum;
	if(vic411_global_param.dt == 24) atmos[rec].prec[vic411_NR] = forcing_data[PREC][day];
	rec++;
      }
    }
  }
  else {
    /* sub-daily prec provided */
    idx = 0;
    for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
      sum = 0;
      for(i = 0; i < vic411_NF; i++) {
	atmos[rec].prec[i] = forcing_data[PREC][idx];
	sum += atmos[rec].prec[i];
	idx++;
      }
      if(vic411_NF>1) atmos[rec].prec[vic411_NR] = sum;
    }
  }

  /*************************************************
    Air Temperature
  *************************************************/

  /************************************************
    Set maximum daily air temperature if provided 
  ************************************************/

  if(vic411_param_set.TYPE[TMAX].SUPPLIED) {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[TMAX].SUPPLIED-1] == 24) {
      /* daily tmax provided */
      for (day = 0; day < Ndays; day++) {
	tmax[day] = forcing_data[TMAX][day];
      }
    }
    else {
      /* sub-daily tmax provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	tmax[rec/stepspday] = forcing_data[TMAX][idx];
	for(i = 0; i < vic411_NF; i++) idx++;
      }
    }
  }

  /************************************************
    Set minimum daily air temperature if provided 
  ************************************************/

  if(vic411_param_set.TYPE[TMIN].SUPPLIED) {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[TMIN].SUPPLIED-1] == 24) {
      /* daily tmin provided */
      for (day = 0; day < Ndays; day++) {
	tmin[day] = forcing_data[TMIN][day];
      }
    }
    else {
      /* sub-daily tmin provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	tmin[rec/stepspday] = forcing_data[TMIN][idx];
	for(i = 0; i < vic411_NF; i++) idx++;
      }
    }
  }

  /*************************************************
    Store sub-daily air temperature if provided
  *************************************************/

  if(vic411_param_set.TYPE[AIR_TEMP].SUPPLIED) {
    /* forcing data defined as equal to or less than SNOW_STEP */
    idx = 0;
    for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
      sum = 0;
      for (i = 0; i < vic411_NF; i++, step++) {
	atmos[rec].air_temp[i] = forcing_data[AIR_TEMP][idx];
	sum += atmos[rec].air_temp[i];
	idx++;
      }
      if(vic411_NF > 1) atmos[rec].air_temp[vic411_NR] = sum / (float)vic411_NF;
    }
  }

  /******************************************************
    Determine Tmax and Tmin from sub-daily temperatures
  ******************************************************/

  if(!(vic411_param_set.TYPE[TMAX].SUPPLIED && vic411_param_set.TYPE[TMIN].SUPPLIED)) {
    rec = 0;
    while(rec < vic411_global_param.nrecs) {
      min = max = atmos[rec].air_temp[0];
      for (j = 0; j < stepspday; j++) {
	for (i = 0; i < vic411_NF; i++, step++) {
	  if ( atmos[rec].air_temp[i] > max ) max = atmos[rec].air_temp[i];
	  if ( atmos[rec].air_temp[i] < min ) min = atmos[rec].air_temp[i];
	}
	rec++;
      }
      tmax[(rec-1)/stepspday] = max;
      tmin[(rec-1)/stepspday] = min;
    }
  }


  /*************************************************
    Pressures
  *************************************************/

  /*************************************************
    If provided, translate specific humidity and atm. pressure
    into vapor pressure
    NOTE: this overwrites any VP data that was supplied
  *************************************************/

  if(vic411_param_set.TYPE[QAIR].SUPPLIED && vic411_param_set.TYPE[PRESSURE].SUPPLIED) {
    /* specific humidity and atm. pressure supplied */
    if (forcing_data[VP] == NULL) {
      forcing_data[VP] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),sizeof(double));
    }
    for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
      forcing_data[VP][idx] = forcing_data[QAIR][idx] * forcing_data[PRESSURE][idx] / EPS;
    }
    vic411_param_set.TYPE[VP].SUPPLIED = vic411_param_set.TYPE[QAIR].SUPPLIED;
  }

  /*************************************************
    If provided, translate relative humidity and atm. pressure
    into vapor pressure
    NOTE: this overwrites any VP data that was supplied
  *************************************************/

  if(vic411_param_set.TYPE[REL_HUMID].SUPPLIED && vic411_param_set.TYPE[PRESSURE].SUPPLIED) {
    /* relative humidity and atm. pressure supplied */
    if (forcing_data[VP] == NULL) {
      forcing_data[VP] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),sizeof(double));
    }
    for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
      forcing_data[VP][idx] = forcing_data[REL_HUMID][idx] * vic411_svp(forcing_data[AIR_TEMP][idx]) / 100.;
    }
    vic411_param_set.TYPE[VP].SUPPLIED = vic411_param_set.TYPE[REL_HUMID].SUPPLIED;
  }


  /*************************************************
    Shortwave, vp, air temp, pressure, and density
  *************************************************/

  /**************************************************
    use the mtclim code to get the hourly shortwave 
    and the daily dew point temperature 

    requires prec, tmax, and tmin
  **************************************************/
  for (i = 0; i < Ndays; i++)
    prec[i] = 0;
  for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
    prec[rec/stepspday] += atmos[rec].prec[vic411_NR];
  }
  vic411_mtclim42_wrapper(0, 0, (theta_l-theta_s)*24./360., elevation, annual_prec,
		   phi, &vic411_global_param, dmy, prec, tmax, tmin, tskc, vp,
		   hourlyrad);

  /***********************************************************
    reaggregate the hourly shortwave to the larger timesteps 
  ***********************************************************/

  if(!vic411_param_set.TYPE[SHORTWAVE].SUPPLIED) {
    for (rec = 0, hour = 0; rec < vic411_global_param.nrecs; rec++) {
      for (i = 0; i < vic411_NF; i++) {
	atmos[rec].shortwave[i] = 0;
	for (j = 0; j < vic411_options.SNOW_STEP; j++, hour++) {
	  atmos[rec].shortwave[i] += hourlyrad[hour];
	}
	atmos[rec].shortwave[i] /= vic411_options.SNOW_STEP;
      }
      if (vic411_NF > 1) {
	atmos[rec].shortwave[vic411_NR] = 0;
	for (i = 0; i < vic411_NF; i++) {
	  atmos[rec].shortwave[vic411_NR] += atmos[rec].shortwave[i];
	}
	atmos[rec].shortwave[vic411_NR] /= vic411_NF;
      }
    }
  }
  else {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[SHORTWAVE].SUPPLIED-1] == 24) {
      /* daily shortwave provided; to get sub-daily, we will take the
         mtclim estimates and scale them to match the supplied daily totals */
      for (rec = 0, hour = 0; rec < vic411_global_param.nrecs; rec++) {
        for (i = 0; i < vic411_NF; i++) {
	  atmos[rec].shortwave[i] = 0;
	  for (j = 0; j < vic411_options.SNOW_STEP; j++, hour++) {
	    atmos[rec].shortwave[i] += hourlyrad[hour];
	  }
	  atmos[rec].shortwave[i] /= vic411_options.SNOW_STEP;
        }
        if (vic411_NF > 1) {
	  atmos[rec].shortwave[vic411_NR] = 0;
	  for (i = 0; i < vic411_NF; i++) {
	    atmos[rec].shortwave[vic411_NR] += atmos[rec].shortwave[i];
	  }
	  atmos[rec].shortwave[vic411_NR] /= vic411_NF;
        }
      }
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	if (forcing_data[SHORTWAVE][day] > 0 && atmos[rec].shortwave[vic411_NR] > 0)
	  factor = forcing_data[SHORTWAVE][day]/atmos[rec].shortwave[vic411_NR];
	else
	  factor = 0;
	for (i = 0; i < stepspday; i++) {
	  for (j = 0; j < vic411_NF; j++)
	    atmos[rec].shortwave[j] *= factor;
	  if(vic411_NF > 1)
	    atmos[rec].shortwave[vic411_NR] *= factor;
	  rec++;
        }
      }
    }
    else {
      /* sub-daily shortwave provided, so it will be used instead
         of the mtclim estimates */
      idx = 0;
      for (rec = 0, hour = 0; rec < vic411_global_param.nrecs; rec++) {
        sum = 0;
        for (i = 0; i < vic411_NF; i++, step++) {
	  atmos[rec].shortwave[i] = ( forcing_data[SHORTWAVE][idx] < 0 ) ? 0 : forcing_data[SHORTWAVE][idx];
	  sum += atmos[rec].shortwave[i];
	  idx++;
        }
        if (vic411_NF > 1) atmos[rec].shortwave[vic411_NR] = sum / (float)vic411_NF;
      }
    }
  }

  /**************************************************************************
    Calculate the hours at which the minimum and maximum temperatures occur
  **************************************************************************/
  if(!vic411_param_set.TYPE[AIR_TEMP].SUPPLIED) {
    vic411_set_max_min_hour(hourlyrad, Ndays, tmaxhour, tminhour);

  /**********************************************************************
    Calculate the subdaily and daily temperature based on tmax and tmin 
  **********************************************************************/

    vic411_HourlyT(vic411_options.SNOW_STEP, Ndays, tmaxhour, tmax, tminhour, tmin, tair);
    for (rec = 0, step = 0; rec < vic411_global_param.nrecs; rec++) {
      for (i = 0; i < vic411_NF; i++, step++) {
	atmos[rec].air_temp[i] = tair[step];
      }
      if (vic411_NF > 1) {
	atmos[rec].air_temp[vic411_NR] = 0;
	for (i = 0; i < vic411_NF; i++) {
	  atmos[rec].air_temp[vic411_NR] += atmos[rec].air_temp[i];
	}
	atmos[rec].air_temp[vic411_NR] /= vic411_NF;
      }
    }
  }

  /**************************************************
    calculate the subdaily and daily vapor pressure 
    and vapor pressure deficit
  **************************************************/

  if(!vic411_param_set.TYPE[VP].SUPPLIED) {
    for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
      atmos[rec].vp[vic411_NR] = vp[rec/stepspday];
      atmos[rec].vpd[vic411_NR] = vic411_svp(atmos[rec].air_temp[vic411_NR]) - atmos[rec].vp[vic411_NR];

      if(atmos[rec].vpd[vic411_NR]<0) {
	atmos[rec].vpd[vic411_NR]=0;
	atmos[rec].vp[vic411_NR]=vic411_svp(atmos[rec].air_temp[vic411_NR]);
      }
      
      for (i = 0; i < vic411_NF; i++) {
	atmos[rec].vp[i]  = atmos[rec].vp[vic411_NR];
	atmos[rec].vpd[i]  = (vic411_svp(atmos[rec].air_temp[i]) - atmos[rec].vp[i]);

	if(atmos[rec].vpd[i]<0) {
	  atmos[rec].vpd[i]=0;
	  atmos[rec].vp[i]=vic411_svp(atmos[rec].air_temp[i]);
	}
      
      }
    }
  }
  else {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[VP].SUPPLIED-1] == 24) {
      /* daily vp provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < vic411_NF; j++) {
	    atmos[rec].vp[j] = forcing_data[VP][day];
	    atmos[rec].vpd[j] = (vic411_svp(atmos[rec].air_temp[j]) 
				 - atmos[rec].vp[j]);
	    sum += atmos[rec].vp[j];
	  }
	  if(vic411_NF > 1) {
	    atmos[rec].vp[vic411_NR] = sum / (float)vic411_NF;
	    atmos[rec].vpd[vic411_NR] = (vic411_svp(atmos[rec].air_temp[vic411_NR]) 
				  - atmos[rec].vp[vic411_NR]);
	  }
	  rec++;
	}
      }
    }
    else {
      /* sub-daily vp provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < vic411_NF; i++) {
	  atmos[rec].vp[i] = forcing_data[VP][idx];
	  atmos[rec].vpd[i] = (vic411_svp(atmos[rec].air_temp[i]) 
			       - atmos[rec].vp[i]);
	  sum += atmos[rec].vp[i];
	  idx++;
	}
	if(vic411_NF > 1) {
	  atmos[rec].vp[vic411_NR] = sum / (float)vic411_NF;
	  atmos[rec].vpd[vic411_NR] = (vic411_svp(atmos[rec].air_temp[vic411_NR]) 
				    - atmos[rec].vp[vic411_NR]);
	}
      }
    }
  }

  /*************************************************
    Longwave
  *************************************************/

  /****************************************************************************
    calculate the daily and sub-daily longwave.  There is a separate case for
    the full energy and the water balance modes.  For water balance mode we 
    need to calculate the net longwave for the daily timestep and the incoming
    longwave for the SNOW_STEPs, for the full energy balance mode we always
    want the incoming longwave. 
  ****************************************************************************/

  if ( !vic411_param_set.TYPE[LONGWAVE].SUPPLIED ) {
    /** Incoming longwave radiation not supplied **/
    for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
      sum = 0;
      for (i = 0; i < vic411_NF; i++) {
	vic411_calc_longwave(&(atmos[rec].longwave[i]), tskc[rec/stepspday],
		      atmos[rec].air_temp[i], atmos[rec].vp[i]);
        sum += atmos[rec].longwave[i];
      }
      if(vic411_NF>1) atmos[rec].longwave[vic411_NR] = sum / (float)vic411_NF;
    }
  }
  else if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[LONGWAVE].SUPPLIED-1] == 24) {
    /* daily incoming longwave radiation provided */
    rec = 0;
    for (day = 0; day < Ndays; day++) {
      for (i = 0; i < stepspday; i++) {
	sum = 0;
	for (j = 0; j < vic411_NF; j++) {
	  atmos[rec].longwave[j] = forcing_data[LONGWAVE][day];
	  sum += atmos[rec].longwave[j];
	}
	if(vic411_NF>1) atmos[rec].longwave[vic411_NR] = sum / (float)vic411_NF; 
	rec++;
      }
    }
  }
  else {
    /* sub-daily incoming longwave radiation provided */
    idx = 0;
    for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
      sum = 0;
      for(i = 0; i < vic411_NF; i++) {
	atmos[rec].longwave[i] = forcing_data[LONGWAVE][idx];
	sum += atmos[rec].longwave[i];
	idx++;
      }
      if(vic411_NF>1) atmos[rec].longwave[vic411_NR] = sum / (float)vic411_NF; 
    }
  }

  /*************************************************
    Wind Speed
  *************************************************/

  /*************************************************
    If provided, translate WIND_E and WIND_N into WIND
    NOTE: this overwrites any WIND data that was supplied
  *************************************************/

  if(vic411_param_set.TYPE[WIND_E].SUPPLIED && vic411_param_set.TYPE[WIND_N].SUPPLIED) {
    /* specific humidity and atm. pressure supplied */
    if (forcing_data[WIND] == NULL) {
      forcing_data[WIND] = (double *)calloc((vic411_global_param.nrecs * vic411_NF),sizeof(double));
    }
    for (idx=0; idx<(vic411_global_param.nrecs*vic411_NF); idx++) {
      forcing_data[WIND][idx] = sqrt( forcing_data[WIND_E][idx]*forcing_data[WIND_E][idx]
                                    + forcing_data[WIND_N][idx]*forcing_data[WIND_N][idx] );
    }
    vic411_param_set.TYPE[WIND].SUPPLIED = vic411_param_set.TYPE[WIND_E].SUPPLIED;
  }

  /********************
    set the windspeed 
  ********************/

  if (vic411_param_set.TYPE[WIND].SUPPLIED) {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[WIND].SUPPLIED-1] == 24) {
      /* daily wind provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < vic411_NF; j++) {
	    if(forcing_data[WIND][day] < vic411_options.MIN_WIND_SPEED)
	      atmos[rec].wind[j] = vic411_options.MIN_WIND_SPEED;
	    else 
	      atmos[rec].wind[j] = forcing_data[WIND][day];
	    sum += atmos[rec].wind[j];
	  }
	  if(vic411_NF>1) atmos[rec].wind[vic411_NR] = sum / (float)vic411_NF;
	  if(vic411_global_param.dt == 24) {
	    if(forcing_data[WIND][day] < vic411_options.MIN_WIND_SPEED)
	      atmos[rec].wind[j] = vic411_options.MIN_WIND_SPEED;
	    else 
	      atmos[rec].wind[vic411_NR] = forcing_data[WIND][day];
	  }
	  rec++;
	}
      }
    }
    else {
      /* sub-daily wind speed provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < vic411_NF; i++) {
	  if(forcing_data[WIND][idx] <vic411_options.MIN_WIND_SPEED)
	    atmos[rec].wind[i] = vic411_options.MIN_WIND_SPEED;
	  else
	    atmos[rec].wind[i] = forcing_data[WIND][idx];
	  sum += atmos[rec].wind[i];
	  idx++;
	}
	if(vic411_NF>1) atmos[rec].wind[vic411_NR] = sum / (float)vic411_NF;
      }
    }
  }
  else {
    /* no wind data provided, use default constant */
    for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
      for (i = 0; i < vic411_NF; i++) {
	atmos[rec].wind[i] = 1.5;
      }
      atmos[rec].wind[vic411_NR] = 1.5;	
    }
  }

  /*************************************************
    Store atmospheric density if provided (kg/m^3)
  *************************************************/

  if(vic411_param_set.TYPE[DENSITY].SUPPLIED) {
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[DENSITY].SUPPLIED-1] == 24) {
      /* daily density provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < vic411_NF; j++) {
	    atmos[rec].density[j] = forcing_data[DENSITY][day];
	    sum += atmos[rec].density[j];
	  }
	  if(vic411_NF>1) atmos[rec].density[vic411_NR] = sum / (float)vic411_NF;
	  rec++;
	}
      }
    }
    else {
      /* sub-daily density provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < vic411_NF; i++) {
	  atmos[rec].density[i] = forcing_data[DENSITY][idx];
	  sum += atmos[rec].density[i];
	  idx++;
	}
	if(vic411_NF>1) atmos[rec].density[vic411_NR] = sum / (float)vic411_NF;
      }
    }
  }

  /**************************************
    Estimate Atmospheric Pressure (Pa) 
  **************************************/

  if(!vic411_param_set.TYPE[PRESSURE].SUPPLIED) {
    if(!vic411_param_set.TYPE[DENSITY].SUPPLIED) {
      /* Estimate pressure */
      if (vic411_options.PLAPSE) {
        /* Assume average virtual temperature in air column
           between ground and sea level = KELVIN+atmos[rec].air_temp[vic411_NR] + 0.5*elevation*LAPSE_PM */
        for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
          atmos[rec].pressure[vic411_NR] = PS_PM*exp(-elevation*G/(Rd*(KELVIN+atmos[rec].air_temp[vic411_NR]+0.5*elevation*LAPSE_PM)));
          for (i = 0; i < vic411_NF; i++) {
            atmos[rec].pressure[i] = PS_PM*exp(-elevation*G/(Rd*(KELVIN+atmos[rec].air_temp[i]+0.5*elevation*LAPSE_PM)));
          }
        }
      }
      else {
        /* set pressure to constant value */
        for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
	  atmos[rec].pressure[vic411_NR] = 95500.;
	  for (i = 0; i < vic411_NF; i++) {
	    atmos[rec].pressure[i] = atmos[rec].pressure[vic411_NR];
	  }
        }
      }
    }
    else {
      /* use observed densities to estimate pressure */
      if (vic411_options.PLAPSE) {
        for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
          atmos[rec].pressure[vic411_NR] = (KELVIN+atmos[rec].air_temp[vic411_NR])*atmos[rec].density[vic411_NR]*Rd;
          for (i = 0; i < vic411_NF; i++) {
            atmos[rec].pressure[i] = (KELVIN+atmos[rec].air_temp[i])*atmos[rec].density[i]*Rd;
          }
        }
      }
      else {
        for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
	  atmos[rec].pressure[vic411_NR] = (275.0 + atmos[rec].air_temp[vic411_NR])
	    *atmos[rec].density[vic411_NR]/0.003486;
	  for (i = 0; i < vic411_NF; i++) {
	    atmos[rec].pressure[i] = (275.0 + atmos[rec].air_temp[i])
	      *atmos[rec].density[i]/0.003486;
	  }
        }
      }
    }
  }
  else {
    /* observed atmospheric pressure supplied */
    if(vic411_param_set.FORCE_DT[vic411_param_set.TYPE[PRESSURE].SUPPLIED-1] == 24) {
      /* daily pressure provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < vic411_NF; j++) {
	    atmos[rec].pressure[j] = forcing_data[PRESSURE][day];
	    sum += atmos[rec].pressure[j];
	  }
	  if(vic411_NF>1) atmos[rec].pressure[vic411_NR] = sum / (float)vic411_NF;
	  rec++;
	}
      }
    }
    else {
      /* sub-daily pressure provided */
      idx = 0;
      for(rec = 0; rec < vic411_global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < vic411_NF; i++) {
	  atmos[rec].pressure[i] = forcing_data[PRESSURE][idx];
	  sum += atmos[rec].pressure[i];
	  idx++;
	}
	if(vic411_NF>1) atmos[rec].pressure[vic411_NR] = sum / (float)vic411_NF;
      }
    }
  }

  /********************************************************
    Estimate Atmospheric Density if not provided (kg/m^3)
  ********************************************************/

  if(!vic411_param_set.TYPE[DENSITY].SUPPLIED) {
    /* use pressure to estimate density */
    if (vic411_options.PLAPSE) {
      for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
        atmos[rec].density[vic411_NR] = atmos[rec].pressure[vic411_NR]/(Rd*(KELVIN+atmos[rec].air_temp[vic411_NR]));
        for (i = 0; i < vic411_NF; i++) {
          atmos[rec].density[i] = atmos[rec].pressure[i]/(Rd*(KELVIN+atmos[rec].air_temp[i]));
        }
      }
    }
    else {
      for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
        atmos[rec].density[vic411_NR] = 0.003486*atmos[rec].pressure[vic411_NR]/
	  (275.0 + atmos[rec].air_temp[vic411_NR]);
        for (i = 0; i < vic411_NF; i++) {
	  atmos[rec].density[i] = 0.003486*atmos[rec].pressure[i]/
	    (275.0 + atmos[rec].air_temp[i]);
        }
      }
    }
  }

  /****************************************************
    Determine if Snow will Fall During Each Time Step
  ****************************************************/

#if !OUTPUT_FORCE
  min_Tfactor = Tfactor[0];
  for (band = 1; band < vic411_options.SNOW_BAND; band++) {
    if (Tfactor[band] < min_Tfactor)
      min_Tfactor = Tfactor[band];
  }
  for (rec = 0; rec < vic411_global_param.nrecs; rec++) {
    atmos[rec].snowflag[vic411_NR] = FALSE;
    for (i = 0; i < vic411_NF; i++) {
      if ((atmos[rec].air_temp[i] + min_Tfactor) < vic411_global_param.MAX_SNOW_TEMP
	  &&  atmos[rec].prec[i] > 0) {
	atmos[rec].snowflag[i] = TRUE;
	atmos[rec].snowflag[vic411_NR] = TRUE;
      }
      else
	atmos[rec].snowflag[i] = FALSE;
    }
  }
#endif // OUTPUT_FORCE
 
  // Free temporary parameters
  free(hourlyrad);
  free(prec);
  free(tair);
  free(tmax);
  free(tmaxhour);
  free(tmin);
  free(tminhour);
  free(tskc);
  free(vp);

  for(i=0;i<N_FORCING_TYPES;i++) 
    if (vic411_param_set.TYPE[i].SUPPLIED) 
      free((char *)forcing_data[i]);
  free((char *)forcing_data);

#if OUTPUT_FORCE_STATS
  calc_forcing_stats(vic411_global_param.nrecs, atmos);
#endif // OUTPUT_FORCE_STATS

#if !OUTPUT_FORCE

  // If COMPUTE_TREELINE is TRUE and the treeline computation hasn't
  // specifically been turned off for this cell (by supplying avgJulyAirTemp
  // and setting it to -999), calculate which snowbands are above the
  // treeline, based on average July air temperature.
  if (vic411_options.COMPUTE_TREELINE) {
    if ( !(vic411_options.JULY_TAVG_SUPPLIED && avgJulyAirTemp == -999) ) {
      vic411_compute_treeline( atmos, dmy, avgJulyAirTemp, Tfactor, AboveTreeLine );
    }
  }

#else

  // If OUTPUT_FORCE is set to TRUE in vic411_user_def.h then the full
  // forcing data array is dumped into a new set of files.
  write_forcing_file(atmos, vic411_global_param.nrecs, out_data_files, out_data);

#endif // OUTPUT_FORCE 

}
