#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_initialize_global.c,v 5.12.2.15 2009/09/20 02:32:07 vicadmin Exp $";

void vic411_initialize_global() {
/*********************************************************************
  vic411_initialize_global              Keith Cherkauer       March 1998

  This subroutine initalizes all global parameters before they are 
  called by the model.

  option_strudt:           structure containing all global model vic411_options
  vic411_options.FULL_ENERGY    = TRUE - compute full energy balance
  vic411_options.FROZEN_SOIL    = TRUE - compute frozen soils
  vic411_options.DIST_PRCP      = TRUE - use distributed precipitation
  vic411_options.INIT_SOIL      = TRUE - use file to initialize soil column
  vic411_options.CORRPREC       = TRUE - correct precipitation measurements to
                           account for gauge loss due to wind
  vic411_options.MOISTFRACT     = TRUE - output moisture as fractional moisture
                           content instead of total water in mm
  vic411_options.BINARY_OUTPUT  = TRUE - create binary output fles
  vic411_options.BLOWING        = TRUE - calculate sublimation from blowing snow
  vic411_options.Nlayer         = Number of soil layers to use in model (at least
                           3 must be used for the energy balance model
  vic411_options.GRID_DECIMAL   = Number of decimal places used in the gridded
                           input and output file names
  vic411_options.SNOW_BAND      = Number of elevation bands over which to solve the
                           enery balance snow model
 
  vic411_debug_struct:            Structure cantains all debugging flags
  debug.DEBUG            = TRUE - turn on all debugging
  debug.PRT_SOIL         = TRUE - print soil parameter debugging files
  debug.PRT_VEGE         = TRUE - print vegetation parameter debugging files
  debug.PRT_GLOBAL       = TRUE - print global parameter debugging files
  debug.PRT_ATMOS        = TRUE - print forcing data debugging files
  debug.PRT_SNOW         = TRUE - print snow pack debugging files
  debug.PRT_FLUX         = TRUE - print energy flux debugging files
  debug.PRT_VAR          = TRUE - print variable debugging files
  debug.PRT_TEMP         = TRUE - 
  debug.PRT_MOIST        = TRUE - 
  debug.PRT_KAPPA        = TRUE - 
  debug.PRT_BALANCE      = TRUE - 
  debug.PRT_GRID         = TRUE - 
  debug.debug_dir        = 

  vic411_param_set.ALBEDO       = FALSE;
  vic411_param_set.AIR_TEMP     = FALSE;
  vic411_param_set.CRAINF       = FALSE;
  vic411_param_set.CSNOWF       = FALSE;
  vic411_param_set.DENSITY      = FALSE;
  vic411_param_set.LONGWAVE     = FALSE;
  vic411_param_set.LSRAINF      = FALSE;
  vic411_param_set.LSSNOWF      = FALSE;
  vic411_param_set.PREC         = FALSE;
  vic411_param_set.PRESSURE     = FALSE;
  vic411_param_set.QAIR         = FALSE;
  vic411_param_set.RAINF        = FALSE;
  vic411_param_set.REL_HUMID    = FALSE;
  vic411_param_set.SHORTWAVE    = FALSE;
  vic411_param_set.SNOWF        = FALSE;
  vic411_param_set.TMAX         = FALSE;
  vic411_param_set.TMIN         = FALSE;
  vic411_param_set.TSKC         = FALSE;
  vic411_param_set.VP           = FALSE;
  vic411_param_set.WIND         = FALSE;
  vic411_param_set.WIND_E       = FALSE;
  vic411_param_set.WIND_N       = FALSE;

  Modifications:
  11-18-02 Added BLOWING_SNOW and PRT_LAKE to global file 
           initialization.                                        KAC
  04-22-03 Removed LAKS and LAKE_PROFILE from pre-processor
           statements that removed them is LAKE_MODEL was FALSE.
           This causes problems with various checks and saves no
           appreciable time of memory.                            KAC
  10-May-04 Initialize ARC_SOIL, COMPRESS, and ARNO_PARAMS to FALSE.
	    Also changed limit on loop over forcing types from
	    hard-coded 17 to variable N_FORCING_TYPES.				TJB
  2005-Mar-24 Modified to handle ALMA-specific global vic411_options.			TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW.			TJB
  2005-May-01 Added the ALMA vars CRainf, CSSnowf, LSRainf, and LSSnowf.	TJB
  2005-May-02 Added the ALMA vars Wind_E, Wind_N.				TJB
  2005-11-29  Added SAVE_STATE (option is now set in global param file) 	GCT
  2006-0913  Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option.		TJB/GCT
  2006-Sep-23 Implemented flexible output configuration; added new
              vic411_options.Noutfiles and organized the vic411_options according
              to function.							TJB
  2006-Dec-29 Added REL_HUMID to the list of supported met input variables.	TJB
  2007-Jan-02 Added ALMA_INPUT option; removed TAIR and PSURF from list
	      of supported met input variables.					TJB
  2007-Jan-15 Added PRT_HEADER option.						TJB
  2007-Apr-25 Added IMPLICIT option.						JCA
  2007-Apr-25 Added EXP_TRANS option.						JCA
  2007-Sep-14 Replaced initialization of vic411_param_set.TYPE[i].SUPPLIED with a loop,
	      and added initialization of vic411_param_set.TYPE[i].SIGNED and
	      vic411_param_set.TYPE[i].multiplier.                                     TJB
  2008-Apr-21 Added SNOW_ALBEDO option.						KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.					TJB
  2009-Jan-12 Added COMPUTE_TREELINE option.					TJB
  2009-Jan-16 Added vic411_options.AERO_RESIST_CANSNOW.				TJB
  2009-May-17 Added vic411_options.MIN_LIQ.						TJB
  2009-May-18 Added vic411_options.PLAPSE.						TJB
  2009-May-20 Set default of vic411_options.AERO_RESIST_CANSNOW to AR_406_FULL.	TJB
  2009-May-20 Added vic411_options.GRND_FLUX_TYPE.					TJB
  2009-Aug-25 Changed default of vic411_options.BINARY_STATE_FILE to FALSE.		TJB
  2009-Sep-19 Moved TFALLBACK to its own separate option.			TJB

*********************************************************************/

  extern vic411_option_struct vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif
  extern vic411_param_set_struct vic411_param_set;

  int i, j;

  /** Initialize model option flags **/

  // simulation modes
  vic411_options.AboveTreelineVeg      = -1;
  vic411_options.AERO_RESIST_CANSNOW   = AR_406_FULL;
  vic411_options.BLOWING               = FALSE;
  vic411_options.COMPUTE_TREELINE      = FALSE;
  vic411_options.CONTINUEONERROR       = TRUE;
  vic411_options.CORRPREC              = FALSE;
  vic411_options.DIST_PRCP             = FALSE;
  vic411_options.EQUAL_AREA            = FALSE;
  vic411_options.EXP_TRANS             = FALSE;
  vic411_options.FROZEN_SOIL           = FALSE;
  vic411_options.FULL_ENERGY           = FALSE;
  vic411_options.GRND_FLUX             = FALSE;
  vic411_options.GRND_FLUX_TYPE        = GF_FULL;
  vic411_options.IMPLICIT              = FALSE;
  vic411_options.LAKES                 = FALSE;
  vic411_options.LAKE_PROFILE          = FALSE;
  vic411_options.MIN_LIQ               = FALSE;
  vic411_options.MIN_WIND_SPEED        = 0.0;
  vic411_options.Nlayer                = 2;
  vic411_options.Nnode                 = 3;
  vic411_options.NOFLUX                = FALSE;
  vic411_options.PLAPSE                = TRUE;
  vic411_options.PREC_EXPT             = 0.6;
  vic411_options.QUICK_FLUX            = TRUE;
  vic411_options.QUICK_SOLVE           = FALSE;
  vic411_options.ROOT_ZONES            = MISSING;
  vic411_options.SNOW_ALBEDO           = USACE;
  vic411_options.SNOW_BAND             = 1;
  vic411_options.SNOW_DENSITY          = DENS_BRAS;
  vic411_options.SNOW_STEP             = 1;
  vic411_options.TFALLBACK             = TRUE;
  // input vic411_options
  vic411_options.ARC_SOIL              = FALSE;
  vic411_options.BASEFLOW              = ARNO;
  vic411_options.GLOBAL_LAI            = FALSE;
  vic411_options.GRID_DECIMAL          = 2;
  vic411_options.JULY_TAVG_SUPPLIED    = FALSE;
  // state vic411_options
  vic411_options.BINARY_STATE_FILE     = FALSE;
  vic411_options.INIT_STATE            = FALSE;
  vic411_options.SAVE_STATE            = FALSE;
  // output vic411_options
  vic411_options.ALMA_OUTPUT           = FALSE;
  vic411_options.BINARY_OUTPUT         = FALSE;
  vic411_options.COMPRESS              = FALSE;
  vic411_options.MOISTFRACT            = FALSE;
  vic411_options.Noutfiles             = 2;
  vic411_options.PRT_HEADER            = FALSE;
  vic411_options.PRT_SNOW_BAND         = FALSE;

#if LINK_DEBUG 

  /** Initialize debugging control flags **/

  debug.DEBUG       = FALSE;
  debug.PRT_SOIL    = FALSE;
  debug.PRT_VEGE    = FALSE;
  debug.PRT_GLOBAL  = FALSE;
  debug.PRT_ATMOS   = FALSE;
  debug.PRT_SNOW    = FALSE;
  debug.PRT_FLUX    = FALSE;
  debug.PRT_VAR     = FALSE;
  debug.PRT_TEMP    = FALSE;
  debug.PRT_MOIST   = FALSE;
  debug.PRT_LAKE    = FALSE;
  debug.PRT_KAPPA   = FALSE;
  debug.PRT_BALANCE = FALSE;
  debug.PRT_GRID    = FALSE;
  strcpy(debug.debug_dir,"./");

#endif // LINK_DEBUG

  /** Initialize forcing file input controls **/

  for(j=0;j<N_FORCING_TYPES;j++) {
    vic411_param_set.TYPE[j].SUPPLIED = FALSE;
    vic411_param_set.TYPE[j].SIGNED   = 1;
    vic411_param_set.TYPE[j].multiplier = 1;
  }
  for(i=0;i<2;i++) {
    vic411_param_set.FORCE_DT[i] = MISSING;
    vic411_param_set.N_TYPES[i] = MISSING;
    vic411_param_set.FORCE_FORMAT[i] = MISSING;
    for(j=0;j<N_FORCING_TYPES;j++) vic411_param_set.FORCE_INDEX[i][j] = MISSING;
  }

}
