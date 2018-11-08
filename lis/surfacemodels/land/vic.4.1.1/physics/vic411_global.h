/**********************************************************************
                        Global Variables

  NOTE: This file exists because global variables that are shared among
        files via the "extern" statement must be initially declared
        (without the word "extern") ONLY once.  Currently, vic411_vicNl_def.h
        is included (via vic411_vicNl.h) in every .c file, meaning that any
        declarations in vic411_vicNl_def.h end up happening multiple times
        (once per .c file).  Thus, these "extern" variables cannot be
        declared in vic411_vicNl_def.h.  This is not a problem for #define
        statements and typedef statements, which is what vic411_vicNl_def.h
        is primarily composed of.

  $Id: vic411_global.h,v 4.9.2.8 2009/10/20 18:42:07 vicadmin Exp $

  29-Oct-03 Added vic411_version string and removed unused vic411_options from
	    vic411_optstring.							TJB
  2009-Jun-09 Added definitions of reference landcover types, used
	      mainly for pot_evap computations but also defines the
	      characteristics of bare soil.				TJB
**********************************************************************/
char *vic411_version = "4.1.1";

char *vic411_optstring = "g:vo";

#if QUICK_FS
double   temps[] = { -1.e-5, -0.075, -0.20, -0.50, -1.00, -2.50, -5, -10 };
#endif

int vic411_flag;

vic411_global_param_struct vic411_global_param;
vic411_veg_lib_struct *vic411_veg_lib;
vic411_option_struct vic411_options;
#if LINK_DEBUG
vic411_debug_struct debug;
#endif
vic411_Error_struct vic411_Error;
vic411_param_set_struct vic411_param_set;

  /**************************************************************************
    Define some reference landcover types that always exist regardless
    of the contents of the library (mainly for potential evap calculations):
    Non-natural:
      satsoil = saturated bare soil
      h2osurf = open water surface (deep enough to have albedo of 0.08)
      short   = short reference crop (grass)
      tall    = tall reference crop (alfalfa)
    Natural:
      natveg  = current vegetation
      vegnocr = current vegetation with canopy resistance set to 0
    NOTE: these are external variables, declared in vic411_vicNl_def.h.
    NOTE2: bare soil roughness and displacement will be overwritten by the
           values found in the soil parameter file; bare soil wind_h will
	   be overwritten by the value specified in the global param file.
  **************************************************************************/

  /* One element for each non-natural PET type */
  char   vic411_ref_veg_over[]        = { 0, 0, 0, 0 };
  double vic411_ref_veg_rarc[]        = { 0.0, 0.0, 25, 25 };
  double vic411_ref_veg_rmin[]        = { 0.0, 0.0, 100, 100 };
  double vic411_ref_veg_lai[]         = { 1.0, 1.0, 2.88, 4.45 };
  double vic411_ref_veg_albedo[]      = { BARE_SOIL_ALBEDO, H2O_SURF_ALBEDO, 0.23, 0.23 };
  double vic411_ref_veg_rough[]       = { 0.001, 0.001, 0.0148, 0.0615 };
  double vic411_ref_veg_displ[]       = { 0.0054, 0.0054, 0.08, 0.3333 };
  double vic411_ref_veg_wind_h[]      = { 10.0, 10.0, 10.0, 10.0 };
  double vic411_ref_veg_RGL[]         = { 0.0, 0.0, 100, 100 };
  double vic411_ref_veg_rad_atten[]   = { 0.0, 0.0, 0.0, 0.0 };
  double vic411_ref_veg_wind_atten[]  = { 0.0, 0.0, 0.0, 0.0 };
  double vic411_ref_veg_trunk_ratio[] = { 0.0, 0.0, 0.0, 0.0 };
  /* One element for each PET type (non-natural or natural) */
  char vic411_ref_veg_ref_crop[] = { FALSE, FALSE, TRUE, TRUE, FALSE, FALSE };

