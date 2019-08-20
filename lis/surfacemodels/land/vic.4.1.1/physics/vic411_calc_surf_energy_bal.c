#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "vic411_vicNl.h"

static char vcid[] = "$Id: vic411_calc_surf_energy_bal.c,v 5.9.2.13 2009/09/20 02:32:07 vicadmin Exp $";

double vic411_calc_surf_energy_bal(double             Le,
			    double             LongUnderIn,
			    double             NetLongSnow, // net LW at snow surface
			    double             NetShortGrnd, // net SW transmitted thru snow
			    double             NetShortSnow, // net SW at snow surface
			    double             OldTSurf,
			    double             ShortUnderIn,
			    double             SnowAlbedo,
			    double             SnowLatent,
			    double             SnowLatentSub,
			    double             SnowSensible,
			    double             Tair, // T of canopy or air
			    double             VPDcanopy,
			    double             VPcanopy,
			    double             advection,
			    double             coldcontent,
			    double             delta_coverage, // change in coverage fraction
			    double             dp,
			    double             ice0,
			    double             melt_energy,
			    double             moist,
			    double             mu,
			    double             snow_coverage,
			    double             snow_depth,
			    double             BareAlbedo,
			    double             surf_atten,
			    double             vapor_flux,
			    double            *aero_resist,
			    double            *aero_resist_used,
			    double            *displacement,
			    double            *melt,
			    double            *ppt,
			    double            *rainfall,
			    double            *ref_height,
			    double            *roughness,
			    double            *snowfall,
			    double            *wind,
			    float             *root,
			    int                INCLUDE_SNOW,
			    int                UnderStory,
			    int                Nnodes,
			    int                Nveg,
			    int                band,
			    int                dt,
			    int                hour,
			    int                iveg,
			    int                nlayer,
			    int                overstory,
			    int                rec,
			    int                veg_class,
			    vic411_atmos_data_struct *atmos,
			    vic411_dmy_struct        *dmy,
			    vic411_energy_bal_struct *energy,
			    vic411_layer_data_struct *layer_dry,
			    vic411_layer_data_struct *layer_wet,
			    vic411_snow_data_struct  *snow,
			    vic411_soil_con_struct   *soil_con,
			    vic411_veg_var_struct    *veg_var_dry,
			    vic411_veg_var_struct    *veg_var_wet,
			    int                nrecs)
/**************************************************************
  vic411_calc_surf_energy_bal.c  Greg O'Donnell and Keith Cherkauer  Sept 9 1997
  
  This function calculates the surface temperature, in the
  case of no snow cover.  Evaporation is computed using the
  previous ground heat flux, and then used to comput latent heat
  in the energy balance routine.  Surface temperature is found
  using the Root Brent method (Numerical Recipies).
  

  modifications:
    02-29-00  Included variables needed to compute energy flux
              through the snow pack.  The ground surface energy
              balance will now be a mixture of snow covered
	      and bare ground, controlled by the snow cover 
	      fraction set in vic411_solve_snow.c                 KAC
    6-8-2000  Modified to make use of spatially distributed 
              soil frost                                   KAC
    03-09-01  Added QUICK_SOLVE vic411_options for finite difference
              soil thermal solution.  By iterating on only a
              few choice nodes near the soil surface the 
              simulation time can be significantly reduced
              with minimal additional energy balance errors.  KAC
    11-18-02  Modified to include the effects of blowing snow
              on the surface energy balance.                 LCB
    07-30-03  Made sure that local NOFLUX variable is always set
              to the vic411_options vic411_flag value.                      KAC
    04-Jun-04 Placed "ERROR" at beginning of screen dump in
	      vic411_error_print_surf_energy_bal.				TJB
    16-Jul-04 Cast the last 6 variables in the parameter list
	      passed to vic411_root_brent, vic411_error_calc_surf_energy_bal
	      and vic411_solve_surf_energy_bal as double, since for
	      some reason letting them remain ints or floats
	      caused them to become garbage in the child
	      functions.						TJB
    16-Jul-04 Modified the cap on vapor_flux to re-scale
	      blowing_flux and surface_flux proportionally
	      so that vapor_flux still = their sum.			TJB
    05-Aug-04 Removed lag_one, sigma_slope, and fetch from
	      parameter list, since these were only used in
	      the call to vic411_root_brent/vic411_func_surf_energy_bal(),
	      which no longer needs them.				TJB
    24-Aug-04 Modified the re-scaling of surface_flux to reduce
	      round-off errors.						TJB
    21-Sep-04 Added ErrorString to store error messages from
	      vic411_root_brent.						TJB
    28-Sep-04 Added aero_resist_used to store the aerodynamic
	      resistance used in flux calculations.			TJB
  2007-Apr-06 Modified to handle grid cell errors by returning to the
	      vic411_main subroutine, rather than ending the simulation.	GCT/KAC
  2007-Apr-24 Changed the read-in order of iveg, and VEG in
	      vic411_error_print_surf_energy_bal to be consistent with the
	      call order, also added year, day, and hour to the
	      argument list.						JCA
  2007-Apr-24 Features included for IMPLICIT frozen soils option.	JCA 
	      including:
	        passing in nrecs
	        passing nrec, nrecs, and iveg to vic411_func_surf_energy_bal
	        passing bulk_density, soil_density, and quartz to
	          vic411_func_surf_energy_bal
  2007-Apr-24 Features included for EXP_TRANS option for frozen soils
	      algorithm.						JCA
  2007-Apr-24 Passing in Zsum_node rather than recalculating.		JCA
  2007-Aug-08 Features included for EXCESS_ICE option for frozen soils
	      algorithm.						JCA
	      including:
	        passing in entire soil_con structure.
  2007-Aug-31 Checked vic411_root_brent return value against -998 rather
	      than -9998.						JCA
  2007-Sep-01 Checked for return value of ERROR from
	      vic411_solve_surf_energy_bal.					JCA
  2007-Nov-09 Modified code to reset NOFLUX boundary to global option
	      value before computing final soil column temperatures.
	      Previously NOFLUX was set to FALSE for initial QUICK_SOLVE
	      estimates, but never reset to reflect actual bottom
	      boundary conditions for final pass through solver.	KAC
  2009-Feb-09 Removed dz_node from variables passed to
	      vic411_func_surf_energy_bal.					KAC via TJB
  2009-May-22 Added TFALLBACK value to vic411_options.CONTINUEONERROR.  This
	      allows simulation to continue when energy balance fails
	      to converge by using previous T value.			TJB
  2009-Jun-19 Added T vic411_flag to indicate whether TFALLBACK occurred.	TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
***************************************************************/
{
  extern vic411_veg_lib_struct *vic411_veg_lib;
  extern vic411_option_struct   vic411_options;

  int      FIRST_SOLN[2];
  int      NOFLUX;
  int      EXP_TRANS;
  int      VEG;
  int      i;
  int      nidx;
  int      tmpNnodes;

  double   Cs1;
  double   Cs2;
  double   D1;
  double   D2;
  double   LongBareIn;
  double   NetLongBare;
  double   NetShortBare;
  double   T1;
  double   T1_old;
  double   T2;
  double   Ts_old;
  double   Tsnow_surf;
  double   Tsurf; 
  char     Tsurf_fbflag; 
  int      Tsurf_fbcount; 
  double   albedo;
  double   atmos_density;
  double   atmos_pressure;
  double   bubble;
  double   delta_t;
  double   emissivity;
  double   error;
  double   expt;
  double   kappa1;
  double   kappa2;
  double   kappa_snow;
  double   max_moist;
  double   ra;
  double   refrozen_water;

  double   Wdew[2];
  double  *T_node;
  double   Tnew_node[MAX_NODES];
  char     Tnew_fbflag[MAX_NODES];
  int      Tnew_fbcount[MAX_NODES];
  double  *Zsum_node;
  double  *kappa_node;
  double  *Cs_node;
  double  *moist_node;
  double  *bubble_node;
  double  *expt_node;
  double  *max_moist_node;
  double  *ice_node;
  double  *alpha;
  double  *beta;
  double  *gamma;
  vic411_layer_data_struct layer[MAX_LAYERS];

  double   T_lower, T_upper;
  double   LongSnowIn;
  double   TmpNetLongSnow;
  double   TmpNetShortSnow;
  double   old_swq, old_depth;
  char ErrorString[MAXSTRING];

  /**************************************************
    Set All Variables For Use
  **************************************************/
  /* Initialize T_fbflag */
  Tsurf_fbflag = 0;
  Tsurf_fbcount = energy->Tsurf_fbcount;
  for (nidx=0; nidx<Nnodes; nidx++) {
    Tnew_fbflag[nidx] = 0;
    Tnew_fbcount[nidx] = 0;
  }

  if(iveg!=Nveg) {
    if(vic411_veg_lib[veg_class].LAI[dmy->month-1] > 0.0) VEG = TRUE;
    else VEG = FALSE;
  }
  else VEG = FALSE;

  T2                  = energy->T[Nnodes-1]; // soil column bottom temp
  Ts_old              = energy->T[0]; // previous surface temperature
  T1_old              = energy->T[1]; // previous first node temperature
  atmos_density       = atmos->density[hour]; // atmospheric density
  atmos_pressure      = atmos->pressure[hour]; // atmospheric pressure
  emissivity          = 1.; // longwave emissivity
  kappa1              = energy->kappa[0]; // top node conductivity
  kappa2              = energy->kappa[1]; // second node conductivity
  Cs1                 = energy->Cs[0]; // top node heat capacity
  Cs2                 = energy->Cs[1]; // second node heat capacity
  D1                  = soil_con->depth[0]; // top node thickness
  D2                  = soil_con->depth[0]; // second node thickness
  delta_t             = (double)dt * 3600.;
  max_moist           = soil_con->max_moist[0] / (soil_con->depth[0]*1000.);
  bubble              = soil_con->bubble[0];
  expt                = soil_con->expt[0];
  Tsnow_surf          = snow->surf_temp;
  Wdew[WET]           = veg_var_wet->Wdew;
  if(vic411_options.DIST_PRCP) Wdew[DRY] = veg_var_dry->Wdew;
  FIRST_SOLN[0] = TRUE;
  FIRST_SOLN[1] = TRUE;
  if ( snow->depth > 0. ) 
    kappa_snow = K_SNOW * (snow->density) * (snow->density) / snow_depth;
  else 
    kappa_snow = 0;

  /** compute incoming and net radiation **/
  NetShortBare  = ( ShortUnderIn * (1. - ( snow_coverage + delta_coverage ) ) 
		    * (1. - BareAlbedo) + ShortUnderIn * ( delta_coverage ) 
		    * ( 1. - SnowAlbedo ) );
  LongBareIn    = (1. - snow_coverage ) * LongUnderIn;
  if ( INCLUDE_SNOW || snow->swq == 0 ) { 
    TmpNetLongSnow  = NetLongSnow;
    TmpNetShortSnow = NetShortSnow;
    LongSnowIn      = snow_coverage * LongUnderIn;
  }
  else {
    TmpNetShortSnow = 0.; 
    TmpNetLongSnow  = 0.; 
    LongSnowIn      = 0.;
  }

  /*************************************************************
    Prepare soil node variables for finite difference solution
  *************************************************************/

  if(!vic411_options.QUICK_FLUX) {

    bubble_node    = soil_con->bubble_node; 
    expt_node      = soil_con->expt_node; 
    max_moist_node = soil_con->max_moist_node;  
    alpha          = soil_con->alpha; 
    beta           = soil_con->beta; 
    gamma          = soil_con->gamma; 
    moist_node     = energy->moist;
    kappa_node     = energy->kappa_node;
    Cs_node        = energy->Cs_node;
    T_node         = energy->T;
    Zsum_node      = soil_con->Zsum_node;   
    ice_node       = energy->ice;

  }
  else {

    bubble_node    = NULL; 
    expt_node      = NULL; 
    max_moist_node = NULL;  
    alpha          = NULL; 
    beta           = NULL; 
    gamma          = NULL; 
    moist_node     = NULL;
    kappa_node     = NULL;
    Cs_node        = NULL;
    T_node         = NULL;
    Zsum_node      = NULL;
    ice_node       = NULL;

  }

  /**************************************************
    Find Surface Temperature Using Root Brent Method
  **************************************************/
  if(vic411_options.FULL_ENERGY) {

    /** If snow included in solution, temperature cannot exceed 0C  **/
    if ( INCLUDE_SNOW ) {
      T_lower = energy->T[0]-SURF_DT;
      T_upper = 0.;
    }
    else {
      T_lower = 0.5*(energy->T[0]+Tair)-SURF_DT;
      T_upper = 0.5*(energy->T[0]+Tair)+SURF_DT;
    }

    if ( vic411_options.QUICK_SOLVE && !vic411_options.QUICK_FLUX ) {
      // Set iterative Nnodes using the depth of the thaw layer
      tmpNnodes = 0;
      for ( nidx = Nnodes-5; nidx >= 0; nidx-- ) 
	if ( T_node[nidx] >= 0 && T_node[nidx+1] < 0 ) tmpNnodes = nidx+1;
      if ( tmpNnodes == 0 ) { 
	if ( T_node[0] <= 0 && T_node[1] >= 0 ) tmpNnodes = Nnodes;
	else tmpNnodes = 3;
      }
      else tmpNnodes += 4;
      NOFLUX = FALSE;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
      EXP_TRANS = FALSE;  // Why would we do this???
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

    }
    else { 
      tmpNnodes = Nnodes;
      NOFLUX = vic411_options.NOFLUX;
      EXP_TRANS = vic411_options.EXP_TRANS;
    }

    Tsurf = vic411_root_brent(T_lower, T_upper, ErrorString, vic411_func_surf_energy_bal,
		       rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t, Cs1, Cs2, D1, D2,
		       T1_old, T2, Ts_old, bubble, dp, 
		       expt, ice0, kappa1, kappa2, 
		       max_moist, moist, root, 
		       UnderStory, overstory, NetShortBare, NetShortGrnd, 
		       TmpNetShortSnow, Tair, atmos_density, 
		       atmos_pressure,
		       emissivity, LongBareIn, LongSnowIn, mu, surf_atten, 
		       VPcanopy, VPDcanopy, 
		       Wdew, displacement, aero_resist, aero_resist_used, 
		       rainfall, ref_height, roughness, wind, Le, 
		       energy->advection, OldTSurf, snow->pack_temp, 
		       Tsnow_surf, kappa_snow, melt_energy, 
		       snow_coverage, 
		       snow->density, snow->swq, snow->surf_water,
		       &energy->deltaCC, &energy->refreeze_energy, 
		       &snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
		       tmpNnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
		       alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
		       ice_node, kappa_node, max_moist_node, moist_node, 
		       soil_con, layer_wet, layer_dry, veg_var_wet, veg_var_dry, 
		       INCLUDE_SNOW, NOFLUX, EXP_TRANS, snow->snow, 
		       FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
		       &energy->deltaH, &energy->fusion, &energy->grnd_flux, 
		       &energy->latent, &energy->latent_sub, 
		       &energy->sensible, &energy->snow_flux, &energy->error);
 
    if(Tsurf <= -998 ) {  
      if (vic411_options.TFALLBACK) {
        Tsurf = Ts_old;
        Tsurf_fbflag = 1;
        Tsurf_fbcount++;
      }
      else {
        fprintf(stderr, "SURF_DT = %.2f\n", SURF_DT);
        error = vic411_error_calc_surf_energy_bal(Tsurf, dmy->year, dmy->month, dmy->day, dmy->hour, VEG, iveg,
					   veg_class, delta_t, Cs1, Cs2, D1, D2, 
					   T1_old, T2, Ts_old, 
					   soil_con->b_infilt, bubble, dp, 
					   expt, ice0, kappa1, kappa2, 
					   soil_con->max_infil, max_moist, 
					   moist, soil_con->Wcr, soil_con->Wpwp, 
					   soil_con->depth, 
					   soil_con->resid_moist, root, 
					   UnderStory, overstory, NetShortBare, 
					   NetShortGrnd, TmpNetShortSnow, Tair, 
					   atmos_density, atmos_pressure, 
					   (double)soil_con->elevation, 
					   emissivity, LongBareIn, LongSnowIn, 
					   mu, surf_atten, VPcanopy, VPDcanopy, 
					   Wdew, displacement, aero_resist, aero_resist_used, 
					   rainfall, ref_height, roughness, 
					   wind, Le, energy->advection, 
					   OldTSurf, snow->pack_temp, 
					   Tsnow_surf, 
					   kappa_snow, melt_energy, 
					   snow_coverage, snow->density, 
					   snow->swq, snow->surf_water, 
					   &energy->deltaCC, 
					   &energy->refreeze_energy, 
					   &snow->vapor_flux, Nnodes, Cs_node, 
					   T_node, Tnew_node, alpha, beta, 
					   bubble_node, Zsum_node, expt_node, 
					   gamma, ice_node, kappa_node, 
					   max_moist_node, moist_node, 
#if SPATIAL_FROST
					   soil_con->frost_fract, 
#endif // SPATIAL_FROST
#if QUICK_FS
					   soil_con->ufwc_table_layer[0], 
					   soil_con->ufwc_table_node, 
#endif // QUICK_FS
					   layer_wet, layer_dry, veg_var_wet, 
					   veg_var_dry, 
					   INCLUDE_SNOW, soil_con->FS_ACTIVE, 
					   NOFLUX, EXP_TRANS,
					   snow->snow, FIRST_SOLN, &NetLongBare, 
					   &TmpNetLongSnow, &T1, &energy->deltaH, 
					   &energy->fusion, &energy->grnd_flux, 
					   &energy->latent, 
					   &energy->latent_sub, 
					   &energy->sensible, 
					   &energy->snow_flux, &energy->error, ErrorString);
        return ( ERROR );
      }
    }

    /**************************************************
      Recalculate Energy Balance Terms for Final Surface Temperature
    **************************************************/

    if ( Ts_old * Tsurf < 0 && vic411_options.QUICK_SOLVE ) {
      tmpNnodes = Nnodes;
      NOFLUX = vic411_options.NOFLUX;
      FIRST_SOLN[0] = TRUE;
      
      Tsurf = vic411_root_brent(T_lower, T_upper, ErrorString, vic411_func_surf_energy_bal,
			 rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t, Cs1, Cs2, D1, D2, 
			 T1_old, T2, Ts_old, bubble, dp, 
			 expt, ice0, kappa1, kappa2, 
			 max_moist, moist, root, 
			 UnderStory, overstory, NetShortBare, NetShortGrnd, 
			 TmpNetShortSnow, Tair, atmos_density, atmos_pressure, 
			 emissivity, LongBareIn, LongSnowIn, mu, surf_atten, 
			 VPcanopy, VPDcanopy, 
			 Wdew, displacement, aero_resist, aero_resist_used, 
			 rainfall, ref_height, roughness, wind, Le, 
			 energy->advection, OldTSurf, snow->pack_temp, 
			 Tsnow_surf, kappa_snow, melt_energy, 
			 snow_coverage, 
			 snow->density, snow->swq, snow->surf_water,
			 &energy->deltaCC, &energy->refreeze_energy, 
			 &snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
			 tmpNnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
			 alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
			 ice_node, kappa_node, max_moist_node, moist_node, 
			 soil_con, layer_wet, layer_dry, veg_var_wet, veg_var_dry, 
			 INCLUDE_SNOW, NOFLUX, EXP_TRANS, snow->snow, 
			 FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
			 &energy->deltaH, &energy->fusion, &energy->grnd_flux, 
			 &energy->latent, &energy->latent_sub, 
			 &energy->sensible, &energy->snow_flux, &energy->error);
      
      if(Tsurf <=  -998 ) {  
        if (vic411_options.TFALLBACK) {
          Tsurf = Ts_old;
          Tsurf_fbflag = 1;
          Tsurf_fbcount++;
        }
        else {
	  error = vic411_error_calc_surf_energy_bal(Tsurf, dmy->year, dmy->month, dmy->day, dmy->hour, VEG, iveg,
					     veg_class, delta_t, Cs1, Cs2, D1, 
					     D2, T1_old, T2, Ts_old, 
					     soil_con->b_infilt, bubble, dp, 
					     expt, ice0, kappa1, kappa2, 
					     soil_con->max_infil, max_moist, 
					     moist, soil_con->Wcr, 
					     soil_con->Wpwp, soil_con->depth, 
					     soil_con->resid_moist, root, 
					     UnderStory, overstory, 
					     NetShortBare, NetShortGrnd, 
					     TmpNetShortSnow, Tair, 
					     atmos_density, atmos_pressure, 
					     (double)soil_con->elevation, 
					     emissivity, LongBareIn, LongSnowIn, 
					     mu, surf_atten, VPcanopy, 
					     VPDcanopy, Wdew, displacement, 
					     aero_resist, aero_resist_used, rainfall, ref_height, 
					     roughness, wind, Le, 
					     energy->advection, 
					     OldTSurf, snow->pack_temp, 
					     Tsnow_surf, 
					     kappa_snow, melt_energy, 
					     snow_coverage, snow->density, 
					     snow->swq, snow->surf_water, 
					     &energy->deltaCC, 
					     &energy->refreeze_energy, 
					     &snow->vapor_flux, Nnodes, Cs_node, 
					     T_node, Tnew_node, alpha, beta, 
					     bubble_node, Zsum_node, expt_node, 
					     gamma, ice_node, kappa_node, 
					     max_moist_node, moist_node, 
#if SPATIAL_FROST
					     soil_con->frost_fract, 
#endif // SPTAIL_FROST
#if QUICK_FS
					     soil_con->ufwc_table_layer[0], 
					     soil_con->ufwc_table_node, 
#endif // QUICK_FS
					     layer_wet, layer_dry, veg_var_wet, 
					     veg_var_dry, INCLUDE_SNOW, 
					     soil_con->FS_ACTIVE, NOFLUX, EXP_TRANS,
					     snow->snow, FIRST_SOLN, &NetLongBare, 
					     &TmpNetLongSnow, &T1, 
					     &energy->deltaH, &energy->fusion, 
					     &energy->grnd_flux, 
					     &energy->latent, 
					     &energy->latent_sub, 
					     &energy->sensible, 
					     &energy->snow_flux, &energy->error, ErrorString);
          return ( ERROR );
        }
      }
    }
  }
  else {

    /** Frozen soil model run with no surface energy balance **/
    Tsurf  = Tair;
    NOFLUX = vic411_options.NOFLUX;
    EXP_TRANS = vic411_options.EXP_TRANS;

  }
  
  if ( vic411_options.QUICK_SOLVE && !vic411_options.QUICK_FLUX ) 
    // Reset model so that it solves thermal fluxes for full soil column
    FIRST_SOLN[0] = TRUE;
  
  error = vic411_solve_surf_energy_bal(Tsurf, 
				rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t, Cs1, Cs2, D1, D2, 
				T1_old, T2, Ts_old, bubble, dp, 
				expt, ice0, kappa1, kappa2, 
				max_moist, moist, root, 
				UnderStory, overstory, NetShortBare, NetShortGrnd, 
				TmpNetShortSnow, Tair, atmos_density, atmos_pressure,
				emissivity, LongBareIn, LongSnowIn, mu, surf_atten, 
				VPcanopy, VPDcanopy, 
				Wdew, displacement, aero_resist, aero_resist_used, 
				rainfall, ref_height, roughness, wind, Le, 
				energy->advection, OldTSurf, snow->pack_temp, 
				Tsnow_surf, kappa_snow, melt_energy, 
				snow_coverage, 
				snow->density, snow->swq, snow->surf_water,
				&energy->deltaCC, &energy->refreeze_energy, 
				&snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
				Nnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
				alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
				ice_node, kappa_node, max_moist_node, moist_node, 
				soil_con, layer_wet, layer_dry, veg_var_wet, veg_var_dry, 
				INCLUDE_SNOW, NOFLUX, EXP_TRANS, snow->snow, 
				FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
				&energy->deltaH, &energy->fusion, &energy->grnd_flux, 
				&energy->latent, &energy->latent_sub, 
				&energy->sensible, &energy->snow_flux, &energy->error);
  if(error == ERROR)
    return(ERROR);
  else
    energy->error = error;
  
  /***************************************************
    Recalculate Soil Moisture and Thermal Properties
  ***************************************************/
  if(vic411_options.GRND_FLUX) {
    if(vic411_options.QUICK_FLUX) {
      
      energy->T[0] = Tsurf;
      energy->T[1] = T1;
      
    }
    else {
      
      vic411_finish_frozen_soil_calcs(energy, layer_wet, layer_dry, layer, soil_con, 
			       Nnodes, iveg, mu, Tnew_node, kappa_node, 
			       Cs_node, moist_node);
      
    }
    
  }
  else {

    energy->T[0] = Tsurf;

  }

  /** Store precipitation that reaches the surface */
  if ( !snow->snow && !INCLUDE_SNOW ) {
    if ( iveg != Nveg ) {
      if ( vic411_veg_lib[veg_class].LAI[dmy->month-1] <= 0.0 ) { 
	veg_var_wet->throughfall = rainfall[WET];
	ppt[WET] = veg_var_wet->throughfall;
	if ( vic411_options.DIST_PRCP ) {
	  veg_var_dry->throughfall = rainfall[DRY];
	  ppt[DRY] = veg_var_dry->throughfall;
	}
      }
      else {
	ppt[WET] = veg_var_wet->throughfall;
	if(vic411_options.DIST_PRCP) ppt[DRY] = veg_var_dry->throughfall;
      }
    }
    else {
      ppt[WET] = rainfall[WET];
      if ( vic411_options.DIST_PRCP ) ppt[DRY] = rainfall[DRY];
    }
  }

  /****************************************
    Store understory energy balance terms 
  ****************************************/

// energy->sensible + energy->latent + energy->latent_sub + NetShortBare + NetLongBare + energy->grnd_flux + energy->deltaH + energy->fusion + energy->snow_flux

  energy->NetShortGrnd = NetShortGrnd;
  if ( INCLUDE_SNOW ) {
    energy->NetLongUnder  = NetLongBare + TmpNetLongSnow;
    energy->NetShortUnder = NetShortBare + TmpNetShortSnow + NetShortGrnd;
    energy->latent        = (energy->latent);
    energy->latent_sub    = (energy->latent_sub);
    energy->sensible      = (energy->sensible);
  }
  else {
    energy->NetLongUnder  = NetLongBare + NetLongSnow;
    energy->NetShortUnder = NetShortBare + NetShortSnow + NetShortGrnd;
/*     energy->latent        = (SnowLatent + (1. - snow_coverage)  */
/* 			     * energy->latent); */
/*     energy->latent_sub    = (SnowLatentSub  */
/* 			     + (1. - snow_coverage) * energy->latent_sub); */
/*     energy->sensible      = (SnowSensible  */
/* 			     + (1. - snow_coverage) * energy->sensible); */
    energy->latent        = (SnowLatent + energy->latent);
    energy->latent_sub    = (SnowLatentSub + energy->latent_sub);
    energy->sensible      = (SnowSensible + energy->sensible);
  }
  energy->LongUnderOut  = LongUnderIn - energy->NetLongUnder;
  energy->AlbedoUnder   = ((1. - ( snow_coverage + delta_coverage ) ) 
			   * BareAlbedo + ( snow_coverage + delta_coverage ) 
			   * SnowAlbedo );
  energy->melt_energy   = melt_energy;
  energy->Tsurf         = (snow->coverage * snow->surf_temp 
			   + (1. - snow->coverage) * Tsurf);

  /*********************************************************************
    adjust snow water balance for vapor mass flux if snowpack included 
  *********************************************************************/

//NEED TO ADJUST SNOW COVERAGE FRACTION - AND STORAGE

  if ( INCLUDE_SNOW ) {

    // don't allow vapor_flux to exceed snow pack
    if (-(snow->vapor_flux) > snow->swq) {
      // if vapor_flux exceeds snow pack, we not only need to
      // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
      snow->blowing_flux *= -( snow->swq / snow->vapor_flux );
//      snow->surface_flux *= -( snow->swq / snow->vapor_flux );
      snow->vapor_flux = -(snow->swq);
      snow->surface_flux = snow->vapor_flux - snow->blowing_flux;
    }

    /* adjust snowpack for vapor flux */
    old_swq           = snow->swq;
    snow->swq        += snow->vapor_flux;
    snow->surf_water += snow->vapor_flux;
    snow->surf_water  = ( snow->surf_water < 0 ) ? 0. : snow->surf_water;

    /* compute snowpack melt or refreeze */
    if (energy->refreeze_energy >= 0.0) {
      refrozen_water = energy->refreeze_energy / ( Lf * RHO_W ) * delta_t; 
      if ( refrozen_water > snow->surf_water) {
        refrozen_water = snow->surf_water;
        energy->refreeze_energy = refrozen_water * Lf * RHO_W / delta_t;
      } 
      snow->surf_water -= refrozen_water;
      assert(snow->surf_water >= 0.0);
      (*melt)           = 0.0;

    }
    else {
      
      /* Calculate snow melt */      
      (*melt) = fabs(energy->refreeze_energy) / (Lf * RHO_W) * delta_t;
      snow->swq -= *melt;
      if ( snow->swq < 0 ) { 
	(*melt) += snow->swq;
	snow->swq = 0;
      }
    }

    if ( snow->swq > 0 ) {

      // set snow energy terms
      snow->surf_temp   = ( Tsurf > 0 ) ? 0 : Tsurf;
      snow->coldcontent = CH_ICE * snow->surf_temp * snow->swq;

      // recompute snow depth
      old_depth   = snow->depth;
      snow->depth = 1000. * snow->swq / snow->density; 
      
      /** Check for Thin Snowpack which only Partially Covers Grid Cell
	  exists only if not snowing and snowpack has started to melt **/
#if SPATIAL_SNOW
      snow->coverage = vic411_calc_snow_coverage(&snow->store_snow, 
					  soil_con->depth_full_snow_cover, 
					  snow_coverage, snow->swq,
					  old_swq, snow->depth, old_depth, 
					  (*melt) - snow->vapor_flux, 
					  &snow->max_swq, snowfall, 
					  &snow->store_swq, 
					  &snow->swq_slope,
					  &snow->store_coverage);
      
#else
      if ( snow->swq > 0 ) snow->coverage = 1.;
      else snow->coverage = 0.;
#endif // SPATIAL_SNOW

      if ( snow->surf_temp > 0 ) 
	energy->snow_flux = ( energy->grnd_flux + energy->deltaH 
			      + energy->fusion );

    }
    else {
      /* snowpack melts completely */
      snow->density    = 0.;
      snow->depth      = 0.;
      snow->surf_water = 0;
      snow->pack_water = 0;
      snow->surf_temp  = 0;
      snow->pack_temp  = 0;
      snow->coverage   = 0;
#if SPATIAL_SNOW
      snow->store_swq = 0.;
#endif // SPATIAL_SNOW
    }
    snow->vapor_flux *= -1;
  }

  /** record T vic411_flag values **/
  energy->Tsurf_fbflag = Tsurf_fbflag;
  energy->Tsurf_fbcount += Tsurf_fbcount;
  for (nidx=0; nidx<Nnodes; nidx++) {
    energy->T_fbflag[nidx] = Tnew_fbflag[nidx];
    energy->T_fbcount[nidx] += Tnew_fbcount[nidx];
  }

  /** Return soil surface temperature **/
  return (Tsurf);
    
}

double vic411_solve_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap, Tsurf);
  error = vic411_func_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double vic411_error_calc_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap, Tsurf);
  error = vic411_error_print_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double vic411_error_print_surf_energy_bal(double Ts, va_list ap) {
/**********************************************************************
  Modifications:
  2009-Mar-03 Fixed format string for print statement, eliminates
	      compiler WARNING.						KAC via TJB
**********************************************************************/

  extern vic411_option_struct vic411_options;

  /* Define imported variables */

  /* general model terms */
  int year,month,day,hour;
  int iveg;
  int VEG;
  int veg_class;

  double delta_t;

  /* soil layer terms */
  double Cs1;
  double Cs2;
  double D1;
  double D2;
  double T1_old;
  double T2;
  double Ts_old;
  double b_infilt;
  double bubble;
  double dp;
  double expt;
  double ice0;
  double kappa1;
  double kappa2;
  double max_infil;
  double max_moist;
  double moist;

  double *Wcr;
  double *Wpwp;
  double *depth;
  double *resid_moist;

  float *root;

  /* meteorological forcing terms */
  int UnderStory;
  int overstory;

  double NetShortBare;  // net SW that reaches bare ground
  double NetShortGrnd;  // net SW that penetrates snowpack
  double NetShortSnow;  // net SW that reaches snow surface
  double Tair;
  double atmos_density;
  double atmos_pressure;
  double elevation;
  double emissivity;
  double LongBareIn; 
  double LongSnowIn; 
  double mu;
  double surf_atten;
  double vp;
  double vpd;

  double *Wdew;
  double *displacement;
  double *ra;
  double *ra_used;
  double *rainfall;
  double *ref_height;
  double *roughness;
  double *wind;
 
  /* latent heat terms */
  double  Le;

  /* snowpack terms */
  double Advection;
  double OldTSurf;
  double TPack;
  double Tsnow_surf;
  double kappa_snow;
  double melt_energy;
  double snow_coverage;
  double vic411_snow_density;
  double snow_swq;
  double snow_water;

  double *deltaCC;
  double *refreeze_energy;
  double *VaporMassFlux;

  /* soil node terms */
  int Nnodes;

  double *Cs_node;
  double *T_node;
  double *Tnew_node;
  double *alpha;
  double *beta;
  double *bubble_node;
  double *Zsum_node;
  double *expt_node;
  double *gamma;
  double *ice_node;
  double *kappa_node;
  double *max_moist_node;
  double *moist_node;

  /* spatial frost terms */
#if SPATIAL_FROST    
  double *frost_fract;
#endif

  /* quick solution frozen soils terms */
#if QUICK_FS
  double **ufwc_table_layer;
  double ***ufwc_table_node;
#endif

  /* model structures */
  vic411_layer_data_struct *layer_wet;
  vic411_layer_data_struct *layer_dry;
  vic411_veg_var_struct *veg_var_wet;
  vic411_veg_var_struct *veg_var_dry;

  /* control flags */
  int INCLUDE_SNOW;
  int FS_ACTIVE;
  int NOFLUX;
  int EXP_TRANS;

  int SNOWING;

  int *FIRST_SOLN;

  /* returned energy balance terms */
  double *NetLongBare; // net LW from snow-free ground
  double *NetLongSnow; // net longwave from snow surface - if INCLUDE_SNOW
  double *T1;
  double *deltaH;
  double *fusion;
  double *grnd_flux;
  double *latent_heat;
  double *latent_heat_sub;
  double *sensible_heat;
  double *snow_flux;
  double *store_error;

  char *ErrorString;

  /* Define internal routine variables */
  int                i;

  /***************************
    Read Variables from List
  ***************************/

  /* general model terms */
  year                   = (int) va_arg(ap, int);
  month                   = (int) va_arg(ap, int);
  day                   = (int) va_arg(ap, int);
  hour                   = (int) va_arg(ap, int);
  VEG                    = (int) va_arg(ap, int);
  iveg                     = (int) va_arg(ap, int);
  veg_class               = (int) va_arg(ap, int);

  delta_t                 = (double) va_arg(ap, double);

  /* soil layer terms */
  Cs1                     = (double) va_arg(ap, double);
  Cs2                     = (double) va_arg(ap, double);
  D1                      = (double) va_arg(ap, double);
  D2                      = (double) va_arg(ap, double);
  T1_old                  = (double) va_arg(ap, double);
  T2                      = (double) va_arg(ap, double);
  Ts_old                  = (double) va_arg(ap, double);
  b_infilt                = (double) va_arg(ap, double);
  bubble                  = (double) va_arg(ap, double);
  dp                      = (double) va_arg(ap, double);
  expt                    = (double) va_arg(ap, double);
  ice0                    = (double) va_arg(ap, double);
  kappa1                  = (double) va_arg(ap, double);
  kappa2                  = (double) va_arg(ap, double);
  max_infil               = (double) va_arg(ap, double);
  max_moist               = (double) va_arg(ap, double);
  moist                   = (double) va_arg(ap, double);

  Wcr                     = (double *) va_arg(ap, double *);
  Wpwp                    = (double *) va_arg(ap, double *);
  depth                   = (double *) va_arg(ap, double *);
  resid_moist             = (double *) va_arg(ap, double *);

  root                    = (float  *) va_arg(ap, float  *);

  /* meteorological forcing terms */
  UnderStory              = (int) va_arg(ap, int);
  overstory               = (int) va_arg(ap, int);

  NetShortBare            = (double) va_arg(ap, double);
  NetShortGrnd            = (double) va_arg(ap, double);
  NetShortSnow            = (double) va_arg(ap, double);
  Tair                    = (double) va_arg(ap, double);
  atmos_density           = (double) va_arg(ap, double);
  atmos_pressure          = (double) va_arg(ap, double);
  elevation               = (double) va_arg(ap, double);
  emissivity              = (double) va_arg(ap, double);
  LongBareIn              = (double) va_arg(ap, double);
  LongSnowIn              = (double) va_arg(ap, double);
  mu                      = (double) va_arg(ap, double);
  surf_atten              = (double) va_arg(ap, double);
  vp                      = (double) va_arg(ap, double);
  vpd                     = (double) va_arg(ap, double);

  Wdew                    = (double *) va_arg(ap, double *);
  displacement            = (double *) va_arg(ap, double *);
  ra                      = (double *) va_arg(ap, double *);
  ra_used                 = (double *) va_arg(ap, double *);
  rainfall                = (double *) va_arg(ap, double *);
  ref_height              = (double *) va_arg(ap, double *);
  roughness               = (double *) va_arg(ap, double *);
  wind                    = (double *) va_arg(ap, double *);

  /* latent heat terms */
  Le                      = (double) va_arg(ap, double);

  /* snowpack terms */
  Advection               = (double) va_arg(ap, double);
  OldTSurf                = (double) va_arg(ap, double);
  TPack                   = (double) va_arg(ap, double);
  Tsnow_surf              = (double) va_arg(ap, double);
  kappa_snow              = (double) va_arg(ap, double);
  melt_energy             = (double) va_arg(ap, double);
  snow_coverage           = (double) va_arg(ap, double);
  vic411_snow_density            = (double) va_arg(ap, double);
  snow_swq                = (double) va_arg(ap, double);
  snow_water              = (double) va_arg(ap, double);

  deltaCC                 = (double *) va_arg(ap, double *);
  refreeze_energy         = (double *) va_arg(ap, double *);
  VaporMassFlux           = (double *) va_arg(ap, double *);

  /* soil node terms */
  Nnodes                  = (int) va_arg(ap, int);

  Cs_node                 = (double *) va_arg(ap, double *);
  T_node                  = (double *) va_arg(ap, double *);
  Tnew_node               = (double *) va_arg(ap, double *);
  alpha                   = (double *) va_arg(ap, double *);
  beta                    = (double *) va_arg(ap, double *);
  bubble_node             = (double *) va_arg(ap, double *);
  Zsum_node               = (double *) va_arg(ap, double *);
  expt_node               = (double *) va_arg(ap, double *);
  gamma                   = (double *) va_arg(ap, double *);
  ice_node                = (double *) va_arg(ap, double *);
  kappa_node              = (double *) va_arg(ap, double *);
  max_moist_node          = (double *) va_arg(ap, double *);
  moist_node              = (double *) va_arg(ap, double *);

#if SPATIAL_FROST    
  frost_fract             = (double *) va_arg(ap, double *);
#endif

#if QUICK_FS
  ufwc_table_layer        = (double **) va_arg(ap, double **);
  ufwc_table_node         = (double ***) va_arg(ap, double ***);
#endif

  /* model structures */
  layer_wet               = (vic411_layer_data_struct *) va_arg(ap, vic411_layer_data_struct *);
  layer_dry               = (vic411_layer_data_struct *) va_arg(ap, vic411_layer_data_struct *);
  veg_var_wet             = (vic411_veg_var_struct *) va_arg(ap, vic411_veg_var_struct *);
  veg_var_dry             = (vic411_veg_var_struct *) va_arg(ap, vic411_veg_var_struct *);

  /* control flags */
  INCLUDE_SNOW            = (int) va_arg(ap, int);
  FS_ACTIVE               = (int) va_arg(ap, int);
  NOFLUX                  = (int) va_arg(ap, int);
  EXP_TRANS               = (int) va_arg(ap, int);
  SNOWING                 = (int) va_arg(ap, int);

  FIRST_SOLN              = (int *) va_arg(ap, int *);

  /* returned energy balance terms */
  NetLongBare             = (double *) va_arg(ap, double *);
  NetLongSnow             = (double *) va_arg(ap, double *);
  T1                      = (double *) va_arg(ap, double *);
  deltaH                  = (double *) va_arg(ap, double *);
  fusion                  = (double *) va_arg(ap, double *);
  grnd_flux               = (double *) va_arg(ap, double *);
  latent_heat             = (double *) va_arg(ap, double *);
  latent_heat_sub         = (double *) va_arg(ap, double *);
  sensible_heat           = (double *) va_arg(ap, double *);
  snow_flux               = (double *) va_arg(ap, double *);
  store_error             = (double *) va_arg(ap, double *);

  ErrorString             = (char *) va_arg(ap, char *);

  /***************
    Main Routine
  ***************/

  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: vic411_calc_surf_energy_bal failed to converge to a solution in vic411_root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  /* Print Variables */
  /* general model terms */
  fprintf(stderr, "iveg = %i\n", iveg);
  fprintf(stderr, "year = %i\n", year);
  fprintf(stderr, "month = %i\n", month);
  fprintf(stderr, "day = %i\n", day);
  fprintf(stderr, "hour = %i\n", hour);
  fprintf(stderr, "VEG = %i\n", VEG);
  fprintf(stderr, "veg_class = %i\n", veg_class);
  fprintf(stderr, "delta_t = %f\n",  delta_t);

  /* soil layer terms */
  fprintf(stderr, "Cs1 = %f\n",  Cs1);
  fprintf(stderr, "Cs2 = %f\n",  Cs2);
  fprintf(stderr, "D1 = %f\n",  D1);
  fprintf(stderr, "D2 = %f\n",  D2);
  fprintf(stderr, "T1_old = %f\n",  T1_old);
  fprintf(stderr, "T2 = %f\n",  T2);
  fprintf(stderr, "Ts_old = %f\n",  Ts_old);
  fprintf(stderr, "b_infilt = %f\n",  b_infilt);
  fprintf(stderr, "bubble = %f\n",  bubble);
  fprintf(stderr, "dp = %f\n",  dp);
  fprintf(stderr, "expt = %f\n",  expt);
  fprintf(stderr, "ice0 = %f\n",  ice0);
  fprintf(stderr, "kappa1 = %f\n",  kappa1);
  fprintf(stderr, "kappa2 = %f\n",  kappa2);
  fprintf(stderr, "max_infil = %f\n",  max_infil);
  fprintf(stderr, "max_moist = %f\n",  max_moist);
  fprintf(stderr, "moist = %f\n",  moist);

  fprintf(stderr, "*Wcr = %f\n",  *Wcr);
  fprintf(stderr, "*Wpwp = %f\n",  *Wpwp);
  fprintf(stderr, "*depth = %f\n",  *depth);
  fprintf(stderr, "*resid_moist = %f\n",  *resid_moist);

  fprintf(stderr, "*root = %f\n",  *root);

  /* meteorological forcing terms */
  fprintf(stderr, "UnderStory = %i\n", UnderStory);
  fprintf(stderr, "overstory = %i\n", overstory);

  fprintf(stderr, "NetShortBare = %f\n",  NetShortBare); 
  fprintf(stderr, "NetShortGrnd = %f\n",  NetShortGrnd); 
  fprintf(stderr, "NetShortSnow = %f\n",  NetShortSnow); 
  fprintf(stderr, "Tair = %f\n",  Tair);
  fprintf(stderr, "atmos_density = %f\n",  atmos_density);
  fprintf(stderr, "atmos_pressure = %f\n",  atmos_pressure);
  fprintf(stderr, "elevation = %f\n",  elevation);
  fprintf(stderr, "emissivity = %f\n",  emissivity);
  fprintf(stderr, "LongBareIn = %f\n",  LongBareIn); 
  fprintf(stderr, "LongSnowIn = %f\n",  LongSnowIn); 
  fprintf(stderr, "mu = %f\n",  mu);
  fprintf(stderr, "surf_atten = %f\n",  surf_atten);
  fprintf(stderr, "vp = %f\n",  vp);
  fprintf(stderr, "vpd = %f\n",  vpd);

  fprintf(stderr, "*Wdew = %f\n",  *Wdew);
  fprintf(stderr, "*displacement = %f\n",  *displacement);
  fprintf(stderr, "*ra = %f\n",  *ra);
  fprintf(stderr, "*ra_used = %f\n",  *ra_used);
  fprintf(stderr, "*rainfall = %f\n",  *rainfall);
  fprintf(stderr, "*ref_height = %f\n",  *ref_height);
  fprintf(stderr, "*roughness = %f\n",  *roughness);
  fprintf(stderr, "*wind = %f\n",  *wind);
 
  /* latent heat terms */
  fprintf(stderr, "Le = %f\n",   Le);

  /* snowpack terms */
  fprintf(stderr, "Advection = %f\n",  Advection);
  fprintf(stderr, "OldTSurf = %f\n",  OldTSurf);
  fprintf(stderr, "TPack = %f\n",  TPack);
  fprintf(stderr, "Tsnow_surf = %f\n",  Tsnow_surf);
  fprintf(stderr, "kappa_snow = %f\n",  kappa_snow);
  fprintf(stderr, "melt_energy = %f\n",  melt_energy);
  fprintf(stderr, "snow_coverage = %f\n",  snow_coverage);
  fprintf(stderr, "vic411_snow_density = %f\n",  vic411_snow_density);
  fprintf(stderr, "snow_swq = %f\n",  snow_swq);
  fprintf(stderr, "snow_water = %f\n",  snow_water);

  fprintf(stderr, "*deltaCC = %f\n",  *deltaCC);
  fprintf(stderr, "*refreeze_energy = %f\n",  *refreeze_energy);
  fprintf(stderr, "*VaporMassFlux = %f\n",  *VaporMassFlux);

  /* soil node terms */
  fprintf(stderr, "Nnodes = %i\n", Nnodes);

  /* spatial frost terms */
#if SPATIAL_FROST    
  fprintf(stderr, "*frost_fract = %f\n",  *frost_fract);
#endif

  /* control flags */
  fprintf(stderr, "INCLUDE_SNOW = %i\n", INCLUDE_SNOW);
  fprintf(stderr, "FS_ACTIVE = %i\n", FS_ACTIVE);
  fprintf(stderr, "NOFLUX = %i\n", NOFLUX);
  fprintf(stderr, "EXP_TRANS = %i\n", EXP_TRANS);
  fprintf(stderr, "SNOWING = %i\n", SNOWING);

  fprintf(stderr, "*FIRST_SOLN = %i\n", *FIRST_SOLN);

  /* returned energy balance terms */
  fprintf(stderr, "*NetLongBare = %f\n",  *NetLongBare); 
  fprintf(stderr, "*NetLongSnow = %f\n",  *NetLongSnow); 
  fprintf(stderr, "*T1 = %f\n",  *T1);
  fprintf(stderr, "*deltaH = %f\n",  *deltaH);
  fprintf(stderr, "*fusion = %f\n",  *fusion);
  fprintf(stderr, "*grnd_flux = %f\n",  *grnd_flux);
  fprintf(stderr, "*latent_heat = %f\n",  *latent_heat);
  fprintf(stderr, "*latent_heat_sub = %f\n",  *latent_heat_sub);
  fprintf(stderr, "*sensible_heat = %f\n",  *sensible_heat);
  fprintf(stderr, "*snow_flux = %f\n",  *snow_flux);
  fprintf(stderr, "*store_error = %f\n",  *store_error);

#if SPATIAL_FROST
  vic411_write_layer(layer_wet, iveg, vic411_options.Nlayer, frost_fract, depth);
#else
  vic411_write_layer(layer_wet, iveg, vic411_options.Nlayer, depth);
#endif
  if(vic411_options.DIST_PRCP) 
#if SPATIAL_FROST
    vic411_write_layer(layer_dry, iveg, vic411_options.Nlayer, frost_fract, depth);
#else
    vic411_write_layer(layer_dry, iveg, vic411_options.Nlayer, depth);
#endif
  vic411_write_vegvar(&(veg_var_wet[0]),iveg);
  if(vic411_options.DIST_PRCP) 
    vic411_write_vegvar(&(veg_var_dry[0]),iveg);

  if(!vic411_options.QUICK_FLUX) {
    fprintf(stderr,"Node\tT\tTnew\tZsum\tkappa\tCs\tmoist\tbubble\texpt\tmax_moist\tice\n");
    for(i=0;i<Nnodes;i++) 
      fprintf(stderr,"%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	      i, T_node[i], Tnew_node[i], Zsum_node[i], kappa_node[i], 
	      Cs_node[i], moist_node[i], bubble_node[i], expt_node[i], 
	      max_moist_node[i], ice_node[i]);
  }
  
  fprintf(stderr,"**********\n**********\nFinished writing vic411_calc_surf_energy_bal variables.\nTry increasing SURF_DT to get model to complete cell.\nThen check output for instabilities.\n**********\n**********\n");

  return (ERROR);
    
}

