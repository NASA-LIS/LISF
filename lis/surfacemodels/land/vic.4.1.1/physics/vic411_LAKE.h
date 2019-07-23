/******************************************************************************
// $Id: vic411_LAKE.h,v 5.9.2.12 2009/10/08 02:03:06 vicadmin Exp $
  Modifications:
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Aug-16 Made return value of initialize_prcp an int.		JCA
  2007-Aug-21 Added features for EXCESS_ICE option.			JCA
  2007-Oct-24 Changed vic411_get_sarea, vic411_get_volume, and vic411_get_depth to return exit
	      status so that errors can be trapped and communicated up the
	      chain of function calls.					KAC via TJB
  2007-Oct-24 Changed vic411_lakeice() to return exit status.			KAC via TJB
  2007-Nov-06 New lake physics parameters.  Modified argument lists for
	      various functions.  Moved vic411_get_dist() to vic411_vicNl.h.		LCB via TJB
  2008-Apr-21 Added argument to vic411_alblake.				LCB via TJB
  2008-Sep-10 Updated values of CONDS and lamwlw to match Laura
	      Bowling's lake work.					LCB via TJB
  2009-Jul-31 Removed lakemain() and wetland_energy(); vic411_initialize_lake
	      no longer takes a snow structure as input.		TJB
  2009-Sep-28 Removed initialize_prcp and update_prcp.  Modified
	      argument list of vic411_initialize_lake.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-05 Added functions for updating/rescaling lake and wetland
	      fluxes and storages when lake area changes.		TJB
******************************************************************************/

//#ifndef LAKE_SET

#define LAKE_SET
#define TMELT 0.0
#define EMICE 0.97      /* Ice emissivity */
#define EMH2O .98
#define RHOSNOW   250.  /* densities of water and snow */
#define RHOICE   917.   /* ice density*/
#define rhosurf 1.275   /* surface air density */
#define MAX_SURFACE_LAKE   .6  /* max. surface layer thickness for E-B (m) */
#define BETA 0.001       /* Curve shape parameter for lake profile. */
#define FRACMIN  0.10   /* min ice thickness in meters */
#define FRACLIM   0.02  /* lower limit on fractional ice cover */
#define CPW_ICE   4200. /* specific heat of ice */
#define DM 1.38889E-07  /* molecular diffusivity of water */
#define SNOWCRIT   0.05  /* for albedo, in m */
//#define G 9.80616
#define ZWATER 0.0045    // 0.004 - original value
#define ZSNOW 0.005
#define CONDI 2.3        /* thermal conductivity of ice */
#define CONDS 0.7       /* thermal conductivity of snow */ 

// attenuation of short and longwave radiation through ice (1/m)
#define lamisw 1.5 // 1.5 in Patterson & Hamblin
#define lamilw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through snow (1/m)
#define lamssw 6.0 // 6.0 in Patterson & Hamblin
#define lamslw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through water (1/m)
#define lamwsw .3  // San Fran Bay data: 0.31 - 29.9 1/m (visible)
#define lamwlw 1.4 // Hostetler and Bartlein assume 0.85 1/m (total)
#define  a1 0.7        /* Percent of radiation in visible band. */
#define  a2 0.3        /* Percent of radiation in infrared band. */
#define QWTAU 86400./2.   /* D. Pollard sub-ice time constant. */
#define RADIUS 6371.228 /* Earth radius in km. */

//#endif // LAKE_SET

/*** Subroutine prototypes ***/

double adjflux(double, double, double ,double, double, double, double,
	       double, double, double, double *, double *);
void vic411_advect_soil_veg_storage(double, double, double, double *, vic411_soil_con_struct *, vic411_veg_con_struct *, vic411_cell_data_struct *, vic411_veg_var_struct *);
void vic411_advect_snow_storage(double, double, double, vic411_snow_data_struct *);
void vic411_alblake(double, double, double *, double *, float *, float *, double, double, 
	     int, int *, double, double, char *, int);
void vic411_alloc_atmos(int, vic411_atmos_data_struct **);
double vic411_calc_density(double);
double vic411_CalcIcePackEnergyBalance(double Tsurf, ...);
void vic411_colavg (double *, double *, double *, float, double *, int, double, double);
float dragcoeff(float, double, double);
void vic411_eddy (int, double, double * , double *, double *, double, int, double, double);
void vic411_energycalc(double *, double *, int, double, double,double *, double *, double *);
double vic411_ErrorIcePackEnergyBalance(double Tsurf, ...);
double vic411_ErrorPrintIcePackEnergyBalance(double, va_list);
int vic411_get_depth(vic411_lake_con_struct, double, double *);
int vic411_get_depth_from_sarea(vic411_lake_con_struct, double, double *);
int vic411_get_sarea(vic411_lake_con_struct, double, double *);
int vic411_get_volume(vic411_lake_con_struct, double, double *);
void vic411_iceform (double *,double *,double ,double,double *,int, int, double, double, double *, double *, double *, double *, double *, double);
void vic411_icerad(double,double ,double,double *, double *,double *);
int vic411_ice_depth(vic411_lake_con_struct, double, double, double *);
int vic411_ice_melt(double, double, double *, double, vic411_snow_data_struct *, vic411_lake_var_struct *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double);
double vic411_IceEnergyBalance(double, va_list);
int vic411_initialize_lake(vic411_lake_var_struct *, vic411_lake_con_struct, vic411_soil_con_struct *, double);
int vic411_lakeice(double *, double, double, double, double, int, 
	    double, double, double *, double, double, int, vic411_dmy_struct, double *, double *, double, double);
void vic411_latsens(double,double, double, double, double, double, double, double,
	     double *, double *, double);
float vic411_lkdrag(float, double, double, double, double);
vic411_lake_con_struct vic411_read_lakeparam(FILE *, vic411_soil_con_struct, vic411_veg_con_struct *);
void vic411_rescale_lake_fluxes(double, double, vic411_lake_var_struct *);
void vic411_rescale_soil_veg_fluxes(double, double, vic411_cell_data_struct *, vic411_veg_var_struct *);
void vic411_rescale_snow_energy_fluxes(double, double, vic411_snow_data_struct *, vic411_energy_bal_struct *);
void vic411_rhoinit(double *, double);
int vic411_solve_lake(double, double, double, double, double, double, double, double, 
		double, double, vic411_lake_var_struct *, vic411_lake_con_struct, 
		vic411_soil_con_struct, int, int, double, vic411_dmy_struct, double);
double vic411_specheat (double);
void vic411_temp_area(double, double, double, double *, double *, double *, double *, int, double *, int, double, double, double*, double *, double *);
void vic411_tracer_mixer(double *, int *, int, double*, int, double, double, double *);
void vic411_tridia(int, double *, double *, double *, double *, double *);
int vic411_water_balance (vic411_lake_var_struct *, vic411_lake_con_struct, int, vic411_dist_prcp_struct *, int, int, int, double, vic411_soil_con_struct, vic411_veg_con_struct,
#if EXCESS_ICE
		    int, double,
#endif		    
		    double, double);
int  vic411_water_energy_balance(int, double*, double*, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double*, double *, double *, double *, double, double *, double *, double *, double *, double *, double);
int vic411_water_under_ice(int, double,  double, double *, double *, double, int, double, double, double, double *, double *, double *, double *, int, double, double, double, double *);
