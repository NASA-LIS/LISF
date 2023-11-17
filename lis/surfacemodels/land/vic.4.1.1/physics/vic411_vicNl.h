/* RCS Id String
 * $Id: vic411_vicNl.h,v 5.20.2.25 2009/10/08 21:30:59 vicadmin Exp $
 */
/************************************************************************
  Modifications:
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data, out_data_files, and save_data structures.	TJB
	      Removed the following functions:
		conv_force_vic2alma
		conv_results_vic2alma
	      Added the following new functions:
		vic411_create_output_list
		vic411_free_out_data_files
		vic411_init_output_list
		vic411_parse_output_info
		vic411_set_output_defaults
		vic411_set_output_var
		vic411_zero_output_list
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Jan-15 Added PRT_HEADER option.					TJB
  2007-Apr-03 Modified the data types of the following functions for
	      CONTINUE_ON_ERROR:					KAC/GTC
	      vic411_CalcAerodynamic
	      vic411_dist_prec
	      vic411_distribute_node_moisture_properties
	      vic411_full_energy
	      vic411_initialize_new_storm
	      vic411_redistribute_during_storm
	      vic411_runoff
	      vic411_snow_intercept
	      vic411_snow_melt
	      vic411_solve_T_profile
	      vic411_surface_fluxes
  2007-Apr-24 Added Ming Pan's new functions for IMPLICIT option.       JCA
              vic411_fda_heat_eqn
              vic411_newt_raph
              vic411_tridiag
              vic411_fdjac3
              vic411_solve_T_profile_implicit
  2007-Apr-21 Added functions:						TJB
	      vic411_free_dmy
	      vic411_free_out_data
	      vic411_free_veglib
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2007-Aug-22 Made vic411_calc_water_balance_error  type double.		JCA
  2007-Nov-06 Moved vic411_get_dist() from vic411_LAKE.h to this file.		JCA
  2008-Feb-17 Changed argument list for vic411_snow_density().			KMA via TJB
  2008-Apr-21 Added snow depth and albedo to vic411_snow_albedo() argument
	      list.							KAC via TJB
  2008-Oct-23 Modified vic411_put_data() to be type int, so that it can
	      return an error status.					TJB
  2009-Jan-16 Added avgJulyAirTemp to argument list of
	      vic411_compute_treeline().					TJB
  2009-Feb-09 Removed dz_node from several functions.			KAC via TJB
  2009-Mar-16 Added resid_moist to argument list of
	      vic411_estimate_layer_ice_content().				TJB
  2009-May-17 Added asat to argument list of vic411_surface_fluxes(),
	      vic411_full_energy(), and wetland_energy().			TJB
  2009-Jun-09 Modified argument lists of some functions that were
	      modified to use extension of vic411_veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-09 Added vic411_compute_pot_evap().					TJB
  2009-Jun-09 Removed unnecessary functions quick_penman() and
	      compute_penman_constants().				TJB
  2009-Jun-19 Added T vic411_flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jun-26 Simplified argument list of vic411_runoff() by passing all cell_data
	      variables via a single reference to the cell data structure.	TJB
  2009-Jul-07 Added soil_con.BandElev[] to vic411_read_snowband() arg list.	TJB
  2009-Jul-31 Removed unused layer_node_fract array from
	      vic411_estimate_layer_ice_content().				TJB 
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-28 Added vic411_collect_wb_terms() and vic411_collect_eb_terms(). Changed
	      argument list of vic411_read_snowband().				TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
************************************************************************/

#include <math.h>
#include <vic411_vicNl_def.h>
#include <vic411_LAKE.h>

/*** SubRoutine Prototypes ***/

double vic411_advected_sensible_heat(double, double, double, double, double);
void vic411_alloc_atmos(int, vic411_atmos_data_struct **);
double vic411_arno_evap(vic411_layer_data_struct *, vic411_layer_data_struct *, double, double, 
		 double, double, double, double, double, double, double, 
		 double, double, double, double, double, 
#if SPATIAL_FROST
		 double, double *);
#else
		 double);
#endif // SPATIAL_FROST
unsigned char vic411_average_moisture_for_storm(double *, double *, double, double);

int   vic411_CalcAerodynamic(char, double, double, double, double, double,
	  	       double *, double *, double *, double *, double *);
void   vic411_calc_cloud_cover_fraction(vic411_atmos_data_struct *, vic411_dmy_struct *, int,
				 int, int, double *);
void   vic411_calc_energy_balance_error(int, double, double, double, double, double);
#if OUTPUT_FORCE_STATS
void   calc_forcing_stats(int, vic411_atmos_data_struct *);
#endif // OUTPUT_FORCE_STATS
void   vic411_calc_longwave(double *, double, double, double);
void   vic411_calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
double vic411_calc_rainonly(double,double,double,double,double);
double vic411_calc_rc(double,double,float,double,double,double,double,char);
void   vic411_calc_root_fractions(vic411_veg_con_struct *, vic411_soil_con_struct *);
double vic411_calc_snow_coverage(int *, double, double, double, double, double, 
                          double, double, double *, double *, double *, 
                          double *, double *);
double calc_snow_ground_flux(int, int, int, int, double, double, double, 
			     double, double, double *, double *, double *, 
			     double *, vic411_energy_bal_struct *, 
			     vic411_snow_data_struct *, vic411_layer_data_struct *,
                             vic411_layer_data_struct *, vic411_soil_con_struct *, char *);
#if QUICK_FS
int    vic411_calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *, 
				double *, double *, double *,double *, 
				double *, double *, double *, 
				double *, double *, double *, double ***, int, int, int, int);
#else
int    vic411_calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *, 
				double *, double *, double *,double *, 
				double *, double *, double *, 
				double *, double *, double *, 
#if EXCESS_ICE
				double *, double *,
#endif // EXCESS_ICE
				int, int, int, int);
#endif // QUICK_FS
double vic411_CalcSnowPackEnergyBalance(double Tsurf, ...);
double vic411_CalcBlowingSnow(double, double, int, double, double, double, double, 
                       double, double, double, double, double, float, 
                       float, double, int, int, float, double, double, double *); 
double vic411_calc_atmos_energy_bal(double, double, double, double, double, double, 
                             double, double, double, double, double, double, 
                             double, double, double, double, 
                             double *, double *, double *, double *, 
                             double *, double *, double *, double *, char *, int *);
double vic411_calc_surf_energy_bal(double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double,
                            double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *,
                            float *, int, int,
                            int, int, int, int, int, int, int, int, int, int,
                            vic411_atmos_data_struct *, vic411_dmy_struct *,
                            vic411_energy_bal_struct *, vic411_layer_data_struct *,
                            vic411_layer_data_struct *,
                            vic411_snow_data_struct *, vic411_soil_con_struct *,
                            vic411_veg_var_struct *, vic411_veg_var_struct *, int);
double calc_trans(double, double);
double vic411_calc_veg_displacement(double);
double vic411_calc_veg_height(double);
double vic411_calc_veg_roughness(double);
double vic411_calc_water_balance_error(int, double, double, double);
double vic411_canopy_evap(vic411_layer_data_struct *, vic411_layer_data_struct *,
		   vic411_veg_var_struct *, vic411_veg_var_struct *, char, int, int, 
		   double, double *, double, double, double, double, 
		   double, double, double, double, double, double, 
		   double *, double *, double *, double *, 
#if SPATIAL_FROST
                   double *, float *);
#else
                   float *);
#endif
void   vic411_check_files(vic411_filep_struct *, vic411_filenames_struct *);
FILE  *vic411_check_state_file(char *, vic411_dmy_struct *, vic411_global_param_struct *, int, int, 
                        int *);
void   vic411_close_files(vic411_filep_struct *, vic411_out_data_file_struct *, vic411_filenames_struct *);
vic411_filenames_struct vic411_cmd_proc(int argc, char *argv[]);
void   vic411_collect_eb_terms(vic411_energy_bal_struct, vic411_snow_data_struct, vic411_cell_data_struct,
                        int *, int *, int *, int *, int *, double, double, double,
                        int, int, int, int, double *, double *,
#if SPATIAL_FROST
                        double *, double,
#endif
                        vic411_out_data_struct *);
void   vic411_collect_wb_terms(vic411_cell_data_struct, vic411_veg_var_struct, vic411_snow_data_struct, vic411_lake_var_struct,
                        double, double, double, double, int, int, int, double *, vic411_out_data_struct *);
void   vic411_compress_files(char string[]);
void   vic411_compute_dz(double *, double *, int, double);
void   vic411_correct_precip(double *, double, double, double, double);
void   vic411_compute_pot_evap(int, vic411_dmy_struct *, int, int, double, double , double, double, double, double **, double *);
void   vic411_compute_soil_layer_thermal_properties(vic411_layer_data_struct *, double *,
					     double *, double *, 
					     double *, 
#if SPATIAL_FROST
                                             double *,
#endif
					     int);
void   vic411_compute_treeline(vic411_atmos_data_struct *, vic411_dmy_struct *, double, double *, char *);
vic411_out_data_struct *vic411_create_output_list();

void   vic411_display_current_settings(int, vic411_filenames_struct *, vic411_global_param_struct *);
int    vic411_dist_prec(vic411_atmos_data_struct *,vic411_dist_prcp_struct *,vic411_soil_con_struct *,
		 vic411_veg_con_struct *, vic411_lake_con_struct *,
		 vic411_dmy_struct *,vic411_global_param_struct *,
		 vic411_filep_struct *, vic411_out_data_file_struct *,
		 vic411_out_data_struct *, vic411_save_data_struct *,
		 int, int, char, char, char *, int *);
#if QUICK_FS
int  vic411_distribute_node_moisture_properties(double *, double *, double *, double *,
					 double *, double *, double *, double ***, 
					 double *, double *, double *,
					 double *, double *, int, int, char);
#else
#if EXCESS_ICE
int  vic411_distribute_node_moisture_properties(double *, double *, double *, double *,
					 double *, double *, double *, double *, 
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, int, int, char);
#else
int  vic411_distribute_node_moisture_properties(double *, double *, double *, 
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, int, int, char);
#endif
#endif
void   distribute_soil_property(double *,double,double,
				double **l_param,
				int, int, double *, double *);

double vic411_error_calc_atmos_energy_bal(double Tcanopy, ...);
double vic411_error_calc_atmos_moist_bal(double , ...);
double vic411_error_calc_canopy_energy_bal(double Tsurf, ...);
double error_calc_snow_ground_flux(double Tsurf, ...);
double vic411_error_calc_surf_energy_bal(double Tsurf, ...);
double vic411_ErrorSnowPackEnergyBalance(double Tsurf, ...);
double vic411_error_print_atmos_energy_bal(double, va_list);
double vic411_error_print_atmos_moist_bal(double, va_list);
double vic411_error_print_canopy_energy_bal(double, va_list);
double error_print_snow_ground_flux(double, va_list);
double vic411_ErrorPrintSnowPackEnergyBalance(double, va_list);
double vic411_error_print_solve_T_profile(double, va_list);
double vic411_error_print_surf_energy_bal(double, va_list);
double vic411_error_solve_T_profile(double Tsurf, ...);
double estimate_dew_point(double, double, double, double, double);
#if QUICK_FS
int vic411_estimate_layer_ice_content(vic411_layer_data_struct *, double *, double *,
			       double *, double ***, double *,
			       double *, double ***, 
#if SPATIAL_FROST
			       double *, double,
#endif // SPATIAL_FROST
			       double *, double *, double *, double *, 
			       int, int, char);
#else
int vic411_estimate_layer_ice_content(vic411_layer_data_struct *, double *, double *,
			       double *, double *, double *, double *,
			       double *, double *, double *, 
#if SPATIAL_FROST
			       double *, double, 
#endif // SPATIAL_FROST
#if EXCESS_ICE
			       double *, double *,
#endif // EXCESS_ICE
			       double *, double *, double *, double *, 
			       int, int, char);
#endif
double vic411_estimate_T1(double, double, double, double, double, double, double, 
		   double, double, double, double);
double vic411_exp_interp(double,double,double,double,double);

double f(double, double, double, double, double, double, double, double,
         double, double, int, double *, double, double, double, double *,
         double *, double *, double *, double *, double *);
void   vic411_fda_heat_eqn(double *, double *, int, int, ...);
void   vic411_fdjac3(double *, double *, double *, double *, double *,
            void (*vecfunc)(double *, double *, int, int, ...), 
            int);
void   vic411_find_0_degree_fronts(vic411_energy_bal_struct *, double *, double *, int);
vic411_layer_data_struct vic411_find_average_layer(vic411_layer_data_struct *, vic411_layer_data_struct *,
				     double, double);
void   find_sublayer_temperatures(vic411_layer_data_struct *, double *, double *,
				  double *, double, double, int, int);
int    vic411_finish_frozen_soil_calcs(vic411_energy_bal_struct *, vic411_layer_data_struct *,
				vic411_layer_data_struct *, vic411_layer_data_struct *,
				vic411_soil_con_struct *, int, int, double, 
				double *, double *, double *, double *);
void   vic411_free_atmos(int nrecs, vic411_atmos_data_struct **atmos);
void   vic411_free_dist_prcp(vic411_dist_prcp_struct *, int);
void   vic411_free_dmy(vic411_dmy_struct **dmy);
void   vic411_free_vegcon(vic411_veg_con_struct **);
void   vic411_free_veglib(vic411_veg_lib_struct **);
void   vic411_free_out_data_files(vic411_out_data_file_struct **);
void   vic411_free_out_data(vic411_out_data_struct **);
int    vic411_full_energy(char, int, int, vic411_atmos_data_struct *, vic411_dist_prcp_struct *,
		   vic411_dmy_struct *, vic411_global_param_struct *, vic411_lake_con_struct *,
                   vic411_soil_con_struct *, vic411_veg_con_struct *);
double func_aero_resist(double,double,double,double,double);
double vic411_func_atmos_energy_bal(double, va_list);
double vic411_func_atmos_moist_bal(double, va_list);
double vic411_func_canopy_energy_bal(double, va_list);
double func_snow_ground_flux(double, va_list);
double vic411_func_surf_energy_bal(double, va_list);

double get_avg_temp(double, double, double *, double *, int);
double vic411_get_dist(double, double, double, double);
void   vic411_get_force_type(char *, int, int *);
vic411_global_param_struct vic411_get_global_param(vic411_filenames_struct *, FILE *);
void   vic411_get_next_time_step(int *, int *, int *, int *, int *, int);

double vic411_hermint(double, int, double *, double *, double *, double *, double *);
void   vic411_hermite(int, double *, double *, double *, double *, double *);
void   vic411_HourlyT(int, int, int *, double *, int *, double *, double *);

void   vic411_init_output_list(vic411_out_data_struct *, int, char *, int, float);
void   vic411_initialize_atmos(vic411_atmos_data_struct *, vic411_dmy_struct *, FILE **, double, 
			double, double, double, double, double, double, 
                        double, double *, 
#if OUTPUT_FORCE
			char *, vic411_out_data_file_struct *, vic411_out_data_struct *);
#else
			char *);
#endif
void   vic411_initialize_global();
int   vic411_initialize_model_state(vic411_dist_prcp_struct *, vic411_dmy_struct,
			      vic411_global_param_struct *, vic411_filep_struct, 
			      int, int, int, int, 
			      double, vic411_soil_con_struct *,
                              vic411_veg_con_struct *, vic411_lake_con_struct,
			      char **, int **, vic411_save_data_struct *);
int    vic411_initialize_new_storm(vic411_cell_data_struct ***, vic411_veg_var_struct ***,
			    int, int, int, double, double);
void   vic411_initialize_snow(vic411_snow_data_struct **, int, int);
void   vic411_initialize_soil(vic411_cell_data_struct **, vic411_soil_con_struct *, vic411_veg_con_struct *, int);
void   vic411_initialize_veg( vic411_veg_var_struct **, vic411_veg_con_struct *,
		       vic411_global_param_struct *, int);

void   vic411_latent_heat_from_snow(double, double, double, double, double, 
                             double, double, double *, double *, 
                             double *, double *, double *);
double vic411_linear_interp(double,double,double,double,double);

vic411_cell_data_struct **vic411_make_cell_data(int, int);
vic411_dist_prcp_struct vic411_make_dist_prcp(int);
vic411_dmy_struct *vic411_make_dmy(vic411_global_param_struct *);
vic411_energy_bal_struct **vic411_make_energy_bal(int);
void vic411_make_in_and_outfiles(vic411_filep_struct *, vic411_filenames_struct *, 
			  vic411_soil_con_struct *, vic411_out_data_file_struct *);
vic411_out_data_struct *make_out_data(int);
vic411_snow_data_struct **vic411_make_snow_data(int);
vic411_veg_var_struct **vic411_make_veg_var(int);
void   vic411_MassRelease(double *,double *,double *,double *);
#if EXCESS_ICE
double vic411_maximum_unfrozen_water(double, double, double, double, double, double);
#else
double vic411_maximum_unfrozen_water(double, double, double, double);
#endif
#if QUICK_FS
double maximum_unfrozen_water_quick(double, double, double **);
#endif
double vic411_modify_Ksat(double);
void vic411_mtclim42_wrapper(int, int, double, double, double, double, 
		      vic411_global_param_struct *, vic411_dmy_struct *, double *, 
		      double *, double *, double *, double *, double *);

double vic411_new_snow_density(double);
int    vic411_newt_raph(void (*vecfunc)(double *, double *, int, int, ...), 
               double *, int);
void   vic411_nrerror(char *);

void   open_debug();
FILE  *vic411_open_file(char string[], char type[]);
FILE  *vic411_open_state_file(vic411_global_param_struct *, vic411_filenames_struct, int, int);

void vic411_parse_output_info(vic411_filenames_struct *, FILE *, vic411_out_data_file_struct **, vic411_out_data_struct *);
double vic411_penman(double, double, double, double, double, double, double);
void   vic411_prepare_full_energy(int, int, int, vic411_dist_prcp_struct *, 
			   vic411_soil_con_struct *, double *, double *); 
double priestley(double, double);
int    vic411_put_data(vic411_dist_prcp_struct *, vic411_atmos_data_struct *,
		vic411_soil_con_struct *, vic411_veg_con_struct *,
                vic411_lake_con_struct *, vic411_out_data_file_struct *,
		vic411_out_data_struct *, vic411_save_data_struct *,
 	        vic411_dmy_struct *, int); 

double vic411_read_arcinfo_value(char *, double, double);
int    vic411_read_arcinfo_info(char *, double **, double **, int **);
void   vic411_read_atmos_data(FILE *, vic411_global_param_struct, int, int, double **);
double **vic411_read_forcing_data(FILE **, vic411_global_param_struct);
void   vic411_read_initial_model_state(FILE *, vic411_dist_prcp_struct *, 
				vic411_global_param_struct *, int, int, int, 
				vic411_soil_con_struct *, int, char *,
				int *, vic411_lake_con_struct);
void   vic411_read_snowband(FILE *, vic411_soil_con_struct *);
void   read_snowmodel(vic411_atmos_data_struct *, FILE *, int, int, int, int);
vic411_soil_con_struct vic411_read_soilparam(FILE *, int);
//vic411_soil_con_struct vic411_read_soilparam_arc(FILE *, char *, int *, int *, int);
vic411_veg_lib_struct *vic411_read_veglib(FILE *, int *);
vic411_veg_con_struct *vic411_read_vegparam(FILE *, int, int);
int    vic411_redistribute_during_storm(vic411_cell_data_struct ***, vic411_veg_var_struct ***,
				 int, int, int, double, double, double, 
				 double *);
void   redistribute_moisture(vic411_layer_data_struct *, double *, double *,
			     double *, double *, double *, int);
unsigned char vic411_redistribute_moisture_for_storm(double *, double *, double, 
					      double, double);
double vic411_root_brent(double, double, char *, double (*Function)(double, va_list), ...);
int    vic411_runoff(vic411_cell_data_struct *, vic411_cell_data_struct *,
              vic411_energy_bal_struct *, vic411_soil_con_struct *, double *,
#if EXCESS_ICE
	      int,
#endif
#if SPATIAL_FROST
              double *, 
#endif
              double, int, int, int, int, int);

void vic411_set_max_min_hour(double *, int, int *, int *);
void vic411_set_node_parameters(double *, double *, double *, double *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *,
#if QUICK_FS
			 double ***,
#endif
#if EXCESS_ICE
			 double *, double *, double *, double *,
#endif
			 int, int, char);
vic411_out_data_file_struct *vic411_set_output_defaults(vic411_out_data_struct *);
int vic411_set_output_var(vic411_out_data_file_struct *, int, int, vic411_out_data_struct *, char *, int, char *, int, float);
double vic411_snow_albedo(double, double, double, double, double, double, int, char);
double vic411_snow_density(vic411_snow_data_struct *, double, double, double, double, double);
int    vic411_snow_intercept(double, double, double, double, double, double, double,
                      double, double, double, double, double, double, double, 
                      double, double, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, char *, int *, double *, double *, double *, 
                      double *, double *, double *, float *, int, int, int, 
                      int, int, int, int, vic411_layer_data_struct *, 
                      vic411_layer_data_struct *, vic411_soil_con_struct *, 
                      vic411_veg_var_struct *, vic411_veg_var_struct *);
int    vic411_snow_melt(double, double, double, double, double *, double, double *, double, 
		 double, double, double, double, double, double, double, 
                 double, double, double, double, double, double, 
                 double *, double *, double *, double *, double *, double *, 
                 double *, double *, double *, double *, double *, double *, 
                 int, int, int, int, vic411_snow_data_struct *, vic411_soil_con_struct *);
double vic411_SnowPackEnergyBalance(double, va_list);
double vic411_soil_conductivity(double, double, double, double, double);
void   soil_thermal_calc(vic411_soil_con_struct *, vic411_layer_data_struct *,
			 vic411_energy_bal_struct, double *, double *, double *,
			 int, int);
double vic411_soil_thermal_eqn(double, va_list);
double vic411_solve_snow(char, double, double, double, double, double, double,
                  double, double, double, double, double, double, double,
                  double, double, double, double, double, double, double,
                  double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  float *,
                  int, int, int, int, int, int, int, int, int,
                  int, int, int, int, int *, vic411_energy_bal_struct *,
                  vic411_layer_data_struct *, vic411_layer_data_struct *,
                  vic411_snow_data_struct *, vic411_soil_con_struct *,
                  vic411_veg_var_struct *, vic411_veg_var_struct *);
double vic411_solve_atmos_energy_bal(double Tcanopy, ...);
double vic411_solve_atmos_moist_bal(double , ...);
double vic411_solve_canopy_energy_bal(double Tfoliage, ...);
double solve_snow_ground_flux(double Tsurf, ...);
double vic411_solve_surf_energy_bal(double Tsurf, ...);
#if QUICK_FS
int    vic411_solve_T_profile(double *, double *, char *, int *, double *, double *,double *, 
		       double *, double, double *, double *, double *,
		       double *, double *, double *, double *, double, double *, double ***,
		       int, int *, int, int, int, int);
#else
int    vic411_solve_T_profile(double *, double *, char *, int *, double *, double *,double *, 
		       double *, double, double *, double *, double *,
		       double *, double *, double *, double *, double, double *,
#if EXCESS_ICE
		       double *, double *,
#endif
		       int, int *, int, int, int, int);

#endif
int   vic411_solve_T_profile_implicit(double *, double *, double *, double *, double *,
			       double *, double, double *, double *, double *,
#if EXCESS_ICE
			       double *, double *,
#endif
			       double *, double *, double *, double *, double, int, int *,
			       int, int, int, int, 
			       double *, double *, double *, double *);
double vic411_StabilityCorrection(double, double, double, double, double, double);
void   store_moisture_for_debug(int,int,double *,vic411_cell_data_struct ***,
				vic411_veg_var_struct ***,vic411_snow_data_struct **,
				vic411_soil_con_struct *);
int    vic411_surface_fluxes(char, double, double, double, double, 
#if EXCESS_ICE
		      int, double *, double *,
#endif
		      double, double, double *, double *, double **,
                      double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *,
		      float *, int, int, int, int, int, 
                      int, int, int, int, vic411_atmos_data_struct *, vic411_dmy_struct *, 
                      vic411_energy_bal_struct *, vic411_global_param_struct *, 
                      vic411_cell_data_struct *, vic411_cell_data_struct *, 
                      vic411_snow_data_struct *, vic411_soil_con_struct *, 
                      vic411_veg_var_struct *, vic411_veg_var_struct *, float, float, float);
double vic411_svp(double);
double vic411_svp_slope(double);

void vic411_transpiration(vic411_layer_data_struct *, int, int, double, double, double, 
		   double, double, double, double, double, double, double, 
		   double *, double *, double *, double *, double *, double *,
#if SPATIAL_FROST
                   double *,
#endif
                   float *);
void tridag(double *,double *,double *,double *,double *,int);
void vic411_tridiag(double *, double *, double *, double *, unsigned);
int vic411_update_thermal_nodes(vic411_dist_prcp_struct *, 
			  int, int, int, vic411_soil_con_struct *, vic411_veg_con_struct *);
void vic411_usage(char *);

void   vic411_vicerror(char *);
double vic411_volumetric_heat_capacity(double,double,double);

void vic411_write_atmosdata(vic411_atmos_data_struct *, int);
void vic411_write_data(vic411_out_data_file_struct *, vic411_out_data_struct *, vic411_dmy_struct *, int);
void write_debug(vic411_atmos_data_struct *, vic411_soil_con_struct *, vic411_cell_data_struct *,
                 vic411_energy_bal_struct *, vic411_snow_data_struct *, vic411_veg_var_struct *,
                 vic411_dmy_struct *, vic411_global_param_struct *,
                 double, double, int, int, int, int, int, char);
void write_dist_prcp(vic411_dist_prcp_struct *);
#if OUTPUT_FORCE
void write_forcing_file(vic411_atmos_data_struct *, int, vic411_out_data_file_struct *, vic411_out_data_struct *);
#endif
void vic411_write_header(vic411_out_data_file_struct *, vic411_out_data_struct *, vic411_dmy_struct *, vic411_global_param_struct);
void vic411_write_layer(vic411_layer_data_struct *, int, int, 
#if SPATIAL_FROST
                 double *,
#endif
                 double *);
void vic411_write_model_state(vic411_dist_prcp_struct *, vic411_global_param_struct *, int, 
		       int, vic411_filep_struct *, vic411_soil_con_struct *, char *,
		       int *, vic411_lake_con_struct);
void vic411_write_snow_data(vic411_snow_data_struct, int, int);
void vic411_write_soilparam(vic411_soil_con_struct *);
void vic411_write_vegparam(vic411_veg_con_struct *);
void vic411_write_vegvar(vic411_veg_var_struct *, int);

void vic411_zero_output_list(vic411_out_data_struct *);
