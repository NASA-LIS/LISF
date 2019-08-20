//<devel -- debugging.patch>
#if 0
//</devel -- debugging.patch>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vic411_vicNl.h>
#include <vic411_global.h>

static char vcid[] = "$Id: vicNl.c,v 5.14.2.13 2009/09/28 21:46:58 vicadmin Exp $";

/** Main Program **/

int vic411_main(int argc, char *argv[])
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
	  model.                                                  KAC
  02-27-01 added controls for lake model                          KAC
  11-18-02 Updated storage of lake water for water balance 
           calculations.                                          LCB
  03-12-03 Modifed to add AboveTreeLine to vic411_soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  04-10-03 Modified to initialize storm parameters using the state
           file.                                                KAC
  04-10-03 Modified to start the model by skipping records until the
           state file date is found.  This replaces the previous method
           of modifying the global file start date, which can change 
           the interpolation of atmospheric forcing data.        KAC
  04-15-03 Modified to store wet and dry fractions when intializing 
           water balance storage.  This accounts for changes in model
           state initialization, which now stores wet and dry fractions
           rather than just averagedvalues.                      KAC
  29-Oct-03 Modified the vic411_version display banner to print the vic411_version
	    string defined in vic411_global.h.					TJB
  01-Nov-04 Updated arglist for vic411_make_dist_prcp(), as part of fix for
	    QUICK_FLUX state file compatibility.			TJB
  02-Nov-04 Updated arglist for vic411_read_lakeparam(), as part of fix for
	    lake fraction readjustment.					TJB
  2005-Apr-13 OUTPUT_FORCE option now calls vic411_close_files().		TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data, out_data_files, and save_data structures.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into vic411_filep_struct;
	      This included merging builtnames into filenames.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-07 Changed statefile to init_state in call to
	      vic411_check_state_file().					TJB
  2007-Jan-15 Added PRT_HEADER option; added call to
	      vic411_write_header().						TJB
  2007-Apr-04 Added option to continue run after a cell fails. 		GCT/KAC
  2007-Apr-21 Added calls to vic411_free_dmy(), vic411_free_out_data_files(),
	      vic411_free_out_data(), and vic411_free_veglib().  Added closing of
	      all parameter files.					TJB
  2007-Aug-21 Return ErrorFlag from vic411_initialize_model_state.		JCA
  2007-Sep-14 Excluded calls to vic411_free_veglib() and closing of parameter
	      files other than the soil param file for the case
	      when OUTPUT_FORCE=TRUE.					TJB
  2007-Nov-06 Moved computation of cell_area from vic411_read_lakeparam() to
	      vic411_read_soilparam() and vic411_read_soilparam_arc().		TJB
  2008-May-05 Added prcp fraction (mu) to initial water storage
	      computation.  This solves water balance errors for the
	      case where DIST_PRCP is TRUE.				TJB
  2009-Jan-16 Added soil_con.avgJulyAirTemp to argument list of
	      vic411_initialize_atmos().					TJB
  2009-Jun-09 Modified to use extension of vic411_veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-07 Added soil_con.BandElev[] to vic411_read_snowband() arg list.	TJB
  2009-Jul-31 Replaced references to N+1st veg tile with references
	      to index of lake/wetland tile.				TJB
  2009-Sep-28 Replaced initial water/energy storage computations and
	      calls to vic411_calc_water_balance_error/vic411_calc_energy_balance_error
	      with an initial call to vic411_put_data.  Modified the call to
	      vic411_read_snowband().						TJB
**********************************************************************/
{

  extern vic411_veg_lib_struct *vic411_veg_lib;
  extern vic411_option_struct vic411_options;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif // LINK_DEBUG
  extern vic411_Error_struct vic411_Error;
  extern vic411_global_param_struct vic411_global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                    *init_STILL_STORM;
  char                     ErrStr[MAXSTRING];
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                     *init_DRY_TIME;
  int                      RUN_MODEL;
  int                      Ncells;
  int                      cell_cnt;
  int                      startrec;
  int                      ErrorFlag;
  float                    mu;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  double                   Clake;
  vic411_dmy_struct              *dmy;
  vic411_atmos_data_struct       *atmos;
  vic411_veg_con_struct          *veg_con;
  vic411_soil_con_struct          soil_con;
  vic411_dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */
  vic411_filenames_struct         filenames;
  vic411_filep_struct             filep;
  vic411_lake_con_struct          lake_con;
  vic411_out_data_file_struct     *out_data_files;
  vic411_out_data_struct          *out_data;
  vic411_save_data_struct         save_data;
  
  /** Read Model Options **/
  vic411_initialize_global();
  filenames = vic411_cmd_proc(argc, argv);

#if VERBOSE
  vic411_display_current_settings(DISP_VERSION,(vic411_filenames_struct*)NULL,(vic411_global_param_struct*)NULL);
#endif

  /** Read Global Control File **/
  filep.globalparam = vic411_open_file(filenames.global,"r");
  vic411_global_param = vic411_get_global_param(&filenames, filep.globalparam);

  /** Set up output data structures **/
  out_data = vic411_create_output_list();
  out_data_files = vic411_set_output_defaults(out_data);
  fclose(filep.globalparam);
  filep.globalparam = vic411_open_file(filenames.global,"r");
  vic411_parse_output_info(&filenames, filep.globalparam, &out_data_files, out_data);

  /** Check and Open Files **/
  vic411_check_files(&filep, &filenames);

#if !OUTPUT_FORCE

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

  /** Read Vegetation Library File **/
  vic411_veg_lib = vic411_read_veglib(filep.veglib,&Nveg_type);

#endif // !OUTPUT_FORCE

  /** Initialize Parameters **/
  if(vic411_options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /** Make Date Data Structure **/
  dmy      = vic411_make_dmy(&vic411_global_param);

  /** allocate memory for the vic411_atmos_data_struct **/
  vic411_alloc_atmos(vic411_global_param.nrecs, &atmos);

  /** Initial state **/
  startrec = 0;
#if !OUTPUT_FORCE
  if ( vic411_options.INIT_STATE ) 
    filep.init_state = vic411_check_state_file(filenames.init_state, dmy, 
					 &vic411_global_param, vic411_options.Nlayer, 
					 vic411_options.Nnode, &startrec);

  /** open state file if model state is to be saved **/
  if ( vic411_options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 )
    filep.statefile = vic411_open_state_file(&vic411_global_param, filenames, vic411_options.Nlayer,
                                         vic411_options.Nnode);
  else filep.statefile = NULL;

#endif // !OUTPUT_FORCE

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  cell_cnt=0;
  while(!MODEL_DONE) {
    if(!vic411_options.ARC_SOIL) {
      if((fscanf(filep.soilparam, "%d", &vic411_flag))!=EOF) {
	if(vic411_flag) RUN_MODEL=TRUE;
	else     RUN_MODEL=FALSE;
      }
      else {
	MODEL_DONE = TRUE;
	RUN_MODEL = FALSE;
      }
      if(!MODEL_DONE) soil_con = vic411_read_soilparam(filep.soilparam, RUN_MODEL);

    }
    else {
      soil_con = vic411_read_soilparam_arc(filep.soilparam, 
				    filenames.soil_dir, &Ncells, 
				    &RUN_MODEL, cell_cnt);
      cell_cnt++;
      if(cell_cnt==Ncells) MODEL_DONE = TRUE;
    }
    if(RUN_MODEL) {
#if LINK_DEBUG
      if(debug.PRT_SOIL) vic411_write_soilparam(&soil_con); 
#endif

#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(vic411_options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
	}
	for(i=0;i<MAX_NODES;i++) {
	  soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
	}
      }
#endif /* QUICK_FS */

      NEWCELL=TRUE;
      cellnum++;

#if !OUTPUT_FORCE

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = vic411_read_vegparam(filep.vegparam, soil_con.gridcel,
                              Nveg_type);
      vic411_calc_root_fractions(veg_con, &soil_con);
#if LINK_DEBUG
      if(debug.PRT_VEGE) vic411_write_vegparam(veg_con); 
#endif /* LINK_DEBUG*/

      if ( vic411_options.LAKES ) 
	lake_con = vic411_read_lakeparam(filep.lakeparam, soil_con, veg_con);

#endif // !OUTPUT_FORCE

      /** Build Gridded Filenames, and Open **/
      vic411_make_in_and_outfiles(&filep, &filenames, &soil_con, out_data_files);

      if (vic411_options.PRT_HEADER) {
        /** Write output file headers **/
        vic411_write_header(out_data_files, out_data, dmy, vic411_global_param);
      }

#if !OUTPUT_FORCE

      /** Read Elevation Band Data if Used **/
      vic411_read_snowband(filep.snowband, &soil_con);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = vic411_make_dist_prcp(veg_con[0].vegetat_type_num);

#endif // !OUTPUT_FORCE

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

      vic411_initialize_atmos(atmos, dmy, filep.forcing,
		       (double)soil_con.time_zone_lng, (double)soil_con.lng,
		       (double)soil_con.lat, soil_con.elevation,
		       soil_con.annual_prec, vic411_global_param.wind_h, 
		       soil_con.rough, soil_con.avgJulyAirTemp,
		       soil_con.Tfactor, 
#if OUTPUT_FORCE
		       soil_con.AboveTreeLine, out_data_files, out_data); 
#else /* OUTPUT_FORCE */
                       soil_con.AboveTreeLine); 
#endif /* OUTPUT_FORCE */

#if !OUTPUT_FORCE
#if LINK_DEBUG
      if(debug.PRT_ATMOS) vic411_write_atmosdata(atmos, vic411_global_param.nrecs);
#endif

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
      ErrorFlag = vic411_initialize_model_state(&prcp, dmy[0], &vic411_global_param, filep, 
			     soil_con.gridcel, veg_con[0].vegetat_type_num,
			     vic411_options.Nnode, Ndist, 
			     atmos[0].air_temp[vic411_NR],
			     &soil_con, veg_con, lake_con,
			     &init_STILL_STORM, &init_DRY_TIME, &save_data);
      if ( ErrorFlag == ERROR ) {
	if ( vic411_options.CONTINUEONERROR == TRUE ) {
	  // Handle grid cell solution error
	  fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	  break;
	} else {
	  // Else exit program on cell solution error as in previous versions
	  sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	  vic411_vicerror(ErrStr);
	}
      }
      
#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

      /** Update vic411_Error Handling Structure **/
      vic411_Error.filep = filep;
      vic411_Error.out_data_files = out_data_files;

      /** Initialize the storage terms in the water and energy balances **/
      /** Sending a negative record number (-vic411_global_param.nrecs) to vic411_dist_prec() will accomplish this **/
      ErrorFlag = vic411_dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
		  &lake_con, dmy, &vic411_global_param, &filep, out_data_files,
		  out_data, &save_data, -vic411_global_param.nrecs, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = startrec ; rec < vic411_global_param.nrecs; rec++ ) {

        if ( rec == vic411_global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;

        ErrorFlag = vic411_dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
		  &lake_con, dmy, &vic411_global_param, &filep,
		  out_data_files, out_data, &save_data, rec, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

        if ( ErrorFlag == ERROR ) {
          if ( vic411_options.CONTINUEONERROR == TRUE ) {
            // Handle grid cell solution error
            fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            break;
          } else {
	    // Else exit program on cell solution error as in previous versions
            sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            vic411_vicerror(ErrStr);
	  }
        }

        NEWCELL=FALSE;
	for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
	  init_DRY_TIME[veg] = -999;

      }	/* End Rec Loop */

#endif /* !OUTPUT_FORCE */

      vic411_close_files(&filep,out_data_files,&filenames); 

#if !OUTPUT_FORCE

#if QUICK_FS
      if(vic411_options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_layer[i][j]);
	  free((char *)soil_con.ufwc_table_layer[i]);
	}
	for(i=0;i<MAX_NODES;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_node[i][j]);
	  free((char *)soil_con.ufwc_table_node[i]);
	}
      }
#endif /* QUICK_FS */
      vic411_free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      vic411_free_vegcon(&veg_con);
      free((char *)soil_con.AreaFract);
      free((char *)soil_con.BandElev);
      free((char *)soil_con.Tfactor);
      free((char *)soil_con.Pfactor);
      free((char *)soil_con.AboveTreeLine);
      free((char*)init_STILL_STORM);
      free((char*)init_DRY_TIME);
#endif /* !OUTPUT_FORCE */
    }	/* End Run Model Condition */
  } 	/* End Grid Loop */

  /** cleanup **/
  vic411_free_atmos(vic411_global_param.nrecs, &atmos);
  vic411_free_dmy(&dmy);
  vic411_free_out_data_files(&out_data_files);
  vic411_free_out_data(&out_data);
#if !OUTPUT_FORCE
  vic411_free_veglib(&vic411_veg_lib);
#endif /* !OUTPUT_FORCE */
  fclose(filep.soilparam);
#if !OUTPUT_FORCE
  fclose(filep.vegparam);
  fclose(filep.veglib);
  if (vic411_options.SNOW_BAND>1)
    fclose(filep.snowband);
  if (vic411_options.LAKES)
    fclose(filep.lakeparam);
#endif /* !OUTPUT_FORCE */

  return EXIT_SUCCESS;
}	/* End Main Program */
//<devel -- debugging.patch>
#endif
//</devel -- debugging.patch>
