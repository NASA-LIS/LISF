//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.4
//
// Copyright (c) 2022 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vic411_vicNl.h"
#include "vic411_global.h"
#include "ftn.h"

#define VIC_BASED ( 0 )
#define LIS_BASED ( 1 )

extern void vic411_lis_scan_soilparam_flat_ascii(FILE *, float *, float *);

extern void vic411_compute_treeline_lis(double, double *, char *);

vic411_soil_con_struct vic411_lis_read_soilparam_arc(FILE *, char *, float *, float *);
vic411_veg_con_struct *vic411_set_vegparam(int tile_idx,       // index of current tile 
                             int vegclass,       // vegetation class determined by LIS
                             int Nveg_type,      // number of vegation types in vic411_veg_lib
                             int gridcel,        // index of grid 
                             FILE *fp_vegparam,  // file pointer of vegetation parameter file
                             float *veg_fracs,   // fractions of all land cover determined by LIS
                             int Nlc_type);      // number of types in LIS land cover data

vic411_soil_con_struct         *vic411_lis_soil_con;
vic411_dmy_struct              *vic411_dmy;
vic411_atmos_data_struct      **vic411_lis_atmos;
vic411_filep_struct             vic411_filep;
vic411_filenames_struct         vic411_filenames;
vic411_out_data_file_struct   **vic411_lis_out_data_files;
vic411_out_data_struct        **vic411_lis_out_data;
vic411_dist_prcp_struct        *vic411_lis_prcp; /* stores information about distributed 
                                      precipitation */
vic411_veg_con_struct         **vic411_lis_veg_con;
vic411_lake_con_struct         *vic411_lis_lake_con;
vic411_save_data_struct        *vic411_lis_save_data;
int                    **vic411_lis_init_DRY_TIME;
char                   **vic411_lis_init_STILL_STORM;
vic411_veg_lib_struct         **vic411_lis_veg_lib;

//BOP
//
// !ROUTINE: setup_vic411
// \label{setup_vic411}
// 
// !REVISION HISTORY: 
//  04 Aug 2011  James Geiger; Initial adaptation from VIC 4.1.1 (vicNl.c).
//  06 Jun 2012  Shugong Wang; Add support to VIC COMPUTE_TREELINE
// !INTERFACE:
void FTN(setup_vic411)(int *NGRIDS, int *NTILES, 
                       float *lat_array, float *lon_array,
                       int *t2gindex, int *vegclasses, 
                       int *vtscheme, int *NT,
                       float *veg_fracs_all,
                       char *global_param_filename, int *dummy)
// !DESCRIPTION: 
//  This routine performs variable initialization.
//  It was adapted from VIC 4.1.1.
//
//  For details about variables, input files and subroutines check: \newline
//  http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC\_home.html
//
//  UNITS: unless otherwise marked:          \newline
//  all water balance components are in mm   \newline
//  all energy balance components are in mks \newline
//  depths and lengths are in m              \newline
//
//  The arguments are:
//  \begin{description}
//  \item[NGRIDS] number of grid cells, considered intent(in)
//  \item[NTILES] number of sub-grid tiles, considered intent(in)
//  \item[lat\_array] array of latitude values for each grid cell,
//                    considered intent(in)
//  \item[lon\_array] array of longitude values for each grid cell,
//                    considered intent(in)
//  \item[t2gindex] array mapping tile indices to parent grid index,
//                  considered intent(in)
//  \item[vegclasses] array of vegetation classification values for each tile,
//                    considered intent(in)
//  \item[vtscheme] vegetation tiling scheme, considered intent(in)
//  \item[NT] total number of vegetation types for VIC lookup table,
//            considered intent(in)
//  \item[global\_param\_filename] the name of VIC configuration file
//  \item[dummy] dummy variable to capture length of global\_param\_filename
//  \end{description}
//
//EOP
{
   extern vic411_veg_lib_struct *vic411_veg_lib;
#if LINK_DEBUG
  extern vic411_debug_struct debug;
#endif // LINK_DEBUG
   extern vic411_global_param_struct vic411_global_param;

   int                      Nveg_type;
   int                      RUN_MODEL;
   int                      startrec;

   int t,i,j;
   int first_tile;

   vic411_soil_con_struct         temp_soil_con;
   FILE *fp_vicveg;   
   int nvegs = *NT; 
   float **veg_fracs;
   // file pointer for writing tile id, grid id, and vegclass
   FILE *fp_t2g_mapping; 
   vic411_lis_prcp = (vic411_dist_prcp_struct *) calloc(*NTILES, sizeof(vic411_dist_prcp_struct));
   vic411_lis_soil_con = (vic411_soil_con_struct *) calloc(*NTILES, sizeof(vic411_soil_con_struct));
   vic411_lis_veg_con = (vic411_veg_con_struct **) calloc(*NTILES, sizeof(vic411_veg_con_struct*));
   vic411_lis_lake_con = (vic411_lake_con_struct *) calloc(*NTILES, sizeof(vic411_lake_con_struct));
   vic411_lis_save_data = (vic411_save_data_struct *) calloc(*NTILES, sizeof(vic411_save_data_struct));
   vic411_lis_atmos=(vic411_atmos_data_struct **) calloc(*NTILES, sizeof(vic411_atmos_data_struct*));
   vic411_lis_init_DRY_TIME = (int **) calloc(*NTILES, sizeof(int*));
   vic411_lis_init_STILL_STORM = (char **) calloc(*NTILES, sizeof(char*));
   vic411_lis_out_data_files = (vic411_out_data_file_struct **) calloc(*NTILES, sizeof(vic411_out_data_file_struct*));
   vic411_lis_out_data = (vic411_out_data_struct **) calloc(*NTILES, sizeof(vic411_out_data_struct*));
   vic411_lis_veg_lib = (vic411_veg_lib_struct **) calloc(*NTILES, sizeof(vic411_veg_lib_struct*));
   veg_fracs = (float **) calloc(*NTILES, sizeof(float*));
  
   j = 0;
   for(t=0; t<*NTILES; t++)
    {
        veg_fracs[t] = (float *)calloc(*NT, sizeof(float));
        for(i=0; i<*NT; i++)
        {
            veg_fracs[t][i]=veg_fracs_all[j];
            j++;
        }
    }
   
   /** Read Model Options **/
   vic411_initialize_global();
   /* vic411_filenames = vic411_cmd_proc(argc, argv); */
   strcpy(vic411_filenames.global, global_param_filename);

   /** Read Global Control File **/
   vic411_filep.globalparam = vic411_open_file(vic411_filenames.global,"r");
   vic411_global_param = vic411_get_global_param(&vic411_filenames, vic411_filep.globalparam);

   /** Set up output data structures **/
   for ( t = 0; t < *NTILES; ++t )
   {
      vic411_lis_out_data[t] = vic411_create_output_list();
      vic411_lis_out_data_files[t] = vic411_set_output_defaults(vic411_lis_out_data[t]);
   }
   fclose(vic411_filep.globalparam);
   for ( t = 0; t < *NTILES; ++t )
   {
      vic411_filep.globalparam = vic411_open_file(vic411_filenames.global,"r");
      vic411_parse_output_info(&vic411_filenames, vic411_filep.globalparam, &vic411_lis_out_data_files[t], vic411_lis_out_data[t]);
   }

   /** Check and Open Files **/
   vic411_check_files(&vic411_filep, &vic411_filenames);

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

   for ( t = 0; t < *NTILES; ++t )
   {
      /** Read Vegetation Library File **/
      vic411_lis_veg_lib[t] = vic411_read_veglib(vic411_filep.veglib, &Nveg_type);
   }

   /** Make Date Data Structure **/
   vic411_dmy = vic411_make_dmy(&vic411_global_param);

   if ( vic411_options.INIT_STATE )
      vic411_filep.init_state = vic411_check_state_file(vic411_filenames.init_state, vic411_dmy,
                                          &vic411_global_param, vic411_options.Nlayer,
                                          vic411_options.Nnode, &startrec);

   for ( t = 0; t < *NTILES; ++t )
   {
      /** allocate memory for the vic411_atmos_data_struct **/
      vic411_alloc_atmos(1, &vic411_lis_atmos[t]);
   }

   if ( !vic411_options.ARC_SOIL )
   {
      for ( i = 0; i < *NGRIDS; ++i )
      {
         // find the first tile within the current grid-cell
         for ( t = 0; t < *NTILES; ++t )
         {
            if ( t2gindex[t] == i+1 ) // t2gindex counts from 1
            {
               break;
            }
         }
         first_tile = t;

         vic411_lis_scan_soilparam_flat_ascii(vic411_filep.soilparam,
                                       &lat_array[t],&lon_array[t]);

         if ( (fscanf(vic411_filep.soilparam, "%d", &vic411_flag) ) != EOF)
         {
            if (vic411_flag)
               RUN_MODEL=TRUE;
            else
               RUN_MODEL=FALSE;
         }
         else
         {
            //MODEL_DONE = TRUE;
            RUN_MODEL = FALSE;
         }

         // read current grid-cell
         vic411_veg_lib = vic411_lis_veg_lib[first_tile];
         temp_soil_con = vic411_read_soilparam(vic411_filep.soilparam, RUN_MODEL);

         // assign to all tiles within the current grid-cell
         for ( t = 0; t < *NTILES; ++t )
         {
            if ( t2gindex[t] == i+1 ) // t2gindex counts from 1
            {
               vic411_lis_veg_lib[t] = vic411_veg_lib;
               vic411_lis_soil_con[t] = temp_soil_con;
            }
         }
      }
   }
   else
   { 
      for ( t = 0; t < *NTILES; ++t )
      {
         vic411_veg_lib = vic411_lis_veg_lib[t];
         vic411_lis_soil_con[t] = vic411_lis_read_soilparam_arc(vic411_filep.soilparam, 
                                                  vic411_filenames.soil_dir,
                                                  &lat_array[t], &lon_array[t]);
      }
   }

   for ( t = 0; t < *NTILES; ++t )
   {
#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if ( vic411_options.FROZEN_SOIL )
      {
         for ( i=0; i<MAX_LAYERS; i++ )
         {
            vic411_lis_soil_con[t].ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
            for ( j=0; j<QUICK_FS_TEMPS+1; j++ )
               vic411_lis_soil_con[t].ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
         }
         for ( i=0; i<MAX_NODES; i++ )
         {
            vic411_lis_soil_con[t].ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

         for ( j=0; j<QUICK_FS_TEMPS+1; j++ ) 
            vic411_lis_soil_con[t].ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
         }
      }
#endif /* QUICK_FS */
      //vic411_write_soilparam(&vic411_lis_soil_con[t]);
   }

/* commented by Shugong Wang
   if ( *vtscheme == LIS_BASED )
   {
      initialize_veg_lookup_table(*NT, vic411_filep.vegparam);
   }
*/ 
   /* tile to grid mapping , added by Shugong Wang for LIS restart. 04/30/2012 */
  // fp_t2g_mapping = fopen("tile2grid_mapping.txt", "w");
   for ( t = 0; t < *NTILES; ++t )
   {
      /** Build Gridded Filenames, and Open **/
      vic411_make_in_and_outfiles(&vic411_filep, &vic411_filenames, &vic411_lis_soil_con[t], vic411_lis_out_data_files[t]);

      /** Read Elevation Band Data if Used **/
      vic411_read_snowband(vic411_filep.snowband, &vic411_lis_soil_con[t]);
      
      // added by Shugong Wang to support COMPUTE_TREELINE option
      // Move to the position behand vic411_read_snowband becase COMPUTE_TREELINE option
      // require parameters of elevation band. 06/23/2012. Shugong Wang
      if( vic411_options.SNOW_BAND>1 && vic411_options.COMPUTE_TREELINE>0)
      {
          // routine vic411_compute_treeline will calculate the 
          // above tree line vic411_flag array in soil_con (AboveTreeLine)
          // according to average July air temperature and Tfactor, 
          // the change in temperature (C) due to elevation in each
          // snow band. If using LIS-VIC in southern hemisphere, the 
          // avgJulyAirTemp should be average December air temperature
          // vic411_compute_treeline doesn't support the computation of 
          // average July air temperature using retrospective forcing data.
          // avgJulyAirTemp should be supplied in the soil parameter file
          // of VIC.  Shugong Wang 06/06/2012
          vic411_compute_treeline_lis(vic411_lis_soil_con[t].avgJulyAirTemp,
                                      vic411_lis_soil_con[t].Tfactor,
                                      vic411_lis_soil_con[t].AboveTreeLine);
       } 

      vic411_veg_lib = vic411_lis_veg_lib[t];

      if ( *vtscheme == VIC_BASED )
      {
         /** Read Grid Cell Vegetation Parameters **/
         vic411_lis_veg_con[t] = vic411_read_vegparam(vic411_filep.vegparam, vic411_lis_soil_con[t].gridcel,
                                        Nveg_type);
      }
      else
      {
         /* vic411_set_vegparam has been revised by Shugong Wang. After revision, this subroutine will 
            read root information (depth and fraction) from the vegetation parameter file of VIC. 
            Meanwhile, this subroutine also calls vic411_read_vegparam, which update LAI and related
            values when GLOBAL_LAI option is turned on.                                        */
         // add t (tile index) as one of arguments of vic411_set_vegparam, Shugong Wang, 06/06/2012
         vic411_lis_veg_con[t]=vic411_set_vegparam(t, vegclasses[t], Nveg_type,  vic411_lis_soil_con[t].gridcel, vic411_filep.vegparam,
                                     veg_fracs[t], *NT); // this is modified by Shugong Wang
        /* write tile to grid mapping */
        // fprintf(fp_t2g_mapping, "%8d %8d %2d\n", t+1, vic411_lis_soil_con[t].gridcel, vegclasses[t]); // in LIS, tile id starts from 1 instead 0
      }
      vic411_calc_root_fractions(vic411_lis_veg_con[t], &vic411_lis_soil_con[t]);
      //vic411_write_vegparam(vic411_lis_veg_con[t]);

      if ( vic411_options.LAKES ) 
	      vic411_lis_lake_con[t] = vic411_read_lakeparam(vic411_filep.lakeparam, vic411_lis_soil_con[t], vic411_lis_veg_con[t]);

      /** Make Precipitation Distribution Control Structure **/
      vic411_lis_prcp[t] = vic411_make_dist_prcp(vic411_lis_veg_con[t][0].vegetat_type_num);
   }
  // fclose(fp_t2g_mapping);

   for ( t = 0; t < *NTILES; ++t )
   {
      vic411_veg_lib = vic411_lis_veg_lib[t];
#if LINK_DEBUG
      if (debug.PRT_SOIL) vic411_write_soilparam(&vic411_lis_soil_con[t]);
      if (debug.PRT_VEGE) vic411_write_vegparam(vic411_lis_veg_con[t],
                                       vic411_lis_soil_con[t].lat,vic411_lis_soil_con[t].lng);
      write_prcp(&vic411_lis_prcp[t], vic411_lis_veg_con[t][0].vegetat_type_num);
      write_lakecon(&vic411_lis_lake_con[t]);
#endif
   }
}
