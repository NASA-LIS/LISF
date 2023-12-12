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
#include <vicNl.h>
#include <global.h>
#include "ftn.h"

#define VIC_BASED ( 0 )
#define LIS_BASED ( 1 )

extern void lis_scan_soilparam_flat_ascii(FILE *, float *, float *);

extern void vic412_compute_treeline(double, double *, char *);

soil_con_struct         *lis_soil_con;
dmy_struct              *dmy;
atmos_data_struct      **lis_atmos;
filep_struct             filep;
filenames_struct         filenames;
out_data_file_struct   **lis_out_data_files;
out_data_struct        **lis_out_data;
dist_prcp_struct        *lis_prcp; 
veg_con_struct         **lis_veg_con;
lake_con_struct         *lis_lake_con;
save_data_struct        *lis_save_data;
int                    **lis_init_DRY_TIME;
char                   **lis_init_STILL_STORM;
veg_lib_struct         **lis_veg_lib;

//soil_con_struct lis_read_soilparam_arc(FILE *, char *, float *, float *);
veg_con_struct *set_vegparam(int tile_idx,       // index of current tile 
                             int vegclass,       // vegetation class determined by LIS
                             int Nveg_type,      // number of vegation types in veg_lib
                             int gridcel,        // index of grid 
                             FILE *fp_vegparam,  // file pointer of vegetation parameter file
                             float *veg_fracs,   // fractions of all land cover determined by LIS
                             int Nlc_type);      // number of types in LIS land cover data


void count_model_state_412(dist_prcp_struct    *prcp,
                           global_param_struct *gp,
                           int                  Nveg,
                           int                  cellnum,
                           filep_struct        *filep,
                           soil_con_struct     *soil_con,
                           char                *STILL_STORM,
                           int                 *DRY_TIME,
                           int                 *count);
//BOP
//
// !ROUTINE: setup_vic412
// \label{setup_vic412}
// 
// !REVISION HISTORY: 
//  04 Aug 2011  James Geiger; Initial adaptation from VIC 4.1.1 (vicNl.c).
//  06 Jun 2012  Shugong Wang; Add support to VIC COMPUTE_TREELINE
// !INTERFACE:
void FTN(setup_vic412)(int *NGRIDS, int *NTILES, 
                       float *lat_array, float *lon_array,
                       int *t2gindex, int *vegclasses, 
                       int *vtscheme, int *NT,
                       float *veg_fracs_all,
                       int *syr, int *smn, int *sda, int *shr, 
                       int *eyr, int *emn, int *eda, int *ehr, 
                       float *ts, float *outInterval, 
                       int *vic412_state_chunk_size, 
                       char *global_param_filename,
                       int *dummy)
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
   extern veg_lib_struct *veg_lib;
#if LINK_DEBUG
  extern debug_struct debug;
#endif // LINK_DEBUG
   extern global_param_struct global_param;

   int                      Nveg_type;
   char                     MODEL_DONE;
   char                     RUN_MODEL;
   int                      startrec;

   int t,i,j;
   int first_tile;
  
   dist_prcp_struct *prcp;
   veg_con_struct   *veg_con;
   soil_con_struct  *soil_con;
   char             *STILL_STORM;
   int              *DRY_TIME;
   lake_con_struct  *lake_con;

   soil_con_struct         temp_soil_con;
   FILE *fp_vicveg;   
   int nvegs = *NT; 
   float **veg_fracs;
   // file pointer for writing tile id, grid id, and vegclass
   FILE *fp_t2g_mapping; 
   lis_prcp = (dist_prcp_struct *) calloc(*NTILES, sizeof(dist_prcp_struct));
   lis_soil_con = (soil_con_struct *) calloc(*NTILES, sizeof(soil_con_struct));
   lis_veg_con = (veg_con_struct **) calloc(*NTILES, sizeof(veg_con_struct*));
   lis_lake_con = (lake_con_struct *) calloc(*NTILES, sizeof(lake_con_struct));
   lis_save_data = (save_data_struct *) calloc(*NTILES, sizeof(save_data_struct));
   lis_atmos=(atmos_data_struct **) calloc(*NTILES, sizeof(atmos_data_struct*));
   lis_init_DRY_TIME = (int **) calloc(*NTILES, sizeof(int*));
   lis_init_STILL_STORM = (char **) calloc(*NTILES, sizeof(char*));
   lis_out_data_files = (out_data_file_struct **) calloc(*NTILES, sizeof(out_data_file_struct*));
   lis_out_data = (out_data_struct **) calloc(*NTILES, sizeof(out_data_struct*));
   lis_veg_lib = (veg_lib_struct **) calloc(*NTILES, sizeof(veg_lib_struct*));
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
   initialize_global();
   strcpy(filenames.global, global_param_filename);

   /** Read Global Control File, now all VIC global configurations are put
       in LIS configuration file. 05/07/2014 Shugong Wang **/
   filep.globalparam = open_file(filenames.global,"r");
   
   // set time step: from second to hour
   global_param.dt = (int)(*ts/3600.0);
   global_param.startyear = *syr;
   global_param.startmonth = *smn;
   global_param.startday = *sda;
   global_param.starthour = *shr;
   global_param.endyear = *eyr;
   global_param.endmonth = *emn;
   global_param.endday = *eda;
   /* read VIC configuration information from LIS configuration file*/
   get_global_param(&filenames, filep.globalparam, &global_param) ;
   
   
   /* over write the following settings with LIS settings*/
   global_param.out_dt = (int)(*outInterval/3600.0);
   
   /** Set up output data structures **/
   for ( t = 0; t < *NTILES; ++t )
   {
      lis_out_data[t] = create_output_list();
      lis_out_data_files[t] = set_output_defaults(lis_out_data[t]);
   }
   fclose(filep.globalparam);
   for ( t = 0; t < *NTILES; ++t )
   {
      filep.globalparam = open_file(filenames.global,"r");
      parse_output_info(&filenames, filep.globalparam, &lis_out_data_files[t], lis_out_data[t]);
   }

   /** Check and Open Files **/
   check_files(&filep, &filenames);

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

   for ( t = 0; t < *NTILES; ++t )
   {
      /** Read Vegetation Library File **/
      lis_veg_lib[t] = read_veglib(filep.veglib, &Nveg_type);
   }

   /** Make Date Data Structure **/
   dmy = make_dmy(&global_param);

   if ( options.INIT_STATE )
      filep.init_state = check_state_file(filenames.init_state, dmy,
                                          &global_param, options.Nlayer,
                                          options.Nnode, &startrec);

   for ( t = 0; t < *NTILES; ++t )
   {
      /** allocate memory for the atmos_data_struct **/
      alloc_atmos(1, &lis_atmos[t]);
   }

    /* read parameters */
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

      lis_scan_soilparam_flat_ascii(filep.soilparam,
                                    &lat_array[t],&lon_array[t]);
/*
      if ( (fscanf(filep.soilparam, "%d", &flag) ) != EOF)
      {
         if (flag)
            RUN_MODEL=TRUE;
         else
            RUN_MODEL=FALSE;
      }
      else
      {
         //MODEL_DONE = TRUE;
         RUN_MODEL = FALSE;
      }
*/
      // read current grid-cell
      veg_lib = lis_veg_lib[first_tile];
      //temp_soil_con = read_soilparam(filep.soilparam, RUN_MODEL);
      MODEL_DONE = FALSE; 
      temp_soil_con = read_soilparam(filep.soilparam, &RUN_MODEL, &MODEL_DONE);

      // assign to all tiles within the current grid-cell
      for ( t = 0; t < *NTILES; ++t )
      {
         if ( t2gindex[t] == i+1 ) // t2gindex counts from 1
         {
            lis_veg_lib[t] = veg_lib;
            lis_soil_con[t] = temp_soil_con;
         }
      }
   }

   for ( t = 0; t < *NTILES; ++t )
   {
#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if ( options.FROZEN_SOIL )
      {
         for ( i=0; i<MAX_LAYERS; i++ )
         {
            lis_soil_con[t].ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
            for ( j=0; j<QUICK_FS_TEMPS+1; j++ )
               lis_soil_con[t].ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
         }
         for ( i=0; i<MAX_NODES; i++ )
         {
            lis_soil_con[t].ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

         for ( j=0; j<QUICK_FS_TEMPS+1; j++ ) 
            lis_soil_con[t].ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
         }
      }
#endif /* QUICK_FS */
      //write_soilparam(&lis_soil_con[t]);
   }

/* commented by Shugong Wang
   if ( *vtscheme == LIS_BASED )
   {
      initialize_veg_lookup_table(*NT, filep.vegparam);
   }
*/ 
   /* tile to grid mapping , added by Shugong Wang for LIS restart. 04/30/2012 */
  // fp_t2g_mapping = fopen("tile2grid_mapping.txt", "w");
   for ( t = 0; t < *NTILES; ++t )
   {
      /** Build Gridded Filenames, and Open **/
      make_in_and_outfiles(&filep, &filenames, &lis_soil_con[t], lis_out_data_files[t]);

      /** Read Elevation Band Data if Used **/
      read_snowband(filep.snowband, &lis_soil_con[t]);
      
      // added by Shugong Wang to support COMPUTE_TREELINE option
      // Move to the position behand read_snowband becase COMPUTE_TREELINE option
      // require parameters of elevation band. 06/23/2012. Shugong Wang
      if( options.SNOW_BAND>1 && options.COMPUTE_TREELINE>0)
      {
          // routine vic412_compute_treeline will calculate the 
          // above tree line flag array in soil_con (AboveTreeLine)
          // according to average July air temperature and Tfactor, 
          // the change in temperature (C) due to elevation in each
          // snow band. If using LIS-VIC in southern hemisphere, the 
          // avgJulyAirTemp should be average December air temperature
          // vic412_compute_treeline doesn't support the computation of 
          // average July air temperature using retrospective forcing data.
          // avgJulyAirTemp should be supplied in the soil parameter file
          // of VIC.  Shugong Wang 06/06/2012
          vic412_compute_treeline(lis_soil_con[t].avgJulyAirTemp,
                                  lis_soil_con[t].Tfactor,
                                  lis_soil_con[t].AboveTreeLine);
       } 

      veg_lib = lis_veg_lib[t];

      if ( *vtscheme == VIC_BASED )
      {
         /** Read Grid Cell Vegetation Parameters **/
         lis_veg_con[t] = read_vegparam(filep.vegparam, lis_soil_con[t].gridcel,
                                        Nveg_type);
      }
      else
      {
         /* set_vegparam has been revised by Shugong Wang. After revision, this subroutine will 
            read root information (depth and fraction) from the vegetation parameter file of VIC. 
            Meanwhile, this subroutine also calls vic412_read_vegparam, which update LAI and related
            values when GLOBAL_LAI option is turned on.                                        */
         // add t (tile index) as one of arguments of set_vegparam, Shugong Wang, 06/06/2012
         lis_veg_con[t]=set_vegparam(t, vegclasses[t], Nveg_type,  lis_soil_con[t].gridcel, filep.vegparam,
                                     veg_fracs[t], *NT); // this is modified by Shugong Wang
      }
      calc_root_fractions(lis_veg_con[t], &lis_soil_con[t]);
      //write_vegparam(lis_veg_con[t]);

      if ( options.LAKES ) 
	      lis_lake_con[t] = read_lakeparam(filep.lakeparam, lis_soil_con[t], lis_veg_con[t]);

      /** Make Precipitation Distribution Control Structure **/
      lis_prcp[t] = make_dist_prcp(lis_veg_con[t][0].vegetat_type_num);
   }

   for ( t = 0; t < *NTILES; ++t )
   {
      veg_lib = lis_veg_lib[t];
#if LINK_DEBUG
      if (debug.PRT_SOIL) write_soilparam(&lis_soil_con[t]);
      if (debug.PRT_VEGE) write_vegparam(lis_veg_con[t],
                                       lis_soil_con[t].lat,lis_soil_con[t].lng);
      write_prcp(&lis_prcp[t], lis_veg_con[t][0].vegetat_type_num);
      write_lakecon(&lis_lake_con[t]);
#endif
   }


   /* determine the size of state chunk */
   // *vic412_state_chunk_size = options.SNOW_BAND*(20 + 2*options.Nlayer + 3*options.Nnode) + 10;
   prcp        = &lis_prcp[0];
   veg_con     = lis_veg_con[0];
   soil_con    = &lis_soil_con[0];
   STILL_STORM = lis_init_STILL_STORM[0];
   DRY_TIME    = lis_init_DRY_TIME[0];
   lake_con    = &lis_lake_con[0];
   
   count_model_state_412(prcp, &global_param, veg_con[0].vegetat_type_num,
                         soil_con->gridcel, &filep, soil_con, 
                         STILL_STORM, DRY_TIME, 
                         vic412_state_chunk_size);
}
