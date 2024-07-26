//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include "vic412_atmos_forcing.h"
#include "ftn.h"

/* Added by Shugong Wang on 05/01/2012 to support reading restart file */
#define VIC_BASED 0  /* VIC based tiling */
#define LIS_BASED 1  /* LIS based tiling */


//BOP
//
// !ROUTINE: vic412_initialize
// \label{vic412_initialize}
// 
// !REVISION HISTORY: 
//  02 Aug 2011  James Geiger; Initial implementation
//  01 May 2012  Shugong Wang; add LIS_VT_SCHEME to pass vegetation tiling
//                             scheme to vic412_initialize.
//  07 May 2012  Shugong Wang; add argument vegclass to calculate MPI-safe
//                             tile id 
//  19 Jul 2013  Shugong Wang; add WRSI and SWI output variables for FLDAS 
//  20 May 2014  Shugong Wang; 
// !INTERFACE:
void FTN(vic412_initialize)(int *LIS_TILE, 
                            int *LIS_TS, 
                            int *LIS_VT_SCHEME, 
                            int *LIS_VEGCLASS,
                            int *LIS_NEST,
                            float *LIS_STATES,
                            int * LIS_START_MODE)

// !DESCRIPTION: 
//  This routine calls VIC's physics routines and diagnoses output.
//  It is based off VIC's vicNl.c file.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_TILE] sub-grid tile index, considered intent(in)
//  \item[LIS\_TS] time-step, considered intent(in)
//  \item[LIS\_VT\_SCHEME] vegetation tiling scheme: \newline
//                         0 = VIC based tiling; \newline
//                         1 = LIS based tiling 
//  \item[LIS\_VEGCLASS] land use type from LIS 
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \end{description}
//
//EOP
{
   extern veg_lib_struct       *veg_lib;
   extern veg_lib_struct      **lis_veg_lib;
   extern global_param_struct   global_param;
   extern option_struct         options;
   extern dmy_struct           *dmy;
   extern soil_con_struct      *lis_soil_con;
   extern atmos_data_struct   **lis_atmos;
   extern filep_struct          filep;
   extern out_data_file_struct **lis_out_data_files;
   extern out_data_struct      **lis_out_data;
   extern filep_struct          filep;
   extern save_data_struct     *lis_save_data;
   extern dist_prcp_struct     *lis_prcp;
   extern lake_con_struct      *lis_lake_con;
   extern veg_con_struct      **lis_veg_con;
   extern Error_struct          Error;
   extern int                 **lis_init_DRY_TIME;
   extern char                **lis_init_STILL_STORM;

   char                     NEWCELL;
   char                     LASTREC;
   char                    *init_STILL_STORM;
   char                     ErrStr[MAXSTRING];
   int                      rec;
   int                      veg;
   int                      Ndist;
   int                      cellnum;
   int                     *init_DRY_TIME;
   int                      ErrorFlag;

   atmos_data_struct       *atmos;
   veg_con_struct          *veg_con;
   soil_con_struct          soil_con;
   dist_prcp_struct         prcp; /* stores information about distributed 
                                     precipitation */
   lake_con_struct          lake_con;
   save_data_struct         save_data;
   dmy_struct              *temp_dmy;
   out_data_file_struct    *out_data_files;
   out_data_struct         *out_data;

   /** added to support restart **/
   int mpi_safe_tile_id; 
   
   if(*LIS_START_MODE == 1)
   {
    /*let vic read initial state*/
    options.INIT_STATE=1;
   }
   else
   {
    options.INIT_STATE=0;
   }

   /** Initialize Parameters **/
   if ( options.DIST_PRCP )
      Ndist = 2;
   else
      Ndist = 1;

   cellnum = *LIS_TILE-1;
   rec     = *LIS_TS-1;

   veg_lib = lis_veg_lib[cellnum];
   prcp = lis_prcp[cellnum];
   soil_con = lis_soil_con[cellnum];
   veg_con = lis_veg_con[cellnum];
   atmos = lis_atmos[cellnum];
   lake_con = lis_lake_con[cellnum];
   save_data = lis_save_data[cellnum];
   temp_dmy = &dmy[rec];
   out_data_files = lis_out_data_files[cellnum];
   out_data = lis_out_data[cellnum];

 
   /*set out_prec, out_rain, out_snow to be zero*/
   /*Shugong Wang 08/06/2012 */
   atmos[0].out_prec = 0.0;
   atmos[0].out_rain = 0.0; 
   atmos[0].out_snow = 0.0;


   if ( rec == 0 )
      NEWCELL = TRUE; // only used by the write_debug routine
   else
      NEWCELL = FALSE;

   if ( rec == 0 )
   {
      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

      /** If VIC based tiling is used, read the restart file with GRID ID as key. Shugong Wang 05/01/2012 **/
      if(*LIS_VT_SCHEME == VIC_BASED) 
      {
          // printf("reading restart file with grid id as key : %d\n", soil_con.gridcel);
          ErrorFlag = initialize_model_state(&prcp, dmy[0], &global_param, filep, 
                                soil_con.gridcel, veg_con[0].vegetat_type_num, /* soil_con.gridcel as the key of restart file*/
                                options.Nnode, Ndist, 
                                atmos[0].air_temp[NR],
                                &soil_con, veg_con, lake_con,
                                &init_STILL_STORM, &init_DRY_TIME, 
                                LIS_STATES);
      }
      else /** If LIS based tiling is used, read the restart file with TILE ID as key. Shugong Wang 05/01/2012 **/
      {
          /** MPI-safe tile ID doesn't rely on processor layout of LIS run. It is 
              calculated based on grid id and vegetation class. For current setting,
              mpi_tid = grid_id * 100 + vegclass, the maximum number of land use type
              is 100. If the total number of land use types exceeds 100, then the 
              the following statement should be revised accordingly. 
              Shugong Wang 05/07/2012                                            **/
          mpi_safe_tile_id = soil_con.gridcel * 100 + *LIS_VEGCLASS; 
          ErrorFlag = initialize_model_state(&prcp, dmy[0], &global_param, filep, 
                                mpi_safe_tile_id, veg_con[0].vegetat_type_num, /**mpi-safe tile id as the key of restart file */
                                options.Nnode, Ndist, 
                                atmos[0].air_temp[NR],
                                &soil_con, veg_con, lake_con,
                                &init_STILL_STORM, &init_DRY_TIME,
                                LIS_STATES);
      }

      lis_init_STILL_STORM[cellnum] = init_STILL_STORM;
      lis_init_DRY_TIME[cellnum] = init_DRY_TIME;
      lis_soil_con[cellnum] = soil_con;

      if ( ErrorFlag == ERROR )
      {
         if ( options.CONTINUEONERROR == TRUE )
         {
            // Handle grid cell solution error
            fprintf(stderr, "ERROR: Grid cell %i failed in record %i\
                            so the simulation has not finished.\
                            An incomplete output file has been generated,\
                            check your inputs before rerunning the simulation.\n", 
                            soil_con.gridcel, rec);
            exit(-1);
         }
         else
         {
            // Else exit program on cell solution error as in previous versions
            sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i \
                            so the simulation has ended. Check your inputs \
                            before rerunning the simulation.\n", soil_con.gridcel, rec);
            vicerror(ErrStr);
         }
      }
   }
}
