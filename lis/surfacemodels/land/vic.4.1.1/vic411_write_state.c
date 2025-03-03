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
#include <string.h>
#include "vic411_vicNl.h" 
#include "ftn.h"
#define VIC_BASED 0  /** VIC based vegetation tiling  Added by Shugong Wang on 05/01/2012 **/
#define LIS_BASED 1  /** LIS based vegetation tiling  Added by Shugong Wang on 05/01/2012 **/
//BOP
//
// !ROUTINE: vic411_write_state
// \label{vic411_write_state}
// 
// !REVISION HISTORY: 
//  09 Feb 2012  James Geiger; Initial implementation
//  01 May 2012  Shugong Wang; Add vegetation tiling scheme as an argument of
//               vic411_write_state.  
//  07 May 2012  Shugong Wang; Add vegclasses as an argument of
//               vic411_write_state to support MPI-safe tile id. 
// !INTERFACE:
void FTN(vic411_write_state)(int * LIS_NEST,
                             int * LIS_NTILES,
                             int *LIS_VT_SCHEME,
                             int *LIS_VEGCLASSES,  
                             int * year, int * month, int * day,
                             char * filename, int * lenf)
// !DESCRIPTION: 
//  This routine calls VIC's write\_model\_state routine to write
//  the restart state.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \item[LIS\_NTILES] number of tiles, considered intent(in)
//  \item[LIS\_VT\_SCHEME] vic411\_flag of vegetation tiling (in) \newline
//                         0 = VIC based tiling;     \newline
//                         1 = LIS based tiling 
//  \item[LIS\_VEGCLASSES] land use types 
//  \item[year] current year, considered intent(in)
//  \item[month] current month, considered intent(in)
//  \item[day] current day, considered intent(in)
//  \item[filename] name of the restart file, considered intent(in)
//  \item[lenf] length of the filename string (this is a hidden argument),
//              considered intent(in)
//  \end{description}
//
//EOP
{

   extern vic411_dist_prcp_struct     *vic411_lis_prcp;
   extern vic411_global_param_struct   vic411_global_param;
   extern vic411_veg_con_struct      **vic411_lis_veg_con;
   extern vic411_soil_con_struct      *vic411_lis_soil_con;
   extern vic411_filep_struct          vic411_filep;
   extern int                 **vic411_lis_init_DRY_TIME;
   extern char                **vic411_lis_init_STILL_STORM;
   extern vic411_lake_con_struct      *vic411_lis_lake_con;
   extern vic411_option_struct         vic411_options;
   extern vic411_filenames_struct      vic411_filenames;

   vic411_dist_prcp_struct * prcp;
   vic411_veg_con_struct   * veg_con;
   vic411_soil_con_struct  * soil_con;
   char             * STILL_STORM;
   int              * DRY_TIME;
   vic411_lake_con_struct  * lake_con;

   int i, t;
   int mpi_safe_tile_id; /** added by Shugong Wang 05/07/2012**/
   vic411_global_param.stateyear  = *year;
   vic411_global_param.statemonth = *month;
   vic411_global_param.stateday   = *day;

   strcpy(vic411_filenames.statefile, filename);

   if ( vic411_options.BINARY_STATE_FILE == TRUE )
   {
      // LIS-VIC requires the restart file to be written in ASCII
      printf("MSG: Resetting BINARY_STATE_FILE to FALSE\n");
      vic411_options.BINARY_STATE_FILE = FALSE;
   }

   vic411_filep.statefile = vic411_open_state_file(&vic411_global_param, vic411_filenames,
                                     vic411_options.Nlayer, vic411_options.Nnode);

   for ( i = 0; i < *LIS_NTILES; ++i )
   {
      prcp        = &vic411_lis_prcp[i];
      veg_con     = vic411_lis_veg_con[i];
      soil_con    = &vic411_lis_soil_con[i];
      STILL_STORM = vic411_lis_init_STILL_STORM[i];
      DRY_TIME    = vic411_lis_init_DRY_TIME[i];
      lake_con    = &vic411_lis_lake_con[i];
      if(*LIS_VT_SCHEME == VIC_BASED)
      {
          vic411_write_model_state(prcp, &vic411_global_param, veg_con[0].vegetat_type_num,
                            soil_con->gridcel, &vic411_filep, soil_con,  /** grid id is the key of restart file when using VIC tiling **/
                            STILL_STORM, DRY_TIME, *lake_con);
      }
      else
      {
          // using MPI-safe tile ID as the key of restart file Shugong Wang 05/07/2012 
          mpi_safe_tile_id = soil_con->gridcel * 100 + LIS_VEGCLASSES[i];
          vic411_write_model_state(prcp, &vic411_global_param, veg_con[0].vegetat_type_num,
                            mpi_safe_tile_id, &vic411_filep, soil_con,/**MPI-safe tile id is the key of restart file when using LIS tiling**/
                            STILL_STORM, DRY_TIME, *lake_con);
      }
   }

   fclose(vic411_filep.statefile);
   vic411_filep.statefile = NULL;
}
