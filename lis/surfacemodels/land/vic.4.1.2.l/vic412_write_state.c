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
#include <vicNl.h>
#include "ftn.h"
#define VIC_BASED 0  /** VIC based vegetation tiling  Added by Shugong Wang on 05/01/2012 **/
#define LIS_BASED 1  /** LIS based vegetation tiling  Added by Shugong Wang on 05/01/2012 **/

extern void pack_model_state(dist_prcp_struct *, global_param_struct *, int, int, filep_struct *, soil_con_struct *, char *, int *, float *, int *);

//BOP
//
// !ROUTINE: vic412_write_state
// \label{vic412_write_state}
// 
// !REVISION HISTORY: 
//  09 Feb 2012  James Geiger; Initial implementation
//  01 May 2012  Shugong Wang; Add vegetation tiling scheme as an argument of
//                             vic412_write_state.  
//  07 May 2012  Shugong Wang; Add vegclasses as an argument of
//                             vic412_write_state to support MPI-safe tile id. 
//  15 May 2014  Shugong Wang; Use LIS 7 to deal with restart writing 
// !INTERFACE:
void FTN(vic412_write_state)(int * LIS_NEST,
                             int * LIS_NTILES,
                             int *LIS_VT_SCHEME,
                             int *LIS_VEGCLASSES,  
                             int * year, int * month, int * day,
                             //char * filename, int * lenf, 
                             float *state_chunk, int *tid)
// !DESCRIPTION: 
//  This routine calls VIC's write\_model\_state routine to write
//  the restart state.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \item[LIS\_NTILES] number of tiles, considered intent(in)
//  \item[LIS\_VT\_SCHEME] flag of vegetation tiling (in) \newline
//                         0 = VIC based tiling; \newline
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

   extern dist_prcp_struct     *lis_prcp;
   extern global_param_struct   global_param;
   extern veg_con_struct      **lis_veg_con;
   extern soil_con_struct      *lis_soil_con;
   extern filep_struct          filep;
   extern int                 **lis_init_DRY_TIME;
   extern char                **lis_init_STILL_STORM;
   extern lake_con_struct      *lis_lake_con;
   extern option_struct         options;
   extern filenames_struct      filenames;

   dist_prcp_struct * prcp;
   veg_con_struct   * veg_con;
   soil_con_struct  * soil_con;
   char             * STILL_STORM;
   int              * DRY_TIME;
   lake_con_struct  * lake_con;
   char filename[128]="dummy";
   int i, t;
   int mpi_safe_tile_id; /** added by Shugong Wang 05/07/2012**/
   int count = 0;

   global_param.stateyear  = *year;
   global_param.statemonth = *month;
   global_param.stateday   = *day;

   strcpy(filenames.statefile, filename);

  /* 
   if ( options.BINARY_STATE_FILE == TRUE )
   {
      // LIS-VIC requires the restart file to be written in ASCII
      printf("MSG: Resetting BINARY_STATE_FILE to FALSE\n");
      options.BINARY_STATE_FILE = FALSE;
   }

   filep.statefile = open_state_file(&global_param, filenames,
                                     options.Nlayer, options.Nnode);

   */
   i = *tid-1; 
   prcp        = &lis_prcp[i];
   veg_con     = lis_veg_con[i];
   soil_con    = &lis_soil_con[i];
   STILL_STORM = lis_init_STILL_STORM[i];
   DRY_TIME    = lis_init_DRY_TIME[i];
   lake_con    = &lis_lake_con[i];
   if(*LIS_VT_SCHEME == VIC_BASED)
   {
       pack_model_state(prcp, &global_param, veg_con[0].vegetat_type_num,
                         soil_con->gridcel, &filep, soil_con,  /** grid id is the key of restart file when using VIC tiling **/
                         STILL_STORM, DRY_TIME, /* *lake_con, */
                         state_chunk, &count);
   }
   else
   {
       // using MPI-safe tile ID as the key of restart file Shugong Wang 05/07/2012 
       mpi_safe_tile_id = soil_con->gridcel * 100 + LIS_VEGCLASSES[i];
       pack_model_state(prcp, &global_param, veg_con[0].vegetat_type_num,
                         mpi_safe_tile_id, &filep, soil_con,/**MPI-safe tile id is the key of restart file when using LIS tiling**/
                         STILL_STORM, DRY_TIME, /* *lake_con */
                         state_chunk, &count);
   }

   //fclose(filep.statefile);
   filep.statefile = NULL;
}
