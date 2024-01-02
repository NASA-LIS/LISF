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
#include <string.h>
#include "vic411_vicNl.h"
#include "vic411_atmos_forcing.h"
#include "ftn.h"

/* Added by Shugong Wang on 05/01/2012 to support reading restart file */
#define VIC_BASED 0  /* VIC based tiling */
#define LIS_BASED 1  /* LIS based tiling */

extern void FTN(vic411_diagnoseoutputvar)(int *, int *, int *, int *, double *, char *, int *, char *, int *);

//BOP
//
// !ROUTINE: vic411_run
// \label{vic411_run}
// 
// !REVISION HISTORY: 
//  02 Aug 2011  James Geiger; Initial implementation
//  01 May 2012  Shugong Wang; add LIS_VT_SCHEME to pass vegetation tiling
//                             scheme to vic411_run.
//  07 May 2012  Shugong Wang; add argument vegclass to calculate
//                             MPI-safe tile id 
//  19 Jul 2013  Shugong Wang; add WRSI and SWI output variables for FLDAS 
// !INTERFACE:
void FTN(vic411_run)(int *LIS_TILE, int *LIS_TS, int *LIS_VT_SCHEME,
                     int *LIS_VEGCLASS,
                     int *LIS_NEST,
                     int *LIS_MOC_EVAP,
                     int *LIS_MOC_QS,
                     int *LIS_MOC_QSB,
                     int *LIS_MOC_CANOPINT,
                     int *LIS_MOC_SMLIQFRAC,
                     int *LIS_MOC_SWNET,
                     int *LIS_MOC_ECANOP,
                     int *LIS_MOC_TVEG,
                     int *LIS_MOC_ESOIL,
                     int *LIS_MOC_ARESIST,
                     int *LIS_MOC_AVGSURFT,
                     int *LIS_MOC_ALBEDO,
                     int *LIS_MOC_LWNET,
                     int *LIS_MOC_SUBSNOW,
                     int *LIS_MOC_RADT,
                     int *LIS_MOC_QLE,
                     int *LIS_MOC_QH,
                     int *LIS_MOC_QG,
                     int *LIS_MOC_QF,
                     int *LIS_MOC_ACOND,
                     int *LIS_MOC_BARESOILT,
                     int *LIS_MOC_DELCOLDCONT,
                     int *LIS_MOC_DELINTERCEPT,
                     int *LIS_MOC_DELSOILMOIST,
                     int *LIS_MOC_DELSURFSTOR,
                     int *LIS_MOC_QFZ,
                     int *LIS_MOC_QSM,
                     int *LIS_MOC_QV,
                     int *LIS_MOC_RAINF,
                     int *LIS_MOC_ROOTMOIST,
                     int *LIS_MOC_SMFROZFRAC,
                     int *LIS_MOC_SNFRALBEDO,
                     int *LIS_MOC_SNOWCOVER,
                     int *LIS_MOC_SNOWDEPTH,
                     int *LIS_MOC_SNOWF,
                     int *LIS_MOC_SNOWT,
                     int *LIS_MOC_SOILMOIST,
                     int *LIS_MOC_SWDOWNFORC,
                     int *LIS_MOC_SWE,
                     int *LIS_MOC_TAIRFORC,
                     int *LIS_MOC_TOTALPRECIP,
                     int *LIS_MOC_VIC_PET_SATSOIL, // PET
                     int *LIS_MOC_VIC_PET_H2OSURF,
                     int *LIS_MOC_VIC_PET_SHORT,  
                     int *LIS_MOC_VIC_PET_TALL,   
                     int *LIS_MOC_VIC_PET_NATVEG, 
                     int *LIS_MOC_VIC_PET_VEGNOCR,
                     int *LIS_MOC_SOILTEMP,
                     int *LIS_MOC_WRSI, 
                     int *LIS_MOC_SWI)


// !DESCRIPTION: 
//  This routine calls VIC's physics routines and diagnoses output.
//  It is based off VIC's vicNl.c file.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_TILE] sub-grid tile index, considered intent(in)
//  \item[LIS\_TS] time-step, considered intent(in)
//  \item[LIS\_VT\_SCHEME] vegetation tiling scheme: \newline
//                         0 = VIC based tiling;     \newline
//                         1 = LIS based tiling 
//  \item[LIS\_VEGCLASS] land use type from LIS 
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \item[LIS\_MOC\_EVAP] output variable index for evaporation,
//     considered intent(in)
//  \item[LIS\_MOC\_QS] output variable index for surface vic411\_runoff, 
//     considered intent(in)
//  \item[LIS\_MOC\_QSB] output variable index for sub-surface vic411\_runoff,
//     considered intent(in)
//  \item[LIS\_MOC\_CANOPINT] output variable index for canopy interception,
//     considered intent(in)
//  \item[LIS\_MOC\_SMLIQFRAC] output variable index for soil moisture 
//     liquid fraction, considered intent(in)
//  \item[LIS\_MOC\_SWNET] output variable index for shortwave net radiation,
//     considered intent(in)
//  \item[LIS\_MOC\_ECANOP] output variable index for interception evaporation,
//     considered intent(in)
//  \item[LIS\_MOC\_TVEG] output variable index for vegetation vic411\_transpiration,
//     considered intent(in)
//  \item[LIS\_MOC\_ESOIL] output variable index for bare soil evaporation,
//     considered intent(in)
//  \item[LIS\_MOC\_ARESIST] output variable index for aerodynamic resistance,
//     considered intent(in)
//  \item[LIS\_MOC\_AVGSURFT] output variable index for average surface
//     temperature, considered intent(in)
//  \item[LIS\_MOC\_ALBEDO] output variable index for albedo,
//     considered intent(in)
//  \item[LIS\_MOC\_LWNET] output variable index for longwave net radiation,
//     considered intent(in)
//  \item[LIS\_MOC\_SUBSNOW] output variable index for snow sublimation,
//     considered intent(in)
//  \item[LIS\_MOC\_RADT] output variable index for surface radiative
//  temperature, considered intent(in)
//  \item[LIS\_MOC\_QLE] output variable index for latent heat flux,
//     considered intent(in)
//  \item[LIS\_MOC\_QH] output variable index for sensible heat flux,
//     considered intent(in)
//  \item[LIS\_MOC\_QG] output variable index for ground heat flux,
//     considered intent(in)
//  \item[LIS\_MOC\_QF] output variable index for energy of fusion,
//     considered intent(in)
//  \item[LIS\_MOC\_ACOND]  aerodynamic conductance  
//      considered intent(in)
//  \item[LIS\_MOC\_BARESOILT]  temperature of bare soil  
//      considered intent(in)
//  \item[LIS\_MOC\_DELCOLDCONT]  change in snow water content 
//      considered intent(in)
//  \item[LIS\_MOC\_DELINTERCEPT]  change in interception storage  
//      considered intent(in)
//  \item[LIS\_MOC\_DELSOILMOIST]  change in soil moisture content 
//      considered intent(in)
//  \item[LIS\_MOC\_DELSURFSTOR]  change in surface water storage  
//      considered intent(in)
//  \item[LIS\_MOC\_QFZ]  refreezing of water in the snowpack  
//      considered intent(in)
//  \item[LIS\_MOC\_QSM]  snowmelt (kg m-2 s-1) 
//      considered intent(in)
//  \item[LIS\_MOC\_QV]  energy of sublimation  
//      considered intent(in)
//  \item[LIS\_MOC\_RAINF]  rainfall rate  
//      considered intent(in)
//  \item[LIS\_MOC\_ROOTMOIST]  root zone soil moisture  
//      considered intent(in)
//  \item[LIS\_MOC\_SMFROZFRAC]  average layer fraction of frozen moisture 
//      considered intent(in)
//  \item[LIS\_MOC\_SNFRALBEDO]  albedo of snow-covered fraction of grid 
//      considered intent(in)
//  \item[LIS\_MOC\_SNOWCOVER]  fraction of snow-covered area in grid 
//      considered intent(in)
//  \item[LIS\_MOC\_SNOWDEPTH]  snow depth 
//      considered intent(in)
//  \item[LIS\_MOC\_SNOWF]  snowfall rate  
//      considered intent(in)
//  \item[LIS\_MOC\_SNOWT]  snow surface temperature 
//      considered intent(in)
//  \item[LIS\_MOC\_SOILMOIST]  aoil moisture contents in all soil layers 
//      considered intent(in)
//  \item[LIS\_MOC\_SWDOWNFORC]  downward shortwave radiation 
//      considered intent(in)
//  \item[LIS\_MOC\_SWE]  snow water equivalent  
//      considered intent(in)
//  \item[LIS\_MOC\_TAIRFORC]  air temperature 
//      considered intent(in)
//  \item[LIS\_MOC\_TOTALPRECIP]  total precipitation 
//      considered intent(in)
//  \end{description}
//
//EOP
{
   extern vic411_veg_lib_struct       *vic411_veg_lib;
   extern vic411_veg_lib_struct      **vic411_lis_veg_lib;
   extern vic411_global_param_struct   vic411_global_param;
   extern vic411_option_struct         vic411_options;
   extern vic411_dmy_struct           *vic411_dmy;
   extern vic411_soil_con_struct      *vic411_lis_soil_con;
   extern vic411_atmos_data_struct   **vic411_lis_atmos;
   extern vic411_out_data_file_struct **vic411_lis_out_data_files;
   extern vic411_out_data_struct      **vic411_lis_out_data;
   extern vic411_filep_struct          vic411_filep;
   extern vic411_save_data_struct     *vic411_lis_save_data;
   extern vic411_dist_prcp_struct     *vic411_lis_prcp;
   extern vic411_lake_con_struct      *vic411_lis_lake_con;
   extern vic411_veg_con_struct      **vic411_lis_veg_con;
   extern vic411_Error_struct          vic411_Error;
   extern int                 **vic411_lis_init_DRY_TIME;
   extern char                **vic411_lis_init_STILL_STORM;

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

   vic411_atmos_data_struct       *atmos;
   vic411_veg_con_struct          *veg_con;
   vic411_soil_con_struct          soil_con;
   vic411_dist_prcp_struct         prcp; /* stores information about distributed 
                                     precipitation */
   vic411_lake_con_struct          lake_con;
   vic411_save_data_struct         save_data;
   vic411_dmy_struct              *temp_dmy;
   vic411_out_data_file_struct    *out_data_files;
   vic411_out_data_struct         *out_data;

   int index;
   int level1, vlevel;
   double value;
   char *units;
   char *direction;
   int lenu;
   int lend;

   /* added by Shugong Wang for VIC output */
   float LS; /* latent heat of ice sublimation J kg-1*/
   float LV; /* latent heat of water evaporation J kg-1*/
   float layer_depth_mm[MAX_LAYERS];
   float total_depth_mm = 0.0;
   float out_dt_sec; /* output time interval in second*/
   /* end of code added by Shugong Wang for VIC output */

   /** added to support restart **/
   int mpi_safe_tile_id; 

   double pet;
   double et; 
   double total_smc;
   double wrsi;
   double swi;

   /** Initialize Parameters **/
   if ( vic411_options.DIST_PRCP )
      Ndist = 2;
   else
      Ndist = 1;

   cellnum = *LIS_TILE-1;
   rec     = *LIS_TS-1;

   vic411_veg_lib = vic411_lis_veg_lib[cellnum];
   prcp = vic411_lis_prcp[cellnum];
   soil_con = vic411_lis_soil_con[cellnum];
   veg_con = vic411_lis_veg_con[cellnum];
   atmos = vic411_lis_atmos[cellnum];
   lake_con = vic411_lis_lake_con[cellnum];
   save_data = vic411_lis_save_data[cellnum];
   temp_dmy = &vic411_dmy[rec];
   out_data_files = vic411_lis_out_data_files[cellnum];
   out_data = vic411_lis_out_data[cellnum];

 
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

#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
      /** If VIC based tiling is used, read the restart file with GRID ID as key. Shugong Wang 05/01/2012 **/
      if(*LIS_VT_SCHEME == VIC_BASED) 
      {
          // printf("reading restart file with grid id as key : %d\n", soil_con.gridcel);
          ErrorFlag = vic411_initialize_model_state(&prcp, vic411_dmy[0], &vic411_global_param, vic411_filep, 
                                soil_con.gridcel, veg_con[0].vegetat_type_num, /* soil_con.gridcel as the key of restart file*/
                                vic411_options.Nnode, Ndist, 
                                atmos[0].air_temp[vic411_NR],
                                &soil_con, veg_con, lake_con,
                                &init_STILL_STORM, &init_DRY_TIME, &save_data);
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
          ErrorFlag = vic411_initialize_model_state(&prcp, vic411_dmy[0], &vic411_global_param, vic411_filep, 
                                mpi_safe_tile_id, veg_con[0].vegetat_type_num, /**mpi-safe tile id as the key of restart file */
                                vic411_options.Nnode, Ndist, 
                                atmos[0].air_temp[vic411_NR],
                                &soil_con, veg_con, lake_con,
                                &init_STILL_STORM, &init_DRY_TIME, &save_data);
      }

      vic411_lis_init_STILL_STORM[cellnum] = init_STILL_STORM;
      vic411_lis_init_DRY_TIME[cellnum] = init_DRY_TIME;
      vic411_lis_soil_con[cellnum] = soil_con;

      if ( ErrorFlag == ERROR )
      {
         if ( vic411_options.CONTINUEONERROR == TRUE )
         {
            // Handle grid cell solution error
            fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            exit(-1);
         }
         else
         {
            // Else exit program on cell solution error as in previous versions
            sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
            vic411_vicerror(ErrStr);
         }
      }

#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

      /** Update vic411_Error Handling Structure **/
      vic411_Error.filep = vic411_filep;
      vic411_Error.out_data_files = out_data_files;

      /** Initialize the storage terms in the water and energy balances **/
      /** Sending a negative record number (-vic411_global_param.nrecs) to vic411_dist_prec() will accomplish this **/
      ErrorFlag = vic411_dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
                         &lake_con, vic411_dmy, &vic411_global_param, &vic411_filep, out_data_files,
                            out_data, &save_data, -vic411_global_param.nrecs, cellnum,
                            NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);
   }
   else
   {
      init_STILL_STORM = vic411_lis_init_STILL_STORM[cellnum];
      init_DRY_TIME = vic411_lis_init_DRY_TIME[cellnum];
   }

   LASTREC = FALSE; // not used
//   The following call of vic411_dist_prec keep using 0 as step index, which will cause
//   warning message in vic411_put_data.c for TreeAdjustFactor.  Shugong Wang 06/08/2012
   ErrorFlag = vic411_dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
                         &lake_con, temp_dmy, &vic411_global_param, &vic411_filep,
                         out_data_files, out_data, &save_data, 0, cellnum,
                         NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);
// The following call of vic411_dist_prec uses rec = *LIS_TS - 1
/*
   ErrorFlag = vic411_dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
                         &lake_con, temp_dmy, &vic411_global_param, &vic411_filep,
                         out_data_files, out_data, &save_data, rec, cellnum,
                         NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);
*/
   if ( ErrorFlag == ERROR )
   {
      if ( vic411_options.CONTINUEONERROR == TRUE )
      {
         // Handle grid cell solution error
         fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
         exit(-1);
      }
      else
      {
         // Else exit program on cell solution error as in previous versions
         sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
         vic411_vicerror(ErrStr);
      }
   }

   for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
      init_DRY_TIME[veg] = -999;

   // diagnose output variables
   level1 = 1;

// Forcing
//   units = "kg m-2 s-1"; lenu = strlen(units);
//   direction = "DN"; lend = strlen(direction);
//   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_RAINFFORC,
//                                 &level1,
//                                 &(out_data[OUT_PREC].data[0]),
//                                 units,&lenu,direction,&lend);
/*
   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_EVAP,
                                 &level1,
                                 &(out_data[OUT_EVAP].data[0]),
                                 units,&lenu,direction,&lend);

   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "OUT"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QS,
                                 &level1,
                                 &(out_data[OUT_RUNOFF].data[0]),
                                 units,&lenu,direction,&lend);

   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "OUT"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QSB,
                                 &level1,
                                 &(out_data[OUT_BASEFLOW].data[0]),
                                 units,&lenu,direction,&lend);

   units = "kg m-2"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_CANOPINT,
                                 &level1,
                                 &(out_data[OUT_WDEW].data[0]),
                                 units,&lenu,direction,&lend);

   units = "-"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   for( index=0; index<vic411_options.Nlayer; ++index )
   {
      vlevel = index+1;
      FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SMLIQFRAC,
                                    &vlevel,
                                    &(out_data[OUT_SOIL_LIQ].data[index]),
                                    units,&lenu,direction,&lend);
   }

   units = "W m-2"; lenu = strlen(units);
   direction = "DN"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWNET,
                                 &level1,
                                 &(out_data[OUT_NET_SHORT].data[0]),
                                 units,&lenu,direction,&lend);

   units = "W m-2"; lenu = strlen(units);
   direction = "DN"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_LWNET,
                                 &level1,
                                 &(out_data[OUT_NET_LONG].data[0]),
                                 units,&lenu,direction,&lend);

   // OUT_R_NET = OUT_NET_SHORT + OUT_NET_LONG
   //units = "W m-2"; lenu = strlen(units);
   //direction = "UP"; lend = strlen(direction);
   //FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_?,
   //                              &level1,
   //                              &(out_data[OUT_R_NET].data[0]),
   //                              units,&lenu,direction,&lend);

   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ECANOP,
                                 &level1,
                                 &(out_data[OUT_EVAP_CANOP].data[0]),
                                 units,&lenu,direction,&lend);

   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TVEG,
                                 &level1,
                                 &(out_data[OUT_TRANSP_VEG].data[0]),
                                 units,&lenu,direction,&lend);

   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ESOIL,
                                 &level1,
                                 &(out_data[OUT_EVAP_BARE].data[0]),
                                 units,&lenu,direction,&lend);

   //units = "W m-2"; lenu = strlen(units);
   //direction = "UP"; lend = strlen(direction);
   //FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_?,
   //                              &level1,
   //                              &(out_data[OUT_SUB_CANOP].data[0]),
   //                              units,&lenu,direction,&lend);

   units = "mm"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW,
                                 &level1,
                                 &(out_data[OUT_SUB_SNOW].data[0]),
                                 units,&lenu,direction,&lend);

   value = out_data[OUT_SUB_SNOW].data[0]/vic411_global_param.dt;
   units = "kg m-2 s-1"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW,
                                 &level1,
                                 &value,
                                 units,&lenu,direction,&lend);

   units = "s/m"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ARESIST,
                                 &level1,
                                 &(out_data[OUT_AERO_RESIST].data[0]),
                                 units,&lenu,direction,&lend);

   units = "K"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_AVGSURFT,
                                 &level1,
                                 &(out_data[OUT_SURF_TEMP].data[0]),
                                 units,&lenu,direction,&lend);

   units = "-"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ALBEDO,
                                 &level1,
                                 &(out_data[OUT_ALBEDO].data[0]),
                                 units,&lenu,direction,&lend);

   //units = "W m-2"; lenu = strlen(units);
   //direction = "UP"; lend = strlen(direction);
   //FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_?,
   //                              &level1,
   //                              &(out_data[OUT_REL_HUMID].data[0]),
   //                              units,&lenu,direction,&lend);

   //units = "W m-2"; lenu = strlen(units);
   //direction = "UP"; lend = strlen(direction);
   //FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_?,
   //                              &level1,
   //                              &(out_data[OUT_IN_LONG].data[0]),
   //                              units,&lenu,direction,&lend);

// Forcing
//   units = "K"; lenu = strlen(units);
//   direction = "-"; lend = strlen(direction);
//   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TAIRFORC,
//                                 &level1,
//                                 &(out_data[OUT_AIR_TEMP].data[0]),
//                                 units,&lenu,direction,&lend);

// Forcing
//   units = "m s-1"; lenu = strlen(units);
//   direction = "-"; lend = strlen(direction);
//   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_WINDFORC,
//                                 &level1,
//                                 &(out_data[OUT_WIND].data[0]),
//                                 units,&lenu,direction,&lend);

   units = "K"; lenu = strlen(units);
   direction = "-"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_RADT,
                                 &level1,
                                 &(out_data[OUT_RAD_TEMP].data[0]),
                                 units,&lenu,direction,&lend);

   units = "W m-2"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QLE,
                                 &level1,
                                 &(out_data[OUT_LATENT].data[0]),
                                 units,&lenu,direction,&lend);

   units = "W m-2"; lenu = strlen(units);
   direction = "UP"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QH,
                                 &level1,
                                 &(out_data[OUT_SENSIBLE].data[0]),
                                 units,&lenu,direction,&lend);

   units = "W m-2"; lenu = strlen(units);
   direction = "DN"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QG,
                                 &level1,
                                 &(out_data[OUT_GRND_FLUX].data[0]),
                                 units,&lenu,direction,&lend);

   units = "W m-2"; lenu = strlen(units);
   direction = "S2L"; lend = strlen(direction);
   FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QF,
                                 &level1,
                                 &(out_data[OUT_FUSION].data[0]),
                                 units,&lenu,direction,&lend);
                                 */
                                 
    
    for(index=0; index<vic411_options.Nlayer; index++)
    {
        layer_depth_mm[index] = soil_con.depth[index] * 1000.0; /* m => mm */
        total_depth_mm += layer_depth_mm[index];
    }
    
    if(vic411_global_param.out_dt == 0)
    {
        out_dt_sec = vic411_global_param.dt * 3600.0;      
    }
    else
    {
        out_dt_sec = vic411_global_param.out_dt * 3600.0;
    }
    
    /** VIC description: "Scene" aerodynamic conductance (tiles with overstory contribute
        overstory conductance; others contribute surface conductance)                    **/
    /** LIS description: Aerodynamic conductance (m s-1) **/
    /*** VIC unit: m s-1 => LIS unit: m s-1 ***/
    units = "m s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ACOND, 
                                  &level1, 
                                  &(out_data[OUT_AERO_COND].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Average surface albedo **/
    /** LIS description: Surface Albedo (-) **/
    /*** VIC unit: - => LIS unit: - ***/
    units = "-"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ALBEDO, 
                                  &level1, 
                                  &(out_data[OUT_ALBEDO].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: "Scene"canopy aerodynamic resistance (tiles with overstory contribute 
        over story resistance; others contribute surface resistance)                       **/
    /** LIS description: no comment in LIS **/
    /*** VIC unit: s/m => LIS unit: s/m ***/
/*    units = ""; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ARESIST, 
                                  &level1, 
                                  &(out_data[OUT_AERO_RESIST].data[0] ),
                                  units,&lenu,direction,&lend);
*/
    /** VIC description: Average surface temperature **/
    /** LIS description: Average Surface Temperature (K) **/
    /*** VIC unit: C => LIS unit: K ***/
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SURF_TEMP].data[0] + VIC_KELVIN;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_AVGSURFT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Bare soil surface temperature **/
    /** LIS description: Temperature of bare soil (K) **/
    /*** VIC unit: C => LIS unit: K ***/
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_BARESOILT].data[0] + VIC_KELVIN;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_BARESOILT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Total moisture interception storage in canopy **/
    /** LIS description:  Total canopy water storage (kg m-2) **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_CANOPINT, 
                                  &level1, 
                                  &out_data[OUT_WDEW].data[0] ,
                                  units,&lenu,direction,&lend);

    /** VIC description: Change in snow water equivalent **/
    /** LIS description: Change in snow water content (J m-2) **/
    /*** VIC unit: mm => LIS unit: J m-2 ***/ /* Heat of water of fusion: 334 KJ/Kg*/
    units = "J m-2"; lenu = strlen(units);
    direction = "INC"; lend = strlen(direction);
    value = out_data[OUT_DELSWE].data[0] * 334.0 * 1000 / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_DELCOLDCONT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Change in canopy interception storage **/
    /** LIS description: Change in interception storage (kg m-2) **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "INC"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_DELINTERCEPT, 
                                  &level1, 
                                  &(out_data[OUT_DELINTERCEPT].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Change in soil water content **/
    /** LIS description: DelSoilMoist (kg m-2) **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "INC"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_DELSOILMOIST, 
                                  &level1, 
                                  &(out_data[OUT_DELSOILMOIST].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Change in surface liquid water storage **/
    /** LIS description: Change in surface water storage (kg m-2) **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "INC"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_DELSURFSTOR, 
                                  &level1, 
                                  &(out_data[OUT_DELSURFSTOR].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net evaporation from canopy interception **/
    /** LIS description: Interception evaporation **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP_CANOP].data[0] / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ECANOP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm hr-1 ***/
    units = "mm hr-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP_CANOP].data[0] / (out_dt_sec / 3600.0);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ECANOP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: W m-2 ***/ /* mm=kg m-2; 2257*1000 J kg-1 * kg m-2 = J m-2; J m-2 / second = W m-2 */
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    LV = 2501000 - 2361 * out_data[OUT_AIR_TEMP].data[0];
    value = out_data[OUT_EVAP_CANOP].data[0] * LV /out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ECANOP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Net evaporation from bare soil **/
    /** LIS description: Bare soil evaporation **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP_BARE].data[0] / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ESOIL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm hr-1 ***/
    units = "mm hr-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP_BARE].data[0] / (out_dt_sec/3600.0); 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ESOIL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: W m-2 ***/ 
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    LV = 2501000 - 2361 * out_data[OUT_AIR_TEMP].data[0];
    value = out_data[OUT_EVAP_BARE].data[0] * LV / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ESOIL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Total net evaporation **/
    /** LIS description: Evapotranspiration **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_EVAP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    LV = 2501000 - 2361 * out_data[OUT_AIR_TEMP].data[0];
    value = out_data[OUT_EVAP].data[0] * LV / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_EVAP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm hr-1 ***/
    units = "mm hr-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_EVAP].data[0] / (out_dt_sec / 3600.0);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_EVAP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);


    /*PET outputs, added by Shugong Wang 09/03/212 */

    /** VIC description: Potential evap from saturated bare soil **/
    /** LIS description: Potential evap from saturated bare soil **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_SATSOIL].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_SATSOIL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm) ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_SATSOIL].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_SATSOIL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Potential evap from open water **/
    /** LIS description: Potential evap from open water **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_H2OSURF].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_H2OSURF, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm) ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_H2OSURF].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_H2OSURF, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description:  Potential evap (vic411_transpiration only) from short reference crop (grass)**/
    /** LIS description:  Potential evap (vic411_transpiration only) from short reference crop (grass)**/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_SHORT].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_SHORT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm) ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_SHORT].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_SHORT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);
    /** VIC description:  Potential evap (vic411_transpiration only) from tall reference crop (alfalfa)**/
    /** LIS description:  Potential evap (vic411_transpiration only) from tall reference crop (alfalfa)**/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_TALL].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_TALL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm) ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_TALL].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_TALL, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);
    /** VIC description: Potential evap (vic411_transpiration only) from current vegetation and current canopy resistance **/
    /** LIS description: Potential evap (vic411_transpiration only) from current vegetation and current canopy resistance **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_NATVEG].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_NATVEG, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm)  ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_NATVEG].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_NATVEG, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);
    /** VIC description: Potential evap (vic411_transpiration only) from current vegetation and 0 canopy resistance **/
    /** LIS description: Potential evap (vic411_transpiration only) from current vegetation and 0 canopy resistance **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_VEGNOCR].data[0] / out_dt_sec ;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_VEGNOCR, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 (mm) ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_PET_VEGNOCR].data[0];
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_VIC_PET_VEGNOCR, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Net downward longwave flux **/
    /** LIS description: Net longwave radiation (surface) (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_LWNET, 
                                  &level1, 
                                  &(out_data[OUT_NET_LONG].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net energy used to melt/freeze soil moisture **/
    /** LIS description: Energy of fusion (W m-2) ( Energy consumed or released during liquid/solid phase changes) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "S2L"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QF, 
                                  &level1, 
                                  &(out_data[OUT_FUSION].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Refreezing of water in the snow **/
    /** LIS description: Refreezing of water in the snowpack (kg m-2 s-1) **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "L2S"; lend = strlen(direction);
    value = out_data[OUT_REFREEZE].data[0] / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QFZ, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Net heat flux into ground **/
    /** LIS description: Ground Heat Flux (W m-2) (Heat flux into the ground, averaged over a grid cell ) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QG, 
                                  &level1, 
                                  &(out_data[OUT_GRND_FLUX].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net upward sensible heat flux **/
    /** LIS description: Sensible Heat Flux (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QH, 
                                  &level1, 
                                  &(out_data[OUT_SENSIBLE].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net upward latent heat flux **/
    /** LIS description: Latent Heat Flux (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QLE, 
                                  &level1, 
                                  &(out_data[OUT_LATENT].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Surface vic411_runoff **/
    /** LIS description: Surface Runoff **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "OUT"; lend = strlen(direction);
    value = out_data[OUT_RUNOFF].data[0] / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QS, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "OUT"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QS, 
                                  &level1, 
                                  &(out_data[OUT_RUNOFF].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Baseflow out of the bottom layer **/
    /** LIS description: Subsurface Runoff **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "OUT"; lend = strlen(direction);
    value = out_data[OUT_BASEFLOW].data[0] /out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QSB, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "OUT"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QSB, 
                                  &level1, 
                                  &(out_data[OUT_BASEFLOW].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Snow melt **/
    /** LIS description: Snowmelt  **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "S2L"; lend = strlen(direction);
    value = out_data[OUT_SNOW_MELT].data[0] /out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QSM, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "S2L"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QSM, 
                                  &level1, 
                                  &(out_data[OUT_SNOW_MELT].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net upward latent heat flux from sublimation **/
    /** LIS description: Energy of sublimation (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_QV, 
                                  &level1, 
                                  &(out_data[OUT_LATENT_SUB].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Average radiative surface temperature **/
    /** LIS description: Surface Radiative Tempearture (K) **/
    /*** VIC unit: K => LIS unit: K ***/
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_RADT, 
                                  &level1, 
                                  &(out_data[OUT_RAD_TEMP].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Rainfall **/
    /** LIS description: Rainfall rate **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    value = out_data[OUT_RAINF].data[0] / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_RAINF, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_RAINF, 
                                  &level1, 
                                  &(out_data[OUT_RAINF].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Root zone soil moisture **/
    /** LIS description: Root zone soil moisture **/
    /*** VIC unit: mm => LIS unit: m^3 m-3 ***/
    units = "m^3 m-3"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_ROOTMOIST].data[0] / total_depth_mm;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ROOTMOIST, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_ROOTMOIST, 
                                  &level1, 
                                  &(out_data[OUT_ROOTMOIST].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Fraction of soil moisture (by mass) that is ice, for each soil layer **/
    /** LIS description: Average layer fraction of frozen moisture **/
    /*** VIC unit: - => LIS unit: - ***/
    units = "-"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1; 
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SMFROZFRAC, 
                                     &vlevel, 
                                     &(out_data[OUT_SMFROZFRAC].data[index] ), 
                                     units,&lenu,direction,&lend);
    } 

    /*** VIC unit: - => LIS unit: m^3 m-3 ***/
    units = "m^3 m-3"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1;
       value =  out_data[OUT_SMFROZFRAC].data[index] * out_data[OUT_SOIL_MOIST].data[index] / layer_depth_mm[index];
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SMFROZFRAC, 
                                     &vlevel, 
                                     &value, 
                                     units,&lenu,direction,&lend);
    } 

    /** VIC description: Fraction of soil moisture (by mass) that is liquid, for each soil layer **/
    /** LIS description: Average layer fraction of liquid moisture **/
    /*** VIC unit: - => LIS unit: - ***/
    units = "-"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1; 
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SMLIQFRAC, 
                                     &vlevel, 
                                     &(out_data[OUT_SMLIQFRAC].data[index]), 
                                     units,&lenu,direction,&lend);
    } 

    /*** vic unit: - => lis unit: m^3 m-3 ***/
    units = "m^3 m-3"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1;
        /* BE AWARE: THE UNIT OF VIC SOIL MOISTURE IS mm*/
       value = out_data[OUT_SMLIQFRAC].data[index] * out_data[OUT_SOIL_MOIST].data[index] / layer_depth_mm[index]; 
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SMLIQFRAC, 
                                     &vlevel, 
                                     &value, 
                                     units,&lenu,direction,&lend);
    } 

    /** VIC description: Snow pack albedo **/
    /** LIS description: snow free albedo **/
    /*** VIC unit: - => LIS unit: - ***/
    units = "-"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNFRALBEDO, 
                                  &level1, 
                                  &(out_data[OUT_SALBEDO].data[0] ),
                                  units,&lenu,direction,&lend);

    /** VIC description: Fractional area of snow cover **/
    /** LIS description: fraction grid cell covered by snow **/
    /*** VIC unit: - => LIS unit: - ***/
    units = "-"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWCOVER, 
                                  &level1, 
                                  &(out_data[OUT_SNOW_COVER].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Depth of snow pack **/
    /** LIS description: Snow Depth **/
    /*** VIC unit: cm => LIS unit: m ***/
    units = "m"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SNOW_DEPTH].data[0] * 0.01;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWDEPTH, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: cm => LIS unit: cm ***/
    units = "cm"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SNOW_DEPTH].data[0] * 1.0;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWDEPTH, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: cm => LIS unit: mm ***/
    units = "mm"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SNOW_DEPTH].data[0] * 10.0;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWDEPTH, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Snowfall **/
    /** LIS description: Snowfall rate **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    value = out_data[OUT_SNOWF].data[0] / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWF, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWF, 
                                  &level1, 
                                  &(out_data[OUT_SNOWF].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Snow surface temperature **/
    /** LIS description: Snow surface temperature (K) **/
    /*** VIC unit: C => LIS unit: K ***/
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SNOW_SURF_TEMP].data[0] + VIC_KELVIN;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SNOWT, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Total soil moisture content for each soil layer **/
    /** LIS description: Soil moisture contents in all soil layers (m^3 m-3, kg m-2) **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1; 
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SOILMOIST, 
                                     &vlevel, 
                                     &(out_data[OUT_SOIL_MOIST].data[index]), 
                                     units,&lenu,direction,&lend);
    } 

    /*** VIC unit: mm => LIS unit: m^3 m-3 ***/
    units = "m^3 m-3"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1; 
       value = out_data[OUT_SOIL_MOIST].data[index] / layer_depth_mm[index];
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SOILMOIST, 
                                     &vlevel, 
                                     &value, 
                                     units,&lenu,direction,&lend);
    } 
    
    /*** SWI ***/
    units = "%";     lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    swi = 0.0;
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       value = out_data[OUT_SOIL_MOIST].data[index] / layer_depth_mm[index];
       swi += value / soil_con.porosity[index];
    } 
    swi /= vic411_options.Nlayer;
    swi *= 100.0;  
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWI, 
                                  &level1, 
                                  &swi, 
                                  units,&lenu,direction,&lend);

                                  
    /*** WRSI ***/
    units = "-";     lenu = strlen(units);
    direction  = "-"; lend = strlen(direction);
    et         = out_data[OUT_EVAP].data[0];
    pet        = out_data[OUT_PET_VEGNOCR].data[0];
    
    if(pet>1e-6 && et>1e-8)
    {
        wrsi        = et / pet;
    }
    else
    {
        wrsi = 0.0; 
    }
    value = wrsi * 100.0; 

    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_WRSI, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Total net sublimation from snow pack (surface and blowing) **/
    /** LIS description: Snow sublimation  **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value =  out_data[OUT_SUB_SNOW].data[0] / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm hr-1 ***/
    units = "mm hr-1"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SUB_SNOW].data[0] / (out_dt_sec/3600.0);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    LS = (677 - 0.07 * out_data[OUT_SNOW_SURF_TEMP].data[0]) * 4.1868 * 1000; /* J kg-1 */
    value = out_data[OUT_SUB_SNOW].data[0] * LS / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm ***/
    units = "mm"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SUBSNOW, 
                                  &level1, 
                                  &(out_data[OUT_SUB_SNOW].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Incoming shortwave **/
    /** LIS description: downward shortwave radiation (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    /*
    units = "W m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWDOWNFORC, 
                                  &level1, 
                                  &(out_data[OUT_SHORTWAVE].data[0]),
                                  units,&lenu,direction,&lend);
    */
    /** VIC description: Snow water equivalent in snow pack (including vegetation-intercepted snow) **/
    /** LIS description: Snow water equivalent **/
    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWE, 
                                  &level1, 
                                  &(out_data[OUT_SWE].data[0] ),
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: m ***/
    units = "m"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_SWE].data[0] / 1000.0; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWE, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** VIC description: Net downward shortwave flux **/
    /** LIS description: Net shortwave radiation (surface) (W m-2) **/
    /*** VIC unit: W m-2 => LIS unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SWNET, 
                                  &level1, 
                                  &(out_data[OUT_NET_SHORT].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Air temperature **/
    /** LIS description: Air temperature (K) **/
    /*** VIC unit: C => LIS unit: K ***/
    /*
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    value = out_data[OUT_AIR_TEMP].data[0] + VIC_KELVIN; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TAIRFORC, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);
    */
    /** VIC description: Incoming precipitation **/
    /** LIS description: Total precipitation rate  **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    value = out_data[OUT_PREC].data[0] / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TOTALPRECIP, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: kg m-2 ***/
    units = "kg m-2"; lenu = strlen(units);
    direction = "DN"; lend = strlen(direction);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TOTALPRECIP, 
                                  &level1, 
                                  &(out_data[OUT_PREC].data[0]),
                                  units,&lenu,direction,&lend);

    /** VIC description: Net vic411_transpiration from vegetation **/
    /** LIS description: Vegetation vic411_transpiration  **/
    /*** VIC unit: mm => LIS unit: kg m-2 s-1 ***/
    units = "kg m-2 s-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_TRANSP_VEG].data[0] / out_dt_sec; 
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TVEG, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** VIC unit: mm => LIS unit: mm hr-1 ***/
    units = "mm hr-1"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    value = out_data[OUT_TRANSP_VEG].data[0] / (out_dt_sec/3600.0);
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TVEG, 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /*** vic unit: mm => lis unit: W m-2 ***/
    units = "W m-2"; lenu = strlen(units);
    direction = "UP"; lend = strlen(direction);
    LV = 2501000 - 2361 * out_data[OUT_AIR_TEMP].data[0];
    value = out_data[OUT_TRANSP_VEG].data[0] * LV / out_dt_sec;
    FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_TVEG, /*lis_nest, lis_tile, lis_moc_tveg,*/ 
                                  &level1, 
                                  &value,
                                  units,&lenu,direction,&lend);

    /** soil tempture K **/
    /*** vic unit: C => LIS unit: K ***/
    units = "K"; lenu = strlen(units);
    direction = "-"; lend = strlen(direction);
    for( index=0; index<vic411_options.Nlayer; ++index ) 
    { 
       vlevel = index+1; 
       value = out_data[OUT_SOIL_TEMP].data[index] +  273.15;  
       FTN(vic411_diagnoseoutputvar)(LIS_NEST, LIS_TILE, LIS_MOC_SOILTEMP, 
                                     &vlevel, 
                                     &value, 
                                     units,&lenu,direction,&lend);
    } 
}
