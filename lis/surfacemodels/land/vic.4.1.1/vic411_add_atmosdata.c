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
#include "vic411_vicNl.h"
#include "vic411_atmos_forcing.h"
#include "ftn.h"
//BOP
//
// !ROUTINE: vic411_add_atmosdata
// \label{vic411_add_atmosdata}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic411_add_atmosdata)(int *LIS_tile, int *step, int *VAR,
                               float *forcing)
// !DESCRIPTION: 
//  This routine packs the forcing data extracted from LIS\_FORC\_State
//  into VIC's atmos data structure.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_TILE] tile index, considered intent(in)
//  \item[step] step index, considered intent(in)
//     VIC's forcing data are arrays of sized either by the number of 
//     snow-steps when running in water-balance mode or by 1 when running in
//     energy-balance mode.  This index specifies which snow-step (or 1 in
//     the case of energy-balance mode) is being added to the atmos data 
//     structure.
//  \item[VAR] variable index, considered intent(in)
//  \item[forcing] forcing value, considered intent(in)
//  \end{description}
//
//EOP
{
   extern vic411_atmos_data_struct     **vic411_lis_atmos;

   vic411_atmos_data_struct      *atmos;

   int tile;


   tile = *LIS_tile - 1;

   atmos = &vic411_lis_atmos[tile][0];


   switch ( *VAR )
   {
      case ATMOS_SNOWFLAG:
         atmos->snowflag[*step-1]  = *forcing;
         break;
      case ATMOS_PREC:
         atmos->prec[*step-1]      = *forcing;
         break;
      case ATMOS_AIR_TEMP:
         atmos->air_temp[*step-1]  = *forcing;
         break;
      case ATMOS_WIND:
         atmos->wind[*step-1]      = *forcing;
         break;
      case ATMOS_VPD:
         atmos->vpd[*step-1]       = *forcing;
         break;
      case ATMOS_VP:
         atmos->vp[*step-1]        = *forcing;
         break;
      case ATMOS_PRESSURE:
         atmos->pressure[*step-1]  = *forcing;
         break;
      case ATMOS_DENSITY:
         atmos->density[*step-1]   = *forcing;
         break;
      case ATMOS_SHORTWAVE:
         atmos->shortwave[*step-1] = *forcing;
         break;
      case ATMOS_LONGWAVE:
         atmos->longwave[*step-1]  = *forcing;
         break;
   }
}

//BOP
//
// !ROUTINE: vic411_nr_atmosdata
// \label{vic411_nr_atmosdata}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic411_nr_atmosdata)(int *NTILES)
// !DESCRIPTION: 
//  For water-balance mode (vic411\_NR $>$ 0), this routine iterates over every tile 
//  and computes the model step average or sum value for the forcing fields.
//
//  This routine computes the average or sum based off the most common
//  computations found in VIC's initialize\_atmos routine.  VIC's forcing
//  data processing is quite complex, so the computations performed here
//  may not exactly match the ones performed by VIC's initialize\_atmos
//  routine.  Thus the forcing for some VIC test-cases may not reproduced
//  exactly by this routine.
//
//  The arguments are:
//  \begin{description}
//  \item[NTILES] number of tiles, considered intent(in)
//  \end{description}
//
//EOP
{
   extern vic411_atmos_data_struct **vic411_lis_atmos;
   vic411_atmos_data_struct         *atmos;

   int i, j;

   if ( vic411_NR > 0 )
   {
      for ( i = 0; i < *NTILES; ++i )
      {
         atmos = &vic411_lis_atmos[i][0];

         atmos->snowflag[vic411_NR]  = FALSE;
         atmos->prec[vic411_NR]      = 0.;
         atmos->air_temp[vic411_NR]  = 0.;
         atmos->wind[vic411_NR]      = 0.;
         atmos->vpd[vic411_NR]       = 0.;
         atmos->vp[vic411_NR]        = 0.;
         atmos->pressure[vic411_NR]  = 0.;
         atmos->density[vic411_NR]   = 0.;
         atmos->shortwave[vic411_NR] = 0.;
         atmos->longwave[vic411_NR]  = 0.;

         for ( j = 0; j < vic411_NF; ++j )
         {
            if ( atmos->snowflag[j] == TRUE )
            {
               atmos->snowflag[vic411_NR] = TRUE;
            }
            atmos->prec[vic411_NR]      += atmos->prec[j];
            atmos->air_temp[vic411_NR]  += atmos->air_temp[j];
            atmos->wind[vic411_NR]      += atmos->wind[j];
//            atmos->vpd[vic411_NR]       += atmos->vpd[j];
            atmos->vp[vic411_NR]        += atmos->vp[j];
            atmos->pressure[vic411_NR]  += atmos->pressure[j];
            atmos->density[vic411_NR]   += atmos->density[j];
            atmos->shortwave[vic411_NR] += atmos->shortwave[j];
            atmos->longwave[vic411_NR]  += atmos->longwave[j];
         }
         atmos->air_temp[vic411_NR]  /= vic411_NF;
         atmos->wind[vic411_NR]      /= vic411_NF;
         atmos->vp[vic411_NR]        /= vic411_NF;
         atmos->pressure[vic411_NR]  /= vic411_NF;
         atmos->density[vic411_NR]   /= vic411_NF;
         atmos->shortwave[vic411_NR] /= vic411_NF;
         atmos->longwave[vic411_NR]  /= vic411_NF;
         atmos->vpd[vic411_NR] = vic411_svp(atmos->air_temp[vic411_NR]) - atmos->vp[vic411_NR];
      }
   }
}
