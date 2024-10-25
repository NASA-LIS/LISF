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
#include <vicNl.h>
#include "vic412_atmos_forcing.h"
#include "ftn.h"
//BOP
//
// !ROUTINE: vic412_add_atmosdata
// \label{vic412_add_atmosdata}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic412_add_atmosdata)(int *LIS_tile, int *step, int *VAR,
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
   extern atmos_data_struct     **lis_atmos;

   atmos_data_struct      *atmos;

   int tile;


   tile = *LIS_tile - 1;

   atmos = &lis_atmos[tile][0];


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
// !ROUTINE: vic412_nr_atmosdata
// \label{vic412_nr_atmosdata}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic412_nr_atmosdata)(int *NTILES)
// !DESCRIPTION: 
//  For water-balance mode (NR $>$ 0), this routine iterates over every tile 
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
   extern atmos_data_struct **lis_atmos;
   atmos_data_struct         *atmos;

   int i, j;

   if ( NR > 0 )
   {
      for ( i = 0; i < *NTILES; ++i )
      {
         atmos = &lis_atmos[i][0];

         atmos->snowflag[NR]  = FALSE;
         atmos->prec[NR]      = 0.;
         atmos->air_temp[NR]  = 0.;
         atmos->wind[NR]      = 0.;
         atmos->vpd[NR]       = 0.;
         atmos->vp[NR]        = 0.;
         atmos->pressure[NR]  = 0.;
         atmos->density[NR]   = 0.;
         atmos->shortwave[NR] = 0.;
         atmos->longwave[NR]  = 0.;

         for ( j = 0; j < NF; ++j )
         {
            if ( atmos->snowflag[j] == TRUE )
            {
               atmos->snowflag[NR] = TRUE;
            }
            atmos->prec[NR]      += atmos->prec[j];
            atmos->air_temp[NR]  += atmos->air_temp[j];
            atmos->wind[NR]      += atmos->wind[j];
//            atmos->vpd[NR]       += atmos->vpd[j];
            atmos->vp[NR]        += atmos->vp[j];
            atmos->pressure[NR]  += atmos->pressure[j];
            atmos->density[NR]   += atmos->density[j];
            atmos->shortwave[NR] += atmos->shortwave[j];
            atmos->longwave[NR]  += atmos->longwave[j];
         }
         atmos->air_temp[NR]  /= NF;
         atmos->wind[NR]      /= NF;
         atmos->vp[NR]        /= NF;
         atmos->pressure[NR]  /= NF;
         atmos->density[NR]   /= NF;
         atmos->shortwave[NR] /= NF;
         atmos->longwave[NR]  /= NF;
         atmos->vpd[NR] = svp(atmos->air_temp[NR]) - atmos->vp[NR];
      }
   }
}
