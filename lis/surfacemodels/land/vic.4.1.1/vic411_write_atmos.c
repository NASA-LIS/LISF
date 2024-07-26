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
#include <vic411_vicNl.h>
#include "ftn.h"

//BOP
//
// !ROUTINE: vic411_write_atmos
// \label{vic411_write_atmos}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic411_write_atmos)(int *LIS_NTILES)
// !DESCRIPTION: 
//  This routine iterates over every grid-cell to write out the values
//  of the forcing data contain within the atmos data structure.
//
//  The actual writing is performed by VIC's write\_atmosdata routine.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_NTILES] number of sub-grid tiles, considered intent(in)
//  \end{description}
//
//EOP
{
   extern vic411_atmos_data_struct     **vic411_lis_atmos;

   vic411_atmos_data_struct      *atmos;

   int i;

   for ( i = 0; i < *LIS_NTILES; ++i )
   {
      atmos = &vic411_lis_atmos[i][0];

      vic411_write_atmosdata(atmos, 1);
   }
}
