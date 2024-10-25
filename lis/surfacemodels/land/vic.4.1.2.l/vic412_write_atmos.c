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
#include "ftn.h"

extern void write_atmosdata(atmos_data_struct *, int);

//BOP
//
// !ROUTINE: vic412_write_atmos
// \label{vic412_write_atmos}
// 
// !REVISION HISTORY: 
//  Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(vic412_write_atmos)(int *LIS_NTILES)
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
   extern atmos_data_struct     **lis_atmos;

   atmos_data_struct      *atmos;

   int i;

   for ( i = 0; i < *LIS_NTILES; ++i )
   {
      atmos = &lis_atmos[i][0];

      write_atmosdata(atmos, 1);
   }
}
