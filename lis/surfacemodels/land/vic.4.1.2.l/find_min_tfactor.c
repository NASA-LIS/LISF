//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.3
//
// Copyright (c) 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include "vicNl.h"
#include "ftn.h"
//BOP
//
// !ROUTINE: find_min_tfactor
// \label{find_min_tfactor}
// 
// !REVISION HISTORY: 
//  21 Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
double FTN(find_min_tfactor)(int *LIS_NEST, int *LIS_TILE)
//
// !DESCRIPTION: 
//  This routine returns the minimum Tfactor value from the soil\_con.Tfactor
//  array.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \item[LIS\_TILE] tile index, considered intent(in)
//  \end{description}
//
//EOP
{
   extern soil_con_struct *lis_soil_con;
   extern option_struct       options;

   double min_Tfactor;

   int nest, tile;
   int band;


   nest = *LIS_NEST - 1;
   tile = *LIS_TILE - 1;

   min_Tfactor = lis_soil_con[tile].Tfactor[0];
   for (band = 1; band < options.SNOW_BAND; band++)
   {
      if (lis_soil_con[tile].Tfactor[band] < min_Tfactor)
         min_Tfactor = lis_soil_con[tile].Tfactor[band];
   }

   return min_Tfactor;
}
