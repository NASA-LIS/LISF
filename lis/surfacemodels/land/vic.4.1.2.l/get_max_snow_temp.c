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
#include "vicNl.h"
#include "ftn.h"
//BOP
//
// !ROUTINE: get_max_snow_temp
// \label{get_max_snow_temp}
// 
// !REVISION HISTORY: 
//  21 Dec 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
double FTN(get_max_snow_temp)(void)
//
// !DESCRIPTION: 
//  This routine returns the MAX\_SNOW\_TEMP value from the global\_param
//  data structure.
//EOP
{
   extern global_param_struct global_param;

   return global_param.MAX_SNOW_TEMP;
}
