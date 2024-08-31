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

void FTN(vic412_getrunoff)(int *LIS_TILE, float *runoff1, float *runoff2)

{
   extern global_param_struct   global_param;
   extern out_data_struct      **lis_out_data;

   float dt_sec; /* time step in seconds */
   int                      cellnum;

   out_data_struct         *out_data;

   cellnum = *LIS_TILE-1;

   out_data = lis_out_data[cellnum];

   dt_sec = global_param.dt * 3600.0;
    
   *runoff1 = (float)out_data[OUT_RUNOFF].data[0] / dt_sec;
   *runoff2 = (float)out_data[OUT_BASEFLOW].data[0] / dt_sec;
}
