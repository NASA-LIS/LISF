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
//BOP
//
// !ROUTINE: lis_scan_soilparam_flat_ascii
// \label{lis_scan_soilparam_flat_ascii}
// 
// !REVISION HISTORY: 
//  Feb 2012  James Geiger; Initial implementation
// 
// !INTERFACE:
void lis_scan_soilparam_flat_ascii(FILE  *soilparam, 
                                   float *lat,
                                   float *lng)
//
// !DESCRIPTION: 
//  This routine scans the flat ASCII soil parameters for the given grid cell,
//  and then resets the position within the FILE stream for subsequent reading.
//
//  The arguments are:
//  \begin{description}
//  \item[soilparam] filename of the soil parameter file
//  \item[lat] latitude of the desired grid-cell
//  \item[lng] longitude of the desired grid-cell
//  \end{description}
//
//
//EOP
{

   long pos;

   int status;

   int flag, gridcell;
   float temp_lat, temp_lng;
   float epsilon;

   char line[MAXSTRING];

   epsilon = 10e-6;

   rewind(soilparam);

   pos = ftell(soilparam);
   while ( ( status = fscanf(soilparam,"%d %d %f %f",
                  &flag, &gridcell, &temp_lat, &temp_lng) ) != EOF )
   {
      if ( ( flag == 1 ) && 
           ( temp_lat - epsilon <= *lat && *lat <= temp_lat + epsilon ) &&
           ( temp_lng - epsilon <= *lng && *lng <= temp_lng + epsilon ) 
         )
      {
         fseek(soilparam, pos, SEEK_SET);
         break;
      }
      else
      {
         fgets( line, MAXSTRING, soilparam );
         pos = ftell(soilparam);
      }
   } 

   if ( status == EOF )
   {
      printf("ERR: Error finding (%f, %f) in the flat ascii soil paramter file.\n", *lat, *lng);
   }
} 
