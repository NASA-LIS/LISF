//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"
#include "ftn_drv.h"

/* EMK Only compile of LIBGEOTIFF support is enabled */
#ifdef USE_LIBGEOTIFF

#include <stdlib.h> 
#include <stdio.h>

#include "geotiffio.h"
#include "xtiffio.h"
#include "ztif.h"

/* A general-use function for just pulling the ZTIF data out at native 
   resolution.  Error status is returned in the ierr argument--values
   other than 0 indicate an error occurred.

   Eric Kemp, NASA GSFC/SSAI */

void
FTN(ztif_frac_slice)(int *map_buffer_slice,
                 int *age_buffer_slice,
                 const int *nc_native,
                 const int *nr_native,
                 char *map_path,
                 char *age_path,
                 const int *offset,
                 const int *jslice,
                 int* ierr)
{

    /* Local variables */
    int status;
    int nc;
    int nr;
    uint8 map_pixel;
    uint8 age_pixel;
    int j;
    int i;

    /* Static variables for keeping data in memory--and files open--between 
       invocations */
    static struct ZTIF map;
    static struct ZTIF age;

    /* For convenience */
    nc = (*nc_native);
    nr = (*nr_native);

    /* Only open the files and check the dimensions for the very first
       row (i.e., the first invocation of this function).  Since the
       map and age structures are static, the files will remain open
       and certain attributes will remain in memory between function
       invocations. */
    if (*jslice == 1) {

	/* Open the snow map file and check the dimensions */
        status = ZTIFOpen(&map, map_path, "r");
        if (status != E_SUCCESS) {
            *ierr = 1;
            return;
        }
        if (nr != map.length) {
            ZTIFClose(&map);
            *ierr = 1;
            return;
        }
	if (nc != map.width) {
            ZTIFClose(&map);
            *ierr = 1;
            return;
        }

	/* Open the age file and check the dimensions  */
        status = ZTIFOpen(&age, age_path, "r");
        if (status != E_SUCCESS) {
            ZTIFClose(&map);
            *ierr = 2;
            return;
        }       
        if (nr != age.length) {
            ZTIFClose(&map);
            ZTIFClose(&age);
            *ierr = 2;
            return;
        }
        if (nc != age.width) {
            ZTIFClose(&map);
            ZTIFClose(&age);
            *ierr = 2;
            return;
        }
    }

    /* Copy a slice of snow map data from the file to the output argument,
       converting to int.  We cannot pull the whole dataset due to massive
       size */
    j = (*jslice) - 1; /* Switch from Fortran to C convention */
    status = ZTIFReadline(&map, j);
    if (status != E_SUCCESS) {
        ZTIFClose(&map);
        ZTIFClose(&age);
        *ierr = 1;
        return;
    }
    for (i=0; i < map.width; i++) {
        map_pixel = ((uint8*) map.rbuf)[i];
        map_buffer_slice[i] = (int) map_pixel;
    }   

    /* Copy the age data from the file to the output argument, converting
       to int */
    j = (*jslice) - 1;
    status = ZTIFReadline(&age, j);
    if (status != E_SUCCESS) {
        ZTIFClose(&map);
        ZTIFClose(&age);
        *ierr = 2;
        return;
    }
    for (i=0; i < age.width; i++) {
        age_pixel = ((uint8*) age.rbuf)[i];
        age_buffer_slice[i] = (int) age_pixel + (*offset); 
    }   

    /* If this is the final row, this is also the last function invocation.
       Therefore, we close the files before exiting. */
    if (*jslice == nr) {
        ZTIFClose(&map);
        ZTIFClose(&age);
    }

    *ierr = 0;
    return;
}

#else

/* Dummy version of ztif_frac_slice */
void
ztif_frac_slice_(int *map_buffer_slice,
                 int *age_buffer_slice,
                 const int *nc_native,
                 const int *nr_native,
                 char *map_path,
                 char *age_path,
                 const int *offset,
                 const int *jslice,
                 int* ierr)
{
    *ierr = 3;
    return;
}
#endif
