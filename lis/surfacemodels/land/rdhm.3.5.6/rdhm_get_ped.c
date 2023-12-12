//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.4
//
// Copyright (c) 2022 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
 /* wrapper for calling get_ped 
    Shugong Wang  10/25/2013  */

extern void get_ped(int, int, int, int, int, float *, float *, float *);

 void rdhm_get_ped_(int *year,            /* */
                    int *month,           /* */
                    int *day,             /* */
                    int *dt_min,          /* time step in minutes   */
                    float *pe_mon,        /* monthly PET, array[12] */
                    float *pe_adj,        /* PET adjustment         */
                    float *pe_day)        /* output, day PET        */
 {
    int npix = 1;
    get_ped(*year,        /* input */
            *month,       /* input */
            *day,         /* input */
            npix,         /* input */
            *dt_min,      /* input */
            pe_mon,       /* input */
            pe_adj,       /* input */
            pe_day);      /* output */
 }
