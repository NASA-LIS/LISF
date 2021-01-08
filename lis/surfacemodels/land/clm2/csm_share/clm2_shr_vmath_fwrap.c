//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.3
//
// Copyright (c) 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
/*
** =============================================================================
** Fortran wrappers for shr_vmath calls for systems that only
** provide a C interface to the system vector math routines.
** =============================================================================
*/

#if (defined IRIX64)

void shr_vmath_fwrap_vsqrt_(double *X, double *Y, int *n)
{
   vsqrt(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vexp_(double *X, double *Y, int *n)
{
   vexp(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vlog_(double *X, double *Y, int *n)
{
   vlog(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vsin_(double *X, double *Y, int *n)
{
   vsin(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vcos_(double *X, double *Y, int *n)
{
   vcos(X, Y, *n, 1, 1);
}

#endif
