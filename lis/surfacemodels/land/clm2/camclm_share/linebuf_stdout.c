//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
/* 
*** FORTRAN-callable routine to line buffer stdout.  Should only be used in
*** debug mode due to system overhead, and the fact that in theory stdout
*** should not have already been written to when this routine is called.
***
*** Original version: Jim Rosinski
***
*/

#include "LIS_misc.h"
#include <cfort.h>

#include <stdio.h>

#if ( defined FORTRANCAPS )

#define linebuf_stdout LINEBUF_STDOUT

#elif ( defined FORTRANUNDERSCORE )

#define linebuf_stdout linebuf_stdout_

#elif ( defined FORTRANDOUBLEUNDERSCORE )

#define linebuf_stdout linebuf_stdout__

#endif

void linebuf_stdout ()
{
  setvbuf (stdout, NULL, _IOLBF, 0);
  return;
}
