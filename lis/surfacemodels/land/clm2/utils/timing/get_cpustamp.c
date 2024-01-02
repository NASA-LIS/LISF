//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include <sys/times.h>  /* times */

#include "gpt.h"

/*
** get_cpustamp: Invoke the proper system timer and return stats.
**
** Output arguments:
**   usr: user time (usec if USE_GETRUSAGE is defined, ticks otherwise)
**   sys: system time (usec if USE_GETRUSAGE is defined, ticks otherwise)
**
** Return value: 0 (success)
*/

int get_cpustamp (long *usr, long *sys)
{
  struct tms buf;

  /*
  ** Throw away the wallclock time from times: use gettimeofday instead
  */
  
  (void) times (&buf);
  *usr = buf.tms_utime;
  *sys = buf.tms_stime;

  return 0;
}
