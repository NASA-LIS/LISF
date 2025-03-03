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
** 
** Fortran wrappers for timing library routines
*/

#if ( defined CRAY ) || ( defined T3D )
#include <fortran.h>
#endif

#include "cfort.h"
#include "gpt.h"

#if ( defined ABSOFT )
#undef FORTRANCAPS
#undef FORTRANUNDERSCORE
#define FORTRANDOUBLEUNDERSCORE
#endif

#if ( defined FORTRANCAPS )

#define t_initializef T_INITIALIZEF
#define t_prf T_PRF
#define t_resetf T_RESETF
#define t_stampf T_STAMPF
#define t_startf T_STARTF
#define t_stopf T_STOPF
#define t_setoptionf T_SETOPTIONF

#elif ( defined FORTRANUNDERSCORE )

#define t_initializef t_initializef_
#define t_prf t_prf_
#define t_resetf t_resetf_
#define t_stampf t_stampf_
#define t_startf t_startf_
#define t_stopf t_stopf_
#define t_setoptionf t_setoptionf_

#elif ( defined FORTRANDOUBLEUNDERSCORE )

#define t_initializef t_initializef__
#define t_prf t_prf__
#define t_resetf t_resetf__
#define t_stampf t_stampf__
#define t_startf t_startf__
#define t_stopf t_stopf__
#define t_setoptionf t_setoptionf__

#endif

#if ( defined CRAY ) || ( defined T3D )

int t_startf (_fcd);
int t_stopf (_fcd);

#else

int t_startf (char *, int);
int t_stopf (char *, int);

#endif

int t_initializef ()
{
  return t_initialize ();
}

int t_prf (int *procid)
{
  return t_pr (*procid);
}

void t_resetf ()
{
  t_reset();
  return;
}

int t_setoptionf (int *option, int *val)
{
  return t_setoption ( (OptionName) *option, (Boolean) *val);
}

int t_stampf (double *wall, double *usr, double *sys)
{
  return t_stamp (wall, usr, sys);
}

#if ( defined CRAY ) || ( defined T3D )

int t_startf (_fcd name)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (_fcdlen (name), MAX_CHARS);
  strncpy (cname, _fcdtocp (name), numchars);
  cname[numchars] = '\0';
  return t_start (cname);
}

int t_stopf (_fcd name)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (_fcdlen (name), MAX_CHARS);
  strncpy (cname, _fcdtocp (name), numchars);
  cname[numchars] = '\0';
  return t_stop (cname);
}

#else

int t_startf (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return t_start (cname);
}

int t_stopf (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return t_stop (cname);
}

#endif
