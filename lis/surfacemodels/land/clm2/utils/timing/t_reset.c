//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <string.h>  /* memset */
#include <stdio.h>

#include "gpt.h"

/*
** t_reset: reset all known timers to 0
**
** Return value: 0 (success) or -1 (failure)
*/

int t_reset ()
{
  int n;             /* index over threads */
  struct node *ptr;  /* linked list index */

#if ( ! defined DISABLE_TIMERS )
  if ( ! t_initialized)
    return t_error ("t_reset: t_initialize has not been called\n");

  /*
  ** Only allow the master thread to reset timers
  */

  if (get_thread_num () != 0)
    return 0;

  for (n = 0; n < numthreads; n++) {
    for (ptr = timers[n]; ptr != NULL; ptr = ptr->next) {
      memset (timers[n], 0, sizeof (struct node));
      printf ("Reset accumulators for timer %s to zero\n", ptr->name);
    }
  }
#endif
  return 0;
}
