//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <sys/time.h>      /* gettimeofday */
#include <unistd.h>        /* gettimeofday */
#include <stdlib.h>        /* malloc */
#include <string.h>        /* memset, strcmp (via STRMATCH) */

#include "gpt.h"

/*
** t_start.c: start a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

int t_start (char *name)        /* timer name */
{
  struct timeval tp1, tp2;      /* argument to gettimeofday */
  struct node *ptr;             /* linked list pointer */

  int numchars;                 /* number of characters in timer */
  int mythread;                 /* thread index (of this thread) */
  int indent_level = 0;         /* required indentation level for this timer */
  int ret;                      /* return code */

  PCL_CNT_TYPE i_pcl_result1[PCL_COUNTER_MAX];     /* init. output fm PCLread */
  PCL_CNT_TYPE i_pcl_result2[PCL_COUNTER_MAX];     /* final output fm PCLread */
  PCL_FP_CNT_TYPE fp_pcl_result[PCL_COUNTER_MAX];  /* required by PCLread */

  /*
  ** 1st system timer call is solely for overhead timing
  */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if (wallenabled)
    gettimeofday (&tp1, NULL);

  if ( ! t_initialized)
    return t_error ("t_start: t_initialize has not been called\n");

  if ((mythread = get_thread_num ()) < 0)
    return t_error ("t_start\n");

  if (ncounter > 0) {
    ret = PCLread (descr[mythread], i_pcl_result1, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_start: error from PCLread: %s\n", t_pclstr (ret));
  }

  /*
  ** Look for the requested timer in the current list.  For those which don't
  ** match but are currently active, increase the indentation level by 1
  */

  for (ptr = timers[mythread]; ptr != NULL && ! STRMATCH (name, ptr->name); 
       ptr = ptr->next) {

    if (ptr->onflg) 
      indent_level++;
  }

  if (indent_level > max_indent_level[mythread])
    max_indent_level[mythread] = indent_level;
    
  /* 
  ** If a new thing is being timed, add a new link and initialize 
  */

  if (ptr == NULL) {

    if ((ptr = (struct node *) malloc (sizeof (struct node))) == NULL)
      return (t_error ("t_start: malloc failed\n"));

    memset (ptr, 0, sizeof (struct node));
    ptr->indent_level = indent_level;
    ptr->next = NULL;

    if (timers[mythread] == NULL)
      timers[mythread] = ptr;
    else
      last[mythread]->next = ptr;

    last[mythread] = ptr;

    /* 
    ** Truncate input name if longer than MAX_CHARS characters 
    */

    numchars = MIN (strlen (name), MAX_CHARS);
    strncpy (ptr->name, name, numchars);
    ptr->name[numchars] = '\0';

  } else {

    /*
    ** If computed indentation level is different than before or was
    ** already ambiguous, reset to ambiguous flag value.  This will likely
    ** happen any time the thing being timed is called from more than 1
    ** branch in the call tree.
    */

    if (ptr->indent_level != indent_level) 
      ptr->indent_level = AMBIGUOUS;

    if (ptr->onflg)
/*      return t_error ("t_start thread %d: timer %s was already on: "
		      "not restarting.\n", mythread, ptr->name);
*/
        return t_error (" ");        
  }

  ptr->onflg = true;

  if (usrsysenabled)
    if (get_cpustamp (&ptr->last_utime, &ptr->last_stime) < 0)
      return t_error ("t_start: get_cpustamp error");

  /*
  ** The 2nd system timer call is used both for overhead estimation and
  ** the input timer
  */

  if (wallenabled) {

    gettimeofday (&tp2, NULL);
    ptr->last_wtime_sec  = tp2.tv_sec;
    ptr->last_wtime_usec = tp2.tv_usec;
    overhead[mythread] +=       (tp2.tv_sec  - tp1.tv_sec) + 
                          1.e-6*(tp2.tv_usec - tp1.tv_usec);
  }

  if (ncounter > 0) {
    int n;
    int index;

    ret = PCLread (descr[mythread], i_pcl_result2, fp_pcl_result, ncounter); 
    if (ret != PCL_SUCCESS)
      return t_error ("t_start: error from PCLread: %s\n", t_pclstr (ret));

    for (n = 0; n < ncounter; n++) {
      ptr->last_pcl_result[n] = i_pcl_result2[n];
    }

    if (pcl_cyclesenabled) {
      index = pcl_cyclesindex;
      overhead_pcl[mythread] += i_pcl_result2[index] - i_pcl_result1[index];
    }
  } 

  return (0);
}

/*
** This stub should never actually be called
*/

#ifndef HAVE_PCL
int PCLread (PCL_DESCR_TYPE descr, PCL_CNT_TYPE *i, PCL_CNT_TYPE *j, int k)
{
  return t_error ("PCLread called when library not there\n");
}
#endif

