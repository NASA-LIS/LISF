//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <sys/time.h>       /* gettimeofday */
#include <unistd.h>         /* gettimeofday */
#include <string.h>         /* strcmp (via STRMATCH) */

#include "gpt.h"

/*
** t_stop: stop a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

int t_stop (char *name)
{
  long delta_wtime_sec;     /* wallclock change fm t_start() to t_stop() */    
  long delta_wtime_usec;    /* wallclock change fm t_start() to t_stop() */
  float delta_wtime;        /* floating point wallclock change */
  struct timeval tp1, tp2;  /* argument to gettimeofday() */
  struct node *ptr;         /* linked list pointer */

  int mythread;             /* thread number for this process */
  int ret;                  /* return code */

  long usr;
  long sys;

  PCL_CNT_TYPE i_pcl_result1[PCL_COUNTER_MAX];     /* init. output fm PCLread */
  PCL_CNT_TYPE i_pcl_result2[PCL_COUNTER_MAX];     /* final output fm PCLread */
  PCL_FP_CNT_TYPE fp_pcl_result[PCL_COUNTER_MAX];  /* required by PCLread */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if ( ! t_initialized)
    return t_error ("t_stop: t_initialize has not been called\n");

  /*
  ** The 1st system timer call is used both for overhead estimation and
  ** the input timer
  */

  if (wallenabled)
    gettimeofday (&tp1, NULL);

  if (usrsysenabled && get_cpustamp (&usr, &sys) < 0)
    return t_error (NULL);

  if ((mythread = get_thread_num ()) < 0)
    return t_error ("t_stop\n");

  if (ncounter > 0) {
    ret = PCLread (descr[mythread], i_pcl_result1, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_stop: error from PCLread: %s\n", t_pclstr (ret));
  }
  
  for (ptr = timers[mythread]; ptr != NULL && ! STRMATCH (name, ptr->name); 
       ptr = ptr->next);

  if (ptr == NULL) 
    return t_error ("t_stop: timer for %s had not been started.\n", name);

  if ( ! ptr->onflg )
    return t_error ("t_stop: timer %s was already off.\n",ptr->name);

  ptr->onflg = false;
  ptr->count++;

  /*
  ** 1st timer stoppage: set max and min to computed values.  Otherwise apply
  ** max or min function
  */

  if (wallenabled) {

    delta_wtime_sec  = tp1.tv_sec  - ptr->last_wtime_sec;
    delta_wtime_usec = tp1.tv_usec - ptr->last_wtime_usec;
    delta_wtime = delta_wtime_sec + 1.e-6*delta_wtime_usec;

    if (ptr->count == 1) {
      ptr->max_wtime = delta_wtime;
      ptr->min_wtime = delta_wtime;
      
    } else {
      
      ptr->max_wtime = MAX (ptr->max_wtime, delta_wtime);
      ptr->min_wtime = MIN (ptr->min_wtime, delta_wtime);
    }

    ptr->accum_wtime_sec  += delta_wtime_sec;
    ptr->accum_wtime_usec += delta_wtime_usec;

    /*
    ** Adjust accumulated wallclock values to guard against overflow in the
    ** microsecond accumulator.
    */

    if (ptr->accum_wtime_usec > 1000000) {
      ptr->accum_wtime_usec -= 1000000;
      ptr->accum_wtime_sec  += 1;
      
    } else if (ptr->accum_wtime_usec < -1000000) {
      
      ptr->accum_wtime_usec += 1000000;
      ptr->accum_wtime_sec  -= 1;
    }

    ptr->last_wtime_sec  = tp1.tv_sec;
    ptr->last_wtime_usec = tp1.tv_usec;

    /*
    ** 2nd system timer call is solely for overhead timing
    */

    gettimeofday (&tp2, NULL);
    overhead[mythread] +=       (tp2.tv_sec  - tp1.tv_sec) + 
                          1.e-6*(tp2.tv_usec - tp1.tv_usec);
  }

  if (usrsysenabled) {
    ptr->accum_utime += usr - ptr->last_utime;
    ptr->accum_stime += sys - ptr->last_stime;
    ptr->last_utime   = usr;
    ptr->last_stime   = sys;
  }

  if (ncounter > 0) {
    int n;
    PCL_CNT_TYPE delta;
    int index;

    for (n = 0; n < ncounter; n++) {
      delta = i_pcl_result1[n] - ptr->last_pcl_result[n];

      /*
      ** Accumulate results only for positive delta
      */

      if (delta < 0) 
	printf ("t_stop: negative delta => probable counter overflow. "
		"Skipping accumulation this round\n"
		"%ld - %ld = %ld\n", (long) i_pcl_result1[n], 
		                     (long) ptr->last_pcl_result[n],
		                     (long) delta);
      else
	ptr->accum_pcl_result[n] += delta;

      ptr->last_pcl_result[n] = i_pcl_result1[n];
    }

    /*
    ** Overhead estimate.  Currently no check for negative delta
    */

    ret = PCLread (descr[mythread], i_pcl_result2, fp_pcl_result, ncounter);
    if (ret != PCL_SUCCESS)
      return t_error ("t_stop: error from PCLread: %s\n", t_pclstr (ret));

    if (pcl_cyclesenabled) {
      index = pcl_cyclesindex;
      overhead_pcl[mythread] += i_pcl_result2[index] - i_pcl_result1[index];
    }
  }

  return 0;
}
