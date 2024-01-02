//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdlib.h>  /* malloc */
#include <unistd.h>  /* sysconf */
#include "gpt.h"

/* 
** Array (1 per thread) of linked lists of timers, and last timer in each list
*/

struct node **timers = NULL;
struct node **last = NULL;

long ticks_per_sec;

/*
** Define lock arrays depending upon the type of threading done
*/

#if ( defined THREADED_OMP )

omp_lock_t lock;

#elif ( defined THREADED_PTHREADS )

pthread_mutex_t t_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_t *threadid;

#endif

float *overhead;                   /* wallclock estimate of timer overhead */
int *max_indent_level;             /* maximum indentation level */
int numthreads            = 1;     /* number of threads.  1 is for no threading */
Boolean t_initialized     = false; /* whether t_initialize has been called */
Boolean wallenabled       = false; /* wallclock timer stats enabled */
Boolean usrsysenabled     = false; /* usr & sys timer stats enabled */
Boolean pclenabled        = false; /* enable PCL library */     
Boolean pcl_cyclesenabled = false; /* enable PCL cycle count */
int pcl_cyclesindex       = -1;    /* index for PCL cycle count */

struct PossibleEvent possible_event[] = {
  {usrsys,               true,  "Usr Sys   "},
  {wall,                 true,  "Wallclock "},
#ifdef HAVE_PCL
  {pcl_start,            false, "          "},  /* bracket PCL entries */
  {pcl_l1dcache_miss,    false, "l1 D miss "},
  {pcl_l2cache_miss,     false, "L2 miss   "},
  {pcl_cycles,           false, "Cycles    "},
  {pcl_elapsed_cycles,   false, "E-Cycles  "},
  {pcl_fp_instr,         false, "FP instr  "},
  {pcl_loadstore_instr,  false, "L/S instr "},
  {pcl_instr,            false, "Instruct  "},
  {pcl_stall,            false, "Stall     "},
  {pcl_end,              false, "          "},  /* bracket PCL entries */
#endif
};

struct Event **event = NULL;
int nevent = 0;
int npossible = sizeof (possible_event) / sizeof (struct PossibleEvent);

/*
** Needed by PCL library: otherwise unused
*/

PCL_DESCR_TYPE *descr;
int counter_list[PCL_COUNTER_MAX];
int ncounter = 0;                  /* number of PCL counters */
PCL_CNT_TYPE *overhead_pcl;        /* overhead counter (cycles) */

/*
** t_initialize (): Initialization routine must be called from single-threaded
**   region before any other timing routines may be called.  The need for this
**   routine could be eliminated if not targetting timing library for threaded
**   capability. 
**
** return value: 0 (success) or -1 (failure)
*/

int t_initialize ()
{
  int n;             /* index */
  int nbytes;        /* number of bytes for malloc */
  int ret;           /* return code */

/*
** Determine number of ticks per second for conversion use by other t_pr(), t_stamp()
*/

  if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
    return t_error ("t_initialize: token _SC_CLK_TCK is not defined\n");

#if ( ! defined DISABLE_TIMERS )

  if (t_initialized)
    return t_error ("t_initialize has already been called\n");

#if ( defined THREADED_OMP )

  /*
  ** OMP: must call init_lock before using the lock (get_thread_num())
  */

  omp_init_lock (&lock);

  numthreads = omp_get_max_threads();

#elif ( defined THREADED_PTHREADS )

  numthreads = MAX_THREADS;

#endif

  /*
  ** Allocate space for global arrays
  */

  nbytes = numthreads * sizeof (struct node *);
  if ((timers = (struct node **) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  if ((last = (struct node **) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (float);
  if ((overhead = (float *) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (PCL_CNT_TYPE);
  if ((overhead_pcl = (PCL_CNT_TYPE *) malloc (nbytes)) == 0)
    return t_error ("malloc failure: %d items\n", numthreads);

  nbytes = numthreads * sizeof (int);
  if ((max_indent_level = (int *) malloc (nbytes)) == 0)
    return t_error ("malloc failure for %d items\n", numthreads);

  /*
  ** Initialize array values
  */

  for (n = 0; n < numthreads; n++) {
    timers[n] = 0;
    last[n] = 0;
    overhead[n] = 0.;
    overhead_pcl[n] = 0;
    max_indent_level[n] = 0;
  }

#ifdef THREADED_PTHREADS

  /*
  ** In the pthreads case, we must manage the threadid array which maps
  ** physical thread id's to logical id's
  */

  nbytes = numthreads * sizeof (pthread_t);
  if ((threadid = (pthread_t *) malloc (nbytes)) == 0)
    return t_error ("malloc failure for %d items\n", numthreads);

  /*
  ** Reset numthreads to 1 and define the threadid array now that initialization 
  ** is done.
  */

  threadid[0] = pthread_self ();
  numthreads = 1;

#endif

  if (get_thread_num () > 0) 
    return t_error ("t_initialize: should only be called by master thread\n");

  for (n = 0; n < npossible; n++) {
    if (possible_event[n].enabled) {
      if (possible_event[n].name == usrsys)
	usrsysenabled = true;

      if (possible_event[n].name == wall)
	wallenabled = true;

      if ((event = realloc (event, (nevent+1) * sizeof (struct Event *))) == NULL)
	return t_error ("realloc failure\n");

      if ((event[nevent] = malloc (sizeof (struct Event))) == NULL)
	return t_error ("realloc failure\n");

      event[nevent]->name = possible_event[n].name;
      strcpy (event[nevent]->string, possible_event[n].string);

#ifdef HAVE_PCL

      /*
      ** Set up PCL stuff based on what t_setoption has provided.
      */

      if (event[nevent]->name > pcl_start && event[nevent]->name < pcl_end) {
	pclenabled = true;
	event[nevent]->index = ncounter;

	switch (possible_event[n].name) {

	case pcl_l1dcache_miss:
	  counter_list[ncounter++] = PCL_L1DCACHE_MISS;
	  break;
	  
	case pcl_l2cache_miss: 
	  counter_list[ncounter++] = PCL_L2CACHE_MISS;
	  break;
	  
	case pcl_cycles: 
	  pcl_cyclesindex = ncounter;
	  pcl_cyclesenabled = true;
	  counter_list[ncounter++] = PCL_CYCLES;
	  break;

	case pcl_elapsed_cycles: 
	  counter_list[ncounter++] = PCL_ELAPSED_CYCLES;
	  break;

	case pcl_fp_instr: 
	  counter_list[ncounter++] = PCL_FP_INSTR;
	  break;

	case pcl_loadstore_instr: 
	  counter_list[ncounter++] = PCL_LOADSTORE_INSTR;
	  break;

	case pcl_instr: 
	  counter_list[ncounter++] = PCL_INSTR;
	  break;

	case pcl_stall: 
	  counter_list[ncounter++] = PCL_STALL;
	  break;
	
	default:
	  break;

	}
      }
#endif
      ++nevent;
    }
  }

#ifdef HAVE_PCL

  if (ncounter > 0) {
    int thread;         /* thread number */

    nbytes = numthreads * sizeof (PCL_DESCR_TYPE);
    if ((descr = (PCL_DESCR_TYPE *) malloc (nbytes)) == 0)
      return t_error ("malloc failure: %d items\n", numthreads);

    /*
    ** PCLinit must be called on a per-thread basis.  Therefore must make the call here
    ** rather than in t_initialize.  null timer list flags not initialized.
    ** Also, the critical section is necessary because PCLstart appears not to be
    ** thread-safe.
    */

#pragma omp parallel for
    
    for (thread = 0; thread < numthreads; thread++) {

      unsigned int flags;           /* mode flags needed by PCL */

#pragma omp critical

      {
	if ((ret = PCLinit (&descr[thread])) != PCL_SUCCESS)
	  return t_error ("unable to allocate PCL handle for thread %d. %s\n",
			  thread, t_pclstr (ret));

	/*
	** Always count user mode only
	*/
      
	flags = PCL_MODE_USER;

	if ((ret = PCLquery (descr[thread], counter_list, ncounter, flags)) != PCL_SUCCESS)
	  return t_error ("Bad return from PCLquery thread %d: %s\n", thread, t_pclstr (ret));

	if ((ret = PCLstart (descr[thread], counter_list, ncounter, flags)) != PCL_SUCCESS)
	  return t_error ("PCLstart failed thread=%d: %s\n", thread, t_pclstr (ret));
      }
    }
  }
#endif
  t_initialized = true;
#endif
  return 0;
}
