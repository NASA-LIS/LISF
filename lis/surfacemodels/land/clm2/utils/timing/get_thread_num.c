//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "gpt.h"
static int lock_mutex (void);
static int unlock_mutex (void);

/*
** get_thread_num: Obtain logical thread number of calling thread.  If new
** thread, adjust global variables.
*/

int get_thread_num ()
{
  int mythread = 0 ;  /* return value: default zero for non-threaded case */

#ifdef THREADED_PTHREADS
  pthread_t mythreadid;  /* thread id from pthreads library */
#endif

/*
** If timers disabled, construct ifdef to avoid function call to things that
** might not be available.
*/

#if ( ! defined DISABLE_TIMERS )

#if ( defined THREADED_OMP )

  if ((mythread = omp_get_thread_num ()) >= numthreads)
    return t_error ("get_thread_num: returned id %d exceed numthreads %d\n",
		    mythread, numthreads);

#elif ( defined THREADED_PTHREADS )

  mythreadid = pthread_self ();

  if (lock_mutex () < 0)
    return t_error ("get_thread_num: mutex lock failure\n");

  /*
  ** Loop over known physical thread id's.  When my id is found, map it 
  ** to logical thread id for indexing.  If not found return a negative 
  ** number.
  ** A critical region is necessary because acess to
  ** the array threadid must be by only one thread at a time.
  */

  {
    int n;              /* loop index over number of threads */
    for (n = 0; n < numthreads; n++)
      if (pthread_equal (mythreadid, threadid[n]))
	break;

  /*
  ** If our thread id is not in the known list, add to it after checking that
  ** we do not have too many threads.
  */

    if (n == numthreads) {
      if (numthreads >= MAX_THREADS)
	return t_error ("get_thread_num: numthreads=%d is too big Recompile "
			"with larger value of MAX_THREADS\n", numthreads);

      threadid[n] = mythreadid;
      numthreads++;
    }

    if (unlock_mutex () < 0)
      return t_error ("get_thread_num: mutex unlock failure\n");

    mythread = n;
  }

#endif
#endif

  return mythread;
}

/*
** lock_mutex: lock a mutex for entry into a critical region
*/

#if ( ! defined DISABLE_TIMERS )

static int lock_mutex ()
{
#if ( defined THREADED_OMP )
  omp_set_lock (&lock);
#elif ( defined THREADED_PTHREADS )
  if (pthread_mutex_lock (&t_mutex) != 0)
    return t_error ("pthread_mutex_lock failure\n");
#endif
  return 0;
}

/*
** unlock_mutex: unlock a mutex for exit from a critical region
*/

static int unlock_mutex ()
{
#if ( defined THREADED_OMP )
  omp_unset_lock (&lock);
#elif ( defined THREADED_PTHREADS )
  if (pthread_mutex_unlock (&t_mutex) != 0)
    return t_error ("pthread_mutex_unlock failure\n");
#endif
  return 0;
}

#endif
