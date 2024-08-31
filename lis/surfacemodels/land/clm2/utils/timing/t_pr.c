//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <string.h>     /* memset */
#include <sys/time.h>   /* gettimeofday */
#include <sys/times.h>  /* times */
#include <unistd.h>     /* sysconf */

#include "gpt.h"

struct Stats {		   
  float usr;	   /* user CPU time */
  float sys;	   /* system CPU time */
  float usrsys;	   /* usr + sys */
  float elapse;	   /* elapsed time */
  float max_wtime; /* max elapsed wallclock time per call */
  float min_wtime; /* min elapsed wallclock time per call */

  long count;	   /* number of invocations of this timer */

  PCL_CNT_TYPE pcl_result[PCL_COUNTER_MAX];
};

static void fillstats (struct Stats *, struct node *);
static void print_header (FILE *, int);
static void print_stats_line (FILE *, struct Stats *);

/*
** t_pr: print stats for all known timers to a file
**
** Input arguments:
**   procid: Designed for MPI jobs to give a unique output file name, 
**     normally the MPI logical process id.
**
** Return value: 0 (success) or -1 (failure)
*/

static int mhz;            /* processor clock rate (from pcl library) */

int t_pr (int procid)
{
  char outfile[11];        /* name of output file: timing.xxx */
			   
  int indent;              /* index over number indentation levels */
  int thread;              /* thread index */
  int ilstart;             /* index over indentation level */
  int n;

  double gttdcost;         /* cost of a single gettimeofday call */
  double deltat;           /* time difference between 2 calls to gettimeofday */
  struct Stats stats;      /* per timer stats */
  struct Stats threadsum;  /* timer sum over threads */
			   
  struct node *ptr, *tptr; /* linked list pointers */
			   
  FILE *fp;                /* output file pointer */

  struct timeval tp1, tp2; /* input to gettimeofday() */
  struct tms buf;          /* input to times() */

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if ( ! t_initialized)
    return t_error ("t_pr: t_initialize has not been called\n");

  /*
  ** Only allow the master thread to print stats
  */

  if (get_thread_num () != 0)
    return 0;

  sprintf (outfile, "%s%d\0","timing.", procid);

  if ((fp = fopen (outfile, "w")) == NULL)
    fp = stderr;

  fprintf (fp,"Procid = %d\n", procid);

  /*
  ** Estimate wallclock timer overhead: 4 times the cost of a call to gettimeofday
  ** (since each start/stop pair involves 4 such calls).
  */

  gettimeofday (&tp1, NULL);
  gettimeofday (&tp2, NULL);
  gttdcost = 1.e6*(tp2.tv_sec  - tp1.tv_sec) + (tp2.tv_usec - tp1.tv_usec);

  fprintf (fp, "Wallclock timer cost est.: %8.3g usec per start/stop pair\n", 
	   gttdcost*4.);

  /*
  ** CPU cost estimate: 2 times the cost of a call to times().  Subtract the
  ** cost of a single gettimeofday call to improve the estimate.
  */

  if (usrsysenabled) {

    gettimeofday (&tp1, NULL);
    (void) times (&buf);
    gettimeofday (&tp2, NULL);

    deltat = 1.e6*(tp2.tv_sec  - tp1.tv_sec) + (tp2.tv_usec - tp1.tv_usec);
    fprintf (fp, "CPU timer cost est.:       %8.3g usec per start/stop pair\n", 
	     2.*deltat - gttdcost);
  }
    
  fprintf (fp, "CPU accumulation interval is %g seconds\n",
           1./(float) ticks_per_sec);

#ifdef HAVE_PCL
  mhz = PCL_determine_mhz_rate();
  fprintf (fp, "Clock speed is %d MHz\n", mhz);
#endif

  for (thread = 0; thread < numthreads; thread++) {

    /*
    ** Only print heading for threads that have 1 or more items to report
    */

    if (timers[thread] == NULL) continue;
    fprintf (fp, "\nStats for thread %d:\n", thread);
#ifdef THREADED_PTHREADS
    fprintf (fp, "pthread id %d:\n", (int) threadid[thread]);
#endif

    print_header (fp, max_indent_level[thread]);

    for (ptr = timers[thread]; ptr != NULL; ptr = ptr->next) {
      if (ptr->onflg) {

	fprintf (fp, "Timer %s was not off.  No stats will be printed\n",
		 ptr->name);

      } else {

	fillstats (&stats, ptr);

	/*
	** Print stats indented.  If indent_level = AMBIGUOUS (a negative 
	** number) the name will be printed with no indentation.
	*/

	for (indent = 0; indent < ptr->indent_level; indent++)
	  fprintf (fp, "  ");
	fprintf (fp, "%-15s", ptr->name);

	/*
	** If indent_level = AMBIGUOUS (a negative number) we want to loop 
	** from 0
	*/

	ilstart = MAX (0, ptr->indent_level);
	for (indent = ilstart; indent < max_indent_level[thread]; indent++)
	  fprintf (fp, "  ");

	print_stats_line (fp, &stats);
      }
    }

    if (usrsysenabled || wallenabled)
      fprintf (fp, "\nTIMER OVERHEAD (wallclock seconds) = %12.6f\n", 
	       overhead[thread]);

    if (pcl_cyclesenabled)
      fprintf (fp, "TIMER OVERHEAD (cycles) = %12.6e\n", 
	       (double) overhead_pcl[thread]);
  }

  /*
  ** Print a vertical summary if data exist for more than 1 thread.  The "2"
  ** passed to print_header is so we'll get an extra 4 spaces of indentation
  ** due to the thread number appended to the timer name.
  */

  if (numthreads > 0 && timers[1] != NULL) {
    fprintf (fp, "\nSame stats sorted by timer with thread number appended:\n");
    print_header (fp, 2);

    /*
    ** Print stats for slave threads that match master
    */

    for (ptr = timers[0]; ptr != NULL; ptr = ptr->next) {

      char name[20];
      Boolean found = false;

      /*
      ** Don't bother printing summation stats when only the master thread
      ** invoked the timer
      */

      for (thread = 1; thread < numthreads; thread++)
	for (tptr = timers[thread]; tptr != NULL; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name))
	    found = true;
	}
      if ( ! found) continue;

      /*
      ** Initialize stats which sum over threads
      */

      memset (&threadsum, 0, sizeof (threadsum));

      if ( ! ptr->onflg) {
	fillstats (&stats, ptr);
	strcpy (name, ptr->name);
	strcat (name, ".0");
	fprintf (fp, "%-19s", name);
	print_stats_line (fp, &stats);
	threadsum = stats;
      }

      /*
      ** loop over slave threads, printing stats for each and accumulating
      ** sum over threads when the name matches
      */

      for (thread = 1; thread < numthreads; thread++) {
	for (tptr = timers[thread]; tptr != NULL; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name)) {
	    if ( ! tptr->onflg) {
	      char num[5];

	      fillstats (&stats, tptr);
	      strcpy (name, tptr->name);
	      sprintf (num, ".%-3d", thread);
	      strcat (name, num);
	      fprintf (fp, "%-19s", name);
	      print_stats_line (fp, &stats);

	      threadsum.usr      += stats.usr;
	      threadsum.sys      += stats.sys;
	      threadsum.usrsys   += stats.usrsys;
	      threadsum.elapse   += stats.elapse;
	      threadsum.max_wtime = MAX (threadsum.max_wtime, stats.max_wtime);
	      threadsum.min_wtime = MIN (threadsum.min_wtime, stats.min_wtime);
	      threadsum.count    += stats.count;

	      for (n = 0; n < ncounter; n++)
		threadsum.pcl_result[n] += stats.pcl_result[n];
	    }
	    break; /* Go to the next thread */
	  }        /* if (STRMATCH (ptr->name, tptr->name) */
	}          /* loop thru linked list of timers for this thread */
      }            /* loop over slave threads */

      strcpy (name, ptr->name);
      strcat (name, ".sum");
      fprintf (fp, "%-19s", name);
      print_stats_line (fp, &threadsum);
      fprintf (fp, "\n");

    } /* loop through master timers */

    for (thread = 0; thread < numthreads; thread++) {
      if (usrsysenabled || wallenabled)
	fprintf (fp, "OVERHEAD.%-3d (wallclock seconds) = %12.6f\n", 
		 thread, overhead[thread]);

      if (pcl_cyclesenabled)
	fprintf (fp, "OVERHEAD.%-3d (cycles) = %12.6e\n", 
		 thread, (double) overhead_pcl[thread]);
    }

  } /* if (numthreads > 0 && timers[1] != NULL */

  return 0;
}

void fillstats (struct Stats *stats, struct node *ptr)
{
  int n;

  stats->usr       = ptr->accum_utime / (float) ticks_per_sec;
  stats->sys       = ptr->accum_stime / (float) ticks_per_sec;
  stats->usrsys    = stats->usr + stats->sys;
  stats->elapse    = ptr->accum_wtime_sec + 1.e-6 * ptr->accum_wtime_usec;
  stats->max_wtime = ptr->max_wtime;
  stats->min_wtime = ptr->min_wtime;
  stats->count     = ptr->count;

  for (n = 0; n < ncounter; n++)
    stats->pcl_result[n] = ptr->accum_pcl_result[n];
}

void print_stats_line (FILE *fp, struct Stats *stats)
{
  int index;
  int n;
  long long cycles;
  long long instr;
  long long flops;
  long long loadstore;
  long long l2cache;
  long long jump;

  float mflops;
  float ipc;
  float memfp;

  fprintf (fp, "%9ld ", stats->count);

  if (usrsysenabled)
    fprintf (fp, "%9.3f %9.3f %9.3f ", stats->usr, stats->sys, stats->usrsys);

  if (wallenabled)
    fprintf (fp, "%9.3f %9.3f %9.3f ", stats->elapse, stats->max_wtime, 
	     stats->min_wtime);
  
  for (n = 0; n < nevent; n++) {
    if (event[n]->name > pcl_start && event[n]->name < pcl_end) {
      index = event[n]->index;
      if (stats->pcl_result[index] > 1.e6)
	fprintf (fp, "%9.3e ", (double) stats->pcl_result[index]);
      else
	fprintf (fp, "%9ld ", (long) stats->pcl_result[index]);
    }
  }

  fprintf (fp, "\n");
}

void print_header (FILE *fp, int indent_level)
{
  int i;
  int n;

  fprintf (fp, "Name           ");
  for (i = 0; i < indent_level; i++)
    fprintf (fp, "  ");
  fprintf (fp, "Called    ");

  if (usrsysenabled)
    fprintf (fp, "Usr       Sys       Usr+Sys   ");

  if (wallenabled)
    fprintf (fp, "Wallclock Max       Min       ");

  for (n = 0; n < nevent; n++)
    if (event[n]->name > pcl_start && event[n]->name <= pcl_end)
      fprintf (fp, event[n]->string);

  fprintf (fp, "\n");
}
