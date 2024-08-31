//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include <string.h>

/*
** Required OMP calls not available on AIX so use PTHREADS
*/

#ifdef AIX
#define THREADED_PTHREADS
#endif

#if ( defined IRIX64 ) || ( defined LINUX ) || ( defined OSF1 )
#undef THREADED_OMP
#endif

#if ( defined USE_GCC)
#undef THREADED_OMP
#endif

/*
** Threading currently doesn't work on SUN so don't enable pthreads or omp
*/

#ifdef SUNOS
#endif

#if ( defined THREADED_OMP )
#include <omp.h>
#elif ( defined THREADED_PTHREADS )
#include <pthread.h>
#endif

#ifdef HAVE_PCL
#include <pcl.h>
#else

/*
** Dummy up pcl stuff if library unavailable
*/

typedef int PCL_CNT_TYPE;
typedef int PCL_FP_CNT_TYPE;
typedef int PCL_DESCR_TYPE;
#define PCL_MODE_USER       -1
#define PCL_L1DCACHE_MISS   -1
#define PCL_L2CACHE_MISS    -1
#define PCL_CYCLES          -1
#define PCL_ELAPSED_CYCLES  -1
#define PCL_FP_INSTR        -1
#define PCL_LOADSTORE_INSTR -1
#define PCL_INSTR           -1
#define PCL_STALL           -1
#define PCL_COUNTER_MAX      1
#define PCL_SUCCESS          0
extern int PCLread (PCL_DESCR_TYPE, PCL_CNT_TYPE *, PCL_CNT_TYPE *, int);
#endif

#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif

#define STRMATCH(X,Y) (strcmp((X),(Y)) == 0)
#define MAX_CHARS 15
#define AMBIGUOUS -1
#define MAX_THREADS 128

typedef enum {false = 0, true = 1} Boolean;

/*
** User specifiable options.  The values must match their counterparts in header.inc
** Also, we must have pcl_start < all valid pcl values < pcl_end.
** To add a new PCL counter: 
** 1) add the new entry to OptionName below.
** 2) add the appropriate array entry for possible_event[] to t_initialize.c.
** 3) add the appropriate code to the "switch" construct in t_initialize.c
*/

typedef enum {
  usrsys               = 1,
  wall                 = 2,
  pcl_start            = 3,   /* bogus entry delimits start of PCL stuff */
#ifdef HAVE_PCL
  pcl_l1dcache_miss    = 4,
  pcl_l2cache_miss     = 5,
  pcl_cycles           = 6,
  pcl_elapsed_cycles   = 7,
  pcl_fp_instr         = 8,
  pcl_loadstore_instr  = 9,
  pcl_instr            = 10,
  pcl_stall            = 11,
#endif
  pcl_end              = 12   /* bogus entry delimits end of PCL stuff */
} OptionName;

struct Event {
  OptionName name;
  char string[9];
  int index;
};

struct PossibleEvent {
  OptionName name;
  Boolean enabled;
  char string[10];
};

struct node {
  char name[MAX_CHARS+1];
  
  int indent_level;        /* indentation level of timer */

  long last_utime;         /* user time from most recent call */
  long last_stime;         /* system time from most recent call */
  long last_wtime_sec;     /* wallclock seconds from most recent call */
  long last_wtime_usec;    /* wallclock microseconds from most recent call */

  long accum_utime;        /* accumulated user time */
  long accum_stime;        /* accumulated system time */
  long accum_wtime_sec;    /* accumulated wallclock seconds */
  long accum_wtime_usec;   /* accumulated wallclock microseconds */

  float max_wtime;         /* maximum wallclock time for each start-stop */
  float min_wtime;         /* minimum wallclock time for each start-stop */

  PCL_CNT_TYPE last_pcl_result[PCL_COUNTER_MAX];
  PCL_CNT_TYPE accum_pcl_result[PCL_COUNTER_MAX];

  Boolean onflg;           /* true => timer is currently on */
  long count;              /* number of calls to t_start for this timer */

  struct node *next;       /* next timer in the linked list */
};

/*
** Globals
*/

extern struct node **timers;
extern struct node **last;
extern struct Options options;
extern long ticks_per_sec;
extern int numthreads;
extern int *max_indent_level;
extern float *overhead;
extern PCL_CNT_TYPE *overhead_pcl;
extern Boolean t_initialized;
extern Boolean wallenabled;
extern Boolean usrsysenabled;
extern struct PossibleEvent possible_event[];

/*
** Needed by PCL library: otherwise unused
*/

extern int counter_list[PCL_COUNTER_MAX];
extern int ncounter;     /* number of counters */
extern int cycles;       /* index of cycle counter */
extern int instr;        /* index of instruction counter */
extern int fp_instr;     /* index of fp instruction counter */
extern int l2cache_miss; /* index of l2cache miss */
extern int jump;
extern PCL_DESCR_TYPE *descr;
extern int nevent;
extern struct Event **event;
extern Boolean pclenabled;
extern Boolean pcl_cyclesenabled;
extern int pcl_cyclesindex;
extern int npossible;

#if ( defined THREADED_OMP )
extern omp_lock_t lock;
#elif ( defined THREADED_PTHREADS )
extern pthread_mutex_t t_mutex;
extern pthread_t *threadid;
#endif

/*
** Function prototypes
*/

extern int add_new_thread (void);
extern int get_cpustamp (long *, long *);
extern int get_thread_num (void);
extern int t_error (const char *, ...);
extern int t_initialize (void);
extern int t_pr (int);
extern int t_reset (void);
extern int t_setoption (OptionName, Boolean);
extern int t_stamp (double *, double *, double *);
extern int t_start (char *);
extern int t_stop (char *);
extern char *t_pclstr (int);
