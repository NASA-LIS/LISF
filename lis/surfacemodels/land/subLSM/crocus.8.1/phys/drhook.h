#ifndef _DRHOOK_H_
#define _DRHOOK_H_

#ifdef _DRHOOK_C_

#if defined(__GNUC__)
#define _GNU_SOURCE 1
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <ctype.h>
#include <signal.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <pthread.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef RS6K
#include <fptrap.h>
#endif

#ifdef VPP
#include <ucontext.h>
#endif

int drhook_lhook = 1;
#else
extern int drhook_lhook;
#endif

#ifndef O_LOCK_DONE
#define O_LOCK_DONE

/* OpenMP/ODB lock type */
/* Keep consistent with "odb/include/privpub.h" */
/* Be ALSO consistent with OML_LOCK_KIND in ifsaux/module/oml_mod.F90 */

typedef long long int o_lock_t; /* i.e. 64-bit integer */

#define INIT_LOCKID_WITH_NAME(mylock, lockname) \
  coml_init_lockid_with_name_(mylock, lockname, strlen(lockname))

extern void coml_set_debug_(const int *konoff, int *kret);
extern void coml_init_lock_();
extern void coml_init_lockid_(o_lock_t *mylock);
extern void coml_init_lockid_with_name_(o_lock_t *mylock, const char *name, int name_len);
extern void coml_set_lock_();
extern void coml_set_lockid_(o_lock_t *mylock);
extern void coml_unset_lock_();
extern void coml_unset_lockid_(o_lock_t *mylock);
extern void coml_test_lock_(int *is_set);
extern void coml_test_lockid_(int *is_set, o_lock_t *mylock);
extern void coml_in_parallel_(int *is_parallel_region);

#endif

/* drhook.c external interfaces */

extern void
c_drhook_getenv_(const char *s,
                 char *value,
                 /* Hidden arguments */
                 int slen,
                 const int valuelen);

extern void
c_drhook_memcounter_(const int *thread_id,
                     const long long int *size,
                     long long int *keyptr_addr);

extern void
c_drhook_raise_(const int *sig);

extern void
c_drhook_print_(const int *ftnunitno,
                const int *thread_id,
                const int *print_option, /* 
                                            1=raw call counts 
                                            2=calling tree
                                            3=profiling info
                                         */
                int *level);

extern void
c_drhook_init_signals_(const int *enforce);

extern void
c_drhook_set_lhook_(const int *lhook);

extern void 
c_drhook_init_(const char *progname,
               const int *num_threads
               /* Hidden length */
               ,int progname_len);

extern void
c_drhook_start_(const char *name, 
                const int *thread_id, 
                double *key,
                const char *filename,
                const int *sizeinfo
                /* Hidden length */
                ,int name_len, int filename_len);

extern void
c_drhook_end_(const char *name,
              const int *thread_id,
              const double *key,
              const char *filename,
              const int *sizeinfo
              /* Hidden length */
              ,int name_len, int filename_len);

extern void
c_drhook_watch_(const int *onoff,
		const char *array_name,
		const void *array_ptr,
		const int *nbytes,
		const int *abort_if_changed,
		const int *printkey,
		const int *nvals,
		const int *print_traceback_when_set
		/* Hidden length */
		,int array_name_len);

extern void
c_drhook_check_watch_(const char *where,
		      const int *allow_abort
		      /* Hidden length */
		      , int where_len);

/**** C-interface to Dr.Hook ****/

extern void
Dr_Hook(const char *name, int option, double *handle, 
        const char *filename, int sizeinfo,
        int name_len, int filename_len);

#define DRHOOK_START_RECUR(name,recur) \
  static const char *drhook_name = #name; \
  static const int drhook_name_len = sizeof(#name) - 1; /* Compile time eval */ \
  static const char *drhook_filename = __FILE__; \
  static const int drhook_filename_len = sizeof(__FILE__) - 1; /* Compile time eval */ \
  double zhook_handle; \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 0, &zhook_handle, \
                                      drhook_filename, 0, \
                                      drhook_name_len, drhook_filename_len); {

#define DRHOOK_START(name) DRHOOK_START_RECUR(name,0)

#define DRHOOK_START_BY_STRING_RECUR(name, recur) \
  static const char *drhook_name = name; \
  static const int drhook_name_len = sizeof(name) - 1; /* Compile time eval */ \
  static const char *drhook_filename = __FILE__; \
  static const int drhook_filename_len = sizeof(__FILE__) - 1; /* Compile time eval */ \
  double zhook_handle; \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 0, &zhook_handle, \
                                      drhook_filename, 0, \
                                      drhook_name_len, drhook_filename_len); {

#define DRHOOK_START_BY_STRING(name) DRHOOK_START_BY_STRING_RECUR(name,0)

#define DRHOOK_RETURN_RECUR(sizeinfo,recur) \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 1, &zhook_handle, \
                                      drhook_filename, sizeinfo, \
                                      drhook_name_len, drhook_filename_len)

#define DRHOOK_RETURN(sizeinfo) DRHOOK_RETURN_RECUR(sizeinfo,0)

#define DRHOOK_END_RECUR(sizeinfo,recur) ; } DRHOOK_RETURN_RECUR(sizeinfo,recur)

#define DRHOOK_END(sizeinfo) DRHOOK_END_RECUR(sizeinfo,0) 

/* Fortran routines */

extern void
dr_hook_prt_(const int *ftnunitno,
             const char *s
             /* Hidden arguments */
             , int s_len);

extern void
dr_hook_procinfo_(int *myproc, int *nproc);

#endif  /* _DRHOOK_H_ */

