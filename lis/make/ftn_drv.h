#ifdef IRIX64
#define FTN(func) func##_
#elif ABSOFT
#define FTN(func) func##_
#elif GFORTRAN
#define FTN(func) func##_
#elif PGI
#define FTN(func) func##_
#elif LAHEY
#define FTN(func) func##_
#elif OSF1
#define FTN(func) func##_
#elif IFC
#define FTN(func) func##_
#else
#define FTN(func) func
#endif


