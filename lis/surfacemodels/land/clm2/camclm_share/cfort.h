#if ( defined SUNOS ) || ( defined IRIX64 ) || ( defined OSF1 ) || ( defined LINUX )
#define FORTRANUNDERSCORE
#elif ( defined CRAY ) || ( defined T3D )
#define FORTRANCAPS
#elif ( defined AIX )
#endif
