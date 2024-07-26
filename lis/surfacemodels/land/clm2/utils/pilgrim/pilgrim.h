!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!-----------------------------------------------------------------------
! MPI Context:
!-----------------------------------------------------------------------
!
!  Now only contains cpp defines.  This file should really be called
!     mpi_defines.h
!

#define  ROOT_PE      0
#define  MAX_PE    1024
#define  MAX_TRF     10
#define  MAX_PAX MAX_PE

!
! Max buffer size for Shared Memory Arenas
!
#define  MAX_BUF  10000000

#if defined(CRAY)
#define CPP_MPI_INTEGER MPI_INTEGER
#define CPP_INTEGER i8
#else
#define CPP_MPI_INTEGER MPI_INTEGER
#define CPP_INTEGER i4
#endif

#if defined(CRAY)
#define CPP_MPI_REAL MPI_REAL
#define CPP_REAL     r8
#else
#define CPP_MPI_REAL MPI_DOUBLE_PRECISION
#define CPP_REAL     r8
#endif
#define CPP_REAL_WIDTH 8
