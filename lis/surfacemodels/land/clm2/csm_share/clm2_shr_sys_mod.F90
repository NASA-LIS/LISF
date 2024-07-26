!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

MODULE clm2_shr_sys_mod

   use clm2_shr_kind_mod, only: CLM2_SHR_KIND_IN, CLM2_SHR_KIND_R8, CLM2_SHR_KIND_I8  ! defines real & integer kinds
   use LIS_logMod, only : LIS_logunit
#if ( defined ABSOFT )
#if (! defined HIDE_MPI)
   use mpi
#endif
#endif
#if ( ! defined ABSOFT )
#if (! defined HIDE_MPI)
   use LIS_mpiMod
#endif
#endif
!
! -- Don't use IMPLICIT none here as it interfere's with the "include "mpif.h" line"
!
! PUBLIC: Public interfaces
!
   private
   public :: clm2_shr_sys_system, clm2_shr_sys_chdir, &
             clm2_shr_sys_getenv, clm2_shr_sys_flush, &
             clm2_shr_sys_abort , clm2_shr_sys_irtc,  &
             clm2_shr_sys_sleep 

CONTAINS

!===============================================================================

SUBROUTINE clm2_shr_sys_system(str,rcode)

   IMPLICIT none

   !----- arguments ---
   character(len=*)    ,intent(in)  :: str    ! system/shell command string
   integer(CLM2_SHR_KIND_IN),intent(out) :: rcode  ! function return error code

#if (defined CRAY)
   !----- functions -----
   integer(CLM2_SHR_KIND_IN),external    :: ishell ! function to envoke shell command
#endif
#if (defined OSF1 || defined SUNOS || (defined LINUX && !defined __GFORTRAN__ && !defined CATAMOUNT))
   !----- functions -----
   integer(CLM2_SHR_KIND_IN),external    :: system ! function to envoke shell command
#endif

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   rcode = 0
#if (defined CRAY)
   rcode=ishell(str)
#endif

#if (defined IRIX64)
   call system(str)
#endif

#if (defined AIX)
   call system(str,rcode)
#endif

#if (defined OSF1 || defined SUNOS || defined __GFORTRAN__ || (defined LINUX && !defined CATAMOUNT))
   rcode = system(str)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX)
   write(*,*) '(clm2_shr_sys_system) ERROR: no implementation for this architecture'
   call clm2_shr_sys_abort('no system routine on this machine')
#endif

END SUBROUTINE clm2_shr_sys_system

!===============================================================================

SUBROUTINE clm2_shr_sys_chdir(path, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: path    ! chdir to this dir
   integer(CLM2_SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(CLM2_SHR_KIND_IN)             :: lenpath ! length of path
#if (defined AIX || defined OSF1 || defined SUNOS || (defined LINUX && !defined __GFORTRAN__) || defined NEC_SX || defined CPRINTEL)
   integer(CLM2_SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenpath=len_trim(path)

#if (defined IRIX64 || defined CRAY)
   call pxfchdir(path, lenpath, rcode)
#endif

#if (defined AIX)
   rcode=chdir(%ref(path(1:lenpath)//char(0)))
#endif

#if (defined OSF1 || defined SUNOS || defined NEC_SX || defined Darwin || (defined LINUX && !defined CPRNAG))
   rcode=chdir(path(1:lenpath))
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX)
   write(*,*) '(clm2_shr_sys_chdir) ERROR: no implementation for this architecture'
   call clm2_shr_sys_abort('no implementation of chdir for this machine')
#endif

END SUBROUTINE clm2_shr_sys_chdir

!===============================================================================

SUBROUTINE clm2_shr_sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(len=*)    ,intent(in)  :: name    ! env var name
#if ( defined ABSOFT )
   character(len=*)    ,intent(inout) :: val     ! env var value
#else
   character(len=*)    ,intent(out) :: val     ! env var value
#endif
   integer(CLM2_SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(CLM2_SHR_KIND_IN)             :: lenname ! length of env var name
   integer(CLM2_SHR_KIND_IN)             :: lenval  ! length of env var value
   character(len=512)               :: tmpval  ! temporary env var value

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenname=len_trim(name)

#if (defined IRIX64 || defined CRAY)
   call pxfgetenv(name, lenname, val, lenval, rcode)
#endif

#if (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0 ) rcode = 1
   if (len_trim(val) > 512 ) rcode = 2
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX)
   write(*,*) '(clm2_shr_sys_getenv) ERROR: no implementation for this architecture'
   call clm2_shr_sys_abort('no implementation of getenv for this machine')
#endif

END SUBROUTINE clm2_shr_sys_getenv

!===============================================================================

SUBROUTINE clm2_shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(CLM2_SHR_KIND_IN) :: unit  ! flush output buffer for this unit

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

#if (defined IRIX64 || defined CRAY || defined OSF1 || defined SUNOS || defined LINUX)
   call flush(unit)
#endif
#if (defined AIX)
   call flush_(unit)
#endif

#if (!defined CRAY && !defined IRIX64 && !defined AIX && !defined OSF1 && !defined SUNOS && !defined LINUX)
   write(*,*) '(clm2_shr_sys_flush) WARNING: no implementation for this architecture'
#endif

END SUBROUTINE clm2_shr_sys_flush

!===============================================================================

SUBROUTINE clm2_shr_sys_abort(string)

   IMPLICIT none

   character*(*),optional :: string    ! error message string

   !----- local -----
   integer(CLM2_SHR_KIND_IN) :: rcode,ierr
   logical              :: flag

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('(clm2_shr_sys_abort) ',a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   call clm2_shr_sys_flush(6)
   if (len_trim(string) > 0) write(6,F00) 'ERROR: '//trim(string)
#if (! defined HIDE_MPI)
   write(6,F00) 'WARNING: calling mpi_abort() and stopping'
   call clm2_shr_sys_flush(6)
   call mpi_initialized(flag,ierr)
   if (flag) call mpi_abort(LIS_mpi_comm,rcode,ierr)
#else
   write(6,F00) 'WARNING: stopping'
#endif
   call clm2_shr_sys_flush(6)
   call abort()
   stop

END SUBROUTINE clm2_shr_sys_abort

!===============================================================================

integer(CLM2_SHR_KIND_I8) FUNCTION clm2_shr_sys_irtc( rate )

   IMPLICIT none
   !----- optional output argument -----
   integer(CLM2_SHR_KIND_I8), optional :: rate

   !----- local -----
   integer(CLM2_SHR_KIND_IN)          :: count
   integer(CLM2_SHR_KIND_IN)          :: count_rate
   integer(CLM2_SHR_KIND_IN)          :: count_max

   integer(CLM2_SHR_KIND_IN),save :: last_count = -1
   integer(CLM2_SHR_KIND_I8),save :: count_offset = 0

!-------------------------------------------------------------------------------
! emulates Cray/SGI irtc function (returns clock tick since last reboot)
!-------------------------------------------------------------------------------
   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   clm2_shr_sys_irtc = count
!
! System clock is a 24-hour clock, so must check each time pass over midnight
!
   if ( last_count /= -1 )then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   clm2_shr_sys_irtc = clm2_shr_sys_irtc + count_offset
   last_count = count

END FUNCTION clm2_shr_sys_irtc

!===============================================================================

SUBROUTINE clm2_shr_sys_sleep(sec)

   IMPLICIT none

   !----- input -----
   real   (clm2_shr_kind_r8),intent(in) :: sec  ! number of seconds to sleep

   !----- local -----
   integer(clm2_shr_kind_in) :: isec   ! integer number of seconds
   integer(clm2_shr_kind_in) :: rcode  ! return code
   character(len=90) :: sleep_var

   !----- i/o formats -----
   character(len=*),parameter :: F00 = "('sleep ',i8 )"

   save

!-------------------------------------------------------------------------------
! PURPOSE: Sleep for approximately sec seconds
!-------------------------------------------------------------------------------

   isec = nint(sec)

   if (isec <= 0.0) then
      write(LIS_logunit,*) 'ERROR: seconds must be > 0, sec=',sec
   else
      write(sleep_var,FMT=F00) isec
      call clm2_shr_sys_system( sleep_var, rcode )
   endif

END SUBROUTINE clm2_shr_sys_sleep

!===============================================================================

END MODULE clm2_shr_sys_mod
