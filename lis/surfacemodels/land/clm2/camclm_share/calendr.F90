!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
   subroutine calendr(nstep   ,dtime   ,ndbase  ,nsbase  ,nbdate  ,  &
  &                   nbsec   ,ndcur   ,nscur   ,ncdate  ,ncsec   ,  &
  &                   calday  )
!-----------------------------------------------------------------------
!
! Purpose:
!
! Compute current date and day information for history tape header.
! Compute current julian day (including fraction) for sun angle
! computations and time interpolation of boundary dataset information.
! One year is defined as 365 days.
!
! Computational notes:
!
! 86400 is the number of seconds in 1 day.
!
! Dividing an integer by 10**n has the effect of right-shifting the
! decimal digits n positions (ex: 861231/100 = 008612).
!
! mod(integer,10**n) has the effect of extracting the rightmost
! n decimal digits of the integer (ex: mod(861231,10000) = 1231).
!
! Author: Jim Rosinski
!
! $Id: calendr.F90,v 1.5 2004/05/07 22:18:33 jim Exp $
! $Author: jim $
!
!-----------------------------------------------------------------------
   use LIS_precisionMod
   use clm2_shr_const_mod, ONLY: SHR_CONST_CDAY
   implicit none
 
! ------------------------ Arguments------------------------------------
   integer , intent(in)  :: nstep   ! current time step (0, 1, ...)
   real(r8), intent(in)  :: dtime   ! length of time step (seconds)
   integer , intent(in)  :: ndbase  ! base day of run (e.g., 0)
   integer , intent(in)  :: nsbase  ! base seconds of base day (e.g., 0)
   integer , intent(in)  :: nbdate  ! base date (yyyymmdd format) of
   integer , intent(in)  :: nbsec   ! base seconds of base date (e.g., 0)
   integer , intent(out) :: ndcur   ! current day (0, 1, ...)
   integer , intent(out) :: nscur   ! current seconds of current day
   integer , intent(out) :: ncdate  ! current date (yyyymmdd format)
   integer , intent(out) :: ncsec   ! current seconds of current date (0, ..., 86400)
   real(r8), intent(out) :: calday  ! current julian day including fraction
!-----------------------------------------------------------------------
 
!---------------------------Local workspace-----------------------------
   integer(i8) nstep8   ! current time step (0, 1, ...)
   integer(i8) nsecs    ! run time in seconds
   integer ndays        ! ndays + nyears = number of day changes since start
   integer nyears       ! ndays + nyears = number of day changes since start
   integer mcyear       ! current year of current date (0-99)
   integer mcmnth       ! current month of current date (1-12)
   integer mcday        ! current day of current date (1-31)
   integer jday         ! Julian day (1-365)
   integer mbmnth       ! base month (1-12) validity check
   integer mbday        ! base day of base date (1-31) validity check
   integer ndm(12)      ! number of days per month
   integer jdcon(12)    ! convert month index to julian day
   integer(i8) nspdy
!   integer(i8), parameter :: nspdy = nint(SHR_CONST_CDAY) ! number of seconds per day
   save ndm,jdcon
   data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
   data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
!-----------------------------------------------------------------------
!
! Check validity of input data
!
   nspdy = nint(SHR_CONST_CDAY)
   mbmnth = mod(nbdate,10000)/100
   mbday = mod(nbdate,100)
   if (mbmnth<1 .or. mbmnth>12) then
     write(6,*)' CALENDR: Invalid base month input:',mbmnth
     call endrun
   end if
   if (mbday<1 .or. mbday>ndm(mbmnth)) then
     write(6,*)' CALENDR: Invalid base day of base date input:',mbday
     call endrun
   end if
   if (nsbase<0 .or. nsbase>=nspdy) then
     write(6,*)' CALENDR: Invalid base seconds(nsbase):',nsbase
     call endrun
   end if
   if (nbsec<0 .or. nbsec>=nspdy) then
     write(6,*)' CALENDR: Invalid base seconds(nbsec):',nbsec
     call endrun
   end if
!
! First: current day, seconds
!
   nstep8 = nstep
   nsecs  = nstep8*nint(dtime)
   ndcur  = ndbase + (nsecs+nsbase)/nspdy
   nscur  = mod((nsecs+nsbase),nspdy)
!
! Next: current date, seconds.
! Uncommenting the next line will have the effect of modifying
! nsecs to include base day and base seconds.  This is the way
! the current date was computed in CCM1.
!
!   nsecs = nsecs + ndbase*nspdy + nsbase
 
   ndays = (nsecs+nbsec)/nspdy
   nyears = ndays/365
   ndays = mod(ndays,365)
!
! Current seconds of current date
!
   ncsec = mod(nsecs+nbsec,nspdy)
!
! Initialize current year, month, day.
!
   mcyear = nbdate/10000 + nyears
   mcmnth = mod(nbdate,10000)/100
   mcday = mod(nbdate,100) + ndays
!
! Now loop through months, converting yyyy, mm, and ddd to yyyymmdd.
! ex: 791235 becomes 800104. 190001370 becomes 19010105.
!
   do while (mcday > ndm(mcmnth) )
     mcday = mcday - ndm(mcmnth)
     mcmnth = mcmnth + 1
     if (mcmnth == 13) then             ! add a year
       mcyear = mcyear + 1
       mcmnth = 1
     end if
   end do
   ncdate = mcyear*10000 + mcmnth*100 + mcday
!
! Convert current month, day, seconds to Julian day + fraction.  Multiply
! nspdy by floating point 1 in order to ensure real*8 arithmetic.  As of
! 9/14/98 this only makes a difference on Linux.
!
   jday = jdcon(mcmnth) + mcday
   calday = float(jday) + float(ncsec)/(1.*nspdy)
 
   return
   end subroutine calendr
 
