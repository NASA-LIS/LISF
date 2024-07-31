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
module LIS_gswpMod
!BOP
!
! !MODULE: LIS_gswpMod
! 
! !DESCRIPTION:
!   This module contains useful routines that generates indices 
!   for reading the GSWP netcdf data. 
! 
! !REVISION HISTORY: 
!  24 Feb 2004    Sujay Kumar  Initial Specification
! 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: getgswp_monindex
  public :: getgswp_timeindex
!EOP
contains
!BOP
! 
! !ROUTINE: getgswp_monindex
! \label{getgswp_monindex}
! 
! !INTERFACE:
  subroutine getgswp_monindex(yr,mo,index)
! !ARGUMENTS:
    integer, intent(out) :: index
    integer, intent(in) :: yr, mo
! !DESCRIPTION: 
!  
!  This routine computes the GSWP monthly index (base time being 1982), 
!  to be used to read monthly data. 
!  
!EOP

    index = 0
    index = index + (yr-1982)*12 + mo
  end subroutine getgswp_monindex

!BOP
! 
! !ROUTINE: getgswp_timeindex
! \label{getgswp_timeindex}
! 
! !INTERFACE: 
  subroutine getgswp_timeindex(yr,mo,da,hr,mn,ss,index)
! 
    implicit none
! !ARGUMENTS: 
    integer, intent(in)  :: yr, mo, da, hr, mn, ss
    integer, intent(out) :: index
! !DESCRIPTION: 
!
!  This routine computes the GSWP timestep index (base time being 1982), 
!  to be used to read the forcing files. 
!  Note:  GSWP forcing data are written into monthly files with
!  a frequency of 3 hours.  The first entry corresponds to 03:00:00
!  of the first day of the given month.  The last entry corresponds
!  to 00:00:00 of the first day of the next month.
!
!  E.g.; for Tair\_cru198207.nc, the data run from 1982-07-01T03:00:00
!  through 1982-08-01T00:00:00, inclusive.
!
!  So, when you are at hour 0 on the first day of the month,
!  reset the day to the last day of the previous month, and reset the
!  hour from 0 to 24.  This will compute the correct index for the last
!  entry in the forcing file.
!
!  E.g.; 1982-08-01T00:00:00 should be re-written as 1982-07-31T24:00:00.
!
!EOP 
    integer :: ryr, rmo, rda, rhr, days1(12)
    integer :: tmp_da, tmp_hr
    data days1 /31,28,31,30,31,30,31,31,30,31,30,31/
    logical :: leap 
    ryr = 1982
    rmo = 7
    rda = 1
    rhr = 3
  
    index = 0
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then
       leap = .true.                   
    else
       leap = .false.
    endif

    if ( leap ) then
       days1(2) = 29
    endif

    tmp_da = da
    tmp_hr = hr

    if ( tmp_da == 1 .and. tmp_hr == 0 .and. mn == 0 .and. ss == 0) then
       if ( mo == 1 ) then
          tmp_da = days1(12)
       else
          tmp_da = days1(mo-1)
       endif
       tmp_hr = 24
    endif

    index = (tmp_da-1)*8 + (tmp_hr+mn/60.0+ss/3600)/3

  end subroutine getgswp_timeindex
end module LIS_gswpMod
