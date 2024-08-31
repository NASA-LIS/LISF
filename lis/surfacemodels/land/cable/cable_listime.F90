!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module cable_listime

      PRIVATE

      public :: calc_localtime

      CONTAINS

!BOP
! !ROUTINE: calc_localtime
! !INTERFACE:
      subroutine calc_localtime(gmt,lon,lhour,change)
! !ARGUMENTS:
      integer :: gmt            ! GMT time (0-23)
      real    :: lon            ! longitude in degrees
      integer :: change         ! the change in number of hours between
      integer :: lhour          ! local hour (0-23) 0= midnight, 23= 11:00 p.m.

      integer :: i              ! working integer
      integer :: zone           ! time zone (1-24)
!EOP
      change = 0 
!----------------------------------------------------------------------
! Determine into which time ZONE (15 degree interval) the
! longitude falls.
!----------------------------------------------------------------------
      do i = 1,25
         if (lon.lt.(-187.5+(15*i))) then
            zone=i
            if (zone.eq.25) zone=1
            goto 60
         endif
      enddo

!----------------------------------------------------------------------
! Calculate change (in number of hours) from GMT time to
! local hour.  Change will be negative for zones < 13 and
! positive for zones > 13.
! There is also a correction for LHOUR < 0 and LHOUR > 23
! to LHOUR between 0 and 23.
!----------------------------------------------------------------------

 60   if (zone.lt.13) then
         change=zone-13
         lhour=gmt+change
      elseif (zone.eq.13) then
         lhour=gmt
      else
         change=zone-13
         lhour=gmt+change
      endif
      if (lhour.lt.0) lhour=lhour+24
      if (lhour.gt.23) lhour=lhour-24
      return
      end subroutine calc_localtime

      end module cable_listime
