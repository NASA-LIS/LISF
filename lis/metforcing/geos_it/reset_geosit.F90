!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: reset_geosit
! \label{reset_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine reset_geosit
! !USES:
      use LIS_coreMod,  only   : LIS_rc
      use LIS_timeMgrMod, only : LIS_date2time
      use geosit_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GEOS-IT forcing.
!
!EOP
      implicit none

      integer :: n

      do n = 1,LIS_rc%nnest
         geosit_struc(n)%startFlag = .true.
         geosit_struc(n)%dayFlag = .true.
         geosit_struc(n)%geosittime1 = 3000.0
         geosit_struc(n)%geosittime2 = 0.0
         geosit_struc(n)%ringtime = 0.0
         geosit_struc(n)%reset_flag = .true.
      enddo

      end subroutine reset_geosit

