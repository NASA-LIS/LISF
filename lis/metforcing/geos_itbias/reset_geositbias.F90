!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: reset_geositbias
! \label{reset_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
!
! !INTERFACE:
      subroutine reset_geositbias
! !USES:
      use LIS_coreMod,  only   : LIS_rc
      use LIS_timeMgrMod, only : LIS_date2time
      use geositbias_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GEOS-ITbias forcing.
!
!EOP
      implicit none

      integer :: n

      do n = 1,LIS_rc%nnest
         geositbias_struc(n)%startFlag = .true.
         geositbias_struc(n)%dayFlag = .true.
         geositbias_struc(n)%geositbiastime1 = 3000.0
         geositbias_struc(n)%geositbiastime2 = 0.0
         geositbias_struc(n)%ringtime = 0.0
         geositbias_struc(n)%reset_flag = .true.
      enddo

      end subroutine reset_geositbias

