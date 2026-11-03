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
! !MODULE: reset_merra2bias
! \label{reset_merra2bias}
!
! !REVISION HISTORY:
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !INTERFACE:
subroutine reset_merra2bias()
! !USES:
  use LIS_coreMod,  only   : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use merra2bias_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for MERRA2bias forcing.
!
!EOP
  implicit none

  integer :: n

  do n = 1,LIS_rc%nnest
     merra2bias_struc(n)%startFlag = .true.
     merra2bias_struc(n)%dayFlag = .true.
     merra2bias_struc(n)%merra2biastime1 = 3000.0
     merra2bias_struc(n)%merra2biastime2 = 0.0
     merra2bias_struc(n)%ringtime = 0.0
     merra2bias_struc(n)%reset_flag = .true.
  enddo

end subroutine reset_merra2bias

