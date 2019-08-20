!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_climatology
! \label{reset_climatology}
! 
! 
! !INTERFACE:
subroutine reset_climatology
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use climatology_forcingMod, only : clim_struc
!
!EOP  
  implicit none

  integer   :: findex

  clim_struc%metforc_time1 = 3000.0
  clim_struc%metforc_time2 = 0.0

end subroutine reset_climatology
