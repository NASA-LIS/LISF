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
!BOP
! !MODULE: reset_metForcGenerated
! \label{reset_metForcGenerated}
! 
! 
! !INTERFACE:
subroutine reset_metForcGenerated
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use metForcGenerated_forcingMod, only : metForcGen_struc
!
!
!EOP  
  implicit none

  integer   :: findex

  metForcGen_struc%metforc_time1 = 3000.0
  metForcGen_struc%metforc_time2 = 0.0


end subroutine reset_metForcGenerated
