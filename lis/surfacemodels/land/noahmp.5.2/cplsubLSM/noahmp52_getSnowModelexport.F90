!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp52_getSnowModelexport
! \label{noahmp52_getSnowModelexport}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_getSnowModelexport(n, LSM2SUBLSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP52_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM2SUBLSM_State
!
! !DESCRIPTION:
!
!
!EOP

end subroutine noahmp52_getSnowModelexport


