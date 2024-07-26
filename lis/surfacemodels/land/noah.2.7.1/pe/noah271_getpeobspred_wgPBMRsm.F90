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
! !ROUTINE: noah271_getpeobspred_wgPBMRsm
!  \label{noah271_getpeobspred_wgPBMRsm}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_getpeobspred_wgPBMRsm(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use noah271_lsmMod, only : noah271_struc
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine defines the obspred for the walnut gulch 
!  PBMR soil moisture observations, which are in mass units. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: smField
  real, allocatable          :: soilm(:)
  integer                :: t
  integer                :: i
  integer                :: status

  n = 1
  call ESMF_StateGet(Obj_Func,"PBMR soil moisture",smField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smField,localDE=0,farrayPtr=soilm,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     soilm(t) = noah271_struc(n)%noah(t)%smc(1)*100.0
  enddo

end subroutine noah271_getpeobspred_wgPBMRsm



