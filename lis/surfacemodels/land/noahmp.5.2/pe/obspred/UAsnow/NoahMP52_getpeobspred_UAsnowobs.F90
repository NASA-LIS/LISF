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
! !ROUTINE: NoahMP52_getpeobspred_UAsnowobs
!  \label{NoahMP52_getpeobspred_UAsnowobs}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine NoahMP52_getpeobspred_UAsnowobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use NoahMP52_lsmMod, only : NoahMP52_struc
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS:
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!
!  This routine retrieves the observation prediction, which is the
!  model's estimate of snow depth.
!
!EOP
  integer                :: n
  type(ESMF_Field)       :: snodField
  real, pointer          :: snod(:)
  type(ESMF_Field)       :: sweField
  real, pointer          :: swe(:)
  integer                :: t
  integer                :: status


  n = 1

  call ESMF_StateGet(Obj_Func,"UA_SNOD",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(Obj_Func,"UA_SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snod(t) = NoahMP52_struc(n)%noahmp52(t)%snowh*1000.0 !mm
     swe(t) = NoahMP52_struc(n)%noahmp52(t)%sneqv
  enddo


end subroutine NoahMP52_getpeobspred_UAsnowobs



