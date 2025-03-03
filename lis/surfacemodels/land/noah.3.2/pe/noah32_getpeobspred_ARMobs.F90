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
! !ROUTINE: noah32_getpeobspred_ARMobs
!  \label{noah32_getpeobspred_ARMobs}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah32_getpeobspred_ARMobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use noah32_lsmMod, only : noah32_struc
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: qleField
  real, allocatable          :: qle(:)
  type(ESMF_Field)       :: qhField
  real, allocatable          :: qh(:)
  type(ESMF_Field)       :: qgField
  real, allocatable          :: qg(:)
  type(ESMF_Field)       :: sfsmField
  real, allocatable          :: sfsm(:)
  type(ESMF_Field)       :: sfstField
  real, allocatable          :: sfst(:)
  integer                :: t
  integer                :: i
  integer                :: status

  n = 1
  call ESMF_StateGet(Obj_Func,"ARM Latent Heat Flux",qleField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qleField,localDE=0,farrayPtr=qle,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(Obj_Func,"ARM Sensible Heat Flux",qhField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qhField,localDE=0,farrayPtr=qh,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(Obj_Func,"ARM Ground Heat Flux",qgField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qgField,localDE=0,farrayPtr=qg,rc=status)
  call LIS_verify(status)
#if 0 
  call ESMF_StateGet(Obj_Func,"ARM Surface Soil Moisture",sfsmField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sfsmField,localDE=0,farrayPtr=sfsm,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(Obj_Func,"ARM Surface Soil Temperature",sfstField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sfstField,localDE=0,farrayPtr=sfst,rc=status)
  call LIS_verify(status)
#endif
  do t=1,LIS_rc%ntiles(n)
     qle(t) = noah32_struc(n)%noah(t)%qle
     qh(t) = noah32_struc(n)%noah(t)%qh
     qg(t) = noah32_struc(n)%noah(t)%qg
!     sfsm(t) = noah32_struc(n)%noah(t)%smc(1)
!     sfst(t) = noah32_struc(n)%noah(t)%stc(1)
  enddo

end subroutine noah32_getpeobspred_ARMobs



