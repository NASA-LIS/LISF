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
! !ROUTINE: noah33_getpeobspred_LPRM_AMSREsmObs
!  \label{noah33_getpeobspred_LPRM_AMSREsmObs}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah33_getpeobspred_LPRM_AMSREsmObs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use noah33_lsmMod, only : noah33_struc
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
  type(ESMF_Field)       :: sfsmField
  real, pointer          :: sfsm(:)
  integer                :: t
  integer                :: status

  n = 1

  call ESMF_StateGet(Obj_Func,"LPRM AMSRE Surface Soil Moisture",sfsmField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sfsmField,localDE=0,farrayPtr=sfsm,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     sfsm(t) = noah33_struc(n)%noah(t)%smc(1)
  enddo

end subroutine noah33_getpeobspred_LPRM_AMSREsmObs



