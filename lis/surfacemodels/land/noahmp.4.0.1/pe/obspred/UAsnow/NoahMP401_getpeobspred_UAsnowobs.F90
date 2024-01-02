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
! !ROUTINE: NoahMP401_getpeobspred_UAsnowobs
!  \label{NoahMP401_getpeobspred_UAsnowobs}
!
! !REVISION HISTORY:
! 02 May 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine NoahMP401_getpeobspred_UAsnowobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use NoahMP401_lsmMod, only : NoahMP401_struc
  use LIS_logMod,       only : LIS_verify, LIS_logunit

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
  integer                :: i
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
     snod(t) = NoahMP401_struc(n)%noahmp401(t)%snowh*1000.0 !mm
     swe(t) = NoahMP401_struc(n)%noahmp401(t)%sneqv
  enddo  


end subroutine NoahMP401_getpeobspred_UAsnowobs



