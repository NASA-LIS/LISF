!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getCROCUSexport
! \label{noahmp401_getCROCUSexport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_getCROCUSexport(n, LSM2SUBLSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM2SUBLSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP



  type(ESMF_Field)   :: gtField
  real, pointer      :: gt(:)
  integer            :: t
  integer            :: status

  call ESMF_StateGet(LSM2SUBLSM_State,"Ground temperature",gtField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gtField,localDE=0,farrayPtr=gt,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gt(t) = NOAHMP401_struc(n)%noahmp401(t)%tgb
  enddo


end subroutine noahmp401_getCROCUSexport


