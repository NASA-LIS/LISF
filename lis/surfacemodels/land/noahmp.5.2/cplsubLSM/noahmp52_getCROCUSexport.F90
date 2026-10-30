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
! !ROUTINE: noahmp52_getCROCUSexport
! \label{noahmp52_getCROCUSexport}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified to work with refactored Noah-MP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_getCROCUSexport(n, LSM2SUBLSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp52_lsmMod

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
  type(ESMF_Field)   :: XWGIField
  type(ESMF_Field)   :: XWGField
  real, pointer      :: gt(:)
  real, pointer      :: XWGI(:)
  real, pointer      :: XWG(:)
  integer            :: t
  integer            :: status

  call ESMF_StateGet(LSM2SUBLSM_State,"Ground temperature",gtField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric liquid water content",&
       XWGField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric frozen water content",&
       XWGIField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gtField,localDE=0,farrayPtr=gt,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(XWGField,localDE=0,farrayPtr=XWG,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(XWGIField,localDE=0,farrayPtr=XWGI,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gt(t) = NoahMP52_struc(n)%noahmp52(t)%tslb(1)
     ! volumetric frozen soil moisture [m3/m3]
     XWGI(t) = NoahMP52_struc(n)%noahmp52(t)%smc(1) - &
          NoahMP52_struc(n)%noahmp52(t)%sh2o(1)
     ! volumetric liquid soil moisture [m3/m3]
     XWG(t)  = NoahMP52_struc(n)%noahmp52(t)%sh2o(1)
  enddo

end subroutine noahmp52_getCROCUSexport
