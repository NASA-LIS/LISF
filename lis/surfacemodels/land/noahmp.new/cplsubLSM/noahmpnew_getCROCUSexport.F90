!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmpnew_getCROCUSexport
! \label{noahmpnew_getCROCUSexport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
! 17 Nov 2020: Mahdi Navari; In analogous to ISBA-Crocus %tgb replaced with %tslb
!  2 Dec 2020: Mahdi Navari; Edited to add soil volumetric liquid and frozen water content
! May 2023: Cenlin He; modified to work with refactored Noah-MP v5 and later
!
! !INTERFACE:
subroutine noahmpnew_getCROCUSexport(n, LSM2SUBLSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmpnew_lsmMod

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
  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric liquid water content",XWGField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric frozen water content",XWGIField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gtField,localDE=0,farrayPtr=gt,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(XWGField,localDE=0,farrayPtr=XWG,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(XWGIField,localDE=0,farrayPtr=XWGI,rc=status)
  call LIS_verify(status)



  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gt(t) = NoahMPnew_struc(n)%noahmpnew(t)%tslb(1)
     !gt(t) = NoahMPnew_struc(n)%noahmpnew(t)%tgb
     XWGI(t) = NoahMPnew_struc(n)%noahmpnew(t)%smc(1) - NoahMPnew_struc(n)%noahmpnew(t)%sh2o(1) ! volumetric frozen soil moisture [m3/m3]
     XWG(t)  = NoahMPnew_struc(n)%noahmpnew(t)%sh2o(1) ! volumetric liquid soil moisture [m3/m3]
  enddo


end subroutine noahmpnew_getCROCUSexport
