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
! !ROUTINE: snowmodel_setLSMimport
! \label{snowmodel_setLSMimport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!  2 Dec 2020: Mahdi Navari; Edited to add soil volumetric liquid and frozen water content
!  2 Aug 2021: Kristi Arsenault; Edited to support SnowModel implementation

! !INTERFACE:
subroutine snowmodel_setLSMimport(n, LSM2SubLSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use snowmodel_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM2SubLSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP

#if 0
  type(ESMF_Field)   :: gtField
  real, pointer      :: gt(:)
  type(ESMF_Field)   :: XWGIField
  real, pointer      :: XWGI(:)
  type(ESMF_Field)   :: XWGField
  real, pointer      :: XWG(:)
  integer            :: t
  integer            :: status

  call ESMF_StateGet(LSM2SubLSM_State,"Ground temperature",gtField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gtField,localDE=0,farrayPtr=gt,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric frozen water content",XWGIField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(XWGIField,localDE=0,farrayPtr=XWGI,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric liquid water content",XWGField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(XWGField,localDE=0,farrayPtr=XWG,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snowmodel_struc(n)%crocus81(t)%TG = gt(t)
     snowmodel_struc(n)%crocus81(t)%XWGI = XWGI(t)
     snowmodel_struc(n)%crocus81(t)%XWG  = XWG(t)
  enddo

#endif

end subroutine snowmodel_setLSMimport


