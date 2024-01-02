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
! !ROUTINE: noah271_qcmultism
! \label{noah271_qcmultism}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_qcmultism(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
!  type(ESMF_Field)       :: sm2Field
!  type(ESMF_Field)       :: sm3Field
!  type(ESMF_Field)       :: sm4Field
  integer                :: t
  integer                :: status
  real, allocatable          :: soilm1(:)
!  real, allocatable          :: soilm2(:)
!  real, allocatable          :: soilm3(:)
!  real, allocatable          :: soilm4(:)
  real                   :: smmax1!,smmax2,smmax3,smmax4
  real                   :: smmin1!,smmin2,smmin3,smmin4
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
!  call LIS_verify(status)
 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
!  call ESMF_FieldGetData(sm2Field,soilm2,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm3Field,soilm3,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm4Field,soilm4,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm2Field,"Max Value",smmax2,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm2Field,"Min Value",smmin2,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm3Field,"Max Value",smmax3,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm3Field,"Min Value",smmin3,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm4Field,"Max Value",smmax4,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetAttribute(sm4Field,"Min Value",smmin4,rc=status)
!  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)

     if(soilm1(t).gt.smmax1) soilm1(t) = smmax1
!     if(soilm2(t).gt.smmax2) soilm2(t) = smmax2
!     if(soilm3(t).gt.smmax3) soilm3(t) = smmax3
!     if(soilm4(t).gt.smmax4) soilm4(t) = smmax4

     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
!     if(soilm2(t).lt.smmin2) soilm2(t) = smmin2
!     if(soilm3(t).lt.smmin3) soilm3(t) = smmin3
!     if(soilm4(t).lt.smmin4) soilm4(t) = smmin4
  enddo

end subroutine noah271_qcmultism

