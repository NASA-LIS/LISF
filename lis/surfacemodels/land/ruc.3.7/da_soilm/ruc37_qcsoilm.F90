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
! !ROUTINE: RUC37_qcsoilm
! \label{RUC37_qcsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine RUC37_qcsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use RUC37_lsmMod

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
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: sm5Field
  type(ESMF_Field)       :: sm6Field
  type(ESMF_Field)       :: sm7Field
  type(ESMF_Field)       :: sm8Field
  type(ESMF_Field)       :: sm9Field


  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: soilm5(:)
  real, pointer          :: soilm6(:)
  real, pointer          :: soilm7(:)
  real, pointer          :: soilm8(:)
  real, pointer          :: soilm9(:)

  real                   :: smmax1,smmax2,smmax3,smmax4,smmax5
  real                   :: smmax6,smmax7,smmax8,smmax9
  real                   :: smmin1,smmin2,smmin3,smmin4,smmin5
  real                   :: smmin6,smmin7,smmin8,smmin9
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 1",sm1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 2",sm2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 3",sm3Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 4",sm4Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 5",sm5Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 6",sm6Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 7",sm7Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 8",sm8Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Level 9",sm9Field,rc=status)
  call LIS_verify(status)

 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm5Field,localDE=0,farrayPtr=soilm5,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm6Field,localDE=0,farrayPtr=soilm6,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm7Field,localDE=0,farrayPtr=soilm7,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm8Field,localDE=0,farrayPtr=soilm8,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm9Field,localDE=0,farrayPtr=soilm9,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm2Field,"Max Value",smmax2,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm2Field,"Min Value",smmin2,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm3Field,"Max Value",smmax3,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm3Field,"Min Value",smmin3,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm4Field,"Max Value",smmax4,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm4Field,"Min Value",smmin4,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm5Field,"Max Value",smmax5,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm5Field,"Min Value",smmin5,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm6Field,"Max Value",smmax6,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm6Field,"Min Value",smmin6,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm7Field,"Max Value",smmax7,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm7Field,"Min Value",smmin7,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm8Field,"Max Value",smmax8,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm8Field,"Min Value",smmin8,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm9Field,"Max Value",smmax9,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm9Field,"Min Value",smmin9,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(soilm1(t).gt.smmax1) soilm1(t) = smmax1
     if(soilm2(t).gt.smmax2) soilm2(t) = smmax2
     if(soilm3(t).gt.smmax3) soilm3(t) = smmax3
     if(soilm4(t).gt.smmax4) soilm4(t) = smmax4
     if(soilm5(t).gt.smmax5) soilm5(t) = smmax5
     if(soilm6(t).gt.smmax6) soilm6(t) = smmax6
     if(soilm7(t).gt.smmax7) soilm7(t) = smmax7
     if(soilm8(t).gt.smmax8) soilm8(t) = smmax8
     if(soilm9(t).gt.smmax9) soilm9(t) = smmax9


     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
     if(soilm2(t).lt.smmin2) soilm2(t) = smmin2
     if(soilm3(t).lt.smmin3) soilm3(t) = smmin3
     if(soilm4(t).lt.smmin4) soilm4(t) = smmin4
     if(soilm5(t).lt.smmin5) soilm5(t) = smmin5
     if(soilm6(t).lt.smmin6) soilm6(t) = smmin6
     if(soilm7(t).lt.smmin7) soilm7(t) = smmin7
     if(soilm8(t).lt.smmin8) soilm8(t) = smmin8
     if(soilm9(t).lt.smmin9) soilm9(t) = smmin9

  enddo

end subroutine RUC37_qcsoilm

