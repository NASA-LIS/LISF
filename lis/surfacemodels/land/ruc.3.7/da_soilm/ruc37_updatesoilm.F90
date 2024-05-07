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
! !ROUTINE: RUC37_updatesoilm
!  \label{RUC37_updatesoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine RUC37_updatesoilm(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use RUC37_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to RUC's
!  model space. 
! 
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
  type(ESMF_Field)       :: sm1IncrField
  type(ESMF_Field)       :: sm2IncrField
  type(ESMF_Field)       :: sm3IncrField
  type(ESMF_Field)       :: sm4IncrField
  type(ESMF_Field)       :: sm5IncrField
  type(ESMF_Field)       :: sm6IncrField
  type(ESMF_Field)       :: sm7IncrField
  type(ESMF_Field)       :: sm8IncrField
  type(ESMF_Field)       :: sm9IncrField

  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: soilm5(:)
  real, pointer          :: soilm6(:)
  real, pointer          :: soilm7(:)
  real, pointer          :: soilm8(:)
  real, pointer          :: soilm9(:)
  real, pointer          :: soilm1Incr(:)
  real, pointer          :: soilm2Incr(:)
  real, pointer          :: soilm3Incr(:)
  real, pointer          :: soilm4Incr(:)
  real, pointer          :: soilm5Incr(:)
  real, pointer          :: soilm6Incr(:)
  real, pointer          :: soilm7Incr(:)
  real, pointer          :: soilm8Incr(:)
  real, pointer          :: soilm9Incr(:)

  integer                :: t
  integer                :: status

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

  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 1",sm1IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 2",sm2IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 3",sm3IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 4",sm4IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 5",sm5IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 6",sm6IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 7",sm7IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 8",sm8IncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Level 9",sm9IncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
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

  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilm1Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2IncrField,localDE=0,farrayPtr=soilm2Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3IncrField,localDE=0,farrayPtr=soilm3Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4IncrField,localDE=0,farrayPtr=soilm4Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4IncrField,localDE=0,farrayPtr=soilm4Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm5IncrField,localDE=0,farrayPtr=soilm5Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm6IncrField,localDE=0,farrayPtr=soilm6Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm7IncrField,localDE=0,farrayPtr=soilm7Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm8IncrField,localDE=0,farrayPtr=soilm8Incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm9IncrField,localDE=0,farrayPtr=soilm9Incr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     soilm1(t) = soilm1(t) + soilm1Incr(t)
     soilm2(t) = soilm2(t) + soilm2Incr(t)
     soilm3(t) = soilm3(t) + soilm3Incr(t)
     soilm4(t) = soilm4(t) + soilm4Incr(t)
     soilm5(t) = soilm5(t) + soilm5Incr(t)
     soilm6(t) = soilm6(t) + soilm6Incr(t)
     soilm7(t) = soilm7(t) + soilm7Incr(t)
     soilm8(t) = soilm8(t) + soilm8Incr(t)
     soilm9(t) = soilm9(t) + soilm9Incr(t)
  enddo
end subroutine RUC37_updatesoilm

