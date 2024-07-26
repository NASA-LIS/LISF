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
! !ROUTINE: clm2_getsoilm
! \label{clm2_getsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_getsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod
  use clm2_varcon,    only : denh2o

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
  type(ESMF_Field)       :: sm10Field

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
  real, pointer          :: soilm10(:)
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 5",sm5Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 6",sm6Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 7",sm7Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 8",sm8Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 9",sm9Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 10",sm10Field,rc=status)
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
  call ESMF_FieldGet(sm10Field,localDE=0,farrayPtr=soilm10,rc=status)
  call LIS_verify(status)



  do t=1,LIS_rc%ntiles(n)
     soilm1(t) = clm2_struc(n)%clm(t)%h2osoi_liq(1)/&
          (clm2_struc(n)%clm(n)%dz(1)*denh2o)
     soilm2(t) = clm2_struc(n)%clm(t)%h2osoi_liq(2)/&
          (clm2_struc(n)%clm(n)%dz(2)*denh2o)
     soilm3(t) = clm2_struc(n)%clm(t)%h2osoi_liq(3)/&
          (clm2_struc(n)%clm(n)%dz(3)*denh2o)
     soilm4(t) = clm2_struc(n)%clm(t)%h2osoi_liq(4)/&
          (clm2_struc(n)%clm(n)%dz(4)*denh2o)
     soilm5(t) = clm2_struc(n)%clm(t)%h2osoi_liq(5)/&
          (clm2_struc(n)%clm(n)%dz(5)*denh2o)
     soilm6(t) = clm2_struc(n)%clm(t)%h2osoi_liq(6)/&
          (clm2_struc(n)%clm(n)%dz(6)*denh2o)
     soilm7(t) = clm2_struc(n)%clm(t)%h2osoi_liq(7)/&
          (clm2_struc(n)%clm(n)%dz(7)*denh2o)
     soilm8(t) = clm2_struc(n)%clm(t)%h2osoi_liq(8)/&
          (clm2_struc(n)%clm(n)%dz(8)*denh2o)
     soilm9(t) = clm2_struc(n)%clm(t)%h2osoi_liq(9)/&
          (clm2_struc(n)%clm(n)%dz(9)*denh2o)
     soilm10(t) = clm2_struc(n)%clm(t)%h2osoi_liq(10)/&
          (clm2_struc(n)%clm(n)%dz(10)*denh2o)
  enddo

end subroutine clm2_getsoilm

