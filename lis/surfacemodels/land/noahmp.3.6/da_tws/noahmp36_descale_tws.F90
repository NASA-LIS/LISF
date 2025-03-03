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
! !ROUTINE: noahmp36_descale_tws
! \label{noahmp36_descale_tws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_descale_tws(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp36_lsmMod
  use LIS_constantsMod,  only : LIS_CONST_RHOFW ! Natt

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Descales twsoisture related state prognostic variables for
!  data assimilation
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

!Wanshu

  type(ESMF_Field)     :: gwField
  real, pointer        :: gws(:)
  type(ESMF_Field)     :: gwIncrField
  real, pointer        :: gwsIncr(:)
  integer              :: t
  integer              :: status

! Natt
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  type(ESMF_Field)       :: sm1IncrField
  type(ESMF_Field)       :: sm2IncrField
  type(ESMF_Field)       :: sm3IncrField
  type(ESMF_Field)       :: sm4IncrField
  type(ESMF_Field)       :: sweIncrField
  type(ESMF_Field)       :: snodIncrField
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real, pointer          :: soilmIncr1(:)
  real, pointer          :: soilmIncr2(:)
  real, pointer          :: soilmIncr3(:)
  real, pointer          :: soilmIncr4(:)
  real, pointer          :: sweincr(:)
  real, pointer          :: snodincr(:)

#if 0
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_settws")

  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp36_settws")

  call ESMF_StateGet(LSM_Incr_State, "Groundwater Storage",gwIncrField, &
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_descale")

  call ESMF_FieldGet(gwIncrField, localDE=0, farrayPtr=gwsIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_descale")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gws(t)   = gws(t)*10.0
     gwsIncr(t) = gwsIncr(t)*10.0
  enddo
#endif

  ! Natt
  ! Descale TWS states
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for sm1 in noahmp36_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for sm2 in noahmp36_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for sm3 in noahmp36_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for sm4 in noahmp36_descale_tws')
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for SWE in noahmp36_descale_tws')
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status, &
       'ESMF_StateGet failed for Snowdepth in noahmp36_descale_tws')

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for sm1 in noahmp36_descale_tws')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for sm2 in noahmp36_descale_tws')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for sm3 in noahmp36_descale_tws')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for sm4 in noahmp36_descale_tws')
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for SWE in noahmp36_descale_tws')
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for Snowdepth in noahmp36_descale_tws')

  ! Natt
  ! Descale TWS state increment
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField, &
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in noahmp36_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 2",sm2IncrField, &
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in noahmp36_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 3",sm3IncrField, &
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in noahmp36_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 4",sm4IncrField, &
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 4 failed in noahmp36_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp36_descale_tws")
  call ESMF_FieldGet(sm2IncrField,localDE=0,farrayPtr=soilmIncr2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp36_descale_tws")
  call ESMF_FieldGet(sm3IncrField,localDE=0,farrayPtr=soilmIncr3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp36_descale_tws")
  call ESMF_FieldGet(sm4IncrField,localDE=0,farrayPtr=soilmIncr4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp36_descale_tws")
  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) / &
          (NOAHMP36_struc(n)%sldpth(1)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilm2(t) = soilm2(t) / &
          (NOAHMP36_struc(n)%sldpth(2)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilm3(t) = soilm3(t) / &
          (NOAHMP36_struc(n)%sldpth(3)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilm4(t) = soilm4(t) / &
          (NOAHMP36_struc(n)%sldpth(4)*LIS_CONST_RHOFW) ! mm -> m3m-3
     swe(t)    = swe(t) / 1000.0  ! mm -> m
     snod(t)   = snod(t) / 1000.0 ! mm -> m

     soilmIncr1(t) = soilmIncr1(t) / &
          (NOAHMP36_struc(n)%sldpth(1)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilmIncr2(t) = soilmIncr2(t) / &
          (NOAHMP36_struc(n)%sldpth(2)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilmIncr3(t) = soilmIncr3(t) / &
          (NOAHMP36_struc(n)%sldpth(3)*LIS_CONST_RHOFW) ! mm -> m3m-3
     soilmIncr4(t) = soilmIncr4(t) / &
          (NOAHMP36_struc(n)%sldpth(4)*LIS_CONST_RHOFW) ! mm -> m3m-3
     sweincr(t)    = sweincr(t) / 1000.0  ! mm -> m
     snodincr(t)   = snodincr(t) / 1000.0 ! mm -> m

  enddo


end subroutine noahmp36_descale_tws

