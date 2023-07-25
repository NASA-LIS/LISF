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
! !ROUTINE: noahmpnew_descale_tws
! \label{noahmpnew_descale_tws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!
! !INTERFACE:
subroutine noahmpnew_descale_tws(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmpnew_lsmMod
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
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: MAX_THRESHOLD , MIN_THRESHOLD

#if 0 
  ! Natt
  ! Descale TWS states
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in noahmpnew_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in noahmpnew_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in noahmpnew_descale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in noahmpnew_descale_tws')
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for SWE in noahmpnew_descale_tws')
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Snowdepth in noahmpnew_descale_tws')
  
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in noahmpnew_descale_tws')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in noahmpnew_descale_tws')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in noahmpnew_descale_tws')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm4 in noahmpnew_descale_tws')
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for SWE in noahmpnew_descale_tws')
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Snowdepth in noahmpnew_descale_tws')
  
  ! Natt
  ! Descale TWS state increment
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in noahmpnew_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 2",sm2IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in noahmpnew_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 3",sm3IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in noahmpnew_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 4",sm4IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 4 failed in noahmpnew_descale_tws")
  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmpnew_descale_tws")
  call ESMF_FieldGet(sm2IncrField,localDE=0,farrayPtr=soilmIncr2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmpnew_descale_tws")
  call ESMF_FieldGet(sm3IncrField,localDE=0,farrayPtr=soilmIncr3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmpnew_descale_tws")
  call ESMF_FieldGet(sm4IncrField,localDE=0,farrayPtr=soilmIncr4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmpnew_descale_tws")
  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)
  
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     !SOILTYP = NoahMPnew_struc(n)%noahmpnew(t)%soiltype
     !MAX_THRESHOLD = NoahMPnew_struc(n)%noahmpnew(t)%param%SMCMAX(1)
     !MIN_THRESHOLD = NoahMPnew_struc(n)%noahmpnew(t)%param%SMCWLT(1)
     soilm1(t) = soilm1(t) / (NoahMPnew_struc(n)%sldpth(1)*1000.0) ! mm -> m3m-3
     soilm2(t) = soilm2(t) / (NoahMPnew_struc(n)%sldpth(2)*1000.0) ! mm -> m3m-3
     soilm3(t) = soilm3(t) / (NoahMPnew_struc(n)%sldpth(3)*1000.0) ! mm -> m3m-3
     soilm4(t) = soilm4(t) / (NoahMPnew_struc(n)%sldpth(4)*1000.0) ! mm -> m3m-3
     !if (soilm1(t) .lt. MIN_THRESHOLD) soilm1(t) = MIN_THRESHOLD
     !if (soilm2(t) .lt. MIN_THRESHOLD) soilm2(t) = MIN_THRESHOLD
     !if (soilm3(t) .lt. MIN_THRESHOLD) soilm3(t) = MIN_THRESHOLD
     !if (soilm4(t) .lt. MIN_THRESHOLD) soilm4(t) = MIN_THRESHOLD
     !if (soilm1(t) .gt. MAX_THRESHOLD) soilm1(t) = MAX_THRESHOLD
     !if (soilm2(t) .gt. MAX_THRESHOLD) soilm2(t) = MAX_THRESHOLD
     !if (soilm3(t) .gt. MAX_THRESHOLD) soilm3(t) = MAX_THRESHOLD
     !if (soilm4(t) .gt. MAX_THRESHOLD) soilm4(t) = MAX_THRESHOLD
     swe(t)    = swe(t) / 1000.0  ! mm -> m
     snod(t)   = snod(t) / 1000.0 ! mm -> m
     
     soilmIncr1(t) = soilmIncr1(t) / (NoahMPnew_struc(n)%sldpth(1)*1000.0) ! mm -> m3m-3
     soilmIncr2(t) = soilmIncr2(t) / (NoahMPnew_struc(n)%sldpth(2)*1000.0) ! mm -> m3m-3
     soilmIncr3(t) = soilmIncr3(t) / (NoahMPnew_struc(n)%sldpth(3)*1000.0) ! mm -> m3m-3
     soilmIncr4(t) = soilmIncr4(t) / (NoahMPnew_struc(n)%sldpth(4)*1000.0) ! mm -> m3m-3
     sweincr(t)    = sweincr(t) / 1000.0  ! mm -> m
     snodincr(t)   = snodincr(t) / 1000.0 ! mm -> m

  enddo
#endif
  
end subroutine noahmpnew_descale_tws
