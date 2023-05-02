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
! !ROUTINE: noahmp401_qctws
! \label{noahmp401_qctws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; Created for Noah-MP4.0.1
!
! !INTERFACE:
subroutine noahmp401_qctws(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_verify
  use noahmp401_lsmMod
  !use module_sf_noahmplsm_401
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE

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
  integer                :: t,gid
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real                   :: smmax
  real                   :: smmin

  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field

!Wanshu
  type(ESMF_Field)       :: gwField, sweField, snodField
  real, pointer          :: gws(:)
  real                   :: gwsmax, gwsmin
  real                   :: MIN_THRESHOLD,MAX_threshold,sm_threshold
  integer                :: SOILTYP

  real, pointer          :: swe(:)
  real, pointer          :: snod(:)

  real                   :: swemax,snodmax
  real                   :: swemin,snodmin

  real                   :: sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
!-------

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 2 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 2 failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 3 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 3 failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 4 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 4 failed in noahmp401_qctws")

  !Wanshu
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for gw in noahmp401_qctws')

  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for gw in noahmp401_qctws')

  call ESMF_AttributeGet(gwField,"Max Value",gwsmax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in noahmp401_qctws")

  call ESMF_AttributeGet(gwField,"Min Value",gwsmin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)


  call ESMF_AttributeGet(sweField,"Max Value",swemax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sweField,"Min Value",swemin,rc=status)
  call LIS_verify(status)

  !-------


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  !Bailing Li: max min soil moisture should be retrieved based on soil type
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)
     MIN_THRESHOLD = SMCWLT_TABLE(SOILTYP)
     sm_threshold = MAX_THRESHOLD - 0.02


     if(soilm1(t).gt.sm_threshold) then
        soilm1(t) = sm_threshold
     endif

     if(soilm1(t).lt.MIN_THRESHOLD) then
        soilm1(t) = MIN_THRESHOLD
     endif

     if(soilm2(t).gt.sm_threshold) then
        soilm2(t) = sm_threshold
     endif
     if(soilm2(t).lt.MIN_THRESHOLD) then
        soilm2(t) = MIN_THRESHOLD
     endif

     if(soilm3(t).gt.sm_threshold) then
        soilm3(t) = sm_threshold
     endif
     if(soilm3(t).lt.MIN_THRESHOLD) then
        soilm3(t) = MIN_THRESHOLD
     endif

     if(soilm4(t).gt.sm_threshold) then
        soilm4(t) = sm_threshold
     endif
     if(soilm4(t).lt.MIN_THRESHOLD) then
        soilm4(t) = MIN_THRESHOLD
     endif
     !Wanshu
     if(gws(t).gt.gwsmax) then
        gws(t) = gwsmax
     endif
     if(gws(t).lt.gwsmin) then
        gws(t) = gwsmin
     endif

     if(swe(t).lt.swemin) then
        swe(t) = swemin
     endif
     !------
  enddo

end subroutine noahmp401_qctws
