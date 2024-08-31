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
! !ROUTINE: noahmp401_gettws
! \label{noahmp401_gettws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 18 Aug 2017: Wanshu Nie; Add groundwater component
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1

! !INTERFACE:
subroutine noahmp401_gettws(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture and groundwater related state prognostic variables for
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
  type(ESMF_Field)       :: gwField
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: gws(:)
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  character*100          :: lsm_state_objs(4)

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in noahmp401_gettws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in noahmp401_gettws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in noahmp401_gettws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in noahmp401_gettws')
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for gw in noahmp401_gettws')
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for SWE in noahmp401_gettws')

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in noahmp401_gettws')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in noahmp401_gettws')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in noahmp401_gettws')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm4 in noahmp401_gettws')
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for gw in noahmp401_gettws')
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for SWE in noahmp401_gettws')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index) !to mm
     soilm1(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(1)
     soilm2(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(2)
     soilm3(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(3)
     soilm4(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(4)
     gws(t)    = NOAHMP401_struc(n)%noahmp401(t)%wa
     swe(t) = noahmp401_struc(n)%noahmp401(t)%sneqv
  enddo

end subroutine noahmp401_gettws
