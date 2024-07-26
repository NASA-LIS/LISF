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
! !ROUTINE: clm2_descale_soilm
! \label{clm2_descale_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_descale_soilm(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
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
!  type(ESMF_Field)       :: sm1Arr
!  type(ESMF_Field)       :: sm2Arr
!  type(ESMF_Field)       :: sm3Arr
!  type(ESMF_Field)       :: sm4Arr
!  integer                :: t
!  integer                :: status
!  real, allocatable          :: soilm1(:)
!  real, allocatable          :: soilm2(:)
!  real, allocatable          :: soilm3(:)
!  real, allocatable          :: soilm4(:)
 
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Arr,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Arr,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Arr,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Arr,rc=status)
!  call LIS_verify(status)
 
!  call ESMF_FieldGet(sm1Arr,soilm1,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm2Arr,soilm2,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm3Arr,soilm3,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm4Arr,soilm4,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)


!  do t=1,LIS_rc%ntiles(n)
!     soilm1(t) = soilm1(t)/1.0
!     soilm2(t) = soilm2(t)/1.0
!     soilm3(t) = soilm3(t)/1.0
!     soilm4(t) = soilm4(t)/1.0
!  enddo

end subroutine clm2_descale_soilm

