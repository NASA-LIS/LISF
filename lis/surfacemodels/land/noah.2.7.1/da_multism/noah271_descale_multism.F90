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
! !ROUTINE: noah271_descale_multism
! \label{noah271_descale_multism}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_descale_multism(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Descales soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
!  type(ESMF_Field)       :: sm1Field
!  type(ESMF_Field)       :: sm2Field
!  type(ESMF_Field)       :: sm3Field
!  type(ESMF_Field)       :: sm4Field
!  integer                :: t
!  integer                :: status
!  real, allocatable          :: soilm1(:)
!  real, allocatable          :: soilm2(:)
!  real, allocatable          :: soilm3(:)
!  real, allocatable          :: soilm4(:)
 
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
!  call LIS_verify(status)
 
!  call ESMF_FieldGet(sm1Field,soilm1,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm2Field,soilm2,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm3Field,soilm3,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGetData(sm4Field,soilm4,ESMF_DATA_REF,rc=status)
!  call LIS_verify(status)


!  do t=1,LIS_rc%ntiles(n)
!     soilm1(t) = soilm1(t)/1.0
!     soilm2(t) = soilm2(t)/1.0
!     soilm3(t) = soilm3(t)/1.0
!     soilm4(t) = soilm4(t)/1.0
!  enddo

end subroutine noah271_descale_multism

