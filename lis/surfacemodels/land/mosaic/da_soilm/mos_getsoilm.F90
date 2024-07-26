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
! !ROUTINE: mos_getsoilm
! \label{mos_getsoilm}
!
! !REVISION HISTORY:
! 06 Oct 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine mos_getsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use mos_lsmMod

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

  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = mos_struc(n)%mos(t)%sowet(1)
     soilm2(t) = mos_struc(n)%mos(t)%sowet(2)
     soilm3(t) = mos_struc(n)%mos(t)%sowet(3)
  enddo

end subroutine mos_getsoilm

