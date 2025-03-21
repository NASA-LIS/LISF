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
! !ROUTINE: Noah39_updatesoilm
!  \label{Noah39_updatesoilm}
!
! !REVISION HISTORY:
! 21Oct2018: Mahdi Navari; Sujay Kumar ; Initial Specification
! 20 Mar 2025: Eric Kemp; Changed to use top soil layer only.  Based on
!   NoahMP401 code from Sujay Kumar and Wanshu Nie used for, e.g.,
!   HydroGlobe.
!
! !INTERFACE:
subroutine Noah39_updatesoilm(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_verify

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State), intent(in) :: LSM_State
  type(ESMF_State), intent(in) :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space.
!
!EOP

  type(ESMF_Field)       :: sm1Field
  real, pointer          :: soilm1(:)
  type(ESMF_Field)       :: sm1IncrField
  real, pointer          :: soilmIncr1(:)
  integer                :: t
  integer                :: status

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in Noah39_updatesoilm")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in Noah39_updatesoilm")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in Noah39_updatesoilm")
  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in Noah39_updatesoilm")
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo
end subroutine Noah39_updatesoilm
