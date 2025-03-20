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
! !ROUTINE: NoahMP401_updatesoilm
!  \label{NoahMP401_updatesoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 20 Mar 2025: Eric Kemp; Changed to use top soil layer only.  Based on
!   code from Sujay Kumar and Wanshu Nie used for, e.g., HydroGlobe.
!
! !INTERFACE:
subroutine NoahMP401_updatesoilm(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_verify
  !use NoahMP401_lsmMod

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
       "ESMF_StateGet: Soil Moisture Layer 1 failed in NoahMP401_updatesoilm")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in NoahMP401_updatesoilm")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in NoahMP401_updatesoilm")
  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in NoahMP401_updatesoilm")
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo
end subroutine NoahMP401_updatesoilm
