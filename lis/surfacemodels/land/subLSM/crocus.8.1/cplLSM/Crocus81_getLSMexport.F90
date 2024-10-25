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
! !ROUTINE: Crocus81_getLSMexport
! \label{Crocus81_getLSMexport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!  9 Dec 2020: Mahdi Navari; edited to take into account the Crocus slope correction 
!
! !INTERFACE:
subroutine Crocus81_getLSMexport(n, SubLSM2LSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use Crocus81_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: SubLSM2LSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP

  type(ESMF_Field)   :: snwdField, sweField
  real, pointer      :: swe(:), snwd(:)
  integer            :: t
  integer            :: status
  real               :: tmp_SLOPE

 
  call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp_SLOPE = CROCUS81_struc(n)%crocus81(t)%SLOPE
     swe(t) = CROCUS81_struc(n)%crocus81(t)%SWE_1D / COS(tmp_SLOPE)
     snwd(t) = CROCUS81_struc(n)%crocus81(t)%SD_1D / COS(tmp_SLOPE)
  enddo

end subroutine Crocus81_getLSMexport


