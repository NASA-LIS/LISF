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
! !ROUTINE: noahmp50_setCROCUSimport
! \label{noahmp50_setCROCUSimport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
! May 2023: Cenlin He; Modified to work with refactored Noah-MP v5 and later
!
! !INTERFACE:
subroutine noahmp50_setCROCUSimport(n, SubLSM2LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp50_lsmMod

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
  real               :: dsneqv,dsnowh
  integer            :: t
  integer            :: status

  call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     dsneqv = swe(t) - NoahMP50_struc(n)%noahmp50(t)%sneqv   !in mm
     dsnowh = snwd(t) - NoahMP50_struc(n)%noahmp50(t)%snowh  !in m

     ! update
     call noahmp50_snow_update(n, t, dsneqv, dsnowh)

  enddo

end subroutine noahmp50_setCROCUSimport


