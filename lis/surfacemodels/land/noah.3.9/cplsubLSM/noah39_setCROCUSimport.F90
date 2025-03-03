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
! !ROUTINE: noah39_setCROCUSimport
! \label{noah39_setCROCUSimport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!  7 Jun 2021: Mahdi Navari; Modified for Noah39
!
! !INTERFACE:
subroutine noah39_setCROCUSimport(n, SubLSM2LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noah39_lsmMod

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
     dsneqv = swe(t) - noah39_struc(n)%noah(t)%sneqv*1000   !swe from Crocus is in mm, the sneqv from Noa39 is in m. 
     dsnowh = snwd(t) - noah39_struc(n)%noah(t)%snowh  !in m
     ! update
     if (dsneqv.gt.0.0 .and. dsnowh.gt.0.0 .and. ((dsneqv/1000)/dsnowh).lt.1.0 )then
     call noah39_snow_update(n, t, dsneqv/1000, dsnowh)
     ! or we can update here
     !noah39_struc(n)%noah(t)%sneqv = dsneqv
     !noah39_struc(n)%noah(t)%snowh = dsnowh
     endif
  enddo

end subroutine noah39_setCROCUSimport


