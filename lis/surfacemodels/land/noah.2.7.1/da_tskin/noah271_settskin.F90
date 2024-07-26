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
! !ROUTINE: noah271_settskin
!  \label{noah271_settskin}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_settskin(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use noah271_lsmMod
  use LIS_logMod,  only  : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the tskin prognostic variables to noah's
!  model space. 
! 
!EOP

!  type(ESMF_Field)       :: avgsurftField
  type(ESMF_Field)       :: soiltField
!  real, allocatable          :: avgsurft(:)
  real, allocatable          :: soilt(:)

  integer                :: t
  integer                :: status

!  call ESMF_StateGet(LSM_State,"Skin Temperature",avgsurftField,rc=status)
!  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Soil Temperature",soiltField,rc=status)
  call LIS_verify(status)

!  call ESMF_FieldGet(avgsurftField,localDE=0,farrayPtr=avgsurft,rc=status)
!  call LIS_verify(status)

  call ESMF_FieldGet(soiltField,localDE=0,farrayPtr=soilt,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
!     noah271_struc(n)%noah(t)%t1 = avgsurft(t)
     noah271_struc(n)%noah(t)%stc(1) = soilt(t)
!     if(t==1) write(LIS_logunit,*) 'set ',soilt(t)
  enddo
end subroutine noah271_settskin

