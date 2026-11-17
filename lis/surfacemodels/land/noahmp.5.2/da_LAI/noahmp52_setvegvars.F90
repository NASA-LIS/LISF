!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp52_setvegvars
!  \label{noahmp52_setvegvars}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp52_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  This routine assigns the LAI prognostic variables to NoahMP's
!  model space.
!
!EOP
  real, parameter        :: MIN_THRESHOLD = 0.02
  type(ESMF_Field)       :: laiField
  integer                :: t
  integer                :: status
  real, pointer          :: lai(:)
  real                   :: lfmass

  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(NoahMP52_struc(n)%noahmp52(t)%param%sla.ne.0) then
        NoahMP52_struc(n)%noahmp52(t)%lai = lai(t)
        lfmass = lai(t)*1000.0/(NoahMP52_struc(n)%noahmp52(t)%param%sla)
        NoahMP52_struc(n)%noahmp52(t)%lfmass = lfmass
     endif
  enddo

end subroutine noahmp52_setvegvars

