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
! !ROUTINE: noahmp36_setvegvars
! \label{noahmp36_setvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp36_lsmMod
  use NOAHMP_VEG_PARAMETERS_36
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to noah's
!  model space. The state vector consists of veg
! 
!EOP

  type(ESMF_Field)       :: laiField,lfmassField

  integer                :: t
  integer                :: status
  real, pointer          :: lai(:)
  real                   :: lfmass
 
  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

!  call ESMF_StateGet(LSM_State,"LeafMass",lfmassField,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGet(lfmassField,localDE=0,farrayPtr=lfmass,rc=status)
!  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     XLAI    = max(LFMASS*LAPM,laimin)

     if(sla(NOAHMP36_struc(n)%noahmp36(t)%vegetype).ne.0) then 
        NOAHMP36_struc(n)%noahmp36(t)%lai = lai(t)
        lfmass = lai(t)/(sla(NOAHMP36_struc(n)%noahmp36(t)%vegetype)/1000.0)
        NOAHMP36_struc(n)%noahmp36(t)%lfmass = lfmass
     endif
  enddo
  
end subroutine noahmp36_setvegvars


