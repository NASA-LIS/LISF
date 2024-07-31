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
! !ROUTINE: noah33_setsnodepvars
! \label{noah33_setsnodepvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!
! !INTERFACE:
subroutine noah33_setsnodepvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain,LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use noah33_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_State)       :: LSM_State
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: npts(LIS_rc%ngrid(n))
  integer                :: t,tid,gid
  integer                :: status 
!  real                   :: temp(LIS_rc%lnc(n), LIS_rc%lnr(n))

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  LIS_snow_struc(n)%sneqv = 0.0
  LIS_snow_struc(n)%snowdepth = 0.0
  npts = 0 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     !transform t to the patch
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
     noah33_struc(n)%noah(t)%sneqv = swe(t)
     
     LIS_snow_struc(n)%sneqv(tid)   = LIS_snow_struc(n)%sneqv(tid)+swe(t)
  enddo
  npts = 0 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     !transform t to the patch
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
     gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     noah33_struc(n)%noah(t)%snowh = snod(t)
     npts(gid) = npts(gid) + 1
     LIS_snow_struc(n)%snowdepth(gid) = LIS_snow_struc(n)%snowdepth(gid)+snod(t)

  enddo

!  temp = -9999.0
  do t=1,LIS_rc%ngrid(n)
     if(npts(t).gt.0) then 
        LIS_snow_struc(n)%snowdepth(t) = LIS_snow_struc(n)%snowdepth(t)/npts(t)
     else
        LIS_snow_struc(n)%snowdepth(t) = 0.0
     endif

  enddo
  
end subroutine noah33_setsnodepvars

