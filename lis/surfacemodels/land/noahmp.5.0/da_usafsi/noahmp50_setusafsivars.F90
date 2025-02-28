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
! !ROUTINE: noahmp50_setusafsivars
! \label{noahmp50_setusafsivars}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:
subroutine noahmp50_setusafsivars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_logMod, only : LIS_verify, LIS_endrun
  use LIS_snowMod, only : LIS_snow_struc
  use NoahMP50_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  This routine assigns the snow progognostic variables to noah's
!  model space. The state vector consists of total SWE and snow depth.
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update.
!
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh
  integer                :: t
  integer                :: status
  integer                :: ncount(LIS_rc%ngrid(n))
  integer                :: tid, gid
  
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     dsneqv = swe(t) - Noahmp50_struc(n)%noahmp50(t)%sneqv   !in mm
     dsnowh = snod(t) - Noahmp50_struc(n)%noahmp50(t)%snowh  !in m

     ! update
     call noahmp50_usafsi_update(n, t, dsneqv, dsnowh)

  enddo

  if (LIS_rc%snowsrc(n) .gt. 0) then

     ncount = 0 ! Number of tiles per grid id (over land)
     LIS_snow_struc(n)%snowdepth = 0 ! At grid points
     LIS_snow_struc(n)%sneqv = 0     ! At tiles

     ! Collect SWE at tiles
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid) = LIS_snow_struc(n)%sneqv(tid) + &
             Noahmp50_struc(n)%noahmp50(t)%sneqv
     end do

     ! Collect mean snow depth at grid points
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        LIS_snow_struc(n)%snowdepth(gid) = &
             LIS_snow_struc(n)%snowdepth(gid) + &
             Noahmp50_struc(n)%noahmp50(t)%snowh
        ncount(gid) = ncount(gid) + 1
     end do
     do t = 1, LIS_rc%ngrid(n)
        if (ncount(t).gt.0) then
           LIS_snow_struc(n)%snowdepth(t) = &
                LIS_snow_struc(n)%snowdepth(t) / ncount(t)
        else
           LIS_snow_struc(n)%snowdepth(t) = 0.0
        endif
     end do
  end if

end subroutine noahmp50_setusafsivars
