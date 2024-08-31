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
! !ROUTINE: jules50_setsnodepvars
! \label{jules50_setsnodepvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 05 Nov 2018: Yeosang Yoon; Modified for Jules 5.0
! 10 Nov 2020: Eric Kemp; Added support for LIS_snow_struc
!
! !INTERFACE:
subroutine jules50_setsnodepvars(n, LSM_State)
! !USES:
  use ESMF
  use jules50_lsmMod
  use LIS_coreMod, only : LIS_rc, LIS_domain,LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_endrun


  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
!
! !DESCRIPTION:
!
!  This routine assigns the snow progognostic variables to JULES
!  model space.
!
!EOP
  type(ESMF_State)       :: LSM_State
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh
  integer                :: t, pft
  integer                :: status
  real                   :: ncount(LIS_rc%ngrid(n))
  integer                :: tid, gid
  integer                :: m, k, start_k, end_k

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     pft = jules50_struc(n)%jules50(t)%pft

     dsneqv = swe(t) - jules50_struc(n)%jules50(t)%snow_mass_ij     ! unit: mm
     dsnowh = snod(t) - jules50_struc(n)%jules50(t)%snowdepth(pft)  ! unit: m

     ! update snow prognostic variables using jules's snow physics
     call jules50_snodep_update(n, t, dsneqv, dsnowh)

  enddo

  if (LIS_rc%snowsrc(n) .gt. 0) then
     ncount = 0 ! Tiles per grid id (land only)
     LIS_snow_struc(n)%snowdepth = 0.0 ! Mean snow depth by grid id
     LIS_snow_struc(n)%sneqv = 0.0     ! SWE by tile

     ! Store SWE by tile
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid) = LIS_snow_struc(n)%sneqv(tid) + &
             jules50_struc(n)%jules50(t)%snow_mass_ij
     end do

     ! Store mean snow depth per grid box
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        pft = jules50_struc(n)%jules50(t)%pft
        gid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%index
        LIS_snow_struc(n)%snowdepth(gid) = &
             LIS_snow_struc(n)%snowdepth(gid) + &
             jules50_struc(n)%jules50(t)%snowdepth(pft)
        ncount(gid) = ncount(gid) + 1
     end do
     do t = 1, LIS_rc%ngrid(n)
        if (ncount(t) .gt. 0) then
           LIS_snow_struc(n)%snowdepth(t) = &
                LIS_snow_struc(n)%snowdepth(t) / ncount(t)
        else
           LIS_snow_struc(n)%snowdepth(t) = 0.0
        endif
     end do

  end if

end subroutine jules50_setsnodepvars

