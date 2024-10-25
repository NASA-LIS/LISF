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
! !ROUTINE: noah39_setusafsivars
! \label{noah39_setusafsivars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  09 Apr 2019: Eric Kemp; Modified for Noah 3.9 and LDT-SI
!  13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
!  16 Nov 2020: Eric Kemp; Added check to LIS_rc%snowsrc(n), and reformatted
!
! !INTERFACE:
subroutine noah39_setusafsivars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_logMod, only : LIS_verify
  use LIS_snowMod, only : LIS_snow_struc
  use noah39_lsmMod

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
  integer                :: t, tid, gid
  integer                :: status

  call ESMF_StateGet(LSM_State, "SWE", sweField, rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State, "Snowdepth", snodField, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField, localDE=0, farrayPtr=swe, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField, localDE=0, farrayPtr=snod, rc=status)
  call LIS_verify(status)

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
    if (noah39_struc(n)%noah(t)%t1 >= 273.15+2.0) then    ! Yeosang Yoon
       if (snod(t) <= 1.e-6 .or. swe(t) <= 1.E-3) THEN
          snod(t) = 0.0
          swe(t) = 0.0
       end if
    end if
  end do

  if (LIS_rc%snowsrc(n) .gt. 0) then
     LIS_snow_struc(n)%sneqv = 0.0
     LIS_snow_struc(n)%snowdepth = 0.0
     npts = 0

     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        !transform t to the patch
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        noah39_struc(n)%noah(t)%sneqv = swe(t)
        LIS_snow_struc(n)%sneqv(tid) = LIS_snow_struc(n)%sneqv(tid) + swe(t)
     end do
     npts = 0
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        !transform t to the patch
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        gid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%index
        noah39_struc(n)%noah(t)%snowh = snod(t)
        npts(gid) = npts(gid) + 1
        LIS_snow_struc(n)%snowdepth(gid) = &
             LIS_snow_struc(n)%snowdepth(gid) + snod(t)
     end do

     do t = 1, LIS_rc%ngrid(n)
        if (npts(t) .gt. 0) then
           LIS_snow_struc(n)%snowdepth(t) = &
                LIS_snow_struc(n)%snowdepth(t) / npts(t)
        else
           LIS_snow_struc(n)%snowdepth(t) = 0.0
        end if
     end do
  end if
end subroutine noah39_setusafsivars

