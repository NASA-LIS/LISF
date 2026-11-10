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
!
! !ROUTINE: NoahMP52_dynsetup
! \label{NoahMP52_dynsetup}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored noahmp v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine NoahMP52_dynsetup(n)
! !USES:
  use LIS_coreMod, only    : LIS_rc, LIS_surface
  use LIS_fileIOMod, only  : LIS_read_param
  use LIS_snowMod, only    : LIS_snow_struc
  use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
  use NoahMP52_lsmMod, only : NoahMP52_struc

!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in NoahMP
!
!EOP
  implicit none
  integer, intent(in) :: n

  integer   :: tid
  integer   :: t, gid
  integer   :: ncount(LIS_rc%ngrid(n))

  if (LIS_rc%snowsrc(n) .gt. 0) then

     ncount = 0 ! Number of tiles per grid id (over land)
     LIS_snow_struc(n)%snowdepth = 0 ! At grid points
     LIS_snow_struc(n)%sneqv = 0     ! At tiles

     ! Collect SWE at tiles
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid) = LIS_snow_struc(n)%sneqv(tid) + &
             NoahMP52_struc(n)%noahmp52(t)%sneqv
     end do

     ! Collect mean snow depth at grid points
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        LIS_snow_struc(n)%snowdepth(gid) = &
             LIS_snow_struc(n)%snowdepth(gid) + &
             NoahMP52_struc(n)%noahmp52(t)%snowh
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

end subroutine NoahMP52_dynsetup

