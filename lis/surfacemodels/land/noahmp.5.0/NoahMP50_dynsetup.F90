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
!
! !ROUTINE: NoahMP50_dynsetup
! \label{NoahMP50_dynsetup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and NoahMP401
!   11 Nov 2020: Eric Kemp, added updates to LIS_snow_struc
!   May 2023: Cenlin He, modified for refactored noahmp v5 and later

! !INTERFACE:
subroutine NoahMP50_dynsetup(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_surface
    use LIS_fileIOMod, only  : LIS_read_param
    use LIS_snowMod, only    : LIS_snow_struc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use NoahMP50_lsmMod, only : NoahMP50_struc

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
               NoahMP50_struc(n)%noahmp50(t)%sneqv
       end do

       ! Collect mean snow depth at grid points
       do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
          LIS_snow_struc(n)%snowdepth(gid) = &
               LIS_snow_struc(n)%snowdepth(gid) + &
                NoahMP50_struc(n)%noahmp50(t)%snowh
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

    !TODO: add code here if needed.
end subroutine NoahMP50_dynsetup

! generate date/time string for reading time-dependent variables
