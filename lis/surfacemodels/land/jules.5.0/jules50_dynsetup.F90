!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: jules50_dynsetup
! \label{jules50_dynsetup}
!
! !INTERFACE:
subroutine jules50_dynsetup(n)
  ! !USES:
  use jules50_lsmMod, only: jules50_struc
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_snowMod, only : LIS_snow_struc

  implicit none
  ! !ARGUMENTS:
  integer, intent(in)  :: n
  !
  ! !DESCRIPTION:
  !
  ! Copies snow depth and snow water equivalent into LIS_snow_struc (used
  ! by AGRMET Ops runmode).
  !
  !  The arguments are:
  !  \begin{description}
  !  \item[n]
  !   index of the nest
  !  \end{description}
  !EOP

  integer :: gid
  integer :: ncount(LIS_rc%ngrid(n))
  integer :: t
  integer :: tid
  integer :: m, k, start_k, end_k, pft

  if (LIS_rc%snowsrc(n) .gt. 0) then

     ncount = 0 ! Tiles per grid id (land only)
     LIS_snow_struc(n)%snowdepth = 0.0 ! Snow depth by grid id
     LIS_snow_struc(n)%sneqv     = 0.0 ! SWE by tile

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

end subroutine jules50_dynsetup
