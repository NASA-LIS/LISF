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
!
! !ROUTINE: cable_dynsetup
! \label{cable_dynsetup}
!
! !REVISION HISTORY:
!  26 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  13 Sep 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_dynsetup(n)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_vegDataMod,         only : LIS_lai
  use cable_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine sets up the time-dependent variables in CABLE
!  The arguments are:
!
!  \begin{description}
!   \item[n]
!    index of the nest
!  \end{description}
!EOP
  integer :: t,tid
  
  if (LIS_rc%uselaimap(n).ne."none") then
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        cable_struc(n)%cable(t)%lai = LIS_lai(n)%tlai(tid)
        !KLUDGE: CABLE doesn't work with a LAI of 0. Problem in cable_canopy with gbhu=0. (ccc)
        if (abs(cable_struc(n)%cable(t)%lai) < 0.009) then
           cable_struc(n)%cable(t)%lai = 0.009
        endif
     enddo
  endif
  
end subroutine cable_dynsetup
