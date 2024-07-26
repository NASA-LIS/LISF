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
! !ROUTINE: vic412_finalize
! \label{vic412_finalize}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
! !INTERFACE:
subroutine vic412_finalize()
! !USES:
  use LIS_coreMod,   only : LIS_rc
  use vic412_lsmMod, only : vic412_struc

! !ARGUMENTS: 

!
! !DESCRIPTION:
!  This routine cleans up the allocated memory structures in 
!  the VIC 4.1.1 lsm.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  implicit none
  integer :: n, t
  integer :: lis_nest, lis_npatch

  do n = 1, LIS_rc%nnest
     lis_nest = n
     lis_npatch = LIS_rc%npatch(n,LIS_rc%lsm_index)
     call finalize_vic412(lis_nest, lis_npatch)
     ! added by Shugong Wang 05/30/2014 
     do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(vic412_struc(n)%vic(t)%state_chunk)
     enddo
     deallocate(vic412_struc(n)%vic)
  enddo
  deallocate(vic412_struc)

end subroutine vic412_finalize
