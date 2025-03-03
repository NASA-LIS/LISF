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
! !ROUTINE: vic411_finalize
! \label{vic411_finalize}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
! !INTERFACE:
subroutine vic411_finalize()
! !USES:
  use LIS_coreMod,   only : LIS_rc
  use vic411_lsmMod, only : vic411_struc

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
  integer :: n
  integer :: lis_nest, lis_npatch

  do n = 1, LIS_rc%nnest
     lis_nest = n
     lis_npatch = LIS_rc%npatch(n,LIS_rc%lsm_index)
     call finalize_vic411(lis_nest, lis_npatch)
     deallocate(vic411_struc(n)%vic)
  enddo
  deallocate(vic411_struc)

end subroutine vic411_finalize
