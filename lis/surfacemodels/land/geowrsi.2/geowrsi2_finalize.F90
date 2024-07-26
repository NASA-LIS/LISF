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
! !ROUTINE: geowrsi2_finalize
! \label{geowrsi2_finalize}
! 
! !INTERFACE:
subroutine geowrsi2_finalize()

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use geowrsi2_lsmMod, only : geowrsi2_struc

! !ARGUMENTS: 
!
! !DESCRIPTION:
! 
!  This routine cleans up the allocated memory structures in 
!   the geowrsi2 (forcing-only option)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  implicit none

  integer :: n, t

  do n = 1, LIS_rc%nnest
     do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
        deallocate(geowrsi2_struc(n)%wrsi(t)%sos_write)
        deallocate(geowrsi2_struc(n)%wrsi(t)%sosa_write)
     enddo
     deallocate(geowrsi2_struc(n)%wrsi)
  enddo
  deallocate(geowrsi2_struc)

end subroutine geowrsi2_finalize
