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
! !ROUTINE: summa1_finalize
! \label{summa1_finalize}
! 
! !INTERFACE:
subroutine summa1_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use summa1_lsmMod, only : summa1_struc

! !ARGUMENTS: 

!
! !DESCRIPTION:
! 
!  This routine cleans up the allocated memory structures in 
!   the summa1 (forcing-only option)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  implicit none
  integer :: n

  do n = 1, LIS_rc%nnest
     deallocate(summa1_struc(n)%summa1)
  enddo
  deallocate(summa1_struc)

end subroutine summa1_finalize
