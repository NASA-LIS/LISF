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
! !ROUTINE: templateOpenWater_finalize
! \label{templateOpenWater_finalize}
! 
! !INTERFACE:
subroutine templateOpenWater_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use templateOpenWaterMod, only : templateOpenWater_struc

! !ARGUMENTS: 

!
! !DESCRIPTION:
! 
!  This routine cleans up the allocated memory structures in 
!   the templateOpenWater (forcing-only option)
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
     deallocate(templateOpenWater_struc(n)%templateOpenWater)
  enddo
  deallocate(templateOpenWater_struc)

end subroutine templateOpenWater_finalize
