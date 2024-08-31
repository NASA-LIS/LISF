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
! !ROUTINE: template_finalize
! \label{template_finalize}
! 
! !INTERFACE:
subroutine template_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use template_lsmMod, only : template_struc

! !ARGUMENTS: 

!
! !DESCRIPTION:
! 
!  This routine cleans up the allocated memory structures in 
!   the template (forcing-only option)
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
     deallocate(template_struc(n)%template)
  enddo
  deallocate(template_struc)

end subroutine template_finalize
