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
! !ROUTINE: vic412_dynsetup
! \label{vic412_dynsetup}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
! !INTERFACE:
subroutine vic412_dynsetup(n)
! !USES:

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
! This routine sets up the time-dependent variables in VIC. 
! Currently this routine is a placeholder.

!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

end subroutine vic412_dynsetup
