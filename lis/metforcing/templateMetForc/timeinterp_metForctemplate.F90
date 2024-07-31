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
! !ROUTINE: timeinterp_metForcTemplate
! \label{timeinterp_metForcTemplate}
!
!               forcing data. 
! !INTERFACE:
subroutine timeinterp_metForcTemplate(n, findex)
! !USES:
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION: 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing
!  \end{description}
! 
!EOP

  return

end subroutine timeinterp_metForcTemplate
