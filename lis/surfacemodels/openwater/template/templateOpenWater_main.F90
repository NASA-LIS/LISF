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
#undef TBOT_TESTING
!BOP
!
! !ROUTINE: templateOpenWater_main
! \label{templateOpenWater_main}
!
! !ROUTINE: templateOpenWater_main.F90
!    Apr 2003; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated code for LISv50
! 
! !INTERFACE:
subroutine templateOpenWater_main(n)

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
! 
!  Calls the run routines for the forcing-only option (templateOpenWater)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP



end subroutine templateOpenWater_main
