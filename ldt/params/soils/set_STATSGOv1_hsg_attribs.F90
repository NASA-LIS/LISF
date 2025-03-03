!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_STATSGOv1_hsg_attribs
!  \label{set_STATSGOv1_hsg_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_STATSGOv1_hsg_attribs()

! !USES:
  use LDT_paramDataMod
  use LDT_soilsMod

  implicit none

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!  \end{description}
!EOP      
!
  LDT_soils_struc(:)%hsg%num_bins = 5
  LDT_soils_struc(:)%hsg%vlevels = 1

end subroutine set_STATSGOv1_hsg_attribs
