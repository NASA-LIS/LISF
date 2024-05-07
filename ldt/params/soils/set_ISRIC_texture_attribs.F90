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
! !ROUTINE: set_ISRIC_texture_attribs
!  \label{set_ISRIC_texture_attribs}
!
! !REVISION HISTORY:
!  31 Jul 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_ISRIC_texture_attribs()

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
  LDT_LSMparam_struc(:)%texture%num_bins = 16
  LDT_LSMparam_struc(:)%texture%vlevels = 1

  LDT_soils_struc(:)%texture_nlyrs%num_bins = 1
  LDT_soils_struc(:)%texture_nlyrs%vlevels = 11

end subroutine set_ISRIC_texture_attribs
