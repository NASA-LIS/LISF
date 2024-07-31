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
module jules50_dausafsi_Mod
!BOP
!
! !MODULE: jules50_dausafsi_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: jules50_dausafsi_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: jules50_dausafsi_init
! \label{jules50_dausafsi_init}
! 
! !INTERFACE:
  subroutine jules50_dausafsi_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none

  end subroutine jules50_dausafsi_init
end module jules50_dausafsi_Mod
