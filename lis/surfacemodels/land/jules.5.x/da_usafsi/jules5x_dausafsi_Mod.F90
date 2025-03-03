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
module jules5x_dausafsi_Mod
!BOP
!
! !MODULE: jules5x_dausafsi_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI
! 17 Feb 2020: Yeosang Yoon; Modified for Jules 5.x
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: jules5x_dausafsi_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: jules5x_dausafsi_init
! \label{jules5x_dausafsi_init}
! 
! !INTERFACE:
  subroutine jules5x_dausafsi_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none

  end subroutine jules5x_dausafsi_init
end module jules5x_dausafsi_Mod
