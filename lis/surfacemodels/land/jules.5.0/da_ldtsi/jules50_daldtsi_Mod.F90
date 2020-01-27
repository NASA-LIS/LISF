!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module jules50_daldtsi_Mod
!BOP
!
! !MODULE: jules50_daldtsi_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: jules50_daldtsi_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: jules50_daldtsi_init
! \label{jules50_daldtsi_init}
! 
! !INTERFACE:
  subroutine jules50_daldtsi_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none

  end subroutine jules50_daldtsi_init
end module jules50_daldtsi_Mod
