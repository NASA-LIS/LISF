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
module noah32_dascf_Mod
!BOP
!
! !MODULE: noah32_dascf_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! !USES:        
! none

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah32_dascf_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
contains
!BOP
! 
! !ROUTINE: noah32_dascf_init
! \label{noah32_dascf_init}
! 
! !INTERFACE:
  subroutine noah32_dascf_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
  end subroutine noah32_dascf_init
end module noah32_lsmMod
