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
module noahmp401_daveg_Mod
!BOP
!
! !MODULE: noahmp401_daveg_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
!
! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp401_daveg_init
!EOP

contains
!BOP
! 
! !ROUTINE: noahmp401_daveg_init
! \label{noahmp401_daveg_init}
! 
! !INTERFACE:
  subroutine noahmp401_daveg_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP   

    implicit none
    integer, intent(in) :: k

  end subroutine noahmp401_daveg_init
end module noahmp401_daveg_Mod
