!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noahmp401_daburnarea_Mod
!BOP
!
! !MODULE: noahmp401_daburnarea_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! 24 Jul 2022: Sujay Kumar; Initial Specification
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
  public :: noahmp401_daburnarea_init
!EOP

contains
!BOP
! 
! !ROUTINE: noahmp401_daburnarea_init
! \label{noahmp401_daburnarea_init}
! 
! !INTERFACE:
  subroutine noahmp401_daburnarea_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP   

    implicit none
    integer, intent(in) :: k

  end subroutine noahmp401_daburnarea_init
end module noahmp401_daburnarea_Mod
