!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module noahmp50_daveg_Mod
!BOP
!
! !MODULE: noahmp50_daveg_Mod
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
  public :: noahmp50_daveg_init
!EOP

contains
!BOP
! 
! !ROUTINE: noahmp50_daveg_init
! \label{noahmp50_daveg_init}
! 
! !INTERFACE:
  subroutine noahmp50_daveg_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP   

    implicit none
    integer, intent(in) :: k

  end subroutine noahmp50_daveg_init
end module noahmp50_daveg_Mod
