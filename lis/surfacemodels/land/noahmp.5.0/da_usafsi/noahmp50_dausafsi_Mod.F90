!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module noahmp50_dausafsi_Mod
!BOP
!
! !MODULE: noahmp50_dausafsi_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
! May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp50_dausafsi_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: noahmp50_dausafsi_init
! \label{noahmp50_dausafsi_init}
! 
! !INTERFACE:
  subroutine noahmp50_dausafsi_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
  end subroutine noahmp50_dausafsi_init
end module noahmp50_dausafsi_Mod
