!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module noahmp52_dausafsi_Mod
!BOP
!
! !MODULE: noahmp52_dausafsi_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !USES:

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp52_dausafsi_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
!
! !ROUTINE: noahmp52_dausafsi_init
! \label{noahmp52_dausafsi_init}
!
! !INTERFACE:
  subroutine noahmp52_dausafsi_init()
! !USES:
! !DESCRIPTION:
!
!EOP
    implicit none
  end subroutine noahmp52_dausafsi_init
end module noahmp52_dausafsi_Mod
