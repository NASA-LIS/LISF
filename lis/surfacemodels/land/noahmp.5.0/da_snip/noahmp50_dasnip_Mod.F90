!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module noahmp50_dasnip_Mod
!BOP
!
! !MODULE: noahmp50_dasnip_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
! 19 Oct 2026: Yeosang Yoon; Initial SNIP interface.
!
! !USES:

  implicit none

  PRIVATE
!--------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!--------------------------------------------------------------------------
  public :: noahmp50_dasnip_init
!--------------------------------------------------------------------------
! !PUBLIC TYPES:
!--------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
!
! !ROUTINE: noahmp50_dasnip_init
! \label{noahmp50_dasnip_init}
!
! !INTERFACE:
  subroutine noahmp50_dasnip_init()
! !USES:
! !DESCRIPTION:
!
!EOP
    implicit none
  end subroutine noahmp50_dasnip_init
end module noahmp50_dasnip_Mod
