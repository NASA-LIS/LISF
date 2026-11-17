!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module noahmp52_datws_Mod
!BOP
!
! !MODULE: noahmp52_datws_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !USES:
  use LIS_coreMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp52_datws_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp52_dasm_struc
!EOP

  type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec

  type(dasm_dec), allocatable :: noahmp52_dasm_struc(:)

contains
!BOP
!
! !ROUTINE: noahmp52_datws_init
! \label{noahmp52_datws_init}
!
! !INTERFACE:
  subroutine noahmp52_datws_init(k)
! !USES:
! !DESCRIPTION:
!
!EOP

    implicit none
    integer                :: k

    if(.not.allocated(noahmp52_dasm_struc)) then
       allocate(noahmp52_dasm_struc(LIS_rc%nnest))
    endif

  end subroutine noahmp52_datws_init
end module noahmp52_datws_Mod
