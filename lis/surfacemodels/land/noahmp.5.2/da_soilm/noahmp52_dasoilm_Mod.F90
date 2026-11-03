!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module NoahMP52_dasoilm_Mod
!BOP
!
! !MODULE: NoahMP52_dasoilm_Mod
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
  public :: NoahMP52_dasoilm_init
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
! !ROUTINE: NoahMP52_dasoilm_init
! \label{NoahMP52_dasoilm_init}
!
! !INTERFACE:
  subroutine NoahMP52_dasoilm_init(k)
! !USES:
! !DESCRIPTION:
!
!EOP
    implicit none
    integer                :: k

    if(.not.allocated(noahmp52_dasm_struc)) then
       allocate(noahmp52_dasm_struc(LIS_rc%nnest))
    endif

  end subroutine NoahMP52_dasoilm_init
end module NoahMP52_dasoilm_Mod
