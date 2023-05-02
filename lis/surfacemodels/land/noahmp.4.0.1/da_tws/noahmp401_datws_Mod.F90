!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!29 May 2020: Bailing Li created for Noah-MP.4.0.1
#include "LIS_misc.h"
module noahmp401_datws_Mod
!BOP
!
! !MODULE: noahmp401_datws_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! 14 Mar 2017: Sujay Kumar; Initial Specification

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
  public :: noahmp401_datws_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp401_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec

  type(dasm_dec), allocatable :: noahmp401_dasm_struc(:)

contains
!BOP
!
! !ROUTINE: noahmp401_datws_init
! \label{noahmp401_datws_init}
!
! !INTERFACE:
  subroutine noahmp401_datws_init(k)
! !USES:
! !DESCRIPTION:
!
!EOP

    implicit none
    integer                :: k
    integer                :: n
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(noahmp401_dasm_struc)) then
       allocate(noahmp401_dasm_struc(LIS_rc%nnest))
    endif

  end subroutine noahmp401_datws_init
end module noahmp401_datws_Mod
