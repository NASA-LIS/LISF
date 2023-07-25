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
module noahmpnew_datws_Mod
!BOP
!
! !MODULE: noahmpnew_datws_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li created for Noah-MP.4.0.1
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
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
  public :: noahmpnew_datws_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmpnew_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noahmpnew_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: noahmpnew_datws_init
! \label{noahmpnew_datws_init}
! 
! !INTERFACE:
  subroutine noahmpnew_datws_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    implicit none
    integer                :: k
    integer                :: n 
    character(len=LIS_CONST_PATH_LEN) :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(noahmpnew_dasm_struc)) then 
       allocate(noahmpnew_dasm_struc(LIS_rc%nnest))
    endif
    
  end subroutine noahmpnew_datws_init
end module noahmpnew_datws_Mod
