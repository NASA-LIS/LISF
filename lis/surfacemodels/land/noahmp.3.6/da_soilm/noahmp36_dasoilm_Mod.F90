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
module noahmp36_dasoilm_Mod
!BOP
!
! !MODULE: noahmp36_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! Sujay Kumar; Initial Code
! 9 Sep 2016: Mahdi Navari; Modified for NoahMP36 !

! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp36_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp36_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noahmp36_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: noahmp36_dasoilm_init
! \label{noahmp36_dasoilm_init}
! 
! !INTERFACE:
  subroutine noahmp36_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k

    if(.not.allocated(noahmp36_dasm_struc)) then 
       allocate(noahmp36_dasm_struc(LIS_rc%nnest))
    endif
    

  end subroutine noahmp36_dasoilm_init
end module noahmp36_dasoilm_Mod
