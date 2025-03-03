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
module noah36_dasoilm_Mod
!BOP
!
! !MODULE: noah36_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!  ---------: Sujay Kumar; Initial Specification
!  28Aug2017: Mahdi Navari; Updated to take into account the latest developments in the SM DA 

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
  public :: noah36_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah36_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noah36_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: noah36_dasoilm_init
! \label{noah36_dasoilm_init}
! 
! !INTERFACE:
  subroutine noah36_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k

    if(.not.allocated(noah36_dasm_struc)) then 
       allocate(noah36_dasm_struc(LIS_rc%nnest))
    endif
    

  end subroutine noah36_dasoilm_init
end module noah36_dasoilm_Mod
