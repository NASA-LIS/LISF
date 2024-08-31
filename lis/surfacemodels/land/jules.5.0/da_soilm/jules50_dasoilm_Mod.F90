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
module jules50_dasoilm_Mod
!BOP
!
! !MODULE: jules50_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! Sujay Kumar; Initial Code
! 23Apr2018: Mahdi Navari; Modified for JULES 5.0 
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
  public :: jules50_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: jules50_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer              	:: nbins
     integer             	:: ntimes
     integer           		:: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: jules50_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: jules50_dasoilm_init
! \label{jules50_dasoilm_init}
! 
! !INTERFACE:
  subroutine jules50_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
    integer                :: n 
    character*100          :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(jules50_dasm_struc)) then 
       allocate(jules50_dasm_struc(LIS_rc%nnest))
    endif
    


  end subroutine jules50_dasoilm_init
end module jules50_dasoilm_Mod
