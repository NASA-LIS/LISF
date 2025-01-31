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
module NoahMP401_dasoilm_Mod
!BOP
!
! !MODULE: NoahMP401_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! 15 Dec 2018: Mahdi Navari, Sujay Kumar ; Modified for NoahMP401 !

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
  public :: NoahMP401_dasoilm_init
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
! !ROUTINE: NoahMP401_dasoilm_init
! \label{NoahMP401_dasoilm_init}
! 
! !INTERFACE:
  subroutine NoahMP401_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k

    if(.not.allocated(noahmp401_dasm_struc)) then 
       allocate(noahmp401_dasm_struc(LIS_rc%nnest))
    endif
    

  end subroutine NoahMP401_dasoilm_init
end module NoahMP401_dasoilm_Mod
