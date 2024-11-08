!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noahmp36_dasoilmLAI_Mod
!BOP
!
! !MODULE: noahmp36_dasoilmLAI_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! Sujay Kumar; Initial Code
! 9 Sep 2016: Mahdi Navari; Modified for NoahMP36 !
! 30 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
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
  public :: noahmp36_dasoilmLAI_init
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
! !ROUTINE: noahmp36_dasoilmLAI_init
! \label{noahmp36_dasoilmLAI_init}
! 
! !INTERFACE:
  subroutine noahmp36_dasoilmLAI_init(k)
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

    if(.not.allocated(noahmp36_dasm_struc)) then 
       allocate(noahmp36_dasm_struc(LIS_rc%nnest))
    endif
   
  ! scaling not active for SM and LAI updating with S1 DA

  end subroutine noahmp36_dasoilmLAI_init
end module noahmp36_dasoilmLAI_Mod
