!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module HYMAP2_daWL_Mod
!BOP
!
! !MODULE: HYMAP2_daWL_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
!  07 Nov 19: Sujay Kumar; Initial specification
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
  public :: HYMAP2_daWL_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: HYMAP2_daWL_struc
!EOP

 type, public :: daWL_dec

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type daWL_dec
  
  type(daWL_dec), allocatable :: HYMAP2_daWL_struc(:)

contains
!BOP
! 
! !ROUTINE: HYMAP2_daWL_init
! \label{HYMAP2_daWL_init}
! 
! !INTERFACE:
  subroutine HYMAP2_daWL_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    
    implicit none

    integer                :: k
    integer                :: n 
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(HYMAP2_daWL_struc)) then 
       allocate(HYMAP2_daWL_struc(LIS_rc%nnest))
    endif
    
  end subroutine HYMAP2_daWL_init
end module HYMAP2_daWL_Mod
