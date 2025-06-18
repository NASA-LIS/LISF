!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module NoahMP50_dasoilm_Mod
!BOP
!
! !MODULE: NoahMP50_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later

! !USES:        
  use LIS_coreMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NoahMP50_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp50_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: noahmp50_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: NoahMP50_dasoilm_init
! \label{NoahMP50_dasoilm_init}
! 
! !INTERFACE:
  subroutine NoahMP50_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
    integer                :: k

    if(.not.allocated(noahmp50_dasm_struc)) then 
       allocate(noahmp50_dasm_struc(LIS_rc%nnest))
    endif
    
  end subroutine NoahMP50_dasoilm_init
end module NoahMP50_dasoilm_Mod
