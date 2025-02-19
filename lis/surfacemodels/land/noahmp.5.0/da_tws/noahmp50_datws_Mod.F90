!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module noahmp50_datws_Mod
!BOP
!
! !MODULE: noahmp50_datws_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li created for Noah-MP.4.0.1
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
! !USES:        
  use LIS_coreMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp50_datws_init
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
! !ROUTINE: noahmp50_datws_init
! \label{noahmp50_datws_init}
! 
! !INTERFACE:
  subroutine noahmp50_datws_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    
    implicit none
    integer                :: k
    !integer                :: n 
    !integer                :: status
    !integer                :: ngrid

    if(.not.allocated(noahmp50_dasm_struc)) then 
       allocate(noahmp50_dasm_struc(LIS_rc%nnest))
    endif
    
  end subroutine noahmp50_datws_init
end module noahmp50_datws_Mod
