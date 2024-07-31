!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMMRSNWDsnow_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the 
!  SMMR snow depth data
! 
! !REVISION HISTORY: 
!  15 July 2016: Sujay Kumar, Initial Specification
!
module SMMRSNWDsnow_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMMRSNWDsnow_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

contains
  
!BOP
! 
! !ROUTINE: SMMRSNWDsnow_obsInit
! \label{SMMRSNWDsnow_obsInit}
! 
! !INTERFACE: 
  subroutine SMMRSNWDsnow_obsinit()
! !USES: 
    use LDT_DAobsDataMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SMMR Snow depth data. Currently the snow DA schemes
! in LIS only require the definition of the observation grid. 
! Therefore this routine is empty. 
! 
!EOP
    integer                 :: n 
    real                    :: gridDesci(20)
    
! Note that the following grid definition is only a placeholder. This 
! does not actually impact any computations. 

    n = 1

    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = 1440
    gridDesci(3) = 360
    gridDesci(4) = 0.125
    gridDesci(5) = -179.875
    gridDesci(6) = 128
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64
    
    call LDT_initializeDAobsEntry(LDT_DAobsData(n)%snowdepth_obs, &
         "m",1,1)
    LDT_DAobsData(n)%snowdepth_obs%selectStats = 1

  end subroutine SMMRSNWDsnow_obsinit
     

end module SMMRSNWDsnow_obsMod
