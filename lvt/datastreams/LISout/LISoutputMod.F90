!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: LISoutputMod
! \label(LISoutputMod)
!
! !INTERFACE:
module LISoutputMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation output as 
!  "observations". 
! 
!  NOTE: The plugin does not work if the "LIS LSM" output is written
!  in binary format and if a different LSM is used to create
!  the "LIS LSM" output. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISoutputInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!
!EOP
  
contains

!BOP
! 
! !ROUTINE: LISoutputInit
! \label{LISoutputInit}
!
! !INTERFACE: 
  subroutine LISoutputInit(source)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: source  
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
    !nothing to be done here since LIS output plugin is handled more generically
  end subroutine LISoutputInit
end module LISoutputMod
