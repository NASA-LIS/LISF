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
! !ROUTINE: clm2_dynsetup
! \label{clm2_dynsetup}
! 
! !REVISION HISTORY: 
! 
! 20 Jan 2003  Sujay Kumar Initial Specification
! 
! !INTERFACE:
subroutine clm2_dynsetup(n)

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in CLM
!  The routine invokes methods to read the LAI data
! 
!EOP    
  call clm2_lairead(n)
  
end subroutine clm2_dynsetup

