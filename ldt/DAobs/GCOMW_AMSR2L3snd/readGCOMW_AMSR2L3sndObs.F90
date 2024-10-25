!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readGCOMW_AMSR2L3sndObs
! \label{readGCOMW_AMSR2L3sndObs}
! 
! !REVISION HISTORY: 
!  15 July 2016: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGCOMW_AMSR2L3sndObs(n)
! !USES:   

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine is a placeholder for the routine that reads
! the AMSR2 snow depth observations. Currently the snow DA schemes
! in LIS only require the definition of the observation grid. 
! Therefore this routine is empty. 
!
!EOP


end subroutine readGCOMW_AMSR2L3sndObs

