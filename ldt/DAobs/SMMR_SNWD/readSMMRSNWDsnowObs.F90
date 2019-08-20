!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readSMMRSNWDsnowObs
! \label{readSMMRSNWDsnowObs}
! 
! !REVISION HISTORY: 
!  15 July 2016: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readSMMRSNWDsnowObs(n)
! !USES:   

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine is a placeholder for the routine that reads
! the SMMR snow depth observations. Currently the snow DA schemes
! in LIS only require the definition of the observation grid. 
! Therefore this routine is empty. 
!
!EOP


end subroutine readSMMRSNWDsnowObs

