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
!BOP
!
! !ROUTINE: vic411_readrst
! \label{vic411_readrst}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 
! !INTERFACE:
subroutine vic411_readrst()

! !USES:

  implicit none
! !ARGUMENTS: 

!
! !DESCRIPTION:
!  This program reads restart files for VIC.  This
!  includes all relevant water/energy storages and tile information. 
!
!  Actually, VIC contains support to read its initial state file/restart file
!  in the initialize\_model\_state routine.  This support cannot be easily
!  pulled out into its own routine because the support to process the 
!  initial state file/restart file requires air temperature to be present.
!  (stand-alone VIC completely processes all forcing data before it initializes
!  the model states.)
!
!  Thus reading a restart file is handled in the vic411\_run routine at the
!  first time-step of the run.
!EOP

end subroutine vic411_readrst
