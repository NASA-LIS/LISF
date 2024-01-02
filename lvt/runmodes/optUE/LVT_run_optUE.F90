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
! !ROUTINE: LVT_run_optUE
!  label(LVT_run_optUE)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine LVT_run_optUE

  use LVT_coreMod,         only : LVT_rc, LVT_endofrun, LVT_ticktime
  use LVT_optUEMod,           only : LVT_readoptUEData, LVT_computeoptUEstats, &
       LVT_optuectl
  use LVT_logMod

  implicit none
  integer :: i 

  do i=1,LVT_optuectl%maxIter
     call LVT_readoptUEData(i)
     call LVT_computeoptUEStats(i)
     flush(LVT_logunit)
  enddo
end subroutine LVT_run_optUE
  
