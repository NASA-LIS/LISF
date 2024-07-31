!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !ROUTINE: LVT_run_DAobs
! \label{LVT_run_DAobs}
!
! !INTERFACE: 
subroutine LVT_run_DAobs
! 
! !USES: 

  use LVT_coreMod,         only : LVT_rc, LVT_endofrun, LVT_ticktime
  use LVT_timeMgrMod,      only : LVT_resetClock
  use LVT_LISModelDataMod, only : LVT_readDAobsData
  use LVT_historyMod,      only : LVT_tavgLISModelData
  use LVT_observationsMod, only : LVT_readObsData, LVT_tavgObsData
  use LVT_statsMod,        only : LVT_diagnoseStats, LVT_computeStats, &
       LVT_readDataMask
  use LVT_logMod

  implicit none
  integer :: i 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This routine defines the run phase of the runmode to 
!   analyse the output from a LIS data assimilation integration.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  do i=1,LVT_rc%pass
     call LVT_resetClock(LVT_rc)
     do while (.NOT. LVT_endofrun())
        call LVT_ticktime   
        call LVT_readDataMask
        call LVT_readDAobsData
        call LVT_readObsData
        call LVT_tavgLISModelData
        call LVT_tavgObsData
        call LVT_diagnoseStats(i)
        call LVT_computeStats(i)
        flush(LVT_logunit)
     enddo    
  enddo

end subroutine LVT_run_DAobs
  
