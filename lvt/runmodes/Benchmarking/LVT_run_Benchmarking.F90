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
! !ROUTINE: LVT_run_Benchmarking
! \label{LVT_run_Benchmarking}
!
! !INTERFACE: 
subroutine LVT_run_Benchmarking

! 
! !USES: 
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_DataStreamsMod
  use LVT_statsMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_domainMod
  use LVT_trainingMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine defines the call structure for performing LVT 
!  analysis using model output from LIS. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  implicit none
  integer :: i 

  do i=LVT_rc%curr_pass,LVT_rc%pass
     call LVT_resetClock(LVT_rc)
     do while (.NOT. LVT_endofrun())
        call LVT_ticktime   
        call LVT_readDataMask
        call LVT_readDataStreams
        call LVT_tavgDataStreams
        call LVT_runTraining(i)
        call LVT_writeTrainingOutput(i)
        call LVT_resetDataStreams
        flush(LVT_logunit)
     enddo    
  enddo
  write(LVT_logunit,*) '[INFO] Finished LVT analysis'
  write(LVT_logunit,*) '[INFO] --------------------------------------'

end subroutine LVT_run_Benchmarking
  
