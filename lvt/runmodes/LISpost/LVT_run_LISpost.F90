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
! !ROUTINE: LVT_run_LISpost
! \label{LVT_run_LISpost}
!
! !INTERFACE: 
subroutine LVT_run_LISpost

! 
! !USES: 
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_LISpostMod
  use LVT_logMod
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

  call LVT_resetClock(LVT_rc)
  do while (.NOT. LVT_endofrun())
     call LVT_ticktime
     call LVT_process_LISoutput()
     flush(LVT_logunit)
  enddo
  
  write(LVT_logunit,*) '[INFO] Finished LVT analysis'
  write(LVT_logunit,*) '[INFO] --------------------------------------'

end subroutine LVT_run_LISpost
  
