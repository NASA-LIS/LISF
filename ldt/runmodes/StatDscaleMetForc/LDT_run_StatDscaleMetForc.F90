!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_StatDscaleMetForc

  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_metforcingMod
  use LDT_statDscaleMod
  use LDT_logMod

  implicit none
  integer :: npass,pass, n

  npass = 1
  if(LDT_rc%statDscaleMode.eq."downscale") then 
     do pass=1,npass   ! Two passes on the metforcing data

        do while (.NOT. LDT_endofrun())
           call LDT_ticktime                         ! Advance time
           do n=1,LDT_rc%nnest
              call LDT_get_met_forcing(n)            
              call LDT_diagnoseForcStatDscale(n,pass) 
           enddo
        enddo
        call LDT_computeForcStatDscaleParams(pass)    
        call LDT_outputForcStatDscaleParams(pass)    
!           call LDT_resetMetDscaleTimeWindow(pass)   
!           call LDT_resetDscaleVars(pass)            
        call LDT_metforcing_reset()               
        call LDT_resetForcingVars()               
     
        flush(LDT_logunit)
     enddo

  elseif(LDT_rc%statDscaleMode.eq."forecast") then 
   ! use the parameters and generate new forcing. 

  endif

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_StatDscaleMetForc
  
