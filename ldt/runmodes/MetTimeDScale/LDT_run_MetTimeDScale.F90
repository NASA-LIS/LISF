!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_MetTimeDScale

  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_metforcingMod
  use LDT_metforcDscaleMod
  use LDT_logMod

  implicit none
  integer :: pass, n

   do while (.NOT. LDT_endofrun())
      do pass=1,2   ! Two passes on the metforcing data
         write(LDT_logunit,*) "[INFO] Time Window Pass :: ",pass
         ! Time-window loop:
         do while(.NOT. LDT_endofMetDscaleTimeWindow())
            call LDT_ticktime                         ! Advance time
            do n=1,LDT_rc%nnest
               call LDT_get_met_forcing(n)            ! Retrieve and timeinterp forcing
               ! Time-window pass #1
               if(pass.eq.1) then 
                  call LDT_diagnoseForMetForcDscale(n) ! Sum and count each field
               endif                             
               ! Time-window pass #2
               if(pass.eq.2) then 
                  call LDT_applyDscaleCorrection(n)   ! TDscale option applied: Multiply weight of F1 to F2
                  call LDT_perturb_forcing(n)         ! Diagnose output (average)
                  call LDT_output_met_forcing(n)      ! Write output
               endif
            enddo  ! End nest loop
         enddo  ! End of time window loop

         call LDT_computeMetForcDscaleParams(pass)    ! Assign sum of F1 to SumF1
         call LDT_resetMetDscaleTimeWindow(pass)      ! Reset Met Dscale Time window
         call LDT_resetDscaleVars(pass)               ! Reset Met Dscale variables
         call LDT_metforcing_reset()                  ! Reset metforcing variables
         call LDT_resetForcingVars()                  ! Reset variables
      enddo     ! Pass 1,2 loops
      flush(LDT_logunit)
   enddo        ! LDT runmode completed running

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_MetTimeDScale
  
