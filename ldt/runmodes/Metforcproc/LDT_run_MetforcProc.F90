!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_MetforcProc

  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_metforcingMod
  use LDT_logMod

  implicit none
  integer :: i,n

   do while (.NOT. LDT_endofrun())
! Run each nest separately 
      call LDT_ticktime
      do n=1,LDT_rc%nnest
!         if(LDT_timeToRunNest(n)) then
            call LDT_get_met_forcing(n)
            call LDT_perturb_forcing(n)
            call LDT_output_met_forcing(n)
!        endif
      enddo
      flush(LDT_logunit)
   enddo


  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_MetforcProc
  
