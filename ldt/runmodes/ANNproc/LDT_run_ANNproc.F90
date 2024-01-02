!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_ANNproc

  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_ANNMod
  use LDT_logMod

  implicit none
  integer :: i,n

  do i=1,LDT_rc%pass
     LDT_rc%pass_id = i
     call LDT_resetClock(LDT_rc)
     do n=1,LDT_rc%nnest
        do while (.NOT. LDT_endofrun())
           call LDT_ticktime   
           call LDT_readANNinputData(n)
           call LDT_readANNoutputData(n)
           call LDT_tavgANNdata(n)
           call LDT_outputANN(n)
           flush(LDT_logunit)
        enddo
     enddo
  enddo
  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_ANNproc
  
