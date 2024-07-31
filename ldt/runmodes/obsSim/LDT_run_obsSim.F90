!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LDT_run_obsSim

  use LDT_coreMod,           only : LDT_rc, LDT_endofrun, LDT_ticktime
  use LDT_timeMgrMod,        only : LDT_resetClock
  use LDT_obsSimMod
  use LDT_logMod

  implicit none
  integer :: i,n

  do i=1,LDT_rc%pass
     LDT_rc%pass_id = i
     call LDT_resetClock(LDT_rc)
     do n=1,LDT_rc%nnest
        do while (.NOT. LDT_endofrun())
           call LDT_ticktime   
           call LDT_readNatureRunData(n)
           call LDT_temporalTransformObsSimData(n)
           call LDT_applyObsSimMask(n)
           call LDT_applyObsSimErrorModel(n)
           call LDT_writeObsSim(n)

        enddo
     enddo
  enddo
  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_obsSim
  
