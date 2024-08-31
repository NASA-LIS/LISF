!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_DApreproc

  use LDT_coreMod,           only : LDT_rc, LDT_endofrun, LDT_ticktime
  use LDT_timeMgrMod,        only : LDT_resetClock
  use LDT_DAobservationsMod, only : LDT_readDAobsData
  use LDT_DAobsDataMod,      only : LDT_tavgDAobsData
  use LDT_DAmetricsMod,      only : LDT_readDAdataMask, LDT_diagnoseDAobsMetrics, &
       LDT_computeDAobsMetrics
  use LDT_logMod

  implicit none
  integer :: i,n

  do i=1,LDT_rc%pass
     LDT_rc%pass_id = i
     call LDT_resetClock(LDT_rc)
     do n=1,LDT_rc%nnest
        do while (.NOT. LDT_endofrun())
           call LDT_ticktime   
           call LDT_readDAdataMask(n)
           call LDT_readDAobsData(n)
           call LDT_tavgDAobsData(n)
           call LDT_diagnoseDAobsMetrics(n,i)
           call LDT_computeDAobsMetrics(n,i)
           flush(LDT_logunit)
        enddo
     enddo
  enddo
  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_DApreproc
  
