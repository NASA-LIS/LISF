!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_climoRstproc

  use LDT_coreMod
  use LDT_logMod
  use LDT_climoRstProcMod 

  implicit none
  
  integer       :: n 
  
  do n=1,LDT_rc%nnest
     do while (.NOT. LDT_endofrun())
        call LDT_ticktime
        call LDT_diagnoseClimoRstData(n)
     enddo
  enddo
!  call LDT_output_rstProc

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

  

end subroutine LDT_run_climoRstproc
  
