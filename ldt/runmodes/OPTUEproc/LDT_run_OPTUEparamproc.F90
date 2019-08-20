!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_OPTUEparamproc

  use LDT_coreMod,         only : LDT_rc
  use LDT_paramProcMod
  use LDT_OPTUEMod
  use LDT_logMod

  implicit none

  integer :: n 

  do n=1,LDT_rc%nnest
     call LDT_paramProcWriteHeader(n)
     call LDT_paramProcWrite(n)
     call LDT_paramProcWriteFinalize(n)
  enddo

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"
  
end subroutine LDT_run_OPTUEparamproc
  
