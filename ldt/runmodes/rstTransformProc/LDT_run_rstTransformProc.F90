!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_rstTransformProc

  use LDT_coreMod
  use LDT_logMod

  implicit none
  
  integer       :: n 
  
  do n=1,LDT_rc%nnest
  enddo
!  call LDT_output_rstProc

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

  

end subroutine LDT_run_rstTransformProc
  
