!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_run_EnsRstpreproc

  use LDT_logMod
  use LDT_ensRstMod 

  implicit none

  call LDT_transformEnsRst

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT run "
  write(LDT_logunit,*) "--------------------------------"

  

end subroutine LDT_run_EnsRstpreproc
  
