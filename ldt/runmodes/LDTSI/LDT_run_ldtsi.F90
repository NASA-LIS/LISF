!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LDT_run_ldtsi()

   ! Imports
   use LDT_logMod
   use LDT_ldtsiMod 

   ! Defaults
   implicit none

   call LDT_ldtsiRun()

   write(LDT_logunit,*) "--------------------------------"
   write(LDT_logunit,*) " Finished LDT run "
   write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_ldtsi
