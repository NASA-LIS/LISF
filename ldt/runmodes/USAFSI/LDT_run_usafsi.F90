!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LDT_run_usafsi()

   ! Imports
   use LDT_logMod
   use LDT_usafsiMod

   ! Defaults
   implicit none

   ! Local variables
   integer :: n

   ! Currently USAFSI assumes single nest
   n = 1
   call LDT_usafsiRun(n)

   write(LDT_logunit,*) "--------------------------------"
   write(LDT_logunit,*) " Finished LDT run "
   write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_usafsi
