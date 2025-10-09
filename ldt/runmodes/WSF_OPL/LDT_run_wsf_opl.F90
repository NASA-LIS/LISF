subroutine LDT_run_wsf_opl()
   use LDT_logMod, only: LDT_logunit
   use LDT_wsf_oplMod, only: LDT_wsf_oplRun
   implicit none
   integer :: n

   n = 1
   call LDT_wsf_oplRun(n)
   write(LDT_logunit,*) "Finished LDT WSF Resampling"
end subroutine LDT_run_wsf_opl