subroutine LDT_init_wsf_opl()
   use LDT_domainMod, only: LDT_setDomainSpecs
   use LDT_logMod, only: LDT_logunit
   use LDT_paramProcMod, only: LDT_paramProcConfig
   use LDT_wsf_oplMod, only: LDT_wsf_oplInit
   implicit none

   write(LDT_logunit,*) "Start of WSF Low Resolution Resampling"
   call LDT_setDomainSpecs()
   call LDT_paramProcConfig()
   call LDT_wsf_oplInit()
   flush(LDT_logunit)
end subroutine LDT_init_wsf_opl