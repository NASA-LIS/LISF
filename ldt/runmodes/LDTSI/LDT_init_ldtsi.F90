!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_init_LDTSI()

   ! Imports
   use LDT_domainMod, only: LDT_setDomainSpecs
   use LDT_logMod, only: LDT_logunit, LDT_flush
   use LDT_paramProcMod, only: LDT_paramProcConfig
   use LDT_ldtsiMod, only: LDT_ldtsiInit

   ! Defaults
   implicit none

   write(LDT_logunit,*) "----------------------------------------"
   write(LDT_logunit,*) " Start of LDTSI processing "
   write(LDT_logunit,*) "----------------------------------------"

   call LDT_setDomainSpecs()
   call LDT_paramProcConfig()
   call LDT_ldtsiInit()
   call LDT_flush(LDT_logunit)

end subroutine LDT_init_ldtsi
