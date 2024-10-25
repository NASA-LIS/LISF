!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_init_USAFSI()

   ! Imports
   use LDT_domainMod, only: LDT_setDomainSpecs
   use LDT_logMod, only: LDT_logunit
   use LDT_paramProcMod, only: LDT_paramProcConfig
   use LDT_usafsiMod, only: LDT_usafsiInit

   ! Defaults
   implicit none

   write(LDT_logunit,*) "----------------------------------------"
   write(LDT_logunit,*) " Start of USAFSI processing "
   write(LDT_logunit,*) "----------------------------------------"

   call LDT_setDomainSpecs()
   call LDT_paramProcConfig()
   call LDT_usafsiInit()
   flush(LDT_logunit)

end subroutine LDT_init_USAFSI
