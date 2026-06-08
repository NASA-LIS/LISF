!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: LDT_init_wsf_opl
! \label{LDT_wsf_oplMod}
!
! !REVISION HISTORY:
! 09 Oct 2025: Ehsan Jalilvand; Initial Specification
!
! !INTERFACE:
subroutine LDT_init_wsf_opl()
! !USES:
   use LDT_domainMod, only: LDT_setDomainSpecs
   use LDT_logMod, only: LDT_logunit
   use LDT_paramProcMod, only: LDT_paramProcConfig
   use LDT_wsf_oplMod, only: LDT_wsf_oplInit

   implicit none
!
! !ARGUMENTS:
! none
!
! !DESCRIPTION:
!  This subroutine initializes the runmode for resampling OPL WSF
!  brightness temperatures.
!EOP

   write(LDT_logunit,*) "Start of WSF Low Resolution Resampling"
   call LDT_setDomainSpecs()
   call LDT_paramProcConfig()
   call LDT_wsf_oplInit()
   flush(LDT_logunit)
end subroutine LDT_init_wsf_opl
