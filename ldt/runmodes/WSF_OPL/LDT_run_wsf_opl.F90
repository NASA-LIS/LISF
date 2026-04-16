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
subroutine LDT_run_wsf_opl()
! !USES:
   use LDT_logMod, only: LDT_logunit
   use LDT_wsf_oplMod, only: LDT_wsf_oplRun

   implicit none
!
! !ARGUMENTS:
! none
!
! !DESCRIPTION:
!  This subroutine runs the runmode for resampling OPL WSF
!  brightness temperatures.
!EOP

   integer :: n

   n = 1
   call LDT_wsf_oplRun(n)
   write(LDT_logunit,*) "Finished LDT WSF Resampling"
end subroutine LDT_run_wsf_opl
