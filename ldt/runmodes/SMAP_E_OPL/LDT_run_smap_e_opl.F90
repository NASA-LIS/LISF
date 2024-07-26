!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LDT_run_SMAP_E_OPL()

   ! Imports
   use LDT_logMod
   use LDT_smap_e_oplMod

   ! Defaults
   implicit none

   ! Local variables
   integer :: n

   ! Currently single nest is assumed
   n = 1
   call LDT_smap_e_oplRun(n)       

   write(LDT_logunit,*) "--------------------------------"
   write(LDT_logunit,*) " Finished LDT run "
   write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_run_SMAP_E_OPL
