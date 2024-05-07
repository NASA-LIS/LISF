!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
