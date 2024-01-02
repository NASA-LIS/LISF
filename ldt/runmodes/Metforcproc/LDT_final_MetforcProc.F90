!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine LDT_final_MetforcProc

  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_metforcingMod
  use LDT_logMod

  implicit none
  integer :: i,n

  print *, "calling LDT_metforcing_finalize"  ! KRA 08/17/2015
  call LDT_metforcing_finalize()

  write(LDT_logunit,*) "--------------------------------"
  write(LDT_logunit,*) " Finished LDT final MetforcProc "
  write(LDT_logunit,*) "--------------------------------"

end subroutine LDT_final_MetforcProc
  
