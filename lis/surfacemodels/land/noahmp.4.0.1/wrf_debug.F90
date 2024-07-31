!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine wrf_debug(unit_number, message)
  use LIS_logMod, only     : LIS_logunit
  implicit none
  integer :: unit_number 
  character(len=*) :: message 
  write(LIS_logunit, *) "Noah-MP.4.0.1 (wrf_debug): ", message 
end subroutine wrf_debug
