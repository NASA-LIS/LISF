!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module jules_lis_exchange
  USE datetime_mod, ONLY : DATETIME_STR_LEN
  implicit none 
  CHARACTER(len=DATETIME_STR_LEN) :: main_run_start
  CHARACTER(len=DATETIME_STR_LEN) ::  main_run_end
end module jules_lis_exchange 

