!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

#if (!defined COUPLED)
subroutine wrf_error_fatal (string)
   character (len=*) :: string
   print *,string
   stop
end subroutine wrf_error_fatal

subroutine wrf_message(message)

    implicit none

    character(len=*), intent(in) :: message

    write(0,*) trim(message)
end subroutine wrf_message

subroutine wrf_error_fatal3(file, line, message)

    implicit none

    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    character(len=*), intent(in) :: message

    write(0,*) trim(file), 'line: ', line, ': ', trim(message)
    stop (1)
end subroutine wrf_error_fatal3
#endif
