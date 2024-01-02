!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Temporary subroutine for writing USAFSI fields in binary

#include "LDT_misc.h"

subroutine USAFSI_binary()
   
   ! Imports
   use LDT_logMod, only: LDT_logunit
   use USAFSI_arraysMod, only: USAFSI_arrays

   ! Defaults
   implicit none

   ! Local constants
   integer, parameter :: lunit = 100

   ! Local variables
   character*125 :: file_path
   integer :: istat

   ! Open file
   file_path = "usafsi.bin"
   write(LDT_logunit,*)"Writing usafsi.bin"
   open(unit=lunit,file=file_path, form='unformatted', action='write', &
        iostat=istat, status='unknown')
   
   ! Write USAFSI fields
   write(lunit) USAFSI_arrays%snoanl
   write(lunit) USAFSI_arrays%snoage
   write(lunit) USAFSI_arrays%icecon
   write(lunit) USAFSI_arrays%icemask
   write(lunit) USAFSI_arrays%iceage
   
   ! Close file
   close(lunit)

end subroutine USAFSI_binary
