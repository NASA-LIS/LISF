!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Temporary subroutine for writing SNIP fields in binary

#include "LDT_misc.h"

subroutine SNIP_binary()

  ! Imports
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_logMod, only: LDT_logunit
  use SNIP_arraysMod, only: SNIP_arrays

  ! Defaults
  implicit none

  ! Local constants
  integer, parameter :: lunit = 100

  ! Local variables
  character(len=LDT_CONST_PATH_LEN) :: file_path
  integer :: istat

  ! Open file
  file_path = "SNIP.bin"
  write(LDT_logunit,*)"Writing SNIP.bin"
  open(unit=lunit,file=file_path, form='unformatted', action='write', &
       iostat=istat, status='unknown')

  ! Write SNIP fields
  write(lunit) SNIP_arrays%snoanl
  write(lunit) SNIP_arrays%snoage
  write(lunit) SNIP_arrays%icecon
  write(lunit) SNIP_arrays%icemask
  write(lunit) SNIP_arrays%iceage

  ! Close file
  close(lunit)

end subroutine SNIP_binary
