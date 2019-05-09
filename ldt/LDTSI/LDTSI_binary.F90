! Temporary subroutine for writing LDTSI fields in binary

#include "LDT_misc.h"

subroutine LDTSI_binary()
   
   ! Imports
   use LDT_logMod, only: LDT_logunit
   use LDTSI_arraysMod, only: LDTSI_arrays

   ! Defaults
   implicit none

   ! Local constants
   integer, parameter :: lunit = 100

   ! Local variables
   character*125 :: file_path
   integer :: istat

   ! Open file
   file_path = "ldtsi.bin"
   write(LDT_logunit,*)"Writing ldtsi.bin"
   open(unit=lunit,file=file_path, form='unformatted', action='write', &
        iostat=istat, status='unknown')
   
   ! Write LDTSI fields
   write(lunit) LDTSI_arrays%snoanl
   write(lunit) LDTSI_arrays%snoage
   write(lunit) LDTSI_arrays%icecon
   write(lunit) LDTSI_arrays%icemask
   write(lunit) LDTSI_arrays%iceage
   
   ! Close file
   close(lunit)

end subroutine LDTSI_binary
