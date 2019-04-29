! Temporary subroutine for writing SNODEP fields in binary

#include "LDT_misc.h"

subroutine snodep_binary()
   
   ! Imports
   use LDT_logMod, only: LDT_logunit
   use SNODEP_arraysMod, only: SNODEP_arrays

   ! Defaults
   implicit none

   ! Local constants
   integer, parameter :: lunit = 100

   ! Local variables
   character*125 :: file_path
   integer :: istat

   ! Open file
   file_path = "snodep.bin"
   write(LDT_logunit,*)"Writing snodep.bin"
   open(unit=lunit,file=file_path, form='unformatted', action='write', &
        iostat=istat, status='unknown')
   
   ! Write SNODEP fields
   write(lunit) SNODEP_arrays%snoanl
   write(lunit) SNODEP_arrays%snoage
   write(lunit) SNODEP_arrays%icecon
   write(lunit) SNODEP_arrays%icemask
   write(lunit) SNODEP_arrays%iceage
   
   ! Close file
   close(lunit)

end subroutine snodep_binary
