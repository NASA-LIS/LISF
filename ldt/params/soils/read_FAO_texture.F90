!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_FAO_texture
! \label{read_FAO_texture}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_FAO_texture(n, array)
! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata 

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n !nest index
  real, intent(inout)  :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves FAO texture data and reprojects
!  it to the map projection used by LDT. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved texture data
!  \end{description}
!EOP

  integer :: ftn
  integer :: c,r
  real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  logical :: file_exists
! ______________________________________

  write(LDT_logunit,*) "[INFO] Reading FAO texture file: ",trim(LDT_rc%txtfile(n))

  inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Texture map ",trim(LDT_rc%txtfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn,file=LDT_rc%txtfile(n),form='unformatted',status='old',&
                     access='direct',recl=4)

  call readLISdata(n, ftn, LDT_rc%soiltext_proj, LDT_rc%soiltext_gridtransform(n), &
                   LDT_rc%soiltext_gridDesc(n,:), 1, array)

  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        if( nint(temp(c,r)) == 0 ) then
           temp(c,r) = 6    ! FAO/Zobler -- Clay-loam
        endif
        array(c,r,1) = temp(c,r)
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

  write(LDT_logunit,*) "[INFO] Done reading FAO texture file"

end subroutine read_FAO_texture
