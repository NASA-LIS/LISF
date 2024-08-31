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
! !ROUTINE: read_FAO_color
! \label{read_FAO_color}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_FAO_color(n, array)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata
  use LDT_soilsMod 

  implicit none
! !ARGUMENTS: 
  integer,   intent(in)    :: n
  real,      intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves FAO soil color data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved soil color 
!  \end{description}
!EOP
  integer :: ftn
  
  select case ( LDT_soils_struc(n)%soilclr_gridtransform )
   case( "none", "neighbor", "mode" )  ! Discrete data type
     write(LDT_logunit,*) "[INFO] Reading soil color file: ",trim(LDT_rc%iscfile(n))
  case default
     write(LDT_logunit,*) "[ERR] Since the soil color field involves discrete data values,"
     write(LDT_logunit,*) "  only 'neighbor' or 'mode' are currently supported spatial"
     write(LDT_logunit,*) "  transform types.  Please check your entries for this parameter."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=LDT_rc%iscfile(n), form='unformatted', status='old',&
            access='direct', recl=4)

  call readLISdata(n, ftn, LDT_soils_struc(n)%soilclr_proj, &
       LDT_soils_struc(n)%soilclr_gridtransform, &
       LDT_soils_struc(n)%soilclr_gridDesc(:), 1, array)  ! 1 indicates 2D layer
  
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit,*) "[INFO] Done reading soil color file."

end subroutine read_FAO_color
