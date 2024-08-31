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
! !ROUTINE: set_glacierfraction_attribs
!  \label{set_glacierfraction_attribs}
!
! !REVISION HISTORY:
!  21 Jun 2020: Mahdi Navari; This code is based on the subroutine 
!               set_gfrac_attribs.F90 initially implemented by Sujay Kumar.
!  26 Jan 2021: Mahdi Navari; Modified for glacier fraction
!
! !INTERFACE:
subroutine set_glacierfraction_attribs( n, source )

! !USES:
  use LDT_glacierFractionMod
  use LDT_logMod

  implicit none

  integer,         intent(in) :: n
  character(len=*),intent(in) :: source

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[source]
!     Glacier fraction dataset source
!   \end{description}
!EOP      
!
   select case( source )

    case( "GLIMS" )
      LDT_glacierfrac_struc(n)%glacierfrac%num_bins = 1
      LDT_glacierfrac_struc(n)%glacierfrac%num_times = 1

    case default
      write (LDT_logunit, *) "[ERR] Glacier fraction source not recognized: ",trim(source)
      write (LDT_logunit, *) " Please select:   GLIMS"
      call LDT_endrun()
       
   end select

end subroutine set_glacierfraction_attribs
