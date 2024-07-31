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
! !ROUTINE: set_irrigation_attribs
!  \label{set_irrigation_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_irrigation_attribs( n, source )

! !USES:
  use LDT_irrigationMod

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
!     Irrigation dataset source
!   \end{description}
!EOP      
!
   select case( source )

    case( "GRIPC" )
      LDT_irrig_struc(n)%irrigtype%num_bins = 4
      LDT_irrig_struc(n)%irrigtype%num_times = 1

    case default
      print *, "[ERR] Irrigation type source not recognized: ",trim(source)
      print *, " Please select:   GRIPC"
      print *, " Program stopping ..."
      stop
!      call LDT_endrun

   end select

end subroutine set_irrigation_attribs
