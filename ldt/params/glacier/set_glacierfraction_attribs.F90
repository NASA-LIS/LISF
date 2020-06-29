!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_glacierfraction_attribs
!  \label{set_glacierfraction_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_glacierfraction_attribs( n, source )

! !USES:
  use LDT_glacierFractionMod

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
      print *, "[ERR] Glacier fraction source not recognized: ",trim(source)
      print *, " Please select:   GLIMS"
      print *, " Program stopping ..."
      stop
!      call LDT_endrun

   end select

end subroutine set_glacierfraction_attribs
