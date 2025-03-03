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
! !ROUTINE: set_gfrac_attribs
!  \label{set_gfrac_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_gfrac_attribs(n,source)

! !USES:
  use LDT_gfracMod

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
!     Greenness fraction dataset source
!   \end{description}
!EOP      
!

   select case( source )

    case( "NCEP_LIS", "NCEP_Native", "CLSMF2.5", "SACHTET.3.5.6" )
      LDT_gfrac_struc(n)%gfrac%num_bins = 1
      LDT_gfrac_struc(n)%gfrac%num_times = 12

    case( "CONSTANT" )
      LDT_gfrac_struc(n)%gfrac%num_bins = 1
      LDT_gfrac_struc(n)%gfrac%num_times = 12

    case default
      write(*,*) "[ERR] Greenness fraction source not recognized: ",trim(source)
      write(*,*) " Please select: NCEP_LIS, NCEP_Native, CLSMF2.5, "
      write(*,*) "                SACHTET.3.5.6, or CONSTANT"
      write(*,*) " Program stopping ..."
      stop
!      call LDT_endrun
   end select

end subroutine set_gfrac_attribs
