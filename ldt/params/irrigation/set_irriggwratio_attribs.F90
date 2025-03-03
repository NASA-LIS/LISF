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
! !ROUTINE: set_irriggwratio_attribs
!  \label{set_irriggwratio_attribs}
!
! !REVISION HISTORY:
!  2 Nov 2018: Wanshu Nie; Initial Specification
!
! !INTERFACE:
subroutine set_irriggwratio_attribs( n, source )

! !USES:
  use LDT_logMod, only : LDT_logunit, LDT_endrun
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
!     Irrigation groundwater ratio dataset source
!   \end{description}
!EOP      
!
   select case( source )

    case( "USGS_Native" )
      LDT_irrig_struc(n)%irriggwratio%num_bins = 1
      LDT_irrig_struc(n)%irriggwratio%num_times = 1

    case default
      write(LDT_logunit,*) '[ERR] Groundwater irrigation ratio source not recognized: ',trim(source)
      write(LDT_logunit,*) '[ERR] Please select:   USGS_Native'
      write(LDT_logunit,*) 'Program stopping ...'
      
      call LDT_endrun

   end select

end subroutine set_irriggwratio_attribs
