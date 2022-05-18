!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: setClimateParmsFullnames
!  \label{setClimateParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!  07 Feb 2022: Eric Kemp/SSAI: Added CHELSAV21 precipitation.
!  13 May 2022: Eric Kemp/SSAI: Added NAFPA background precipitation.
! !INTERFACE:
subroutine setClimateParmsFullnames(n, datatype, source)

! !USES:
  use LDT_climateParmsMod
  use LDT_logMod
  use LDT_paramDataMod

  implicit none

  integer,          intent(in) :: n
  character(len=*), intent(in) :: datatype
  character(len=*), intent(in) :: source

! !ARGUMENTS:

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[datatype]
!     Climate data type
!   \item[source]
!     Climate dataset source
!   \end{description}
!EOP
!

  select case( datatype )

  case( "climppt" )
     select case( source )
     case( "PRISM" )
        LDT_climate_struc(n)%climppt%standard_name =&
             "PRISM PPT climatology fields"
     case( "WORLDCLIM" )
        LDT_climate_struc(n)%climppt%standard_name =&
             "WorldCLIM PPT climatology fields"
     case ( "CHELSAV21" )
        LDT_climate_struc(n)%climppt%standard_name =&
             "CHELSAV21 PPT climatology fields"
     case ( "NAFPA_BACK_GFS" )
        LDT_climate_struc(n)%climppt%standard_name = &
             "NAFPA_BACK_GFS PPT climatology fields"
     case( "NAFPA_BACK_GALWEM" )
        LDT_climate_struc(n)%climppt%standard_name = &
             "NAFPA_BACK_GALWEM PPT climatology fields"
     end select

  case( "mintemp" )
     select case( source )
     case( "PRISM" )
        LDT_climate_struc(n)%climtmin%standard_name =&
             "PRISM minimum Temp climatology fields"
     end select

  case( "maxtemp" )
     select case( source )
     case( "PRISM" )
        LDT_climate_struc(n)%climtmax%standard_name =&
             "PRISM maximum Temp climatology fields"
     end select

  case default
     write(ldt_logunit,*)'[ERR] Climate data type not recognized: ', &
          trim(source)
     write(ldt_logunit,*)'[ERR] Program stopping ...'
     call LDT_endrun()
  end select

end subroutine setClimateParmsFullnames
