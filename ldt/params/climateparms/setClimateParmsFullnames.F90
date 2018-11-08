!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: setClimateParmsFullnames
!  \label{setClimateParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setClimateParmsFullnames(n,datatype,source)

! !USES:
  use LDT_climateParmsMod
  use LDT_logMod
  use LDT_paramDataMod

  implicit none
  integer,         intent(in) :: n
  character(len=*),intent(in) :: datatype
  character(len=*),intent(in) :: source

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
      print *, "[ERR] Climate data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setClimateParmsFullnames
