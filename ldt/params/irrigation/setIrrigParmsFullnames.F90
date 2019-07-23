!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: setIrrigParmsFullnames
!  \label{setIrrigParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setIrrigParmsFullnames(n,datatype,source)

! !USES:
  use LDT_irrigationMod
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
!     Irrig data type
!   \item[source]
!     Irrig dataset source
!   \end{description}
!EOP      
!

   select case( datatype )

    case( "irrigtype" )
      select case( source )
        case( "GRIPC" )
          LDT_irrig_struc(n)%irrigtype%standard_name =&
             "GRIPC (Salmon,2013) Irrigation type (tiles)"
      end select

    case( "irrigfrac" )
      select case( source )
        case( "GRIPC" )
          LDT_irrig_struc(n)%irrigfrac%standard_name =&
              "GRIPC (Salmon,2013) Irrig gridcell fraction"
        case( "MODIS_OG" )
          LDT_irrig_struc(n)%irrigfrac%standard_name =&
              "MODIS (Ozdogan+Gutman,2008) Irrig gridcell fraction"
      end select

    case default
      print *, "[ERR] Irrig data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setIrrigParmsFullnames
