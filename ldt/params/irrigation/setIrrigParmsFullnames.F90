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
! !ROUTINE: setIrrigParmsFullnames
!  \label{setIrrigParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!  11 Apr 2021: Wanshu Nie; add support for reading irrigation groundwater ratio
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

        case( "UserDerived" )
          LDT_irrig_struc(n)%irrigfrac%standard_name =&
              "User Derived Irrig gridcell fraction"
      end select

    case( "irriggwratio" )
      select case ( source )
        case( "USGS_Native" )
          LDT_irrig_struc(n)%irriggwratio%standard_name =&
             "USGS groundwater irrigation ratio (0.125 deg gridcell) "
      end select

    case default
      write(LDT_logunit,*) "[ERR] Irrig data type not recognized: ",trim(source)
      write(LDT_logunit,*) " Program stopping ..."
      stop
   end select

end subroutine setIrrigParmsFullnames
