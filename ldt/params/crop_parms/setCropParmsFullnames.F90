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
! !ROUTINE: setCropParmsFullnames
!  \label{setCropParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!  27 Feb 2020: H. Beaudoing; Added optional parameters for MIRCA data
!
! !INTERFACE:
subroutine setCropParmsFullnames(n,datatype,source)

! !USES:
  use LDT_LSMCropModifier_Mod
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
!     Crop data type
!   \item[source]
!     Crop dataset source
!   \end{description}
!EOP      
!

   select case( datatype )

    case( "croptype" )
      select case( source )
        case( "UMDCROPMAP" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "UMD+CROPMAP crop types (based on Leff et al, 2004)"

        case( "Monfreda08" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "Monfreda et al (2008) crop types"

       case( "MIRCA" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "MIRCA-2000 crop types"

       case( "MIRCA52" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "MIRCA-2000 Irrigated+Rainfed crop types"

      end select

    case( "plantday" )
      select case( source )
       case( "MIRCA" )
          LDT_LSMCrop_struc(n)%plantday%standard_name =&
              "MIRCA-2000 crop planting date"

       case( "MIRCA52" )
          LDT_LSMCrop_struc(n)%plantday%standard_name =&
              "MIRCA-2000 Irrigated+Rainfed crop planting date"

      end select

    case( "harvestday" )
      select case( source )
       case( "MIRCA" )
          LDT_LSMCrop_struc(n)%harvestday%standard_name =&
              "MIRCA-2000 crop harvesting date"

       case( "MIRCA52" )
          LDT_LSMCrop_struc(n)%harvestday%standard_name =&
              "MIRCA-2000 Irrigated+Rainfed crop harvesting date"
      end select

    case( "irrigcrop" )
      select case( source )
       case( "MIRCA", "MIRCA52" )
          LDT_LSMCrop_struc(n)%irrigcrop%standard_name =&
              "MIRCA-2000 irrigation crop fraction"

      end select

    case( "rainfedcrop" )
      select case( source )
       case( "MIRCA", "MIRCA52" )
          LDT_LSMCrop_struc(n)%rainfedcrop%standard_name =&
              "MIRCA-2000 rainfed crop fraction"

      end select

    case default
      print *, "[ERR] Crop data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select


end subroutine setCropParmsFullnames
