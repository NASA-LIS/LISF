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
      end select

    case default
      print *, "[ERR] Crop data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setCropParmsFullnames
