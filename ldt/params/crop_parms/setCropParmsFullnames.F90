!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: setCropParmsFullnames
!  \label{setCropParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!  21 May 2019: H. Beaudoing; Added MIRCA2000 crop map
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

       case( "MIRCAIrrig" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "MIRCA-2000 Irrigated crop types"

       case( "MIRCA52" )
          LDT_LSMCrop_struc(n)%croptype%standard_name =&
              "MIRCA-2000 Irrigated+Rainfed crop types"

      end select

    case default
      print *, "[ERR] Crop data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setCropParmsFullnames
