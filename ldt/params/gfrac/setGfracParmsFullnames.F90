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
! !ROUTINE: setGfracParmsFullnames
!  \label{setGfracParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setGfracParmsFullnames(n,datatype,source)

! !USES:
  use LDT_gfracMod
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
!     Gfrac data type
!   \item[source]
!     Gfrac dataset source
!   \end{description}
!EOP      
!

   select case( datatype )

    case( "gfrac" )
      select case( source )
        case( "NCEP_LIS" )
          LDT_gfrac_struc(n)%gfrac%standard_name =&
             "NCEP (LIS-based) greenness fraction climatology"
        case( "NCEP_Native" )
          LDT_gfrac_struc(n)%gfrac%standard_name =&
              "NCEP (Native) monthly greenness fraction climatology"
        case( "SACHTET.3.5.6" )
          LDT_gfrac_struc(n)%gfrac%standard_name =&
             "SACHTET 3.5.6 greenness fraction climatology"
        case( "CLSMF2.5" )
          LDT_gfrac_struc(n)%gfrac%standard_name =&
             "CLSMF2.5 greenness fraction climatology"
        case( "CONSTANT" )
          LDT_gfrac_struc(n)%gfrac%standard_name =&
             "CONSTANT greenness fraction climatology"
      end select

    case( "shdmin" )
      select case( source )
        case( "NCEP_LIS" )
          LDT_gfrac_struc(n)%shdmin%standard_name =&
             "NCEP (LIS-based) min greenness"
        case( "NCEP_Native" )
          LDT_gfrac_struc(n)%shdmin%standard_name =&
             "NCEP (Native) min greenness"
        case( "CONSTANT" )
          LDT_gfrac_struc(n)%shdmin%standard_name =&
             "CONSTANT min greenness"
      end select

    case( "shdmax" )
      select case( source )
        case( "NCEP_LIS" )
          LDT_gfrac_struc(n)%shdmax%standard_name =&
             "NCEP (LIS-based) max greenness"
        case( "NCEP_Native" )
          LDT_gfrac_struc(n)%shdmax%standard_name =&
             "NCEP (Native) max greenness"
        case( "CONSTANT" )
          LDT_gfrac_struc(n)%shdmax%standard_name =&
             "CONSTANT max greenness"
      end select

    case default
      print *, "[ERR] Gfrac data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setGfracParmsFullnames
