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
! !ROUTINE: setLAISAIParmsFullnames
!  \label{setLAISAIParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setLAISAIParmsFullnames(n,datatype,source)

! !USES:
  use LDT_laisaiMod
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
!     LAI/SAI data type
!   \item[source]
!     LAI/SAI dataset source
!   \end{description}
!EOP      
!

   select case( datatype )

    case( "lai" )
      select case( source )
        case( "AVHRR" )
          LDT_laisai_struc(n)%lai%standard_name =&
             "AVHRR LAI climatology"
        case( "MODIS" )
          LDT_laisai_struc(n)%lai%standard_name =&
              "MODIS  monthly LAI climatology"
        case( "CLSMF2.5" )
          LDT_laisai_struc(n)%lai%standard_name =&
             "CLSMF2.5 LAI climatology"
        case( "CONSTANT" )
          LDT_laisai_struc(n)%lai%standard_name =&
             "CONSTANT LAI climatology"
      end select

    case( "sai" )
      select case( source )
        case( "AVHRR" )
          LDT_laisai_struc(n)%sai%standard_name =&
             "AVHRR SAI climatology"
        case( "MODIS" )
          LDT_laisai_struc(n)%lai%standard_name =&
             "MODIS SAI climatology"
        case( "CONSTANT" )
          LDT_laisai_struc(n)%lai%standard_name =&
             "CONSTANT SAI climatology"
      end select

    case( "laimin" )
      select case( source )
        case( "AVHRR" )
          LDT_laisai_struc(n)%laimin%standard_name =&
             "AVHRR min LAI"
        case( "MODIS" )
          LDT_laisai_struc(n)%laimin%standard_name =&
             "MODIS min LAI"
        case( "CLSMF2.5" )
          LDT_laisai_struc(n)%laimin%standard_name =&
             "CLSMF2.5 min LAI"
        case( "CONSTANT" )
          LDT_laisai_struc(n)%laimin%standard_name =&
             "CONSTANT min LAI"
      end select

    case( "laimax" )
      select case( source )
        case( "AVHRR" )
          LDT_laisai_struc(n)%laimax%standard_name =&
             "AVHRR max LAI"
        case( "MODIS" )
          LDT_laisai_struc(n)%laimax%standard_name =&
             "MODIS max LAI"
        case( "CLSMF2.5" )
          LDT_laisai_struc(n)%laimax%standard_name =&
             "CLSMF2.5 max LAI"
        case( "CONSTANT" )
          LDT_laisai_struc(n)%laimax%standard_name =&
             "CONSTANT max LAI"
      end select

    case default
      print *, "[ERR] LAI/SAI data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setLAISAIParmsFullnames
