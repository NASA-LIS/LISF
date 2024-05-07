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
! !ROUTINE: setTopoParmsFullnames
!  \label{setTopoParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setTopoParmsFullnames(n,datatype,source)

! !USES:
  use LDT_topoMod
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
!     Topographic data type
!   \item[source]
!     Topographic dataset source
!   \end{description}
!EOP      
!

   select case( datatype )

    case( "elevation" )
      select case( source )
        case( "GTOPO30_LIS" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
             "GTOPO30 (LIS-based) elevation"
        case( "GTOPO30_Native" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "GTOPO30 (Native) elevation"
        case( "SRTM_LIS" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "SRTM (LIS-based) elevation"
        case( "SRTM_Native" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "SRTM (Native) elevation"
        case( "GTOPO30_GFS" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "GTOPO30 'GFS' elevation"
        case( "CONSTANT" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "CONSTANT elevation"
        case( "MERIT_1K" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "MERIT '1K' elevation"
        case( "NED_SM" )
          LDT_LSMparam_struc(n)%elevation%standard_name =&
              "SnowModel-based NED elevation"
      end select

    case( "slope" )
      select case( source )
        case( "GTOPO30_LIS" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
             "GTOPO30 (LIS-based) slope"
        case( "GTOPO30_Native" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "GTOPO30 (Native) slope"
        case( "SRTM_LIS" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "SRTM (LIS-based) slope"
        case( "SRTM_Native" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "SRTM (Native) slope"
        case( "GTOPO30_GFS" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "GTOPO30 'GFS' slope"
        case( "CONSTANT" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "CONSTANT slope"
        case( "MERIT_1K" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "MERIT '1K' slope"
        case( "NED_SM" )
          LDT_LSMparam_struc(n)%slope%standard_name =&
              "SnowModel-based NED-derived slope"
      end select

    case( "aspect" )
      select case( source )
        case( "GTOPO30_LIS" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
             "GTOPO30 (LIS-based) aspect"
        case( "GTOPO30_Native" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "GTOPO30 (Native) aspect"
        case( "SRTM_LIS" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "SRTM (LIS-based) aspect"
        case( "SRTM_Native" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "SRTM (Native) aspect"
        case( "GTOPO30_GFS" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "GTOPO30 'GFS' aspect"
        case( "CONSTANT" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "CONSTANT aspect"
        case( "MERIT_1K" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "MERIT '1K' aspect"
        case( "NED_SM" )
          LDT_LSMparam_struc(n)%aspect%standard_name =&
              "SnowModel-based NED-derived aspect"
      end select

    case( "curvature" )
      select case( source )
        case( "GTOPO30_LIS" )
          LDT_LSMparam_struc(n)%curvature%standard_name =&
             "GTOPO30 (LIS-based) curvature"
        case( "NED_SM" )
          LDT_LSMparam_struc(n)%curvature%standard_name =&
              "SnowModel-based NED-derived curvature"
      end select

    case default
      print *, "[ERR] Topographic data type not recognized: ",trim(source)
      print *, " Program stopping ..."
      stop
   end select

end subroutine setTopoParmsFullnames
