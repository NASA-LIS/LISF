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
! !ROUTINE: setAlbedoParmsFullnames
!  \label{setAlbedoParmsFullnames}
!
! !REVISION HISTORY:
!  19 Sep 2014: K. Arsenault; Initial Specification
!
! !INTERFACE:
subroutine setAlbedoParmsFullnames(n,datatype,source)

! !USES:
  use LDT_albedoMod
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
!     Albedo data type
!   \item[source]
!     Albedo dataset source
!   \end{description}
!EOP      
!

  select case( datatype )

  case( "albedo" )
     select case( source )
     case( "NCEP_LIS" )
        LDT_albedo_struc(n)%albedo%standard_name =&
             "NCEP (LIS-based) albedo clim"
     case( "NCEP_Native" )
        LDT_albedo_struc(n)%albedo%standard_name =&
             "NCEP (Native) monthly albedo clim"
     case( "NCEP_NativeQtr" )
        LDT_albedo_struc(n)%albedo%standard_name =&
             "NCEP (Native) quarterly albedo clim"
     case( "NCEP_GFS" )
        LDT_albedo_struc(n)%albedo%standard_name =&
             "NCEP 'GFS' albedo clim"
     case( "CONSTANT" )
        LDT_albedo_struc(n)%albedo%standard_name =&
             "CONSTANT albedo clim"
     end select

  case( "maxsnowalb" )
     select case( source )
     case( "NCEP_LIS" )
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "NCEP (LIS-based) max snow albedo"
     case( "NCEP_Native" )
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "NCEP (Native) max snow albedo"
     case( "MODIS" )
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "MODIS (LIS-based) max snow albedo"
     case( "NCEP_GFS" )
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "NCEP 'GFS' max snow albedo"
     case( "Barlage_Native")
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "Barlage (Native) max snow albedo"
     case( "CONSTANT" )
        LDT_albedo_struc(n)%mxsnoalb%standard_name =&
             "CONSTANT max snow albedo"
     end select

  case default
     print *, "[ERR] Albedo data type not recognized: ",trim(source)
     print *, " Program stopping ..."
     stop
  end select

end subroutine setAlbedoParmsFullnames
