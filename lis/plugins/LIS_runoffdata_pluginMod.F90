!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_runoffdata_pluginMod
!BOP
!
! !MODULE: LIS_runoffdata_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   introducing new runoffdatas (Radiative Transfer Models) for use in LIS
!
! !REVISION HISTORY:
!  6 Jan 2016    Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_runoffdata_plugin
contains
!BOP
! !ROUTINE: LIS_runoffdata_plugin
! \label{LIS_runoffdata_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
  subroutine LIS_runoffdata_plugin
!EOP

#if ( ( defined ROUTE_HYMAP2_ROUTER ) || ( defined ROUTE_HYMAP3_ROUTER ) )
    use LIS_pluginIndices
#if ( defined ROUTE_HYMAP2_ROUTER )
    use LISrunoffdataMod
#endif
#if ( defined ROUTE_HYMAP3_ROUTER )
    use HYMAP3_LISrunoffdataMod
#endif

!   use GLDAS1runoffdataMod
!   use GLDAS2runoffdataMod
!   use NLDAS2runoffdataMod
!   use MERRA2runoffdataMod
!   use ERAILandrunoffdataMod
!   use GWBMIPrunoffdataMod

    implicit none

    external :: registerinitrunoffdata
    external :: registerreadrunoffdata

#if ( defined ROUTE_HYMAP2_ROUTER )
    external readLISrunoffdata
#endif

#if ( defined ROUTE_HYMAP3_ROUTER )
    external HYMAP3_readLISrunoffdata
#endif


!   external readGLDAS1runoffdata
!   external readGLDAS2runoffdata
!   external readNLDAS2runoffdata
!   external readMERRA2runoffdata
!   external readERAILandrunoffdata
!   external readGWBMIPrunoffdata

#if ( defined ROUTE_HYMAP2_ROUTER )
    call registerinitrunoffdata(trim(LIS_LISrunoffdataId)//char(0), &
        LISrunoffdata_init)
    call registerreadrunoffdata(trim(LIS_LISrunoffdataId)//char(0), &
         readLISrunoffdata)
#endif

#if ( defined ROUTE_HYMAP3_ROUTER )
    call registerinitrunoffdata( &
         trim(LIS_HYMAP3_LISrunoffdataId)//char(0), &
         HYMAP3_LISrunoffdata_init)
    call registerreadrunoffdata( &
         trim(LIS_HYMAP3_LISrunoffdataId)//char(0), &
         HYMAP3_readLISrunoffdata)
#endif

!   call registerinitrunoffdata(trim(LIS_GLDAS1runoffdataId)//char(0), &
!        GLDAS1runoffdata_init)
!   call registerreadrunoffdata(trim(LIS_GLDAS1runoffdataId)//char(0), &
!        readGLDAS1runoffdata)

!   call registerinitrunoffdata(trim(LIS_GLDAS2runoffdataId)//char(0), &
!        GLDAS2runoffdata_init)
!   call registerreadrunoffdata(trim(LIS_GLDAS2runoffdataId)//char(0), &
!        readGLDAS2runoffdata)

!   call registerinitrunoffdata(trim(LIS_NLDAS2runoffdataId)//char(0), &
!        NLDAS2runoffdata_init)
!   call registerreadrunoffdata(trim(LIS_NLDAS2runoffdataId)//char(0), &
!        readNLDAS2runoffdata)

!   call registerinitrunoffdata(trim(LIS_MERRA2runoffdataId)//char(0), &
!        MERRA2runoffdata_init)
!   call registerreadrunoffdata(trim(LIS_MERRA2runoffdataId)//char(0), &
!        readMERRA2runoffdata)

!   call registerinitrunoffdata(trim(LIS_ERAIlandrunoffdataId)//char(0), &
!        ERAILandrunoffdata_init)
!   call registerreadrunoffdata(trim(LIS_ERAIlandrunoffdataId)//char(0), &
!        readERAILandrunoffdata)

!   call registerinitrunoffdata(trim(LIS_GWBMIPrunoffdataId)//char(0), &
!        GWBMIPrunoffdata_init)
!   call registerreadrunoffdata(trim(LIS_GWBMIPrunoffdataId)//char(0), &
!        readGWBMIPrunoffdata)
#endif
  end subroutine LIS_runoffdata_plugin

end module LIS_runoffdata_pluginMod
