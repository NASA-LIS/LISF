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
module LIS_irrigationmodel_pluginMod
!BOP
!
! !MODULE: LIS_irrigationmodel_pluginMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  11 Nov 2012    Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_irrigationmodel_plugin
contains
!BOP
! !ROUTINE: LIS_irrigationmodel_plugin
! \label{LIS_irrigationmodel_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
subroutine LIS_irrigationmodel_plugin
!EOP

   use LIS_pluginIndices

#if ( defined IRR_SPRINKLER )
   use sprinkler_irrigationMod
#endif

#if ( defined IRR_FLOOD )
   use flood_irrigationMod
#endif

#if ( defined IRR_DRIP )
   use drip_irrigationMod
#endif

#if ( defined IRR_SPRINKLER )
   call registerirrigationschemeinit(trim(LIS_sprinklerIrrigationId)//char(0),&
                                     sprinkler_irrigation_init)
   call registerirrigationupdate(trim(LIS_sprinklerIrrigationId)//char(0),&
                                 sprinkler_irrigation_updates)
#endif

#if ( defined IRR_FLOOD )
   call registerirrigationschemeinit(trim(LIS_floodIrrigationId)//char(0),&
                                     flood_irrigation_init)
   call registerirrigationupdate(trim(LIS_floodIrrigationId)//char(0),&
                                 flood_irrigation_updates)
#endif

#if ( defined IRR_DRIP )
   call registerirrigationschemeinit(trim(LIS_dripIrrigationId)//char(0),&
                                     drip_irrigation_init)
   call registerirrigationupdate(trim(LIS_dripIrrigationId)//char(0),&
                                 drip_irrigation_updates)
#endif

end subroutine LIS_irrigationmodel_plugin
end module LIS_irrigationmodel_pluginMod
