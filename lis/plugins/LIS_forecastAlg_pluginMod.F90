!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
#include "LIS_plugins.h"
module LIS_forecastAlg_pluginMod
!BOP
!
! !MODULE: LIS_forecastAlg_pluginMod
!
! !DESCRIPTION:

!
! !REVISION HISTORY:

!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_forecastAlg_plugin

contains
!BOP
! !ROUTINE: LIS_forecastAlg_plugin
! \label{LIS_forecastAlg_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
  subroutine LIS_forecastAlg_plugin
!EOP
    use LIS_pluginIndices
#if ( defined FA_ESPBOOT )   
   use ESPboot_Mod
#endif
#if ( defined FA_ESPCONV )   
   use ESPconv_Mod
#endif

#if ( defined FA_ESPBOOT )   
   call registerforecastalginit(trim(LIS_ESPbootId)//char(0), &
        ESPboot_initialize)
   call registerforecastsampledate(trim(LIS_ESPbootId)//char(0), &
        ESPboot_sampledate)
#endif
#if ( defined FA_ESPCONV )   
   call registerforecastalginit(trim(LIS_ESPconvId)//char(0), &
        ESPconv_initialize)
   call registerforecastsampledate(trim(LIS_ESPconvId)//char(0), &
        ESPconv_sampledate)
#endif

 end subroutine LIS_forecastAlg_plugin
end module LIS_forecastAlg_pluginMod

