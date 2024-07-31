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
module LIS_lsmirrigation_pluginMod
!BOP
!
! !MODULE: LIS_lsmirrigation_pluginMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  13 Nov 2012    Sujay Kumar  Initial Specification
!
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_lsmirrigation_plugin
!EOP

contains
!BOP
! !ROUTINE: LIS_lsmirrigation_plugin
!  \label{LIS_lsmirrigation_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
subroutine LIS_lsmirrigation_plugin
!EOP
   use LIS_pluginIndices

#if ( defined SM_NOAH_3_3 )
   external noah33_getirrigationstates
#endif

#if ( defined SM_CLSM_F2_5 )
   external clsmf25_getirrigationstates
#endif

#if ( defined SM_RUC_3_7 )
   external ruc37_getirrigationstates
#endif

#if ( defined SM_NOAHMP_3_6 )
   external noahmp36_getirrigationstates
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   external noahmp401_getirrigationstates
#endif

#if ( defined IRR_SPRINKLER )
#if ( defined SM_NOAH_3_3 )
   call registerlsmirrigationgetstates(trim(LIS_noah33Id)//"+"//&
        trim(LIS_sprinklerIrrigationId)//char(0),noah33_getirrigationstates)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsmirrigationgetstates(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_sprinklerIrrigationId)//char(0),clsmf25_getirrigationstates)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsmirrigationgetstates(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_sprinklerIrrigationId)//char(0),ruc37_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_sprinklerIrrigationId)//char(0),NoahMP36_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_sprinklerIrrigationId)//char(0),NoahMP401_getirrigationstates)
#endif
#endif

#if ( defined IRR_FLOOD )
#if ( defined SM_NOAH_3_3 )
   call registerlsmirrigationgetstates(trim(LIS_noah33Id)//"+"//&
        trim(LIS_floodIrrigationId)//char(0),noah33_getirrigationstates)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsmirrigationgetstates(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_floodIrrigationId)//char(0),clsmf25_getirrigationstates)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsmirrigationgetstates(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_floodIrrigationId)//char(0),ruc37_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_floodIrrigationId)//char(0),NoahMP36_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_floodIrrigationId)//char(0),NoahMP401_getirrigationstates)
#endif
#endif

#if ( defined IRR_DRIP )
#if ( defined SM_NOAH_3_3 )
   call registerlsmirrigationgetstates(trim(LIS_noah33Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),noah33_getirrigationstates)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsmirrigationgetstates(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),clsmf25_getirrigationstates)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsmirrigationgetstates(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),ruc37_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),NoahMP36_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),NoahMP36_getirrigationstates)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmirrigationgetstates(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_dripIrrigationId)//char(0),NoahMP401_getirrigationstates)
#endif
#endif

end subroutine LIS_lsmirrigation_plugin
end module LIS_lsmirrigation_pluginMod
