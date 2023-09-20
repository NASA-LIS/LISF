!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_lsmrouting_pluginMod
!BOP
!
! !MODULE: LIS_lsmrouting_pluginMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  16 Jul 09    Sujay Kumar  Initial Specification
!  01 Jun 17    Augusto Getirana: Add HyMAP2
!  12 Sep 19    Augusto Getirana: Add 2-way coupling
!  17 Mar 21    Yeosang Yoon: Add RAPID
! 
  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_lsmrouting_plugin  
!EOP  

contains
!BOP
! !ROUTINE: LIS_lsmrouting_plugin
!  \label{LIS_lsmrouting_plugin}
!
! !DESCRIPTION:
!
!  
! !INTERFACE:
subroutine LIS_lsmrouting_plugin
!EOP

   use LIS_pluginIndices
   use LIS_routingMod, only : LIS_runoff_state

#if ( ( defined ROUTE_HYMAP_ROUTER ) ||  ( defined ROUTE_HYMAP2_ROUTER ) || ( defined ROUTE_NLDAS_ROUTER ) || ( defined ROUTE_RAPID_ROUTER ))

#if ( defined SM_LSM_TEMPLATE )
   external template_getrunoffs
#endif

#if ( defined SM_NOAH_2_7_1 )
   external noah271_getrunoffs
   external noah271_getrunoffs_mm
#endif

#if ( defined SM_NOAH_3_2 )
   external noah32_getrunoffs
   external noah32_getrunoffs_mm
#endif

#if ( defined SM_NOAH_3_3 )
   external noah33_getrunoffs
   external noah33_getrunoffs_et
   external noah33_getrunoffs_mm
#endif

#if ( defined SM_NOAH_3_6 )
   external noah36_getrunoffs
   external noah36_getrunoffs_mm
   external noah36_getrunoffs_hymap2
   external noah36_getsws_hymap2
#endif

#if ( defined SM_NOAH_3_9 )
   external noah39_getrunoffs
   external noah39_getrunoffs_mm
   external noah39_getrunoffs_rapid
   external noah39_getrunoffs_hymap2
   external noah39_getsws_hymap2
#endif

#if ( defined SM_NOAHMP_3_6 )
   external noahmp36_getrunoffs
   external noahmp36_getrunoffs_mm
   external noahmp36_getrunoffs_hymap2
   external noahmp36_getsws_hymap2
   external noahmp36_getrunoffs_rapid
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   external noahmp401_getrunoffs
   external noahmp401_getrunoffs_mm
   external noahmp401_getrunoffs_hymap2
   external noahmp401_getsws_hymap2
   external noahmp401_getrunoffs_rapid   
#endif

#if ( defined SM_RUC_3_7 )
   external ruc37_getrunoffs
   external ruc37_getrunoffs_mm
#endif

#if ( defined SM_CLSM_F2_5 )
   external clsmf25_getrunoffs
   external clsmf25_getrunoffs_mm
   external clsmf25_getrunoffs_hymap2
   external clsmf25_getsws_hymap2
#endif

#if ( defined SM_VIC_4_1_2 )
   external vic412_getrunoffs
   external vic412_getrunoffs_mm
#endif

#if ( defined SM_JULES_5_0 )
   external jules50_getrunoffs_mm
   external jules50_getrunoffs_rapid
   external jules50_getrunoffs_hymap2
   external jules50_getsws_hymap2
#endif

#if ( defined SM_JULES_5_2 )
   external jules52_getrunoffs_mm
#endif

#if ( defined SM_JULES_5_3 )
   external jules53_getrunoffs_mm
#endif

#if ( defined SM_JULES_5_X )
   external jules5x_getrunoffs_mm
   external jules5x_getrunoffs_rapid
#endif

#if ( defined ROUTE_HYMAP_ROUTER )
#if ( defined SM_LSM_TEMPLATE )
   call registerlsmroutinggetrunoff(trim(LIS_templateLSMId)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        template_getrunoffs)
#endif
#if ( defined SM_NOAH_2_7_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noah271Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noah271_getrunoffs_mm)
#endif


#if ( defined SM_NOAH_3_2 )
   call registerlsmroutinggetrunoff(trim(LIS_noah32Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noah32_getrunoffs_mm)
#endif

#if ( defined SM_NOAH_3_3 )
   call registerlsmroutinggetrunoff(trim(LIS_noah33Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noah33_getrunoffs_mm)
#endif

#if ( defined SM_NOAH_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noah36Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noah36_getrunoffs_mm)
#endif

#if ( defined SM_NOAH_3_9 )
   call registerlsmroutinggetrunoff(trim(LIS_noah39Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noah39_getrunoffs_mm)
!   call registerlsmroutinggetsws(trim(LIS_noah39Id)//"+"//&
!        trim(LIS_HYMAP2routerId)//char(0), &
!        noah39_getsws_hymap2)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noahmp36_getrunoffs_mm)
   call registerlsmroutinggetsws(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmp36_getsws_hymap2)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noahmp401_getrunoffs_mm)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsmroutinggetrunoff(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        ruc37_getrunoffs_mm)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsmroutinggetrunoff(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        clsmf25_getrunoffs_mm)
#endif

#if ( defined SM_VIC_4_1_2 )
   call registerlsmroutinggetrunoff(trim(LIS_vic412Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        vic412_getrunoffs_mm)
#endif

#if ( defined SM_JULES_5_0 )
   call registerlsmroutinggetrunoff(trim(LIS_jules50Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        jules50_getrunoffs_mm)
#endif

#if ( defined SM_JULES_5_2 )
   call registerlsmroutinggetrunoff(trim(LIS_jules52Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        jules52_getrunoffs_mm)
#endif

#if ( defined SM_JULES_5_3 )
   call registerlsmroutinggetrunoff(trim(LIS_jules53Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        jules53_getrunoffs_mm)
#endif

#if ( defined SM_JULES_5_X )
   call registerlsmroutinggetrunoff(trim(LIS_jules5xId)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        jules5x_getrunoffs_mm)
#endif

#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
#if ( defined SM_LSM_TEMPLATE )
   call registerlsmroutinggetrunoff(trim(LIS_templateLSMId)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        template_getrunoffs)
#endif

#if ( defined SM_NOAH_3_2 )
   call registerlsmroutinggetrunoff(trim(LIS_noah32Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah32_getrunoffs_mm)
#endif

#if ( defined SM_NOAH_3_3 )
   call registerlsmroutinggetrunoff(trim(LIS_noah33Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah33_getrunoffs_et)
#endif

#if ( defined SM_NOAH_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noah36Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah36_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_noah36Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah36_getsws_hymap2) 
#endif

#if ( defined SM_NOAH_3_9 )
   call registerlsmroutinggetrunoff(trim(LIS_noah39Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah39_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_noah39Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noah39_getsws_hymap2)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmp36_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmp36_getsws_hymap2)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmp401_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmp401_getsws_hymap2)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsmroutinggetrunoff(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        ruc37_getrunoffs_mm)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsmroutinggetrunoff(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        clsmf25_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        clsmf25_getsws_hymap2)
#endif

#if ( defined SM_VIC_4_1_2 )
   call registerlsmroutinggetrunoff(trim(LIS_vic412Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        vic412_getrunoffs_mm)
#endif

#if ( defined SM_JULES_5_0 )
   call registerlsmroutinggetrunoff(trim(LIS_jules50Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        jules50_getrunoffs_hymap2)
   call registerlsmroutinggetsws(trim(LIS_jules50Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        jules50_getsws_hymap2)
#endif
#endif

#if ( defined ROUTE_NLDAS_ROUTER )
#if ( defined SM_NOAH_2_7_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noah271Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noah271_getrunoffs)
#endif
   
#if ( defined SM_NOAH_3_2 )
   call registerlsmroutinggetrunoff(trim(LIS_noah32Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noah32_getrunoffs)
#endif
   
#if ( defined SM_NOAH_3_3 )
   call registerlsmroutinggetrunoff(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noah33_getrunoffs)
#endif
   
#if ( defined SM_NOAH_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noah36_getrunoffs)
#endif

#if ( defined SM_NOAH_3_9 )
   call registerlsmroutinggetrunoff(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noah39_getrunoffs)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noahmp36_getrunoffs)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        noahmp401_getrunoffs)
#endif
   
#if ( defined SM_RUC_3_7 )
   call registerlsmroutinggetrunoff(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        ruc37_getrunoffs)
#endif
   
#if ( defined SM_CLSM_F2_5 )
   call registerlsmroutinggetrunoff(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        clsmf25_getrunoffs)
#endif
   
#if ( defined SM_VIC_4_1_2 )
   call registerlsmroutinggetrunoff(trim(LIS_vic412Id)//"+"//&
        trim(LIS_NLDASrouterId)//char(0), &
        vic412_getrunoffs)
#endif
#endif

#if ( defined ROUTE_RAPID_ROUTER )
#if ( defined SM_NOAH_3_9 )
   call registerlsmroutinggetrunoff(trim(LIS_noah39Id)//"+"//&
        trim(LIS_RAPIDrouterId)//char(0), &
        noah39_getrunoffs_rapid)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_RAPIDrouterId)//char(0), &
        noahmp36_getrunoffs_rapid)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmroutinggetrunoff(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_RAPIDrouterId)//char(0), &
        noahmp401_getrunoffs_rapid)
#endif

#if ( defined SM_JULES_5_0 )
   call registerlsmroutinggetrunoff(trim(LIS_jules50Id)//"+"//&
        trim(LIS_RAPIDrouterId)//char(0), &
        jules50_getrunoffs_rapid)
#endif

#if ( defined SM_JULES_5_X )
   call registerlsmroutinggetrunoff(trim(LIS_jules5xId)//"+"//&
        trim(LIS_RAPIDrouterId)//char(0), &
        jules5x_getrunoffs_rapid)
#endif
#endif
#endif
   
end subroutine LIS_lsmrouting_plugin
end module LIS_lsmrouting_pluginMod
