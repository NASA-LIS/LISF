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
module LIS_routing_pluginMod
!BOP
!
! !MODULE: LIS_routing_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   routing model initialization, execution, reading and
!   writing of restart files and other relevant routing
!   model computations.
!
! !REVISION HISTORY:
!  6 May 11    Sujay Kumar  Initial Specification
! 17 Mar 21    Yeosang Yoon: Add RAPID
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_routing_plugin

contains
!BOP
! !ROUTINE: LIS_routing_plugin
!  \label{LIS_routing_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new routing scheme.
! The interface mandates that the following routines be implemented
! and registered for each of the routing scheme that is included in LIS.
!

! !INTERFACE:
subroutine LIS_routing_plugin
!EOP

   use LIS_pluginIndices

#if ( defined ROUTE_NLDAS_ROUTER )
   use NLDAS_routingMod, only : NLDAS_routingInit
#endif

#if ( defined ROUTE_HYMAP_ROUTER )
   use HYMAP_routingMod, only : HYMAP_routingInit
#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
    use HYMAP2_routingMod, only : HYMAP2_routingInit
#endif

#if ( defined ROUTE_RAPID_ROUTER )
    use RAPID_routingMod, only : RAPID_routingInit
#endif

   implicit none

#if ( defined ROUTE_NLDAS_ROUTER )
   external NLDAS_routing_readrst
   external NLDAS_routing_run
   external NLDAS_routing_output
   external NLDAS_routing_writerst
#endif

#if ( defined ROUTE_HYMAP_ROUTER )
   external HYMAP_routing_readrst
   external HYMAP_routing_run
   external HYMAP_routing_output
   external HYMAP_routing_writerst
#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
   external HYMAP2_routing_readrst
   external HYMAP2_routing_run
   external HYMAP2_routing_output
   external HYMAP2_routing_writerst
#endif

#if ( defined ROUTE_RAPID_ROUTER )
   external RAPID_routing_readrst
   external RAPID_routing_run
   external RAPID_routing_output
   external RAPID_routing_writerst
#endif

#if ( defined ROUTE_NLDAS_ROUTER )
   call registerroutinginit(trim(LIS_NLDASrouterId)//char(0),NLDAS_routingInit)
   call registerroutingreadrestart(trim(LIS_NLDASrouterId)//char(0), &
                                   NLDAS_routing_readrst)
   call registerroutingrun(trim(LIS_NLDASrouterId)//char(0),NLDAS_routing_run)
   call registerroutingoutput(trim(LIS_NLDASrouterId)//char(0), &
                              NLDAS_routing_output)
   call registerroutingwriterestart(trim(LIS_NLDASrouterId)//char(0), &
                                    NLDAS_routing_writerst)
#endif

#if ( defined ROUTE_HYMAP_ROUTER )
   call registerroutinginit(trim(LIS_HYMAProuterId)//char(0),HYMAP_routingInit)
   call registerroutingreadrestart(trim(LIS_HYMAProuterId)//char(0), &
                                   HYMAP_routing_readrst)
   call registerroutingrun(trim(LIS_HYMAProuterId)//char(0),HYMAP_routing_run)
   call registerroutingoutput(trim(LIS_HYMAProuterId)//char(0), &
                              HYMAP_routing_output)
   call registerroutingwriterestart(trim(LIS_HYMAProuterId)//char(0), &
                                    HYMAP_routing_writerst)
#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
   call registerroutinginit(trim(LIS_HYMAP2routerId)//char(0),HYMAP2_routingInit)
   call registerroutingreadrestart(trim(LIS_HYMAP2routerId)//char(0), &
                                   HYMAP2_routing_readrst)
   call registerroutingrun(trim(LIS_HYMAP2routerId)//char(0),HYMAP2_routing_run)
   call registerroutingoutput(trim(LIS_HYMAP2routerId)//char(0), &
                              HYMAP2_routing_output)
   call registerroutingwriterestart(trim(LIS_HYMAP2routerId)//char(0), &
                                    HYMAP2_routing_writerst)
#endif

#if ( defined ROUTE_RAPID_ROUTER )
   call registerroutinginit(trim(LIS_RAPIDrouterId)//char(0),RAPID_routingInit)
   call registerroutingreadrestart(trim(LIS_RAPIDrouterId)//char(0), &
                                   RAPID_routing_readrst)
   call registerroutingrun(trim(LIS_RAPIDrouterId)//char(0),RAPID_routing_run)
   call registerroutingoutput(trim(LIS_RAPIDrouterId)//char(0), &
                              RAPID_routing_output)
   call registerroutingwriterestart(trim(LIS_RAPIDrouterId)//char(0), &
                                    RAPID_routing_writerst)
#endif

end subroutine LIS_routing_plugin

end module LIS_routing_pluginMod
