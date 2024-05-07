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
module LIS_glacierrouting_pluginMod
!BOP
!
! !MODULE: LIS_glacierrouting_pluginMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
! 
  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_glacierrouting_plugin  
!EOP  

contains
!BOP
! !ROUTINE: LIS_glacierrouting_plugin
!  \label{LIS_glacierrouting_plugin}
!
! !DESCRIPTION:
!
!  
! !INTERFACE:
subroutine LIS_glacierrouting_plugin
!EOP

   use LIS_pluginIndices

#if ( ( defined ROUTE_HYMAP_ROUTER ) ||  ( defined ROUTE_HYMAP2_ROUTER ) || ( defined ROUTE_NLDAS_ROUTER ) )

#if ( defined SM_GLACIER_TEMPLATE )
   external templateGL_getrunoffs_mm
#endif

#if ( defined SM_NOAHMP_GLACIER_3_9_1_1 )
   external noahmpglacier3911_getrunoffs_mm
   external noahmpglacier3911_getrunoffs_hymap2
#endif

#if ( defined SM_GLACIER_TEMPLATE )
   call registerglacierroutinggetrunoff(trim(LIS_templateGLId)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        templateGL_getrunoffs_mm)
#endif

#if ( defined SM_NOAHMP_GLACIER_3_9_1_1 )
   call registerglacierroutinggetrunoff(trim(LIS_noahmpglacier3911Id)//"+"//&
        trim(LIS_HYMAProuterId)//char(0), &
        noahmpglacier3911_getrunoffs_mm)

   call registerglacierroutinggetrunoff(trim(LIS_noahmpglacier3911Id)//"+"//&
        trim(LIS_HYMAP2routerId)//char(0), &
        noahmpglacier3911_getrunoffs_hymap2)
#endif

#endif

end subroutine LIS_glacierrouting_plugin
end module LIS_glacierrouting_pluginMod
