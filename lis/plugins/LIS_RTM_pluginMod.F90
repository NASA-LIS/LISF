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
module LIS_RTM_pluginMod
!BOP
!
! !MODULE: LIS_RTM_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   introducing new RTMs (Radiative Transfer Models) for use in LIS
!
! !REVISION HISTORY:
!  21 Mar 09    Sujay Kumar  Initial Specification
!  20 Oct 10    Yudong Tian  Added CRTM2 and CRTM2 EMonly
!  08 Feb 11    Yudong Tian  Added CMEM3
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_RTM_plugin
contains
!BOP
! !ROUTINE: LIS_RTM_plugin
! \label{LIS_RTM_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
subroutine LIS_RTM_plugin
!EOP

#if ( defined RTMS )
   use LIS_pluginIndices

#if ( defined RTMS_CRTM )
!   use CRTM_handlerMod, only : CRTM_Kmatrix_initialize, &
!       CRTM_Kmatrix_f2t, CRTM_Kmatrix_geometry, CRTM_Kmatrix_run, &
!       CRTM_Kmatrix_output
#endif

#if ( defined RTMS_CRTM2 )
!   use CRTM2_handlerMod, only : CRTM2_Forward_initialize, &
!       CRTM2_Forward_f2t, CRTM2_Forward_geometry, CRTM2_Forward_run, &
!       CRTM2_Forward_output
#endif

#if ( defined RTMS_CRTM2EM )
! CRTM 2.x Emissivity only
   use CRTM2_EMMod, only : CRTM2_EMonly_initialize, &
       CRTM2_EMonly_f2t, CRTM2_EMonly_geometry, CRTM2_EMonly_run
#endif

#if ( defined RTMS_CMEM )
! CMEM3 only
   use CMEM3_Mod, only : CMEM3_initialize, CMEM3_f2t, CMEM3_geometry, &
       CMEM3_run
#endif

#if ( defined RTMS_TAU_OMEGA )
   use TauOmegaRTM_Mod, only : TauOmegaRTM_initialize, TauOmegaRTM_f2t,&
       TauOmegaRTM_geometry, &
       TauOmegaRTM_run
#endif

#if ( defined RTMS_CRTM )
! CRTM 1.x
!   call registerrtminit(trim(LIS_crtmId)//char(0), CRTM_Kmatrix_initialize)
!   call registerrtmf2t(trim(LIS_crtmId)//char(0), CRTM_Kmatrix_f2t)
!   call registergeometry2rtm(trim(LIS_crtmId)//char(0), CRTM_Kmatrix_geometry)
!   call registerrtmrun(trim(LIS_crtmId)//char(0), CRTM_Kmatrix_run)
!   call registerrtmoutput(trim(LIS_crtmId)//char(0), CRTM_Kmatrix_output)
#endif

#if ( defined RTMS_CRTM2 )
! CRTM 2.x
!   call registerrtminit(trim(LIS_crtm2Id)//char(0), CRTM2_Forward_initialize)
!   call registerrtmf2t(trim(LIS_crtm2Id)//char(0), CRTM2_Forward_f2t)
!   call registergeometry2rtm(trim(LIS_crtm2Id)//char(0), CRTM2_Forward_geometry)
!   call registerrtmrun(trim(LIS_crtm2Id)//char(0), CRTM2_Forward_run)
!   call registerrtmoutput(trim(LIS_crtm2Id)//char(0), CRTM2_Forward_output)
#endif

#if ( defined RTMS_CRTM2EM )
! CRTM 2.x Emissivity only
   call registerrtminit(trim(LIS_crtm2EMId)//char(0),CRTM2_EMonly_initialize)
   call registerrtmf2t(trim(LIS_crtm2EMId)//char(0),CRTM2_EMonly_f2t)
   call registergeometry2rtm(trim(LIS_crtm2EMId)//char(0),CRTM2_EMonly_geometry)
   call registerrtmrun(trim(LIS_crtm2EMId)//char(0),CRTM2_EMonly_run)
#endif

#if ( defined RTMS_CMEM )
! CMEM3
   call registerrtminit(trim(LIS_cmem3Id)//char(0),CMEM3_initialize)
   call registerrtmf2t(trim(LIS_cmem3Id)//char(0),CMEM3_f2t)
   call registergeometry2rtm(trim(LIS_cmem3Id)//char(0),CMEM3_geometry)
   call registerrtmrun(trim(LIS_cmem3Id)//char(0),CMEM3_run)
#endif

#if ( defined RTMS_TAU_OMEGA )
! TauOmegaRTM
   call registerrtminit(trim(LIS_tauomegartmId)//char(0),TauOmegaRTM_initialize)
   call registerrtmf2t(trim(LIS_tauomegartmId)//char(0),TauOmegaRTM_f2t)
   call registergeometry2rtm(trim(LIS_tauomegartmId)//char(0), &
                             TauOmegaRTM_geometry)
   call registerrtmrun(trim(LIS_tauomegartmId)//char(0),TauOmegaRTM_run)
#endif
#endif
end subroutine LIS_RTM_plugin
end module LIS_RTM_pluginMod
