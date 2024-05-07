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
module LIS_RTMoptue_pluginMod
!BOP
!
! !MODULE: LIS_RTMoptue_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that
!   need to be defined to enable the use of a land surface
!   model in parameter estimation setup.
!
! !REVISION HISTORY:
!  16 Jul 09    Sujay Kumar  Initial Specification
!
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_RTMoptue_plugin
!EOP

contains
!BOP
! !ROUTINE: LIS_RTMoptue_plugin
!  \label{LIS_RTMoptue_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new RTM
! in a parameter estimation mode. The interface mandates that
! a number of routines be implemented and registered for
! each of the RTM that is used in a parameter estimation mode.
!
! !INTERFACE:
subroutine LIS_RTMoptue_plugin
!EOP
#if ( defined RTMS )
#if ( ( defined OPTUE_ALG_ES )        || \
      ( defined OPTUE_ALG_LM )        || \
      ( defined OPTUE_ALG_GA )        || \
      ( defined OPTUE_ALG_SCEUA )     || \
      ( defined OPTUE_ALG_MCSIM )     || \
      ( defined OPTUE_ALG_RWMCMC )    || \
      ( defined OPTUE_ALG_DEMC )      || \
      ( defined OPTUE_ALG_DEMCZ ) )

   use LIS_pluginIndices
#if ( defined PE_OBS_EMPTYOBS )
   use Empty_obsMod,  only : Empty_getpeobspred, Empty_setupobspred  ! used for MCSIM as does not rely on obs
#endif
   use LIS_coreMod,   only : LIS_rc

#if ( defined RTMS_CRTM2EM )
   use crtm2em_peMod, only : crtm2em_setup_pedecvars
#endif
#if ( defined RTMS_CMEM )
   use cmem3_peMod,   only : cmem3_setup_pedecvars
#endif

!    use CRTM2_EMMod
!    use CMEM3_Mod
!    use noah271_peMod, only : noah271_setup_pedecvars
!    use noah32_peMod, only : noah32_setup_pedecvars

#if ( defined RTMS_CRTM2EM )
   external CRTM2EM_qcdec
   external CRTM2EM_set_pedecvars
   external CRTM2EM_param_reset
   external CRTM2EM_setupobspred_CNRS
   external CRTM2EM_setupobspred_CNRS_MPDI
   external CRTM2EM_getpeobspred_CNRS
   external CRTM2EM_getpeobspred_CNRS_MPDI
   external CRTM2EM_setupobspred_AMSRE_SR
   external CRTM2EM_getpeobspred_AMSRE_SR
#endif

#if ( defined RTMS_CMEM )
   external CMEM3_qcdec
   external CMEM3_set_pedecvars
   external CMEM3_param_reset
!    external CMEM3_setupobspred_CNRS
!    external CMEM3_getpeobspred_CNRS
!    external CMEM3_getpeobspred_CNRS_MPDI
   external CMEM3_setupobspred_AMSRE_SR
   external CMEM3_getpeobspred_AMSRE_SR
#endif

#if ( defined RTMS_CRTM2EM )
   call registerrtmpesetupdecisionspace(trim(LIS_CRTM2EMId)//char(0), &
                                        crtm2em_setup_pedecvars)
   call registerrtmpesetdecisionspace(trim(LIS_CRTM2EMId)//char(0), &
                                      crtm2em_set_pedecvars)

   call registerrtmpesetupobspred(trim(LIS_CRTM2EMId)//"+"//    &
                                  trim(LIS_CNRSObsId)//char(0), &
                                  CRTM2EM_setupobspred_CNRS)
   call registerrtmpesetupobspred(trim(LIS_CRTM2EMId)//"+"//         &
                                  trim(LIS_CNRS_MPDIObsId)//char(0), &
                                  CRTM2EM_setupobspred_CNRS) !note same call as non-mpdi
   call registerrtmpesetupobspred(trim(LIS_CRTM2EMId)//"+"// &
                                  trim(LIS_AMSRE_SRObsId)//char(0), &
                                  CRTM2EM_setupobspred_AMSRE_SR)
   call registerrtmpesetupobspred(trim(LIS_CRTM2EMId)//"+"// &
                                  trim(LIS_EmptyObsId)//char(0), &
                                  Empty_setupobspred)

   call registerrtmpegetobspred(trim(LIS_CRTM2EMId)//"+"// &
                                trim(LIS_CNRSObsId)//char(0), &
                                CRTM2EM_getpeobspred_CNRS)
   call registerrtmpegetobspred(trim(LIS_CRTM2EMId)//"+"// &
                                trim(LIS_CNRS_MPDIObsId)//char(0), &
                                CRTM2EM_getpeobspred_CNRS_MPDI)
   call registerrtmpegetobspred(trim(LIS_CRTM2EMId)//"+"// &
                                trim(LIS_AMSRE_SRObsId)//char(0), &
                                CRTM2EM_getpeobspred_AMSRE_SR)
   call registerrtmpegetobspred(trim(LIS_CRTM2EMId)//"+"// &
                                trim(LIS_EmptyObsId)//char(0), &
                                Empty_getpeobspred)
#endif

#if ( defined RTMS_CMEM )
   call registerrtmpesetupdecisionspace(trim(LIS_CMEM3Id)//char(0), &
                                        cmem3_setup_pedecvars)
   call registerrtmpesetdecisionspace(trim(LIS_CMEM3Id)//char(0), &
                                      cmem3_set_pedecvars)

!    call registerrtmpesetupobspred(trim(LIS_CMEM3Id)//"+"// &
!                                   trim(LIS_CNRSObsId)//char(0), &
!                                   CMEM3_setupobspred_CNRS)
!    call registerrtmpesetupobspred(trim(LIS_CMEM3Id)//"+"// &
!                                   trim(LIS_CNRS_MPDIObsId)//char(0), &
!                                   CMEM3_setupobspred_CNRS) !note same call as non-mpdi
   call registerrtmpesetupobspred(trim(LIS_CMEM3Id)//"+"// &
                                  trim(LIS_AMSRE_SRObsId)//char(0), &
                                  CMEM3_setupobspred_AMSRE_SR)
   call registerrtmpesetupobspred(trim(LIS_CMEM3Id)//"+"// &
                                  trim(LIS_EmptyObsId)//char(0), &
                                  Empty_setupobspred)

!    call registerrtmpegetobspred(trim(LIS_CMEM3Id)//"+"// &
!                                 trim(LIS_CNRSObsId)//char(0), &
!                                 CMEM3_getpeobspred_CNRS)
!    call registerrtmpegetobspred(trim(LIS_CMEM3Id)//"+"// &
!                                 trim(LIS_CNRS_MPDIObsId)//char(0), &
!                                 CMEM3_getpeobspred_CNRS_MPDI)
   call registerrtmpegetobspred(trim(LIS_CMEM3Id)//"+"// & 
                                trim(LIS_AMSRE_SRObsId)//char(0), &
                                CMEM3_getpeobspred_AMSRE_SR)
   call registerrtmpegetobspred(trim(LIS_CMEM3Id)//"+"// &
                                trim(LIS_EmptyObsId)//char(0), &
                                Empty_getpeobspred)
#endif
#endif
#endif
end subroutine LIS_RTMoptue_plugin
end module LIS_RTMoptue_pluginMod
