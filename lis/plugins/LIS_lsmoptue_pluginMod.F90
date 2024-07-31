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
module LIS_lsmoptue_pluginMod
!BOP
!
! !MODULE: LIS_lsmoptue_pluginMod
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
  PUBLIC :: LIS_lsmoptue_plugin
!EOP

contains
!BOP
! !ROUTINE: LIS_lsmoptue_plugin
!  \label{LIS_lsmoptue_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new LSM
! in a parameter estimation mode. The interface mandates that
! a number of routines be implemented and registered for
! each of the LSM that is used in a parameter estimation mode.
!
! !INTERFACE:
subroutine LIS_lsmoptue_plugin
!EOP

#if ( ( defined OPTUE_ALG_ES )        || \
      ( defined OPTUE_ALG_LM )        || \
      ( defined OPTUE_ALG_GA )        || \
      ( defined OPTUE_ALG_SCEUA )     || \
      ( defined OPTUE_ALG_MCSIM )     || \
      ( defined OPTUE_ALG_RWMCMC )    || \
      ( defined OPTUE_ALG_DEMC )      || \
      ( defined OPTUE_ALG_DEMCZ ) )

   use LIS_pluginIndices
   use LIS_coreMod,  only : LIS_rc

#if ( defined PE_OBS_EMPTYOBS )
   use Empty_obsMod, only: Empty_getpeobspred, Empty_setupobspred  ! used for MCSIM as does not rely on obs
#endif

!    use noah271_peMod, only : noah271_setup_pedecvars
!    use noah32_peMod, only : noah32_setup_pedecvars
#if ( defined SM_NOAH_3_3 )
   use noah33_peMod, only : noah33_setup_pedecvars
#endif
#if ( defined SM_NOAH_3_6 )
   use noah36_peMod, only : noah36_setup_pedecvars
#endif
#if ( defined SM_NOAHMP_3_6 )
   use NoahMP36_peMod, only : NoahMP36_setup_pedecvars
#endif
#if ( defined SM_NOAHMP_4_0_1 )
   use NoahMP401_peMod, only : NoahMP401_setup_pedecvars
#endif

!    external noah271_f2t
!    external noah271_set_pedecvars
!    external noah271_setupobspred_wgPBMRsm
!    external noah271_getpeobspred_wgPBMRsm
!

!    external noah32_f2t
!    external noah32_set_pedecvars
!    external noah32_setupobspred_ARMObs
!    external noah32_getpeobspred_ARMObs


#if ( defined SM_NOAH_3_3 )
   external noah33_f2t
   external noah33_set_pedecvars

   external noah33_getpeobspred_LPRM_AMSREsmObs
   external noah33_setupobspred_LPRM_AMSREsmObs

   external noah33_getpeobspred_USDA_ARSsmObs
   external noah33_setupobspred_USDA_ARSsmObs

   external noah33_setupobspred_FLUXNETObs
   external noah33_getpeobspred_FLUXNETObs

!    external noah33_setupobspred_ARMObs
!    external noah33_getpeobspred_ARMObs
#endif

#if ( defined SM_NOAH_3_6 )
   external noah36_f2t
   external noah36_set_pedecvars

!    external noah33_setupobspred_ARMObs
!    external noah33_getpeobspred_ARMObs
#endif

#if ( defined SM_NOAHMP_3_6 )

   external NoahMP36_f2t
   external NoahMP36_set_pedecvars

   external NoahMP36_getpeobspred_ARSsmobs
   external NoahMP36_setupobspred_ARSsmobs

   external NoahMP36_getpeobspred_ISMNsmobs
   external NoahMP36_setupobspred_ISMNsmobs

   external NoahMP36_getpeobspred_SMAPsmobs
   external NoahMP36_setupobspred_SMAPsmobs

#endif

#if ( defined SM_NOAHMP_4_0_1)

   external NoahMP401_f2t
   external NoahMP401_set_pedecvars

   external NoahMP401_getpeobspred_UAsnowobs
   external NoahMP401_setupobspred_UAsnowobs

#endif

!    call registerlsmf2t(trim(LIS_noah271Id)//char(0), &
!                        trim(LIS_paramEstimRunId)//char(0),noah271_f2t)
!    call registersetuplsmdecisionspace(trim(LIS_noah271Id)//char(0), &
!         noah271_setup_pedecvars)
!    call registersetlsmdecisionspace(trim(LIS_noah271Id)//char(0), &
!         noah271_set_pedecvars)
!    call registersetuplsmpeobspred(trim(LIS_noah271Id)//char(0),  &
!                                   trim(LIS_wgPBMRsmId)//char(0), &
!                                   noah271_setupobspred_wgPBMRsm)
!    call registergetlsmpeobspred(trim(LIS_noah271Id)//char(0),  &
!                                 trim(LIS_wgPBMRsmId)//char(0), &
!                                 noah271_getpeobspred_wgPBMRsm)
!    call registersetuplsmpeobspred(trim(LIS_noah271Id)//char(0),  &
!                                   trim(LIS_EmptyObsId)//char(0), &
!                                   Empty_setupobspred)
!    call registergetlsmpeobspred(trim(LIS_noah271Id)//char(0),  &
!                                 trim(LIS_EmptyObsId)//char(0), &
!                                 Empty_getpeobspred)
!
!    call registerlsmf2t(trim(LIS_noah32Id)//char(0), &
!                        trim(LIS_paramEstimRunId)//char(0),noah32_f2t)
!    call registersetuplsmdecisionspace(trim(LIS_noah32Id)//char(0), &
!                                       noah32_setup_pedecvars)
!    call registersetlsmdecisionspace(trim(LIS_noah32Id)//char(0), &
!                                     noah32_set_pedecvars)
!    call registersetuplsmpeobspred(trim(LIS_noah32Id)//char(0), &
!                                   trim(LIS_ARMObsId)//char(0), &
!                                   noah32_setupobspred_ARMobs)
!    call registergetlsmpeobspred(trim(LIS_noah32Id)//char(0), &
!                                 trim(LIS_ARMObsId)//char(0), &
!                                 noah32_getpeobspred_ARMobs)
!    call registersetuplsmpeobspred(trim(LIS_noah32Id)//char(0),       &
!                                   trim(LIS_CNRS_MPDIObsId)//char(0), &
!                                   noah32_setupobspred_ARMobs)
!    call registergetlsmpeobspred(trim(LIS_noah32Id)//char(0),       &
!                                 trim(LIS_CNRS_MPDIObsId)//char(0), &
!                                 noah32_getpeobspred_ARMobs)
!    call registersetuplsmpeobspred(trim(LIS_noah32Id)//char(0),   &
!                                   trim(LIS_EmptyObsId)//char(0), &
!                                   Empty_setupobspred)
!    call registergetlsmpeobspred(trim(LIS_noah32Id)//char(0),   &
!                                 trim(LIS_EmptyObsId)//char(0), &
!                                 Empty_getpeobspred)
!

#if ( defined SM_NOAH_3_3 )
   call registerlsmf2t(trim(LIS_noah33Id)//"+"// &
                       trim(LIS_paramEstimRunId)//char(0),noah33_f2t)
   call registerlsmpesetupdecisionspace(trim(LIS_noah33Id)//char(0), &
                                        noah33_setup_pedecvars)
   call registerlsmpesetdecisionspace(trim(LIS_noah33Id)//char(0), &
                                      noah33_set_pedecvars)
!    call registersetuplsmpeobspred(trim(LIS_noah33Id)//char(0), &
!                                   trim(LIS_ARMObsId)//char(0), &
!                                   noah33_setupobspred_ARMobs)
!    call registergetlsmpeobspred(trim(LIS_noah33Id)//char(0), &
!                                 trim(LIS_ARMObsId)//char(0), &
!                                 noah33_getpeobspred_ARMobs)
   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//          &
                                  trim(LIS_CNRS_MPDIObsId)//char(0), &
                                  Empty_setupobspred)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//          &
                                trim(LIS_CNRS_MPDIObsId)//char(0), &
                                Empty_setupobspred)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//     &
                                  trim(LIS_CNRSObsId)//char(0), &
                                  Empty_setupobspred)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//     &
                                trim(LIS_CNRSObsId)//char(0), &
                                Empty_setupobspred)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//         &
                                  trim(LIS_AMSRE_SRObsId)//char(0), &
                                  Empty_setupobspred)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//         &
                                trim(LIS_AMSRE_SRObsId)//char(0), &
                                Empty_setupobspred)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//               &
                                  trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                                  noah33_setupobspred_LPRM_AMSREsmObs)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//               &
                                trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                                noah33_getpeobspred_LPRM_AMSREsmObs)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//             &
                                  trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                                  noah33_setupobspred_USDA_ARSsmObs)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//             &
                                trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                                noah33_getpeobspred_USDA_ARSsmObs)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//      &
                                  trim(LIS_EmptyObsId)//char(0), &
                                  Empty_setupobspred)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//      &
                                trim(LIS_EmptyObsId)//char(0), &
                                Empty_getpeobspred)

   call registerlsmpesetupobspred(trim(LIS_noah33Id)//"+"//          &
                                  trim(LIS_FLUXNETpeObsId)//char(0), &
                                  noah33_setupobspred_FLUXNETObs)
   call registerlsmpegetobspred(trim(LIS_noah33Id)//"+"//          &
                                trim(LIS_FLUXNETpeObsId)//char(0), &
                                noah33_getpeobspred_FLUXNETObs)

#endif

#if ( defined SM_NOAH_3_6)
   call registerlsmf2t(trim(LIS_noah36Id)//"+"// &
                       trim(LIS_paramEstimRunId)//char(0),noah36_f2t)
   call registerlsmpesetupdecisionspace(trim(LIS_noah36Id)//char(0), &
                                        noah36_setup_pedecvars)
   call registerlsmpesetdecisionspace(trim(LIS_noah36Id)//char(0), &
                                      noah36_set_pedecvars)

   call registerlsmpesetupobspred(trim(LIS_noah36Id)//"+"//      &
                                  trim(LIS_EmptyObsId)//char(0), &
                                  Empty_setupobspred)
   call registerlsmpegetobspred(trim(LIS_noah36Id)//"+"//      &
                                trim(LIS_EmptyObsId)//char(0), &
                                Empty_getpeobspred)

#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsmf2t(trim(LIS_noahmp36Id)//"+"// &
                       trim(LIS_paramEstimRunId)//char(0),NoahMP36_f2t)
   call registerlsmpesetupdecisionspace(trim(LIS_noahmp36Id)//char(0), &
                                        NoahMP36_setup_pedecvars)
   call registerlsmpesetdecisionspace(trim(LIS_noahmp36Id)//char(0), &
                                      NoahMP36_set_pedecvars)

   call registerlsmpesetupobspred(trim(LIS_noahmp36Id)//"+"//      &
                                  trim(LIS_ARSsmobsId)//char(0), &
                                  NoahMP36_setupobspred_ARSsmobs)
   call registerlsmpegetobspred(trim(LIS_noahmp36Id)//"+"//      &
                                trim(LIS_ARSsmobsId)//char(0), &
                                NoahMP36_getpeobspred_ARSsmobs)

   call registerlsmpesetupobspred(trim(LIS_noahmp36Id)//"+"//      &
                                  trim(LIS_ISMNsmobsId)//char(0), &
                                  NoahMP36_setupobspred_ISMNsmobs)
   call registerlsmpegetobspred(trim(LIS_noahmp36Id)//"+"//      &
                                trim(LIS_ISMNsmobsId)//char(0), &
                                NoahMP36_getpeobspred_ISMNsmobs)

   call registerlsmpesetupobspred(trim(LIS_noahmp36Id)//"+"//      &
                                  trim(LIS_SMAPsmobsId)//char(0), &
                                  NoahMP36_setupobspred_SMAPsmobs)
   call registerlsmpegetobspred(trim(LIS_noahmp36Id)//"+"//      &
                                trim(LIS_SMAPsmobsId)//char(0), &
                                NoahMP36_getpeobspred_SMAPsmobs)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsmf2t(trim(LIS_noahmp401Id)//"+"// &
                       trim(LIS_paramEstimRunId)//char(0),NoahMP401_f2t)
   call registerlsmpesetupdecisionspace(trim(LIS_noahmp401Id)//char(0), &
                                        NoahMP401_setup_pedecvars)
   call registerlsmpesetdecisionspace(trim(LIS_noahmp401Id)//char(0), &
                                      NoahMP401_set_pedecvars)

   call registerlsmpesetupobspred(trim(LIS_noahmp401Id)//"+"//      &
                                  trim(LIS_UAsnowobsId)//char(0), &
                                  NoahMP401_setupobspred_UAsnowobs)
   call registerlsmpegetobspred(trim(LIS_noahmp401Id)//"+"//      &
                                trim(LIS_UAsnowobsId)//char(0), &
                                NoahMP401_getpeobspred_UAsnowobs)
#endif
#endif
end subroutine LIS_lsmoptue_plugin
end module LIS_lsmoptue_pluginMod
