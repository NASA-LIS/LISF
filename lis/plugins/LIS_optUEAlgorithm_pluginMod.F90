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
module LIS_optUEAlgorithm_pluginMod
!BOP
!
! !MODULE: LIS_optUEAlgorithm_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining optimization /uncertainty estimation algorithm algorithms
!   The user defined functions are incorporated into
!   the appropriate registry to be later invoked through generic calls.
!
! !REVISION HISTORY:
!  2 Feb 2008;   Sujay Kumar  Initial Specification
!
!EOP
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_optUEAlgorithm_plugin

contains
!BOP
! !ROUTINE: LIS_optUEAlgorithm_plugin
! \label{LIS_optUEAlgorithm_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new optimization
!  /uncertainty estimation scheme.
!
! !INTERFACE:
subroutine LIS_optUEAlgorithm_plugin
!EOP
! !USES:

#if ( ( defined OPTUE_ALG_ES )        || \
      ( defined OPTUE_ALG_LM )        || \
      ( defined OPTUE_ALG_GA )        || \
      ( defined OPTUE_ALG_SCEUA )     || \
      ( defined OPTUE_ALG_MCSIM )     || \
      ( defined OPTUE_ALG_RWMCMC )    || \
      ( defined OPTUE_ALG_DEMC )      || \
      ( defined OPTUE_ALG_DEMCZ ) )

   use LIS_pluginIndices

#if ( defined OPTUE_ALG_ES ) 
   use EnumeratedSearch, only : ESOpt_init, ESOpt_setup, ESOpt_run,           &
                                ESOpt_checkConvergence, ES_getdecspacevalues, &
                                ES_setdecSpaceValues, ES_getNparam
#endif

#if ( defined OPTUE_ALG_LM )
   use LevenbergMarquardt, only : LMOpt_init, LMOpt_setup, LMOpt_run, &
                                  LMOpt_checkConvergence,             &
                                  LMOpt_getdecspaceValues,            &
                                  LMOpt_getNparam, LMOpt_readrestart
#endif

#if ( defined OPTUE_ALG_GA )
   use GeneticAlgorithm, only : GAOpt_init, GAOpt_setup, GAOpt_run, &
                                GAOpt_checkConvergence,             &
                                GAOpt_getdecspaceValues,            &
                                GAOpt_getNparam, GAOpt_readrestart
#endif

#if ( defined OPTUE_ALG_SCEUA )
   use ShuffledComplexEvolution, only : SCEUAOpt_init, SCEUAOpt_setup,  &
                                        SCEUAOpt_run,                   &
                                        SCEUAOpt_checkOptStopCriterion, &
                                        SCEUA_getdecspaceValues,        &
                                        SCEUA_setdecSpaceValues, SCEUA_getNparam
#endif

#if ( defined OPTUE_ALG_MCSIM )
   use MCSIMAlgorithm, only : MCSIM_init, MCSIM_run, MCSIM_checkconvergence, &
                              MCSIM_getNparam, MCSIM_getdecspacevalues,      &
                              MCSIM_setup, MCSIM_readrestart
#endif

#if ( defined OPTUE_ALG_RWMCMC )
   use RWMCMCAlgorithm, only  : RWMCMC_init, RWMCMC_run,                    &
                                RWMCMC_checkConvergence,                    &
                                RWMCMC_getNparam, RWMCMC_getdecspacevalues, &
                                RWMCMC_setup
#endif

#if ( defined OPTUE_ALG_DEMC )
   use DEMCAlgorithm, only  : DEMC_init, DEMC_run, DEMC_checkConvergence, &
                              DEMC_getNparam, DEMC_getdecspacevalues,     &
                              DEMC_setup, DEMC_readrestart
#endif

#if ( defined OPTUE_ALG_DEMCZ )
   use DEMCzAlgorithm, only  : DEMCz_init, DEMCz_run, DEMCz_checkConvergence, &
                               DEMCz_getNparam, DEMCz_getdecspacevalues,      &
                               DEMCz_setup, DEMCz_readrestart
#endif

#if ( defined OPTUE_ALG_ES ) 
   call registeroptuealginit(trim(LIS_ESOptId)//char(0),ESOpt_init)
   call registeroptuealgsetup(trim(LIS_ESOptId)//char(0),ESOpt_setup)
   call registeroptuealgrun(trim(LIS_ESOptId)//char(0),ESOpt_run)
   call registeroptueconvergencecheck(trim(LIS_ESOptId)//char(0), &
                                      ESOpt_checkConvergence)
   call registeroptuegetdecisionspace(trim(LIS_ESoptId)//char(0), &
                                      ES_getdecSpaceValues)
   call registeroptuesetdecisionspace(trim(LIS_ESoptId)//char(0), &
                                      ES_setdecSpaceValues)
   call registeroptuegetnparam(trim(LIS_ESoptId)//char(0), ES_getNparam)
#endif

#if ( defined OPTUE_ALG_LM )
   call registeroptuealginit(trim(LIS_LMOptId)//char(0),LMOpt_init)
   call registeroptuealgsetup(trim(LIS_LMOptId)//char(0),LMOpt_setup)
   call registeroptuealgrun(trim(LIS_LMOptId)//char(0),LMOpt_run)
   call registeroptueconvergencecheck(trim(LIS_LMOptId)//char(0), &
                                      LMOpt_checkConvergence)
   call registeroptuegetdecisionspace(trim(LIS_LMoptId)//char(0), &
                                      LMOpt_getdecSpaceValues)
! get rid of setdecspace for all in LIS7    call registeroptuesetdecisionspace(trim(LIS_LMoptId)//char(0), LMOpt_setdecSpaceValues)
   call registeroptuegetnparam(trim(LIS_LMoptId)//char(0), LMOpt_getNparam)
   call registeroptuereadrestart(trim(LIS_LMoptId)//char(0), LMOpt_readrestart)
#endif

#if ( defined OPTUE_ALG_GA )
   call registeroptuealginit(trim(LIS_GAOptId)//char(0),GAOpt_init)
   call registeroptuealgsetup(trim(LIS_GAOptId)//char(0),GAOpt_setup)
   call registeroptuealgrun(trim(LIS_GAOptId)//char(0),GAOpt_run)
   call registeroptueconvergencecheck(trim(LIS_GAOptId)//char(0), &
                                      GAOpt_checkConvergence)
   call registeroptuegetdecisionspace(trim(LIS_GAoptId)//char(0), &
                                      GAOpt_getdecSpaceValues)
   call registeroptuegetnparam(trim(LIS_GAoptId)//char(0), GAOpt_getNparam)
   call registeroptuereadrestart(trim(LIS_GAoptId)//char(0), GAOpt_readrestart)
#endif

#if ( defined OPTUE_ALG_SCEUA )
   call registeroptuealginit(trim(LIS_SCEUAOptId)//char(0),SCEUAOpt_init)
   call registeroptuealgsetup(trim(LIS_SCEUAOptId)//char(0),SCEUAOpt_setup)
   call registeroptuealgrun(trim(LIS_SCEUAOptId)//char(0),SCEUAOpt_run)
   call registeroptueconvergencecheck(trim(LIS_SCEUAOptId)//char(0), &
                                      SCEUAOpt_checkOptStopCriterion)
   call registeroptuegetdecisionspace(trim(LIS_SCEUAoptId)//char(0), &
                                      SCEUA_getdecSpaceValues)
   call registeroptuesetdecisionspace(trim(LIS_SCEUAoptId)//char(0), &
                                      SCEUA_setdecSpaceValues)
   call registeroptuegetnparam(trim(LIS_SCEUAoptId)//char(0), SCEUA_getNparam)
#endif

#if ( defined OPTUE_ALG_MCSIM )
   call registeroptuealginit(trim(LIS_MCSIMId)//char(0), MCSIM_init)
   call registeroptuealgsetup(trim(LIS_MCSIMId)//char(0),MCSIM_setup)
   call registeroptuealgrun(trim(LIS_MCSIMId)//char(0), MCSIM_run)
   call registeroptueconvergencecheck(trim(LIS_MCSIMId)//char(0), &
                                      MCSIM_checkconvergence)
   call registeroptuegetdecisionspace(trim(LIS_MCSIMId)//char(0), &
                                      MCSIM_getdecspacevalues)
   call registeroptuegetnparam(trim(LIS_MCSIMId)//char(0), MCSIM_getNparam)
   call registeroptuereadrestart(trim(LIS_MCSIMId)//char(0), MCSIM_readrestart)
#endif

#if ( defined OPTUE_ALG_RWMCMC )
   call registeroptuealginit(trim(LIS_RWMCMCId)//char(0), RWMCMC_init)
   call registeroptuealgsetup(trim(LIS_RWMCMCId)//char(0),RWMCMC_setup)
   call registeroptuealgrun(trim(LIS_RWMCMCId)//char(0), RWMCMC_run)
   call registeroptueconvergencecheck(trim(LIS_RWMCMCId)//char(0), &
                                      RWMCMC_checkconvergence)
   call registeroptuegetdecisionspace(trim(LIS_RWMCMCId)//char(0), &
                                      RWMCMC_getdecspacevalues)
!    call registeroptuesetdecisionspace(trim(LIS_RWMCMCId)//char(0), RWMCMC_setdecspacevalues)
   call registeroptuegetnparam(trim(LIS_RWMCMCId)//char(0), RWMCMC_getNparam)
#endif

#if ( defined OPTUE_ALG_DEMC )
   call registeroptuealginit(trim(LIS_DEMCId)//char(0), DEMC_init)
   call registeroptuealgsetup(trim(LIS_DEMCId)//char(0),DEMC_setup)
   call registeroptuealgrun(trim(LIS_DEMCId)//char(0), DEMC_run)
   call registeroptueconvergencecheck(trim(LIS_DEMCId)//char(0), &
                                      DEMC_checkconvergence)
   call registeroptuegetdecisionspace(trim(LIS_DEMCId)//char(0), &
                                      DEMC_getdecspacevalues)
   call registeroptuegetnparam(trim(LIS_DEMCId)//char(0), DEMC_getNparam)
   call registeroptuereadrestart(trim(LIS_DEMCId)//char(0), DEMC_readrestart)
#endif

#if ( defined OPTUE_ALG_DEMCZ )
   call registeroptuealginit(trim(LIS_DEMCzId)//char(0), DEMCz_init)
   call registeroptuealgsetup(trim(LIS_DEMCzId)//char(0),DEMCz_setup)
   call registeroptuealgrun(trim(LIS_DEMCzId)//char(0), DEMCz_run)
   call registeroptueconvergencecheck(trim(LIS_DEMCzId)//char(0), &
                                      DEMCz_checkconvergence)
   call registeroptuegetdecisionspace(trim(LIS_DEMCzId)//char(0), &
                                      DEMCz_getdecspacevalues)
   call registeroptuegetnparam(trim(LIS_DEMCzId)//char(0), DEMCz_getNparam)
   call registeroptuereadrestart(trim(LIS_DEMCzId)//char(0), DEMCz_readrestart)
#endif
#endif
end subroutine LIS_optUEAlgorithm_plugin
end module LIS_optUEAlgorithm_pluginMod
