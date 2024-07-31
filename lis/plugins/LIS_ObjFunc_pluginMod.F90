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
module LIS_ObjFunc_pluginMod
!BOP
!
! !MODULE: LIS_ObjFunc_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that specify
!   different methods (metrics) to compute objective function for use
!   in optimization.
!
!   The user defined functions are incorporated into
!   the appropriate registry to be later invoked through generic calls.
!
! !REVISION HISTORY:
!  16 Jul 2009;   Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_ObjFunc_plugin
contains
!BOP
! !ROUTINE: LIS_ObjFunc_plugin
!  \label{LIS_ObjFunc_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing different methods for objective
!  function computations (in parameter estimation)
!
!
! !INTERFACE:
subroutine LIS_ObjFunc_plugin
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

#if ( defined PE_OBJFUNC_LS )
   use LSobjFunc_Mod, only : initializeLSObjFunc
#endif
#if ( defined PE_OBJFUNC_LM )
   use LMobjFunc_Mod, only : initializeLMObjFunc
#endif
#if ( defined PE_OBJFUNC_LL )
   use LLobjFunc_Mod, only : initializeLLObjFunc
#endif
#if ( defined PE_OBJFUNC_P )
   use PobjFunc_Mod,  only : initializePObjFunc
#endif

#if ( defined PE_OBJFUNC_LS )
   external updateLSestimate, computeLSestimate, resetLSestimate
#endif
#if ( defined PE_OBJFUNC_LM )
   external updateLMestimate, computeLMestimate, resetLMestimate
#endif
#if ( defined PE_OBJFUNC_LL )
   external updateLLestimate, computeLLestimate, resetLLestimate
#endif
#if ( defined PE_OBJFUNC_P )
   external updatePestimate, computePestimate, resetPestimate
#endif


! These methods are specified for different type of objective
! function evaluations

#if ( defined PE_OBJFUNC_LM )
   call registerinitobjfunctype(trim(LIS_LMestimateId)//char(0),    &
                                initializeLMObjFunc)
   call registerupdateobjfunctype(trim(LIS_LMestimateId)//char(0),  &
                                  updateLMestimate)
   call registercomputeobjfunctype(trim(LIS_LMestimateId)//char(0), &
                                   computeLMestimate)
   call registerresetobjfunctype(trim(LIS_LMestimateId)//char(0),   &
                                 resetLMestimate)
#endif

#if ( defined PE_OBJFUNC_LS )
   call registerinitobjfunctype(trim(LIS_LSestimateId)//char(0),    &
                                initializeLSObjFunc)
   call registerupdateobjfunctype(trim(LIS_LSestimateId)//char(0),  &
                                  updateLSestimate)
   call registercomputeobjfunctype(trim(LIS_LSestimateId)//char(0), &
                                   computeLSestimate)
   call registerresetobjfunctype(trim(LIS_LSestimateId)//char(0),   &
                                 resetLSestimate)
#endif

#if ( defined PE_OBJFUNC_LL )
   call registerinitobjfunctype(trim(LIS_LLestimateId)//char(0),    &
                                initializeLLObjFunc)
   call registerupdateobjfunctype(trim(LIS_LLestimateId)//char(0),  &
                                  updateLLestimate)
   call registercomputeobjfunctype(trim(LIS_LLestimateId)//char(0), &
                                   computeLLestimate)
   call registerresetobjfunctype(trim(LIS_LLestimateId)//char(0),   &
                                 resetLLestimate)
#endif

#if ( defined PE_OBJFUNC_P )
   call registerinitobjfunctype(trim(LIS_PestimateId)//char(0),     &
                                initializePObjFunc)
   call registerupdateobjfunctype(trim(LIS_PestimateId)//char(0),   &
                                  updatePestimate)
   call registercomputeobjfunctype(trim(LIS_PestimateId)//char(0),  &
                                   computePestimate)
   call registerresetobjfunctype(trim(LIS_PestimateId)//char(0),    &
                                 resetPestimate)

#endif
#endif
end subroutine LIS_ObjFunc_plugin
end module LIS_ObjFunc_pluginMod
