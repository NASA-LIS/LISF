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
module LIS_routingda_pluginMod
!BOP
!
! !MODULE: LIS_routingda_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that
!   need to be defined to enable the use of a routing model
!   in a data assimilation setup.
!
! !REVISION HISTORY:
! 07 Nov 2019: Sujay KUmar, Initial specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_routingda_plugin
contains
!BOP
! !ROUTINE: LIS_routingda_plugin
!  \label{LIS_routingda_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new routing scheme
! in a data assimilation mode. The interface mandates that
! a number of routines be implemented and registered for
! each of the routing model that is used in a data assimilation setup.
!
!  For use with EnKF, the following routines need to be implemented.
!   \begin{description}
!    \item[Return the state prognostic variables]
!      Routine that retrieves the specified state prognostic variable(s)
!      (to be registered using {\tt registerroutingdagetstatevar} and later called
!       through {\tt routingdagetstatevar} method)
!    \item[Set the state prognostic variables]
!      Routine that sets the specified state prognostic variable(s)
!      (to be registered using {\tt registerroutingdasetstatevar} and later called
!       through {\tt routingdasetstatevar} method)
!     \item[Retrieve ``Obspred'']
!      Routine that retrieves the model's estimate of observations.
!      (to be registered using {\tt registerroutingdagetobspred} and later
!      called through {\tt routingdagetobspred} method)
!     \item[QC the Routing state]
!      Routine that QCs the given Routing state for physical consistency.
!      (to be registered using {\tt registerroutingdaqcstate} and later
!      called through {\tt routingdaqcstate} method)
!     \item[QC the OBS state]
!      Routine that QCs the given OBS state based on routing model states
!      (e.g. filter out observations when dense vegetation is present)
!      (to be registered using {\tt registerroutingdaqcobsstate} and later
!      called through {\tt routingdaqcobsstate} method)
!     \item[Scale Routing state]
!      Scale the Routing state variables in order to change the variables
!      to similar scales so that the matrices are well-conditioned
!      (to be registered using {\tt registerscaleroutingstate} and later
!      called through {\tt scaleroutingstate} method)
!     \item[Descale Routing state]
!      Descale the Routing state variables from the scaled state.
!      (to be registered using \newline
!      {\tt registerdescaleroutingstate} and later
!      called through {\tt descaleroutingstate} method)
!     \item[Update routing state]
!      updates the state variables using the values of the Increment
!      object (note that the set method actually sets the variables to the
!      routing state, whereas the update method simply changes the Routing\_State
!      object).
!      (to be registered using {\tt registerroutingdaupdatestate} and later
!      called through {\tt routingdaupdatestate} method)
!    \end{description}
!
!   For use with DI, the following routines are required.
!    \begin{description}
!    \item[Return the state prognostic variables]
!      Routine that retrieves the specified state prognostic variable(s)
!      (to be registered using {\tt registerroutingdagetstatevar} and later called
!       through {\tt routingdagetstatevar} method)
!    \item[Set the state prognostic variables]
!      Routine that sets the specified state prognostic variable(s)
!      (to be registered using {\tt registerroutingdasetstatevar} and later called
!       through {\tt routingdasetstatevar} method)
!     \item[Transform Observations]
!      Routine that transforms observations to the Routing state space
!      (could be as simple as a unit conversion).
!      (to be registered using {\tt registerroutingdaobstransform} and later
!      called through {\tt routingdaobstransform} method)
!     \item[Map observations to routing space]
!      Routine that maps observations to the routing state variables.
!      (to be registered using {\tt registerroutingdamapobstorouting} and later
!      called through {\tt routingdamapobstorouting} method)
!     \end{description}
!
!
! !INTERFACE:
subroutine LIS_routingda_plugin
!EOP

#if ( ( defined DA_DIRECT_INSERTION ) || \
      ( defined DA_ENKS )             || \
      ( defined DA_ENKF ) )

   use LIS_pluginIndices

#if ( defined ROUTE_HYMAP2_ROUTER )
   use HYMAP2_dawL_Mod
#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
   external HYMAP2_getStateSpaceSize
   external HYMAP2_getWL
   external HYMAP2_setWL
   external HYMAP2_getWLpred
   external HYMAP2_qcWL
   external HYMAP2_qc_WLobs
   external HYMAP2_scale_WL
   external HYMAP2_descale_WL
   external HYMAP2_updateWL
   external HYMAP2_setPertStates
#endif

#if ( defined ROUTE_HYMAP2_ROUTER )
   call registerroutingdainit(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_daWL_init)
   call registerroutingdagetstatespacesize(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_getStateSpaceSize)
   call registerroutingdagetstatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_getWL)
   call registerroutingdasetstatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_setWL)
   call registerroutingdagetobspred(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_getWLpred)
   call registerroutingdaqcstate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_qcWL)
   call registerroutingdaqcobsstate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_qc_WLobs)
   call registerroutingdascalestatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_scale_WL)
   call registerroutingdadescalestatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_descale_WL)
   call registerroutingdaupdatestate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_updateWL)
   call registerroutingdasetpertstates(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_synwlId)//char(0),HYMAP2_setPertStates)


   call registerroutingdainit(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_daWL_init)
   call registerroutingdagetstatespacesize(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_getStateSpaceSize)
   call registerroutingdagetstatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_getWL)
   call registerroutingdasetstatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_setWL)
   call registerroutingdagetobspred(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_getWLpred)
   call registerroutingdaqcstate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_qcWL)
   call registerroutingdaqcobsstate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_qc_WLobs)
   call registerroutingdascalestatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_scale_WL)
   call registerroutingdadescalestatevar(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_descale_WL)
   call registerroutingdaupdatestate(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_updateWL)
   call registerroutingdasetpertstates(trim(LIS_HYMAP2routerId)//"+"//&
        trim(LIS_hydrowebwlId)//char(0),HYMAP2_setPertStates)
#endif

#endif
 end subroutine LIS_routingda_plugin
end module LIS_routingda_pluginMod
