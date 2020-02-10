!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
! each of the LSM that is used in a data assimilation setup.
! Currently two algorithms are supported, Direct Insertion (DI),
! and Ensemble Kalman Filter (EnKF).
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
!     \item[QC the LSM state]
!      Routine that QCs the given LSM state for physical consistency.
!      (to be registered using {\tt registerroutingdaqcstate} and later
!      called through {\tt routingdaqcstate} method)
!     \item[QC the OBS state]
!      Routine that QCs the given OBS state based on land model states
!      (e.g. filter out observations when dense vegetation is present)
!      (to be registered using {\tt registerroutingdaqcobsstate} and later
!      called through {\tt routingdaqcobsstate} method)
!     \item[Scale LSM state]
!      Scale the LSM state variables in order to change the variables
!      to similar scales so that the matrices are well-conditioned
!      (to be registered using {\tt registerscalelsmstate} and later
!      called through {\tt scalelsmstate} method)
!     \item[Descale LSM state]
!      Descale the LSM state variables from the scaled state.
!      (to be registered using \newline
!      {\tt registerdescalelsmstate} and later
!      called through {\tt descalelsmstate} method)
!     \item[Update LSM state]
!      updates the state variables using the values of the Increment
!      object (note that the set method actually sets the variables to the
!      lsm state, whereas the update method simply changes the LSM\_State
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
!      Routine that transforms observations to the LSM state space
!      (could be as simple as a unit conversion).
!      (to be registered using {\tt registerroutingdaobstransform} and later
!      called through {\tt routingdaobstransform} method)
!     \item[Map observations to LSM space]
!      Routine that maps observations to the LSM state variables.
!      (to be registered using {\tt registerroutingdamapobstolsm} and later
!      called through {\tt routingdamapobstolsm} method)
!     \end{description}
!
!  The user-defined functions are included in the registry using a
!  user-selected indices. For example, consider the Noah LSM is
!  incorporated in the registry by the following calls. In case of
!  registry functions defined with two indices, the first index
!  refers to Noah LSM and the second index refers to the ``assimilation
!  set'' (in the following example, soil moisture assimilation using
!  synthetic soil moisture). For the definition of plugin indices,
!  please see LIS\_pluginIndices.F90
!
!  \begin{verbatim}
!    call registerroutingdagetstatevar(1,1,noah271_getsoilm)
!    call registerroutingdasetstatevar(1,1,noah271_setsoilm)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call routingdagetstatevar(1,1) - calls noah271_getsoilm
!    call routingdasetstatevar(1,1) - calls noah271_setsoilm
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call routingdagetstatevar(LIS_rc%lsm, LIS_rc%daset)
!    call routingdasetstatevar(LIS_rc%lsm, LIS_rc%daset)
!    call routingdaobstransform(LIS_rc%lsm, LIS_rc%daset)
!   \end{verbatim}
!   where $LIS\_rc\%lsm$ and $LIS\_rc\%daset$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
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

#endif

#endif
 end subroutine LIS_routingda_plugin
end module LIS_routingda_pluginMod
