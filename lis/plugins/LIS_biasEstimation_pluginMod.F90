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
module LIS_biasEstimation_pluginMod
!BOP
!
! !MODULE: LIS_biasEstimation_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining routines that performs bias correction in data assimilation
!   The user defined functions are incorporated into
!   the appropriate registry to be later invoked through generic calls.
!
! !REVISION HISTORY:
!  30 Nov 2007;   Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_biasEstimation_plugin
contains
!BOP
! !ROUTINE: LIS_biasEstimation_plugin
! \label{LIS_biasEstimation_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new bias estimation
!  scheme. This is an optional extention to the overall data assimilation
!  implementation.
!
!  The following interfaces should be implemented in this routine.
!  \begin{description}
!  \item[init]
!      Initialization of global datastructures
!      (to be registered using {\tt registerbiasestimationinit}
!       and later called through {\tt biasestimationinit})
!  \item[setup]
!      setup and initialization of data structures specific to each
!      instance of bias estimation algorithm
!      (to be registered using {\tt registerbiasestimationsetup}
!       and later called through {\tt biasestimationsetup})
!  \item[bias estimation computation]
!      Routines that compute bias correction estimates prior to
!      model propagation
!      (to be registered using {\tt registerbiasestimationcompute} and later called
!       through {\tt computebiascorrection})
!  \item[bias estimation update]
!      Routines that updates and applies the bias correction estimates.
!      (to be registered using {\tt registerbiasestimationupdate} and later called
!       through {\tt applybiascorrection})
!  \item[finalize]
!      Routines that cleans up the allocated memory structures
!   \end{description}
!
!  A single index is used to register a routine,
!  corresponding to the type of bias estimation algorithm used. For example,
!  consider an instance where method 'FOO' is employed, and is denoted
!  by a numerical index y. The methods
!  should be defined in the registry as follows.
!
!  \begin{verbatim}
!    call registerbiasestimationinit(y,FOO_init)
!    call registerbiasestimationsetup(y,FOO_setup)
!    call registerbiasestimationcompute(y,FOO_comp)
!    call registerbiasestimationupdate(y,FOO_update)
!    call registerbiasestimationfinalize(y,FOO_final)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call biasestimationinit(y) -  calls FOO_init
!    call biasestimationsetup(y) - calls FOO_run
!    call computebiascorrection(y) - calls FOO_comp
!    call applybiascorrection(y) - calls FOO_update
!    call biasestimationfinalize(y) - calls FOO_final
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call biasestimationinit(LIS_rc%biasalg(j))
!    call biasestimationsetup(LIS_rc%biasalg(j))
!   \end{verbatim}
!   where $LIS\_rc\%biasalg$ is set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_biasEstimation_plugin
!EOP
! !USES:
#if ( defined BE_bias_estimation )
   use LIS_pluginIndices
   use gmaobias_estimationMod, only : gmaobiasestimation_init,          &
                                      gmaobiasestimation_setup,         &
                                      gmaobiasestimation_calc,          &
                                      gmaobiasestimation_update,        &
                                      gmaobiasestimation_write_restart, &
                                      gmaobiasestimation_finalize

!----------------------------------------------------------------------
! The following convention is used to define the index of registered
! functions. The user may change these options, and the lis.config
! should be changed appropriately to ensure that the correct function
! is called at run time
!----------------------------------------------------------------------

   call registerbiasestimationinit(trim(LIS_gmaobiasId)//char(0), &
                                   gmaobiasestimation_init)
   call registerbiasestimationsetup(trim(LIS_gmaobiasId)//char(0), &
                                    gmaobiasestimation_setup)
!    call registerbiasestimationrun(trim(LIS_gmaobiasId)//char(0), &
!                                   gmaobiasestimation)
   call registerbiasestimationcompute(trim(LIS_gmaobiasId)//char(0), &
                                      gmaobiasestimation_calc)
   call registerbiasestimationupdate(trim(LIS_gmaobiasId)//char(0), &
                                     gmaobiasestimation_update)
   call registerbiasestimationrestart(trim(LIS_gmaobiasId)//char(0), &
                                      gmaobiasestimation_write_restart)
   call registerbiasestimationfinalize(trim(LIS_gmaobiasId)//char(0), &
                                       gmaobiasestimation_finalize)

#endif
end subroutine LIS_biasEstimation_plugin
end module LIS_biasEstimation_pluginMod
