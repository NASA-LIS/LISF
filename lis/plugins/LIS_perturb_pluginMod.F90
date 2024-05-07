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
module LIS_perturb_pluginMod
!BOP
!
! !MODULE: LIS_perturb_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining routines that perform pertubations of model states,
!   forcing, or parameters. The user defined functions are
!   incorporated into the appropriate registry to be later invoked
!   through generic calls.
!
! !REVISION HISTORY:
!  08Jul2005;   Sujay Kumar  Initial Specification
!
! !INTERFACE:
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_perturb_plugin
contains
!BOP
! !ROUTINE: LIS_perturb_plugin
! \label{LIS_perturb_plugin}
!
! !DESCRIPTION:
!
! This is a plugin point for introducing a new forcing perturbation
! scheme. The interface mandates that the following interfaces be
! implemented for each scheme.
!
!  \begin{description}
!  \item[Setup]
!      Initialization of data and memory structures
!      (to be registered using {\tt registerperturbsetup} and later
!       called using the generic call {\tt perturbinit})
!  \item[Forecast]
!      Routines to compute perturbations
!      (to be registered using {\tt registerperturb} and later
!       called using the generic call {\tt perturb})
!   \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider a new scheme called `FOO'
!  The methods should be defined in the registry as follows,
!  if the index of the scheme is defined to be 1.
!
!  \begin{verbatim}
!    call registerperturbsetup(1,foo_setup)
!    call registerperturb(1,foo_method)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call perturbinit(1)  - calls foo_setup
!    call perturb(1)      - calls foo_method
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call perturbinit(lis%perturb)
!    call perturb(lis%perturb)
!   \end{verbatim}
!   where $lis\%perturb$ is set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_perturb_plugin
!EOP
#if ( defined PERT_perturbations )
   use LIS_pluginIndices
   use uniformPert_Mod
   use gmaopert_Mod


   call registerperturbinit(trim(LIS_uniformpertId)//char(0),uniformPert_init)
   call registerperturbsetup(trim(LIS_uniformpertId)//char(0),uniformPert_setup)
   call registerperturbmethod(trim(LIS_uniformpertId)//char(0),uniformperturb)
   call registerperturbreadrst(trim(LIS_uniformpertId)//char(0),&
                               uniformpert_readrestart)
   call registerperturbwriterst(trim(LIS_uniformpertId)//char(0),&
                                uniformpert_writerestart)

   call registerperturbinit(trim(LIS_gmaopertId)//char(0),gmaoPert_init)
   call registerperturbsetup(trim(LIS_gmaopertId)//char(0),gmaoPert_setup)
   call registerperturbmethod(trim(LIS_gmaopertId)//char(0),gmaoperturb)
   call registerperturbreadrst(trim(LIS_gmaopertId)//char(0),&
                               gmaopert_readrestart)
   call registerperturbwriterst(trim(LIS_gmaopertId)//char(0),&
                                gmaopert_writerestart)
#endif
end subroutine LIS_perturb_plugin
end module LIS_perturb_pluginMod
