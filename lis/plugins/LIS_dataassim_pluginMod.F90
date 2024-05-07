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
module LIS_dataassim_pluginMod
!BOP
!
! !MODULE: LIS_dataassim_pluginMod.F90
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining routines that performs data assimilation.
!   The user defined functions are incorporated into
!   the appropriate registry to be later invoked through generic calls.
!
! !REVISION HISTORY:
!  27 Feb 2005;   Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_dataassim_plugin

contains
!BOP
! !ROUTINE: LIS_dataassim_plugin
! \label{LIS_dataassim_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new data assimilation
!  scheme. As explaind in the {\tt LIS\_dataassim\_module}, there are
!  three different abstractions associated with the
!  data assimilation implementation. The implementations defined in this
!  subroutine complete the ``wirings'' required to complete the
!  data assimilation algorithm-related abstractions. Other required
!  interfaces defined in {\tt lsmda\_pluginMod} and {\tt dataobs\_pluginMod}
!  should be completed to complete a successful implementation.
!
!  The following interfaces should be implemented in this routine.
!  \begin{description}
!  \item[init]
!      Initialization of data and memory structures
!      (to be registered using {\tt registerdainit} and later called
!       through {\tt dataassiminit})
!  \item[Assimilate]
!      Routines to perform data assimilation using observations.
!      (to be registered using {\tt registerassim} and later called
!       through {\tt assimilate})
!  \item[Output]
!       Routines to perform data assimilation related output
!       (to be registered using {\tt registerdaoutput} and later called
!       through {\tt daoutput})
!  \item[Finalize]
!       Routines to perform data assimilation related cleanups
!       (to be registered using {\tt registerdafinalize} and later called
!       throuh {\tt dafinalize})
!   \end{description}
!
!  A single index is used to register a routine,
!  corresponding to the type of data assimilation algorithm used. For example,
!  consider an instance where Ensemble
!  Kalman Filter (EnKF) algorithm is employed, and is denoted
!  by an index 3. The methods
!  should be defined in the registry as follows.
!
!  \begin{verbatim}
!    call registerdainit(3,enkf_init)
!    call registerassim(3,enkf_assim)
!    call registerdaoutput(3,enkf_output)
!    call registerdafinalize(3,enkf_final)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call dainit(3) -  calls enkf_init
!    call assimilate(3) - calls enkf_assim
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call dainit(lis%daalg)
!    call assimilate(lis%daalg)
!   \end{verbatim}
!   where $lis\%daalg$ is set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_dataassim_plugin
!EOP
! !USES:

   use LIS_pluginIndices

#if ( defined DA_DIRECT_INSERTION )
   use directInsertion_Mod, only : di_init, di_increments, di_final,       &
                                   di_setup, di_update, di_diagnostics
#endif

#if ( defined DA_EKF )
   use ekf_Mod,            only : ekf_init, ekf_setup,                  &
                                   ekf_increments, ekf_update,           &
                                   ekf_final,ekf_diagnostics
#endif

#if ( defined DA_ENKF )
   use enkf_Mod,            only : enkf_init, enkf_setup,                  &
                                   enkf_increments, enkf_update,           &
                                   enkf_final,enkf_diagnostics
#endif

#if ( defined DA_ENSRF )
   use ensrf_Mod,            only : ensrf_init, ensrf_setup,                  &
                                   ensrf_increments, ensrf_update,           &
                                   ensrf_final,ensrf_diagnostics
#endif

#if ( defined DA_ENKS )
   use enksgrace_Mod,       only : enksgrace_init, enksgrace_setup,        &
                                   enksgrace_increments, enksgrace_update, &
                                   enksgrace_final,enksgrace_diagnostics
#endif

#if ( defined DA_PF )
   use pf_Mod,            only : pf_init, pf_setup,                  &
                                   pf_increments, pf_update,           &
                                   pf_final,pf_diagnostics
#endif
!----------------------------------------------------------------------
! The following convention is used to define the index of registered
! functions. The user may change these options, and the lis.config
! should be changed appropriately to ensure that the correct function
! is called at run time
!----------------------------------------------------------------------

#if ( defined DA_DIRECT_INSERTION )
   call registerdainit(trim(LIS_diId)//char(0),di_init)
   call registerdasetup(trim(LIS_diId)//char(0),di_setup)
   call registercomputeincrements(trim(LIS_diId)//char(0),di_increments)
   call registerapplyincrements(trim(LIS_diId)//char(0),di_update)
   call registerdaoutput(trim(LIS_diId)//char(0),di_diagnostics)
   call registerdafinalize(trim(LIS_diId)//char(0),di_final)
#endif

#if ( defined DA_EKF )
   call registerdainit(trim(LIS_ekfId)//char(0),ekf_init)
   call registerdasetup(trim(LIS_ekfId)//char(0),ekf_setup)
   call registercomputeincrements(trim(LIS_ekfId)//char(0),&
        ekf_increments)
   call registerapplyincrements(trim(LIS_ekfId)//char(0),ekf_update)
   call registerdaoutput(trim(LIS_ekfId)//char(0),ekf_diagnostics)
   call registerdafinalize(trim(LIS_ekfId)//char(0),ekf_final)
#endif

#if ( defined DA_ENKF )
!EnKF-GMAO version
   call registerdainit(trim(LIS_enkfId)//char(0),enkf_init)
   call registerdasetup(trim(LIS_enkfId)//char(0),enkf_setup)
   call registercomputeincrements(trim(LIS_enkfId)//char(0),&
        enkf_increments)
   call registerapplyincrements(trim(LIS_enkfId)//char(0),enkf_update)
   call registerdaoutput(trim(LIS_enkfId)//char(0),enkf_diagnostics)
   call registerdafinalize(trim(LIS_enkfId)//char(0),enkf_final)
#endif

#if ( defined DA_ENSRF )
!Ensrf
   call registerdainit(trim(LIS_ensrfId)//char(0),ensrf_init)
   call registerdasetup(trim(LIS_ensrfId)//char(0),ensrf_setup)
   call registercomputeincrements(trim(LIS_ensrfId)//char(0),&
        ensrf_increments)
   call registerapplyincrements(trim(LIS_ensrfId)//char(0),ensrf_update)
   call registerdaoutput(trim(LIS_ensrfId)//char(0),ensrf_diagnostics)
   call registerdafinalize(trim(LIS_ensrfId)//char(0),ensrf_final)
#endif

#if ( defined DA_ENKS )
!GRACE ENKS version
   call registerdainit(trim(LIS_enksId)//char(0),enksgrace_init)
   call registerdasetup(trim(LIS_enksId)//char(0),enksgrace_setup)
   call registercomputeincrements(trim(LIS_enksId)//char(0), &
                                  enksgrace_increments)
   call registerapplyincrements(trim(LIS_enksId)//char(0),enksgrace_update)
   call registerdaoutput(trim(LIS_enksId)//char(0),enksgrace_diagnostics)
   call registerdafinalize(trim(LIS_enksId)//char(0),enksgrace_final)
#endif

#if ( defined DA_PF )
   call registerdainit(trim(LIS_pfId)//char(0),pf_init)
   call registerdasetup(trim(LIS_pfId)//char(0),pf_setup)
   call registercomputeincrements(trim(LIS_pfId)//char(0),&
        pf_increments)
   call registerapplyincrements(trim(LIS_pfId)//char(0),pf_update)
   call registerdaoutput(trim(LIS_pfId)//char(0),pf_diagnostics)
   call registerdafinalize(trim(LIS_pfId)//char(0),pf_final)
#endif

end subroutine LIS_dataassim_plugin
end module LIS_dataassim_pluginMod
