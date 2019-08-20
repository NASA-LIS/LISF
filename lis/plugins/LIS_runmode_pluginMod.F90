!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_runmode_pluginMod
!BOP
!
! !MODULE: LIS_runmode_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   LIS initialization, execution, and finalization
!   for different running modes in LIS
!
! !REVISION HISTORY:
!  21 Oct 05    Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_runmode_plugin
contains
!BOP
! !ROUTINE: LIS_runmode_plugin
! \label{LIS_runmode_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new running
! mode. The interface mandates that the following routines be implemented
! and registered for each of the running modes in LIS
!
!  \begin{description}
!  \item[Initialization]
!      Defining the init routines needed for each
!      running mode in LIS.
!      (to be registered using {\tt registerlisinit} and later called
!       using {\tt lisinit})
!  \item[Run]
!      Define the execution routines needed for each
!      running mode in LIS
!      (to be registered using {\tt registerlisrun} and later called
!       using {\tt lisrun})
!  \item[Finalize]
!      Define the finalize routines needed for each
!      LIS running mode.
!      (to be registered using {\tt registerlisfinalize} and later called
!       using {\tt lisfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the definition of a
!  retrospective running mode in LIS. The methods
!  should be defined in the registry as follows, if the index of the
!  mode is defined to be 1.
!
!  \begin{verbatim}
!    call registerlisinit(1,lisinit_retrospective)
!    call registerlisrun(1,lisrun_retrospective)
!    call registerlisfinalize(1,lisfinal_retrospective)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call lisinit(1)  - calls lisinit_retrospective
!    call lisrun(1)   - calls lisrun_retrospective
!    call lisfinal(1)   - calls lisfinal_retrospective
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call lisinit(lis%runmode)
!    call lisrun(lis%runmode)
!    call lisfinal(lis%runmode)
!   \end{verbatim}
!   where $lis\%runmode$ is set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_runmode_plugin
!EOP

   use LIS_pluginIndices

#if ( defined RM_RETROSPECTIVE )
   use retrospective_runMod,  only : lis_init_retrospective, &
                                     lis_run_retrospective,  &
                                     lis_final_retrospective
#endif

#if ( defined RM_AGRMETMODE )
   use agrmet_runMod,         only : lis_init_agrmet, &
                                     lis_run_agrmet,  &
                                     lis_final_agrmet
#endif

#if ( defined RM_PARAM_ESTIMATION )
   use paramEstim_runMod,     only : lis_init_paramEstim, &
                                     lis_run_paramEstim,  &
                                     lis_final_paramEstim
#endif

#if ( defined RM_RTM_FORWARD )
   use RTMforward_runMod,     only : lis_init_RTMforward, &
                                     lis_run_RTMforward,  &
                                     lis_final_RTMforward
#endif

#if ( defined RM_ENSEMBLE_SMOOTHER )
   use smootherDA_runMod,     only : lis_init_smootherDA, &
                                     lis_run_smootherDA,  &
                                     lis_final_smootherDA
#endif

#if ( defined RM_FORECAST )
   use forecast_runMod,       only : lis_init_forecast, &
                                     lis_run_forecast,  &
                                     lis_final_forecast
#endif


#if ( defined RM_RETROSPECTIVE )
   call registerlisinit(trim(LIS_retroId)//char(0),lis_init_retrospective)
   call registerlisrun(trim(LIS_retroId)//char(0),lis_run_retrospective)
   call registerlisfinalize(trim(LIS_retroId)//char(0),lis_final_retrospective)
#endif

#if ( defined RM_AGRMETMODE )
   call registerlisinit(trim(LIS_agrmetrunId)//char(0),lis_init_agrmet)
   call registerlisrun(trim(LIS_agrmetrunId)//char(0),lis_run_agrmet)
   call registerlisfinalize(trim(LIS_agrmetrunId)//char(0),lis_final_agrmet)
#endif

#if ( defined RM_PARAM_ESTIMATION )
   call registerlisinit(trim(LIS_paramEstimRunId)//char(0),lis_init_paramEstim)
   call registerlisrun(trim(LIS_paramEstimRunId)//char(0),lis_run_paramEstim)
   call registerlisfinalize(trim(LIS_paramEstimRunId)//char(0), &
                            lis_final_paramEstim)
#endif

#if ( defined RM_RTM_FORWARD )
   call registerlisinit(trim(LIS_RTMforwardId)//char(0),lis_init_RTMforward)
   call registerlisrun(trim(LIS_RTMforwardId)//char(0),lis_run_RTMforward)
   call registerlisfinalize(trim(LIS_RTMforwardId)//char(0), &
                            lis_final_RTMforward)
#endif

#if ( defined RM_ENSEMBLE_SMOOTHER )
   call registerlisinit(trim(LIS_smootherDAId)//char(0),lis_init_smootherDA)
   call registerlisrun(trim(LIS_smootherDAId)//char(0),lis_run_smootherDA)
   call registerlisfinalize(trim(LIS_smootherDAId)//char(0), &
                            lis_final_smootherDA)
#endif

#if ( defined RM_FORECAST )
   call registerlisinit(trim(LIS_forecastrunId)//char(0),lis_init_forecast)
   call registerlisrun(trim(LIS_forecastrunId)//char(0),lis_run_forecast)
   call registerlisfinalize(trim(LIS_forecastrunId)//char(0), &
                            lis_final_forecast)
#endif

end subroutine LIS_runmode_plugin
end module LIS_runmode_pluginMod
