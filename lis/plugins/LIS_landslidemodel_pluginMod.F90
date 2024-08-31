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
module LIS_landslidemodel_pluginMod
!BOP
!
! !MODULE: LIS_landslidemodel_pluginMod
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
  PUBLIC :: LIS_landslidemodel_plugin
contains
!BOP
! !ROUTINE: LIS_landslidemodel_plugin
! \label{LIS_landslidemodel_plugin}
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
!    call lisinit(lis%landslidemodel)
!    call lisrun(lis%landslidemodel)
!    call lisfinal(lis%landslidemodel)
!   \end{verbatim}
!   where $lis\%landslidemodel$ is set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_landslidemodel_plugin
!EOP

   use LIS_pluginIndices
#if ( defined APP_GLS )
   use GLSMod, only : GLS_init, GLS_run, GLS_output, GLS_final
#endif
#if ( defined APP_TRIGRS )
   use TRIGRSMod, only : TRIGRS_init, TRIGRS_run, TRIGRS_output, TRIGRS_final
#endif

#if ( defined APP_GLS )
   call registerlandslidemodelinit(trim(LIS_GLSId)//char(0),GLS_init)
   call registerlandslidemodelrun(trim(LIS_GLSId)//char(0),GLS_run)
   call registerlandslidemodeloutput(trim(LIS_GLSId)//char(0),GLS_output)
   call registerlandslidemodelfinalize(trim(LIS_GLSId)//char(0),GLS_final)
#endif

#if ( defined APP_TRIGRS )
   call registerlandslidemodelinit(trim(LIS_TRIGRSId)//char(0),TRIGRS_init)
   call registerlandslidemodelrun(trim(LIS_TRIGRSId)//char(0),TRIGRS_run)
   call registerlandslidemodeloutput(trim(LIS_TRIGRSId)//char(0),TRIGRS_output)
   call registerlandslidemodelfinalize(trim(LIS_TRIGRSId)//char(0),TRIGRS_final)
#endif

end subroutine LIS_landslidemodel_plugin
end module LIS_landslidemodel_pluginMod
