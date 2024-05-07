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
module LIS_openwatermodel_pluginMod
!BOP
!
! !MODULE: LIS_openwatermodel_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   land surface model initialization, execution, reading and
!   writing of restart files and other relevant land surface
!   model computations, corresponding to each of the OPENWATERs used in LIS.
!
! !REVISION HISTORY:
!  09 Oct 03    Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_openwatermodel_plugin

contains
!BOP
! !ROUTINE: LIS_openwatermodel_plugin
!  \label{LIS_openwatermodel_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new OPENWATER.
! The interface mandates that the following routines be implemented
! and registered for each of the OPENWATER that is included in LIS.
!
!  \begin{description}
!  \item[Initialization]
!      Definition of OPENWATER variables
!      (to be registered using {\tt registeropenwaterinit} and later called
!       using {\tt openwaterinit})
!  \item[Setup]
!      Initialization of parameters
!      (to be registered using {\tt registeropenwatersetup} and later called
!       using {\tt openwatersetup})
!  \item[Run]
!      Routines to execute OPENWATER on a single gridcell for single timestep
!      (to be registered using {\tt registeropenwaterrun} and later called
!       using {\tt openwaterrun})
!  \item[Read restart]
!      Routines to read a restart file for an OPENWATER run
!      (to be registered using {\tt registeropenwaterrestart} and later called
!       using {\tt openwaterrestart})
!  \item[Output]
!      Routines to write output
!      (to be registered using {\tt registeropenwateroutput} and later called
!       using {\tt openwateroutput})
!  \item[Forcing transfer to model tiles]
!      Routines to transfer an array of given forcing to model tiles
!      (to be registered using {\tt registeropenwaterf2t} and later called
!       using {\tt openwaterf2t})
!  \item[Write restart]
!      Routines to write a restart file
!      (to be registered using {\tt registeropenwaterwrst} and later called
!       using {\tt openwaterwrst})
!  \item[Finalize]
!      Routines to cleanup OPENWATER data structures
!      (to be registered using {\tt registeropenwaterfinalize} and later called
!       using {\tt openwaterfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the Noah OPENWATER is
!  incorporated in the registry with an index of 1 and
!  is invoked later by the following calls
!
!  \begin{verbatim}
!    call registeropenwaterinit(1,noah_openwatermodel_ini)
!    call registeropenwatersetup(1,noah_setup)
!    call registeropenwaterf2t(1,noah_f2t)
!    call registeropenwaterrun(1,noah_main)
!    call registeropenwaterrestart(1,noahrst)
!    call registeropenwaterdynsetup(1,noah_dynsetup)
!    call registeropenwateroutput(1,noah_output)
!    call registeropenwaterwrst(1,noah_writerst)
!    call registeropenwaterfinalize(1,noah_finalize)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call openwaterinit(1)      - calls noah_openwatermodel_ini
!    call openwatersetup(1)    - calls noah_setup
!    call openwaterf2t(1)      - calls noah_f2t
!    call openwaterrun(1)      - calls noah_main
!    call openwaterdynsetup(1) - calls noah_dynsetup
!    call openwateroutput(1)   - calls noah_output
!    call openwaterwrst(1)     - calls noah_writerst
!    call openwaterfinalize(1) - calls noah_finalize
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call openwaterinit(lis%openwater)
!    call openwatersetup(lis%openwater)
!    call openwaterf2t(lis%openwater)
!    call openwaterrun(lis%openwater)
!    call openwaterdynsetup(lis%openwater)
!    call openwateroutput(lis%openwater)
!    call openwaterwrst(lis%openwater)
!    call openwaterfinalize(lis%openwater)
!   \end{verbatim}
!   where $lis\%openwater$ is  set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_openwatermodel_plugin
!EOP
   use LIS_pluginIndices

#if ( defined SM_TEMPLATE_OPEN_WATER )
   use templateOpenWaterMod, only : template_openwater_ini
#endif

#if ( defined SM_TEMPLATE_OPEN_WATER )
   external templateOpenWater_main
   external templateOpenWater_setup
   external templateOpenWater_dynsetup
   external templateOpenWater_readrst
   external templateOpenWater_f2t
   external templateOpenWater_writerst
   external templateOpenWater_finalize
#endif

#if ( defined SM_TEMPLATE_OPEN_WATER )
   call registeropenwaterinit(trim(LIS_templateOpenWaterId)//char(0),&
        template_openwater_ini)
   call registeropenwatersetup(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_setup)
   call registeropenwaterrun(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_main)
   call registeropenwaterrestart(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_readrst)
   call registeropenwaterdynsetup(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_dynsetup)
   call registeropenwaterf2t(trim(LIS_templateOpenWaterId)//"+"//&
        trim(LIS_retroId)//char(0),templateOpenWater_f2t)
   call registeropenwaterf2t(trim(LIS_templateOpenWaterId)//"+"//&
        trim(LIS_agrmetrunId)//char(0),templateOpenWater_f2t)
   call registeropenwaterwrst(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_writerst)
   call registeropenwaterfinalize(trim(LIS_templateOpenWaterId)//char(0),&
        templateOpenWater_finalize)
#endif

end subroutine LIS_openwatermodel_plugin
end module LIS_openwatermodel_pluginMod
