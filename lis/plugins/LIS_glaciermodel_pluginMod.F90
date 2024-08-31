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
module LIS_glaciermodel_pluginMod
!BOP
!
! !MODULE: LIS_glaciermodel_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   land surface model initialization, execution, reading and
!   writing of restart files and other relevant land surface
!   model computations, corresponding to each of the GLACIERMODELs used in LIS.
!
! !REVISION HISTORY:
!  06 Apr 2018    Sujay Kumar  Initial Specification
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_glaciermodel_plugin

contains
!BOP
! !ROUTINE: LIS_glaciermodel_plugin
!  \label{LIS_glaciermodel_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new GLACIERMODEL.
! The interface mandates that the following routines be implemented
! and registered for each of the GLACIERMODEL that is included in LIS.
!
!  \begin{description}
!  \item[Initialization]
!      Definition of GLACIERMODEL variables
!      (to be registered using {\tt registerglaciermodelinit} and later called
!       using {\tt glaciermodelinit})
!  \item[Setup]
!      Initialization of parameters
!      (to be registered using {\tt registerglaciermodelsetup} and later called
!       using {\tt glaciermodelsetup})
!  \item[Run]
!      Routines to execute GLACIERMODEL on a single gridcell for single timestep
!      (to be registered using {\tt registerglaciermodelrun} and later called
!       using {\tt glaciermodelrun})
!  \item[Read restart]
!      Routines to read a restart file for an GLACIERMODEL run
!      (to be registered using {\tt registerglaciermodelrestart} and later called
!       using {\tt glaciermodelrestart})
!  \item[Output]
!      Routines to write output
!      (to be registered using {\tt registerglaciermodeloutput} and later called
!       using {\tt glaciermodeloutput})
!  \item[Forcing transfer to model tiles]
!      Routines to transfer an array of given forcing to model tiles
!      (to be registered using {\tt registerglaciermodelf2t} and later called
!       using {\tt glaciermodelf2t})
!  \item[Write restart]
!      Routines to write a restart file
!      (to be registered using {\tt registerglaciermodelwrst} and later called
!       using {\tt glaciermodelwrst})
!  \item[Finalize]
!      Routines to cleanup GLACIERMODEL data structures
!      (to be registered using {\tt registerglaciermodelfinalize} and later called
!       using {\tt glaciermodelfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the Noah GLACIERMODEL is
!  incorporated in the registry with an index of 1 and
!  is invoked later by the following calls
!
!  \begin{verbatim}
!    call registerglaciermodelinit(1,noah_glaciermodel_ini)
!    call registerglaciermodelsetup(1,noah_setup)
!    call registerglaciermodelf2t(1,noah_f2t)
!    call registerglaciermodelrun(1,noah_main)
!    call registerglaciermodelrestart(1,noahrst)
!    call registerglaciermodeldynsetup(1,noah_dynsetup)
!    call registerglaciermodeloutput(1,noah_output)
!    call registerglaciermodelwrst(1,noah_writerst)
!    call registerglaciermodelfinalize(1,noah_finalize)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call glaciermodelinit(1)      - calls noah_glaciermodel_ini
!    call glaciermodelsetup(1)    - calls noah_setup
!    call glaciermodelf2t(1)      - calls noah_f2t
!    call glaciermodelrun(1)      - calls noah_main
!    call glaciermodeldynsetup(1) - calls noah_dynsetup
!    call glaciermodeloutput(1)   - calls noah_output
!    call glaciermodelwrst(1)     - calls noah_writerst
!    call glaciermodelfinalize(1) - calls noah_finalize
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call glaciermodelinit(lis%glaciermodel)
!    call glaciermodelsetup(lis%glaciermodel)
!    call glaciermodelf2t(lis%glaciermodel)
!    call glaciermodelrun(lis%glaciermodel)
!    call glaciermodeldynsetup(lis%glaciermodel)
!    call glaciermodeloutput(lis%glaciermodel)
!    call glaciermodelwrst(lis%glaciermodel)
!    call glaciermodelfinalize(lis%glaciermodel)
!   \end{verbatim}
!   where $lis\%glaciermodel$ is  set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_glaciermodel_plugin
!EOP
   use LIS_pluginIndices
#if ( defined SM_GLACIER_TEMPLATE )
   use templateGL_Mod, only :templateGL_ini
#endif

#if ( defined SM_NOAHMP_GLACIER_3_9_1_1)
   use noahmpglacier3911_Mod, only : noahmpglacier3911_ini
#endif

#if ( defined SM_GLACIER_TEMPLATE )
   external templateGL_main
   external templateGL_setup
   external templateGL_dynsetup
   external templateGL_readrst
   external templateGL_f2t
   external templateGL_writerst
   external templateGL_finalize
#endif

#if ( defined SM_NOAHMP_GLACIER_3_9_1_1 )
   external noahmpglacier3911_main
   external noahmpglacier3911_setup
   external noahmpglacier3911_dynsetup
   external noahmpglacier3911_readrst
   external noahmpglacier3911_f2t
   external noahmpglacier3911_writerst
   external noahmpglacier3911_finalize
#endif

#if ( defined SM_GLACIER_TEMPLATE )
   call registerglaciermodelinit(trim(LIS_templateGLId)//char(0),&
        templateGL_ini)
   call registerglaciermodelsetup(trim(LIS_templateGLId)//char(0),&
        templateGL_setup)
   call registerglaciermodelrun(trim(LIS_templateGLId)//char(0),&
        templateGL_main)
   call registerglaciermodelrestart(trim(LIS_templateGLId)//char(0),&
        templateGL_readrst)
   call registerglaciermodeldynsetup(trim(LIS_templateGLId)//char(0),&
        templateGL_dynsetup)
   call registerglaciermodelf2t(trim(LIS_templateGLId)//"+"//&
        trim(LIS_retroId)//char(0),templateGL_f2t)
   call registerglaciermodelf2t(trim(LIS_templateGLId)//"+"//&
        trim(LIS_agrmetrunId)//char(0),templateGL_f2t)
   call registerglaciermodelwrst(trim(LIS_templateGLId)//char(0),&
        templateGL_writerst)
   call registerglaciermodelfinalize(trim(LIS_templateGLId)//char(0),&
        templateGL_finalize)
#endif

#if ( defined SM_NOAHMP_GLACIER_3_9_1_1 )
   call registerglaciermodelinit(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_ini)
   call registerglaciermodelsetup(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_setup)
   call registerglaciermodelrun(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_main)
   call registerglaciermodelrestart(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_readrst)
   call registerglaciermodeldynsetup(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_dynsetup)
   call registerglaciermodelf2t(trim(LIS_noahmpglacier3911Id)//"+"//&
        trim(LIS_retroId)//char(0),noahmpglacier3911_f2t)
   call registerglaciermodelf2t(trim(LIS_noahmpglacier3911Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noahmpglacier3911_f2t)
   call registerglaciermodelwrst(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_writerst)
   call registerglaciermodelfinalize(trim(LIS_noahmpglacier3911Id)//char(0),&
        noahmpglacier3911_finalize)
#endif
end subroutine LIS_glaciermodel_plugin
end module LIS_glaciermodel_pluginMod
