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
module LIS_lakemodel_pluginMod
!BOP
!
! !MODULE: LIS_lakemodel_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   land surface model initialization, execution, reading and
!   writing of restart files and other relevant land surface
!   model computations, corresponding to each of the LAKEMODELs used in LIS.
!
! !REVISION HISTORY:
!  09 Oct 03    Sujay Kumar  Initial Specification
!  22 May 13    Shugong Wang Add FLake 1.0
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_lakemodel_plugin

contains
!BOP
! !ROUTINE: LIS_lakemodel_plugin
!  \label{LIS_lakemodel_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new LAKEMODEL.
! The interface mandates that the following routines be implemented
! and registered for each of the LAKEMODEL that is included in LIS.
!
!  \begin{description}
!  \item[Initialization]
!      Definition of LAKEMODEL variables
!      (to be registered using {\tt registerlakemodelinit} and later called
!       using {\tt lakemodelinit})
!  \item[Setup]
!      Initialization of parameters
!      (to be registered using {\tt registerlakemodelsetup} and later called
!       using {\tt lakemodelsetup})
!  \item[Run]
!      Routines to execute LAKEMODEL on a single gridcell for single timestep
!      (to be registered using {\tt registerlakemodelrun} and later called
!       using {\tt lakemodelrun})
!  \item[Read restart]
!      Routines to read a restart file for an LAKEMODEL run
!      (to be registered using {\tt registerlakemodelrestart} and later called
!       using {\tt lakemodelrestart})
!  \item[Output]
!      Routines to write output
!      (to be registered using {\tt registerlakemodeloutput} and later called
!       using {\tt lakemodeloutput})
!  \item[Forcing transfer to model tiles]
!      Routines to transfer an array of given forcing to model tiles
!      (to be registered using {\tt registerlakemodelf2t} and later called
!       using {\tt lakemodelf2t})
!  \item[Write restart]
!      Routines to write a restart file
!      (to be registered using {\tt registerlakemodelwrst} and later called
!       using {\tt lakemodelwrst})
!  \item[Finalize]
!      Routines to cleanup LAKEMODEL data structures
!      (to be registered using {\tt registerlakemodelfinalize} and later called
!       using {\tt lakemodelfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the Noah LAKEMODEL is
!  incorporated in the registry with an index of 1 and
!  is invoked later by the following calls
!
!  \begin{verbatim}
!    call registerlakemodelinit(1,noah_lakemodel_ini)
!    call registerlakemodelsetup(1,noah_setup)
!    call registerlakemodelf2t(1,noah_f2t)
!    call registerlakemodelrun(1,noah_main)
!    call registerlakemodelrestart(1,noahrst)
!    call registerlakemodeldynsetup(1,noah_dynsetup)
!    call registerlakemodeloutput(1,noah_output)
!    call registerlakemodelwrst(1,noah_writerst)
!    call registerlakemodelfinalize(1,noah_finalize)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call lakemodelinit(1)      - calls noah_lakemodel_ini
!    call lakemodelsetup(1)    - calls noah_setup
!    call lakemodelf2t(1)      - calls noah_f2t
!    call lakemodelrun(1)      - calls noah_main
!    call lakemodeldynsetup(1) - calls noah_dynsetup
!    call lakemodeloutput(1)   - calls noah_output
!    call lakemodelwrst(1)     - calls noah_writerst
!    call lakemodelfinalize(1) - calls noah_finalize
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call lakemodelinit(lis%lakemodel)
!    call lakemodelsetup(lis%lakemodel)
!    call lakemodelf2t(lis%lakemodel)
!    call lakemodelrun(lis%lakemodel)
!    call lakemodeldynsetup(lis%lakemodel)
!    call lakemodeloutput(lis%lakemodel)
!    call lakemodelwrst(lis%lakemodel)
!    call lakemodelfinalize(lis%lakemodel)
!   \end{verbatim}
!   where $lis\%lakemodel$ is  set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_lakemodel_plugin
!EOP
   use LIS_pluginIndices

#if ( defined SM_FLAKE_1_0 )
   use FLake1_Mod, only : FLake1_ini
#endif

#if ( defined SM_FLAKE_1_0 )
   external FLake1_main
   !external FLake1_ini
   external FLake1_setup
   external FLake1_dynsetup
   external FLake1_readrst
   external FLake1_f2t
   external FLake1_writerst
   external FLake1_finalize
#endif

#if ( defined SM_FLAKE_1_0 )
   call registerlakemodelinit(trim(LIS_FLake1Id)//char(0),FLake1_ini)
   call registerlakemodelsetup(trim(LIS_FLake1Id)//char(0),FLake1_setup)
   call registerlakemodelrun(trim(LIS_FLake1Id)//char(0),FLake1_main)
   call registerlakemodelrestart(trim(LIS_FLake1Id)//char(0),FLake1_readrst)
   call registerlakemodeldynsetup(trim(LIS_FLake1Id)//char(0),FLake1_dynsetup)
   call registerlakemodelf2t(trim(LIS_FLake1Id)//"+"//&
        trim(LIS_retroId)//char(0),FLake1_f2t)
   call registerlakemodelf2t(trim(LIS_FLake1Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),FLake1_f2t)
   call registerlakemodelwrst(trim(LIS_FLake1Id)//char(0),FLake1_writerst)
   call registerlakemodelfinalize(trim(LIS_FLake1Id)//char(0),FLake1_finalize)
#endif
end subroutine LIS_lakemodel_plugin
end module LIS_lakemodel_pluginMod
