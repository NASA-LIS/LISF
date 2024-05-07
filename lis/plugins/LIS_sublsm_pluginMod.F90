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
module LIS_sublsm_pluginMod
!BOP
!
! !MODULE: LIS_sublsm_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   land surface model initialization, execution, reading and
!   writing of restart files and other relevant land surface
!   model computations, corresponding to each of the LSMs used in LIS.
!
! !REVISION HISTORY:
!  09 Oct 2003    Sujay Kumar  Initial Specification
!  04 Jun 2021    Mahdi Navari Modified for Naoh.3.9
!  12 Aug 2021    Kristi Arsenault  Added SnowModel 
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_sublsm_plugin

contains
!BOP
! !ROUTINE: LIS_sublsm_plugin
!  \label{LIS_sublsm_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new LSM.
! The interface mandates that the following routines be implemented
! and registered for each of the LSM that is included in LIS.
!
!  \begin{description}
!  \item[Initialization]
!      Definition of LSM variables
!      (to be registered using {\tt registersublsminit} and later called
!       using {\tt sublsminit})
!  \item[Setup]
!      Initialization of parameters
!      (to be registered using {\tt registersublsmsetup} and later called
!       using {\tt sublsmsetup})
!  \item[Run]
!      Routines to execute LSM on a single gridcell for single timestep
!      (to be registered using {\tt registersublsmrun} and later called
!       using {\tt sublsmrun})
!  \item[Read restart]
!      Routines to read a restart file for an LSM run
!      (to be registered using {\tt registersublsmrestart} and later called
!       using {\tt sublsmrestart})
!  \item[Forcing transfer to model tiles]
!      Routines to transfer an array of given forcing to model tiles
!      (to be registered using {\tt registersublsmf2t} and later called
!       using {\tt sublsmf2t})
!  \item[Write restart]
!      Routines to write a restart file
!      (to be registered using {\tt registersublsmwrst} and later called
!       using {\tt sublsmwrst})
!  \item[Finalize]
!      Routines to cleanup LSM data structures
!      (to be registered using {\tt registersublsmfinalize} and later called
!       using {\tt sublsmfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the Crocus LSM is
!  incorporated in the registry with an index of 1 and
!  is invoked later by the following calls
!
!  \begin{verbatim}
!    call registersublsminit(1,noah_lsm_ini)
!    call registersublsmsetup(1,noah_setup)
!    call registersublsmf2t(1,noah_f2t)
!    call registersublsmrun(1,noah_main)
!    call registersublsmrestart(1,noahrst)
!    call registersublsmdynsetup(1,noah_dynsetup)
!    call registersublsmwrst(1,noah_writerst)
!    call registersublsmfinalize(1,noah_finalize)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call sublsminit(1)     - calls noah_lsm_ini
!    call sublsmsetup(1)    - calls noah_setup
!    call sublsmf2t(1)      - calls noah_f2t
!    call sublsmrun(1)      - calls noah_main
!    call sublsmdynsetup(1) - calls noah_dynsetup
!    call sublsmwrst(1)     - calls noah_writerst
!    call sublsmfinalize(1) - calls noah_finalize
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call sublsminit(lis%lsm)
!    call sublsmsetup(lis%lsm)
!    call sublsmf2t(lis%lsm)
!    call sublsmrun(lis%lsm)
!    call sublsmdynsetup(lis%lsm)
!    call sublsmwrst(lis%lsm)
!    call sublsmfinalize(lis%lsm)
!   \end{verbatim}
!   where $lis\%lsm$ is  set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LIS_sublsm_plugin
!EOP
    use LIS_pluginIndices

#if ( defined SM_Crocus_8_1 )
   use Crocus81_lsmMod, only : Crocus81_ini
#endif

#if ( defined SM_SNOWMODEL )
   use snowmodel_lsmMod, only : snowmodel_init
#endif

   implicit none

#if ( defined SM_Crocus_8_1 )
   external Crocus81_main
   external Crocus81_setup
   external Crocus81_readrst
   external Crocus81_f2t
   external Crocus81_dynsetup
   external Crocus81_writerst
   external Crocus81_finalize

   external Crocus81_setLSMimport
   external Crocus81_getLSMexport

#if ( defined SM_LSM_TEMPLATE )
   external template_getCROCUSexport
   external template_setCROCUSimport
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   external NoahMP401_getCROCUSexport
   external NoahMP401_setCROCUSimport
#endif

#if ( defined SM_NOAH_3_9 )
   external Noah39_getCROCUSexport
   external Noah39_setCROCUSimport
#endif

#endif

#if ( defined SM_SNOWMODEL )
   external snowmodel_main
   external snowmodel_setup
   external snowmodel_readrst
   external snowmodel_dynsetup
   external snowmodel_f2t
   external snowmodel_writerst
   external snowmodel_finalize
   external snowmodel_reset

   external snowmodel_setLSMimport
   external snowmodel_getLSMexport

#if ( defined SM_NOAHMP_4_0_1 )
   external NoahMP401_getSnowModelexport
   external NoahMP401_setSnowModelimport
#endif

#endif

#if ( defined SM_Crocus_8_1 )
   call registersublsminit(trim(LIS_Crocus81Id)//char(0),Crocus81_ini)
   call registersublsmsetup(trim(LIS_Crocus81Id)//char(0),Crocus81_setup)
   call registersublsmf2t(trim(LIS_Crocus81Id)//"+"//&
        trim(LIS_retroId)//char(0),Crocus81_f2t)
   call registersublsmf2t(trim(LIS_Crocus81Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),Crocus81_f2t)
   call registersublsmrun(trim(LIS_Crocus81Id)//char(0),Crocus81_main)
   call registersublsmrestart(trim(LIS_Crocus81Id)//char(0),Crocus81_readrst)
   call registersublsmdynsetup(trim(LIS_Crocus81Id)//char(0),Crocus81_dynsetup)
   call registersublsmwrst(trim(LIS_Crocus81Id)//char(0),Crocus81_writerst)
   call registersublsmfinalize(trim(LIS_Crocus81Id)//char(0),Crocus81_finalize)

   !wirings between NoahMP and CROCUS
   call registersublsmsetLSMimport(trim(LIS_Crocus81Id)//char(0),&
        Crocus81_setLSMimport)

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsm2sublsmgetexport(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_Crocus81Id)//char(0),NoahMP401_getCROCUSexport)
   call registerlsmsetsublsmimport(trim(LIS_noahmp401Id)//char(0),&
        NoahMP401_setCROCUSimport)
   call registersublsm2lsmgetexport(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_Crocus81Id)//char(0),Crocus81_getLSMexport)
#endif

#if ( defined SM_NOAH_3_9 )
   call registerlsm2sublsmgetexport(trim(LIS_noah39Id)//"+"//&
        trim(LIS_Crocus81Id)//char(0),Noah39_getCROCUSexport)
   call registerlsmsetsublsmimport(trim(LIS_noah39Id)//char(0),&
        Noah39_setCROCUSimport)
   call registersublsm2lsmgetexport(trim(LIS_noah39Id)//"+"//&
        trim(LIS_Crocus81Id)//char(0),Crocus81_getLSMexport)
#endif


#endif


#if ( defined SM_SNOWMODEL )
   call registersublsminit(trim(LIS_snowmodelId)//char(0),snowmodel_init)
   call registersublsmsetup(trim(LIS_snowmodelId)//char(0),snowmodel_setup)
   call registersublsmf2t(trim(LIS_snowmodelId)//"+"&
        //trim(LIS_retroId)//char(0),snowmodel_f2t)
   call registersublsmf2t(trim(LIS_snowmodelId)//"+"//&
        trim(LIS_agrmetrunId)//char(0),snowmodel_f2t)
   call registersublsmrun(trim(LIS_snowmodelId)//char(0),snowmodel_main)
   call registersublsmrestart(trim(LIS_snowmodelId)//char(0),snowmodel_readrst)
   call registersublsmdynsetup(trim(LIS_snowmodelId)//char(0),snowmodel_dynsetup)
   call registersublsmwrst(trim(LIS_snowmodelId)//char(0),snowmodel_writerst)
   call registersublsmfinalize(trim(LIS_snowmodelId)//char(0),snowmodel_finalize)
   call registersublsmreset(trim(LIS_snowmodelId)//char(0),snowmodel_reset)

   ! Wirings between NoahMP and SnowModel
   call registersublsmsetLSMimport(trim(LIS_snowmodelId)//char(0),&
        Snowmodel_setLSMimport)

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsm2sublsmgetexport(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snowmodelId)//char(0),NoahMP401_getSnowModelexport)
   call registerlsmsetsublsmimport(trim(LIS_noahmp401Id)//char(0),&
        NoahMP401_setSnowModelimport)
   call registersublsm2lsmgetexport(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snowmodelId)//char(0),SnowModel_getLSMexport)
#endif
#endif

  end subroutine LIS_sublsm_plugin

end module LIS_sublsm_pluginMod
