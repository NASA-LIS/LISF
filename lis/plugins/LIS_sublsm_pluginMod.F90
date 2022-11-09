!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
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
!  09 Oct 03    Sujay Kumar  Initial Specification
!  04 Jun 21    Mahdi Navari Modified for Naoh.3.9
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
!      (to be registered using {\tt registerlsminit} and later called
!       using {\tt lsminit})
!  \item[Setup]
!      Initialization of parameters
!      (to be registered using {\tt registerlsmsetup} and later called
!       using {\tt lsmsetup})
!  \item[Run]
!      Routines to execute LSM on a single gridcell for single timestep
!      (to be registered using {\tt registerlsmrun} and later called
!       using {\tt lsmrun})
!  \item[Read restart]
!      Routines to read a restart file for an LSM run
!      (to be registered using {\tt registerlsmrestart} and later called
!       using {\tt lsmrestart})
!  \item[Forcing transfer to model tiles]
!      Routines to transfer an array of given forcing to model tiles
!      (to be registered using {\tt registerlsmf2t} and later called
!       using {\tt lsmf2t})
!  \item[Write restart]
!      Routines to write a restart file
!      (to be registered using {\tt registerlsmwrst} and later called
!       using {\tt lsmwrst})
!  \item[Finalize]
!      Routines to cleanup LSM data structures
!      (to be registered using {\tt registerlsmfinalize} and later called
!       using {\tt lsmfinalize})
!  \end{description}
!
!  The user-defined functions are included in the registry using a
!  single index. For example, consider the Noah LSM is
!  incorporated in the registry with an index of 1 and
!  is invoked later by the following calls
!
!  \begin{verbatim}
!    call registerlsminit(1,noah_lsm_ini)
!    call registerlsmsetup(1,noah_setup)
!    call registerlsmf2t(1,noah_f2t)
!    call registerlsmrun(1,noah_main)
!    call registerlsmrestart(1,noahrst)
!    call registerlsmdynsetup(1,noah_dynsetup)
!    call registerlsmwrst(1,noah_writerst)
!    call registerlsmfinalize(1,noah_finalize)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call lsminit(1)      - calls noah_lsm_ini
!    call lsmsetup(1)    - calls noah_setup
!    call lsmf2t(1)      - calls noah_f2t
!    call lsmrun(1)      - calls noah_main
!    call lsmdynsetup(1) - calls noah_dynsetup
!    call lsmwrst(1)     - calls noah_writerst
!    call lsmfinalize(1) - calls noah_finalize
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call lsminit(lis%lsm)
!    call lsmsetup(lis%lsm)
!    call lsmf2t(lis%lsm)
!    call lsmrun(lis%lsm)
!    call lsmdynsetup(lis%lsm)
!    call lsmwrst(lis%lsm)
!    call lsmfinalize(lis%lsm)
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
  end subroutine LIS_sublsm_plugin

end module LIS_sublsm_pluginMod
