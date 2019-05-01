!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_lsm_pluginMod
!BOP
!
! !MODULE: LIS_lsm_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   land surface model initialization, execution, reading and
!   writing of restart files and other relevant land surface
!   model computations, corresponding to each of the LSMs used in LIS.
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
  PUBLIC :: LIS_lsm_plugin

contains
!BOP
! !ROUTINE: LIS_lsm_plugin
!  \label{LIS_lsm_plugin}
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
subroutine LIS_lsm_plugin
!EOP
   use LIS_pluginIndices
#if ( defined SM_LSM_TEMPLATE )
   use template_lsmMod, only : template_lsm_ini
#endif

#if ( defined SM_NOAH_2_7_1 )
   use noah271_lsmMod, only : noah271_lsm_ini
#endif

#if ( defined SM_NOAH_3_2 )
   use noah32_lsmMod, only : noah32_lsm_ini
#endif

#if ( defined SM_NOAH_3_3 )
   use noah33_lsmMod, only : noah33_lsm_ini
#endif

#if ( defined SM_NOAH_3_6 )
   use noah36_lsmMod, only : noah36_lsm_ini
#endif

#if ( defined SM_NOAH_3_9 )
   use noah39_lsmMod, only : noah39_lsm_ini
#endif

#if ( defined SM_NOAHMP_3_6 )
   use NoahMP36_lsmMod, only : noahmp36_ini
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   use NoahMP401_lsmMod, only : noahmp401_ini
#endif

#if ( defined SM_RUC_3_7 )
   use RUC37_lsmMod, only : ruc37_ini
#endif

#if ( defined SM_JULES_4_3 )
   use jules43_lsmMod, only : jules43_ini
#endif

#if ( defined SM_JULES_5_0 )
   use jules50_lsmMod, only : jules50_ini
#endif

#if ( defined SM_JULES_5_2 )
   use jules52_lsmMod, only : jules52_ini
#endif

#if ( defined SM_JULES_5_3 )
   use jules53_lsmMod, only : jules53_ini
#endif

#if ( defined SM_CLM_2 )
   use clm2_lsmMod, only : clm2_lsm_ini
   use clm2_atmdrvMod, only : clm2_atmdrv
#endif

#if ( defined SM_VIC_4_1_1 )
   use vic411_lsmMod, only : vic411_lsm_ini
#endif

#if ( defined SM_VIC_4_1_2 )
   use vic412_lsmMod, only : vic412_lsm_ini
#endif

#if ( defined SM_MOSAIC )
   use mos_lsmMod, only : mos_lsm_ini
#endif

#if ( defined SM_HYSSIB )
   use hyssib_lsmMod, only : hyssib_lsm_ini
#endif

#if ( defined SM_CLSM_F2_5 )
   use clsmf25_lsmMod, only : clsmf25_varalloc
#endif

#if ( defined SM_CABLE )
   use cable_lsmMod, only : cable_lsm_ini
#endif

#if ( defined SM_GEOWRSI_2 )
   use geowrsi2_lsmMod, only : geowrsi2_lsm_ini
#endif

#if ( defined SM_HTESSEL )
!    use tess_lsmMod,  only : tess_lsm_ini
#endif

#if ( defined SM_FASST )
!    use fasst_lsmMod, only : fasst_lsm_ini
#endif

#if ( defined SM_RDHM_3_5_6 )
   use RDHM356_lsmMod, only : RDHM356_ini
#endif

#if ( defined SM_SUMMA_1_0 )
   use summa1_lsmMod,  only : summa1_lsm_ini
#endif

#if ( defined SM_LSM_TEMPLATE )
   external template_main
   external template_setup
   external template_readrst
   external template_f2t
   external template_dynsetup
   external template_writerst
   external template_finalize
   external template_reset
#endif

#if ( defined SM_JULES_4_3 )
   external jules43_main
   external jules43_setup
   external jules43_readrst
   external jules43_f2t
   external jules43_dynsetup
   external jules43_writerst
   external jules43_finalize
#endif

#if ( defined SM_JULES_5_0 )
   external jules50_main
   external jules50_setup
   external jules50_readrst
   external jules50_f2t
   external jules50_dynsetup
   external jules50_writerst
   external jules50_finalize
#endif

#if ( defined SM_JULES_5_2 )
   external jules52_main
   external jules52_setup
   external jules52_readrst
   external jules52_f2t
   external jules52_dynsetup
   external jules52_writerst
   external jules52_finalize
#endif

#if ( defined SM_JULES_5_3 )
   external jules53_main
   external jules53_setup
   external jules53_readrst
   external jules53_f2t
   external jules53_dynsetup
   external jules53_writerst
   external jules53_finalize
#endif

#if ( defined SM_MOSAIC )
   external mos_main
   external mos_setup
   external mos_readrestart
   external mosdynp
   external mos_f2t
   external mos_writerst
   external mos_finalize
#endif

#if ( defined SM_HYSSIB )
   external hyssib_main
   external hyssib_setup
   external hyssib_readrst
   external hyssib_f2t
   external hyssib_writerst
   external hyssib_dynsetup
   external hyssib_finalize
#endif

#if ( defined SM_NOAH_3_6 )
   external noah36_main
   external noah36_setup
   external noah36_readrst
   external noah36_dynsetup
   external noah36_f2t
   external noah36_writerst
   external noah36_finalize
   external noah36_reset
#endif

#if ( defined SM_NOAH_3_9 )
   external noah39_main
   external noah39_setup
   external noah39_readrst
   external noah39_dynsetup
   external noah39_f2t
   external noah39_writerst
   external noah39_finalize
   external noah39_reset
#endif

#if ( defined SM_NOAHMP_3_6 )
   external noahmp36_main
   external noahmp36_setup
   external noahmp36_readrst
   external noahmp36_dynsetup
   external noahmp36_f2t
   external noahmp36_writerst
   external noahmp36_finalize
   external noahmp36_reset
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   external noahmp401_main
   external noahmp401_setup
   external noahmp401_readrst
   external noahmp401_dynsetup
   external noahmp401_f2t
   external noahmp401_writerst
   external noahmp401_finalize
   external noahmp401_reset
#endif

#if ( defined SM_RUC_3_7 )
   external ruc37_main
   external ruc37_setup
   external ruc37_readrst
   external ruc37_dynsetup
   external ruc37_f2t
   external ruc37_writerst
   external ruc37_finalize
   external ruc37_reset
#endif

#if ( defined SM_NOAH_3_3 )
   external noah33_main
   external noah33_setup
   external noah33_readrst
   external noah33_dynsetup
   external noah33_f2t
   external noah33_writerst
   external noah33_finalize
   external noah33_reset
#endif

#if ( defined SM_NOAH_3_2 )
   external noah32_main
   external noah32_setup
   external noah32_readrst
   external noah32_dynsetup
   external noah32_f2t
   external noah32_writerst
   external noah32_finalize
#endif

#if ( defined SM_NOAH_2_7_1 )
   external noah271_main
   external noah271_setup
   external noah271_readrst
   external noah271_dynsetup
   external noah271_f2t
   external noah271_writerst
   external noah271_finalize
#endif

#if ( defined SM_CLM_2 )
   external clm2_main
   external clm2_setup
   external clm2_readrestart
   external clm2_dynsetup
   external clm2_writerestart
   external clm2_finalize
#endif

#if ( defined SM_VIC_4_1_1 )
   external vic411_main
   external vic411_setup
   external vic411_readrst
   external vic411_dynsetup
   external vic411_f2t
   external vic411_writerst
   external vic411_finalize
#endif

#if ( defined SM_VIC_4_1_2 )
   external vic412_main
   external vic412_setup
   external vic412_readrst
   external vic412_dynsetup
   external vic412_f2t
   external vic412_writerst
   external vic412_finalize
#endif

#if ( defined SM_CLSM_F2_5 )
   external clsmf25_main
   external clsmf25_setup
   external clsmf25_dynsetup
   external clsmf25_readrst
   external clsmf25_f2t
   external clsmf25_writerst
   external clsmf25_finalize
   external clsmf25_reset
#endif

#if ( defined SM_HTESSEL )
!    external tess_main
!    external tess_setup
!    external tess_dynsetup
!    external tess_readrestart
!    external tess_f2t
!    external tess_writerestart
!    external tess_finalize
#endif

#if ( defined SM_CABLE )
   external cable_driver
   external cable_setup
   external cable_dynsetup
   external cable_readrst
   external cable_f2t
   external cable_writerst
   external cable_finalize
#endif

#if ( defined SM_GEOWRSI_2 )
   external geowrsi2_main
   external geowrsi2_setup
   external geowrsi2_dynsetup
   external geowrsi2_readrst
   external geowrsi2_f2t
   external geowrsi2_writerst
   external geowrsi2_finalize
#endif

#if ( defined SM_RDHM_3_5_6 )
   external RDHM356_main
   external RDHM356_setup
   external RDHM356_dynsetup
   external RDHM356_readrst
   external RDHM356_f2t
   external RDHM356_writerst
   external RDHM356_finalize
#endif

#if ( defined SM_SUMMA_1_0 )
   external summa1_setup
   external summa1_main
   external summa1_readrst
   external summa1_dynsetup
   external summa1_f2t
   external summa1_writerst
   external summa1_finalize
#endif

#if ( defined SM_LSM_TEMPLATE )
   call registerlsminit(trim(LIS_templateLSMId)//char(0),template_lsm_ini)
   call registerlsmsetup(trim(LIS_templateLSMId)//char(0),template_setup)
   call registerlsmf2t(trim(LIS_templateLSMId)//"+"//&
        trim(LIS_retroId)//char(0),template_f2t)
   call registerlsmf2t(trim(LIS_templateLSMId)//"+"//&
        trim(LIS_forecastrunId)//char(0),template_f2t)
   call registerlsmrun(trim(LIS_templateLSMId)//char(0),template_main)
   call registerlsmrestart(trim(LIS_templateLSMId)//char(0),template_readrst)
   call registerlsmdynsetup(trim(LIS_templateLSMId)//char(0),template_dynsetup)
   call registerlsmwrst(trim(LIS_templateLSMId)//char(0),template_writerst)
   call registerlsmfinalize(trim(LIS_templateLSMId)//char(0),template_finalize)
   call registerlsmreset(trim(LIS_templateLSMId)//char(0),template_reset)
#endif

#if ( defined SM_JULES_4_3 )
   call registerlsminit(trim(LIS_jules43Id)//char(0),jules43_ini)
   call registerlsmsetup(trim(LIS_jules43Id)//char(0),jules43_setup)
   call registerlsmf2t(trim(LIS_jules43Id)//"+"//&
        trim(LIS_retroId)//char(0),jules43_f2t)
   call registerlsmf2t(trim(LIS_jules43Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),jules43_f2t)
   call registerlsmrun(trim(LIS_jules43Id)//char(0),jules43_main)
   call registerlsmrestart(trim(LIS_jules43Id)//char(0),jules43_readrst)
   call registerlsmdynsetup(trim(LIS_jules43Id)//char(0),jules43_dynsetup)
   call registerlsmwrst(trim(LIS_jules43Id)//char(0),jules43_writerst)
   call registerlsmfinalize(trim(LIS_jules43Id)//char(0),jules43_finalize)
#endif

#if ( defined SM_JULES_5_0 )
   call registerlsminit(trim(LIS_jules50Id)//char(0),jules50_ini)
   call registerlsmsetup(trim(LIS_jules50Id)//char(0),jules50_setup)
   call registerlsmf2t(trim(LIS_jules50Id)//"+"//&
        trim(LIS_retroId)//char(0),jules50_f2t)
   call registerlsmf2t(trim(LIS_jules50Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),jules50_f2t)
   call registerlsmrun(trim(LIS_jules50Id)//char(0),jules50_main)
   call registerlsmrestart(trim(LIS_jules50Id)//char(0),jules50_readrst)
   call registerlsmdynsetup(trim(LIS_jules50Id)//char(0),jules50_dynsetup)
   call registerlsmwrst(trim(LIS_jules50Id)//char(0),jules50_writerst)
   call registerlsmfinalize(trim(LIS_jules50Id)//char(0),jules50_finalize)
#endif

#if ( defined SM_JULES_5_2 )
   call registerlsminit(trim(LIS_jules52Id)//char(0),jules52_ini)
   call registerlsmsetup(trim(LIS_jules52Id)//char(0),jules52_setup)
   call registerlsmf2t(trim(LIS_jules52Id)//"+"//&
        trim(LIS_retroId)//char(0),jules52_f2t)
   call registerlsmf2t(trim(LIS_jules52Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),jules52_f2t)
   call registerlsmrun(trim(LIS_jules52Id)//char(0),jules52_main)
   call registerlsmrestart(trim(LIS_jules52Id)//char(0),jules52_readrst)
   call registerlsmdynsetup(trim(LIS_jules52Id)//char(0),jules52_dynsetup)
   call registerlsmwrst(trim(LIS_jules52Id)//char(0),jules52_writerst)
   call registerlsmfinalize(trim(LIS_jules52Id)//char(0),jules52_finalize)
#endif

#if ( defined SM_JULES_5_3 )
   call registerlsminit(trim(LIS_jules53Id)//char(0),jules53_ini)
   call registerlsmsetup(trim(LIS_jules53Id)//char(0),jules53_setup)
   call registerlsmf2t(trim(LIS_jules53Id)//"+"//&
        trim(LIS_retroId)//char(0),jules53_f2t)
   call registerlsmf2t(trim(LIS_jules53Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),jules53_f2t)
   call registerlsmrun(trim(LIS_jules53Id)//char(0),jules53_main)
   call registerlsmrestart(trim(LIS_jules53Id)//char(0),jules53_readrst)
   call registerlsmdynsetup(trim(LIS_jules53Id)//char(0),jules53_dynsetup)
   call registerlsmwrst(trim(LIS_jules53Id)//char(0),jules53_writerst)
   call registerlsmfinalize(trim(LIS_jules53Id)//char(0),jules53_finalize)
#endif

#if ( defined SM_NOAH_3_6 )
   call registerlsminit(trim(LIS_noah36Id)//char(0),noah36_lsm_ini)
   call registerlsmsetup(trim(LIS_noah36Id)//char(0),noah36_setup)
   call registerlsmf2t(trim(LIS_noah36Id)//"+"&
        //trim(LIS_retroId)//char(0),noah36_f2t)
   call registerlsmf2t(trim(LIS_noah36Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noah36_f2t)
   call registerlsmrun(trim(LIS_noah36Id)//char(0),noah36_main)
   call registerlsmrestart(trim(LIS_noah36Id)//char(0),noah36_readrst)
   call registerlsmdynsetup(trim(LIS_noah36Id)//char(0),noah36_dynsetup)
   call registerlsmwrst(trim(LIS_noah36Id)//char(0),noah36_writerst)
   call registerlsmfinalize(trim(LIS_noah36Id)//char(0),noah36_finalize)
   call registerlsmreset(trim(LIS_noah36Id)//char(0),noah36_reset)
#endif

#if ( defined SM_NOAH_3_9 )
   call registerlsminit(trim(LIS_noah39Id)//char(0),noah39_lsm_ini)
   call registerlsmsetup(trim(LIS_noah39Id)//char(0),noah39_setup)
   call registerlsmf2t(trim(LIS_noah39Id)//"+"&
        //trim(LIS_retroId)//char(0),noah39_f2t)
   call registerlsmf2t(trim(LIS_noah39Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noah39_f2t)
   call registerlsmrun(trim(LIS_noah39Id)//char(0),noah39_main)
   call registerlsmrestart(trim(LIS_noah39Id)//char(0),noah39_readrst)
   call registerlsmdynsetup(trim(LIS_noah39Id)//char(0),noah39_dynsetup)
   call registerlsmwrst(trim(LIS_noah39Id)//char(0),noah39_writerst)
   call registerlsmfinalize(trim(LIS_noah39Id)//char(0),noah39_finalize)
   call registerlsmreset(trim(LIS_noah39Id)//char(0),noah39_reset)
#endif

#if ( defined SM_NOAHMP_3_6 )
   call registerlsminit(trim(LIS_noahmp36Id)//char(0),noahmp36_ini)
   call registerlsmsetup(trim(LIS_noahmp36Id)//char(0),noahmp36_setup)
   call registerlsmf2t(trim(LIS_noahmp36Id)//"+"//trim(LIS_retroId)//char(0),&
        noahmp36_f2t)
   ! ------------wanshu----add registry for smootherDA for NoahMP-----------------
   call registerlsmf2t(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_smootherDAId)//char(0), noahmp36_f2t)
   ! -----------------------------------------------------------------------------
   call registerlsmf2t(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noahmp36_f2t)
   call registerlsmf2t(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_forecastrunId)//char(0),noahmp36_f2t)
   call registerlsmrun(trim(LIS_noahmp36Id)//char(0),noahmp36_main)
   call registerlsmrestart(trim(LIS_noahmp36Id)//char(0),noahmp36_readrst)
   call registerlsmdynsetup(trim(LIS_noahmp36Id)//char(0),noahmp36_dynsetup)
   call registerlsmwrst(trim(LIS_noahmp36Id)//char(0),noahmp36_writerst)
   call registerlsmfinalize(trim(LIS_noahmp36Id)//char(0),noahmp36_finalize)
   call registerlsmreset(trim(LIS_noahmp36Id)//char(0),noahmp36_reset)
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   call registerlsminit(trim(LIS_noahmp401Id)//char(0),noahmp401_ini)
   call registerlsmsetup(trim(LIS_noahmp401Id)//char(0),noahmp401_setup)
   call registerlsmf2t(trim(LIS_noahmp401Id)//"+"&
        //trim(LIS_retroId)//char(0),noahmp401_f2t)
   call registerlsmf2t(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noahmp401_f2t)
   call registerlsmrun(trim(LIS_noahmp401Id)//char(0),noahmp401_main)
   call registerlsmrestart(trim(LIS_noahmp401Id)//char(0),noahmp401_readrst)
   call registerlsmdynsetup(trim(LIS_noahmp401Id)//char(0),noahmp401_dynsetup)
   call registerlsmwrst(trim(LIS_noahmp401Id)//char(0),noahmp401_writerst)
   call registerlsmfinalize(trim(LIS_noahmp401Id)//char(0),noahmp401_finalize)
   call registerlsmreset(trim(LIS_noahmp401Id)//char(0),noahmp401_reset)
#endif

#if ( defined SM_RUC_3_7 )
   call registerlsminit(trim(LIS_ruc37Id)//char(0),ruc37_ini)
   call registerlsmsetup(trim(LIS_ruc37Id)//char(0),ruc37_setup)
   call registerlsmf2t(trim(LIS_ruc37Id)//"+"//trim(LIS_retroId)//char(0),&
        ruc37_f2t)
   call registerlsmf2t(trim(LIS_ruc37Id)//"+"//trim(LIS_agrmetrunId)//char(0),&
        ruc37_f2t)
   call registerlsmrun(trim(LIS_ruc37Id)//char(0),ruc37_main)
   call registerlsmrestart(trim(LIS_ruc37Id)//char(0),ruc37_readrst)
   call registerlsmdynsetup(trim(LIS_ruc37Id)//char(0),ruc37_dynsetup)
   call registerlsmwrst(trim(LIS_ruc37Id)//char(0),ruc37_writerst)
   call registerlsmfinalize(trim(LIS_ruc37Id)//char(0),ruc37_finalize)
#endif

#if ( defined SM_NOAH_3_3 )
   call registerlsminit(trim(LIS_noah33Id)//char(0),noah33_lsm_ini)
   call registerlsmsetup(trim(LIS_noah33Id)//char(0),noah33_setup)
   call registerlsmf2t(trim(LIS_noah33Id)//"+"&
        //trim(LIS_retroId)//char(0),noah33_f2t)
   call registerlsmf2t(trim(LIS_noah33Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noah33_f2t)
   call registerlsmf2t(trim(LIS_noah33Id)//"+"//&
        trim(LIS_forecastrunId)//char(0),noah33_f2t)
   call registerlsmrun(trim(LIS_noah33Id)//char(0),noah33_main)
   call registerlsmrestart(trim(LIS_noah33Id)//char(0),noah33_readrst)
   call registerlsmdynsetup(trim(LIS_noah33Id)//char(0),noah33_dynsetup)
   call registerlsmwrst(trim(LIS_noah33Id)//char(0),noah33_writerst)
   call registerlsmfinalize(trim(LIS_noah33Id)//char(0),noah33_finalize)
   call registerlsmreset(trim(LIS_noah33Id)//char(0),noah33_reset)
#endif

#if ( defined SM_NOAH_3_2 )
   call registerlsminit(trim(LIS_noah32Id)//char(0),noah32_lsm_ini)
   call registerlsmsetup(trim(LIS_noah32Id)//char(0),noah32_setup)
   call registerlsmf2t(trim(LIS_noah32Id)//"+"&
        //trim(LIS_retroId)//char(0),noah32_f2t)
   call registerlsmf2t(trim(LIS_noah32Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noah32_f2t)
   call registerlsmrun(trim(LIS_noah32Id)//char(0),noah32_main)
   call registerlsmrestart(trim(LIS_noah32Id)//char(0),noah32_readrst)
   call registerlsmdynsetup(trim(LIS_noah32Id)//char(0),noah32_dynsetup)
   call registerlsmwrst(trim(LIS_noah32Id)//char(0),noah32_writerst)
   call registerlsmfinalize(trim(LIS_noah32Id)//char(0),noah32_finalize)
#endif

#if ( defined SM_NOAH_2_7_1 )
   call registerlsminit(trim(LIS_noah271Id)//char(0),noah271_lsm_ini)
   call registerlsmsetup(trim(LIS_noah271Id)//char(0),noah271_setup)
   call registerlsmf2t(trim(LIS_noah271Id)//"+"&
        //trim(LIS_retroId)//char(0),noah271_f2t)
   call registerlsmf2t(trim(LIS_noah271Id)//"+"//&
        trim(LIS_agrmetrunId)//char(0),noah271_f2t)
   call registerlsmrun(trim(LIS_noah271Id)//char(0),noah271_main)
   call registerlsmrestart(trim(LIS_noah271Id)//char(0),noah271_readrst)
   call registerlsmdynsetup(trim(LIS_noah271Id)//char(0),noah271_dynsetup)
   call registerlsmwrst(trim(LIS_noah271Id)//char(0),noah271_writerst)
   call registerlsmfinalize(trim(LIS_noah271Id)//char(0),noah271_finalize)
#endif


#if ( defined SM_CLM_2 )
   call registerlsminit(trim(LIS_clm2Id)//char(0), clm2_lsm_ini)
   call registerlsmsetup(trim(LIS_clm2Id)//char(0), clm2_setup)
   call registerlsmf2t(trim(LIS_clm2Id)//"+"&
        //trim(LIS_retroId)//char(0),clm2_atmdrv)
   call registerlsmrun(trim(LIS_clm2Id)//char(0),clm2_main)
   call registerlsmrestart(trim(LIS_clm2Id)//char(0),clm2_readrestart)
   call registerlsmdynsetup(trim(LIS_clm2Id)//char(0),clm2_dynsetup)
   call registerlsmwrst(trim(LIS_clm2Id)//char(0),clm2_writerestart)
   call registerlsmfinalize(trim(LIS_clm2Id)//char(0),clm2_finalize)
#endif

#if ( defined SM_VIC_4_1_1 )
   call registerlsminit(trim(LIS_vic411Id)//char(0),vic411_lsm_ini)
   call registerlsmsetup(trim(LIS_vic411Id)//char(0),vic411_setup)
   call registerlsmf2t(trim(LIS_vic411Id)//"+"//trim(LIS_retroId)//char(0),&
        vic411_f2t)
   call registerlsmrun(trim(LIS_vic411Id)//char(0),vic411_main)
   call registerlsmdynsetup(trim(LIS_vic411Id)//char(0),vic411_dynsetup)
   call registerlsmrestart(trim(LIS_vic411Id)//char(0),vic411_readrst)
   call registerlsmwrst(trim(LIS_vic411Id)//char(0),vic411_writerst)
   call registerlsmfinalize(trim(LIS_vic411Id)//char(0),vic411_finalize)
#endif

#if ( defined SM_VIC_4_1_2 )
   call registerlsminit(trim(LIS_vic412Id)//char(0),vic412_lsm_ini)
   call registerlsmsetup(trim(LIS_vic412Id)//char(0),vic412_setup)
   call registerlsmf2t(trim(LIS_vic412Id)//"+"//trim(LIS_retroId)//char(0),&
        vic412_f2t)
   call registerlsmrun(trim(LIS_vic412Id)//char(0),vic412_main)
   call registerlsmdynsetup(trim(LIS_vic412Id)//char(0),vic412_dynsetup)
   call registerlsmrestart(trim(LIS_vic412Id)//char(0),vic412_readrst)
   call registerlsmwrst(trim(LIS_vic412Id)//char(0),vic412_writerst)
   call registerlsmfinalize(trim(LIS_vic412Id)//char(0),vic412_finalize)
#endif

#if ( defined SM_MOSAIC )
   call registerlsminit(trim(LIS_mosaicId)//char(0),mos_lsm_ini)
   call registerlsmsetup(trim(LIS_mosaicId)//char(0), mos_setup)
   call registerlsmf2t(trim(LIS_mosaicId)//"+"//&
        trim(LIS_retroId)//char(0),mos_f2t)
   call registerlsmrun(trim(LIS_mosaicId)//char(0), mos_main)
   call registerlsmrestart(trim(LIS_mosaicId)//char(0), mos_readrestart)
   call registerlsmdynsetup(trim(LIS_mosaicId)//char(0),mosdynp)
   call registerlsmwrst(trim(LIS_mosaicId)//char(0), mos_writerst)
   call registerlsmfinalize(trim(LIS_mosaicId)//char(0),mos_finalize)
#endif

#if ( defined SM_HYSSIB )
   call registerlsminit(trim(LIS_hyssibId)//char(0),hyssib_lsm_ini)
   call registerlsmsetup(trim(LIS_hyssibId)//char(0), hyssib_setup)
   call registerlsmf2t(trim(LIS_hyssibId)//"+"//&
        trim(LIS_retroId)//char(0), hyssib_f2t)
   call registerlsmrun(trim(LIS_hyssibId)//char(0), hyssib_main)
   call registerlsmrestart(trim(LIS_hyssibId)//char(0),hyssib_readrst)
   call registerlsmdynsetup(trim(LIS_hyssibId)//char(0), hyssib_dynsetup)
   call registerlsmwrst(trim(LIS_hyssibId)//char(0),hyssib_writerst)
   call registerlsmfinalize(trim(LIS_hyssibId)//char(0), hyssib_finalize)
#endif

#if ( defined SM_HTESSEL )
!    call registerlsminit(trim(LIS_tessId)//char(0),tess_lsm_ini)
!    call registerlsmsetup(trim(LIS_tessId)//char(0),tess_setup)
!    call registerlsmrun(trim(LIS_tessId)//char(0),tess_main)
!    call registerlsmrestart(trim(LIS_tessId)//char(0),tess_readrestart)
!    call registerlsmdynsetup(trim(LIS_tessId)//char(0),tess_dynsetup)
!    call registerlsmf2t(trim(LIS_tessId)//"+"//&
!         trim(LIS_retroId)//char(0),tess_f2t)
!    call registerlsmf2t(trim(LIS_tessId)//"+"//&
!         trim(LIS_agrmetrunId)//char(0),tess_f2t)
!    call registerlsmwrst(trim(LIS_tessId)//char(0),tess_writerestart)
!    call registerlsmfinalize(trim(LIS_tessId)//char(0),tess_finalize)
#endif

#if ( defined SM_CABLE )
   call registerlsminit(trim(LIS_cableId)//char(0),cable_lsm_ini)
   call registerlsmsetup(trim(LIS_cableId)//char(0),cable_setup)
   call registerlsmrun(trim(LIS_cableId)//char(0),cable_driver)
   call registerlsmrestart(trim(LIS_cableId)//char(0),cable_readrst)
   call registerlsmdynsetup(trim(LIS_cableId)//char(0),cable_dynsetup)
   call registerlsmf2t(trim(LIS_cableId)//"+"//&
        trim(LIS_retroId)//char(0),cable_f2t)
   call registerlsmwrst(trim(LIS_cableId)//char(0),cable_writerst)
   call registerlsmfinalize(trim(LIS_cableId)//char(0),cable_finalize)
#endif

#if ( defined SM_GEOWRSI_2 )
   call registerlsminit(trim(LIS_geowrsi2Id)//char(0),geowrsi2_lsm_ini)
   call registerlsmsetup(trim(LIS_geowrsi2Id)//char(0),geowrsi2_setup)
   call registerlsmrun(trim(LIS_geowrsi2Id)//char(0),geowrsi2_main)
   call registerlsmrestart(trim(LIS_geowrsi2Id)//char(0),geowrsi2_readrst)
   call registerlsmdynsetup(trim(LIS_geowrsi2Id)//char(0),geowrsi2_dynsetup)
   call registerlsmf2t(trim(LIS_geowrsi2Id)//"+"//&
        trim(LIS_retroId)//char(0),geowrsi2_f2t)
   call registerlsmwrst(trim(LIS_geowrsi2Id)//char(0),geowrsi2_writerst)
   call registerlsmfinalize(trim(LIS_geowrsi2Id)//char(0),geowrsi2_finalize)
#endif

#if ( defined SM_CLSM_F2_5 )
   call registerlsminit(trim(LIS_clsmf25Id)//char(0),clsmf25_varalloc)
   call registerlsmsetup(trim(LIS_clsmf25Id)//char(0),clsmf25_setup)
   call registerlsmrun(trim(LIS_clsmf25Id)//char(0),clsmf25_main)
   call registerlsmrestart(trim(LIS_clsmf25Id)//char(0),clsmf25_readrst)
   call registerlsmdynsetup(trim(LIS_clsmf25Id)//char(0),clsmf25_dynsetup)
   call registerlsmf2t(trim(LIS_clsmf25Id)//"+"//trim(LIS_retroId)//char(0),&
        clsmf25_f2t)
   call registerlsmf2t(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_smootherDAId)//char(0), clsmf25_f2t)
   call registerlsmf2t(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_forecastrunId)//char(0),clsmf25_f2t)
   call registerlsmwrst(trim(LIS_clsmf25Id)//char(0),clsmf25_writerst)
   call registerlsmfinalize(trim(LIS_clsmf25Id)//char(0),clsmf25_finalize)
   call registerlsmreset(trim(LIS_clsmf25Id)//char(0),clsmf25_reset)
#endif

#if ( defined SM_RDHM_3_5_6 )
   call registerlsminit(trim(LIS_rdhm356lsmId)//char(0),RDHM356_ini)
   call registerlsmsetup(trim(LIS_rdhm356lsmId)//char(0),RDHM356_setup)
   call registerlsmrun(trim(LIS_rdhm356lsmId)//char(0),RDHM356_main)
   call registerlsmrestart(trim(LIS_rdhm356lsmId)//char(0),RDHM356_readrst)
   call registerlsmdynsetup(trim(LIS_rdhm356lsmId)//char(0),RDHM356_dynsetup)
   call registerlsmf2t(trim(LIS_rdhm356lsmId)//"+"//&
        trim(LIS_retroId)//char(0),RDHM356_f2t)
   call registerlsmwrst(trim(LIS_rdhm356lsmId)//char(0),RDHM356_writerst)
   call registerlsmfinalize(trim(LIS_rdhm356lsmId)//char(0),RDHM356_finalize)
#endif

#if ( defined SM_SUMMA_1_0 )
   call registerlsminit(trim(LIS_summa1Id)//char(0),summa1_lsm_ini)
   call registerlsmsetup(trim(LIS_summa1Id)//char(0),summa1_setup)
   call registerlsmrun(trim(LIS_summa1Id)//char(0),summa1_main)
   call registerlsmrestart(trim(LIS_summa1Id)//char(0),summa1_readrst)
   call registerlsmdynsetup(trim(LIS_summa1Id)//char(0),summa1_dynsetup)
   call registerlsmwrst(trim(LIS_summa1Id)//char(0),summa1_writerst)
   call registerlsmfinalize(trim(LIS_summa1Id)//char(0),summa1_finalize)
   call registerlsmf2t(trim(LIS_summa1Id)//"+"//&
        trim(LIS_retroId)//char(0),summa1_f2t)
#endif

end subroutine LIS_lsm_plugin
end module LIS_lsm_pluginMod
