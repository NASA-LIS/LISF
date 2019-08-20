!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_lsmda_pluginMod
!BOP
!
! !MODULE: LIS_lsmda_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that
!   need to be defined to enable the use of a land surface
!   model in a data assimilation setup.
!
! !REVISION HISTORY:
!  29 Mar 07    Sujay Kumar  Separated from the lsm module implementation.
!  30 Oct 2014: David Mocko, re-organized and added Noah-3.6
!  1  Aug 2016: Mahdi Navari, added NoahMP.3.6
!  Sep 2017: Mahdi Navari added JULES 5.0
!  Oct  2018 : Mhdi Navari , added Noah.3.9 
!  21 Oct 2018: Mahdi Navari, added NoahMP.3.9
!  Dec 2018: Mahdi Navari: added Noah-MP.4.0.1
!  13 May 2019: Yeosang Yoon, added SNODEP & LDTSI Assimilation for NoahMP.4.0.1
!  06 Jun 2019: Yeosang Yoon, added SNODEP Assimilation for Jules 5.0
!
!EOP
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_lsmda_plugin
contains
!BOP
! !ROUTINE: LIS_lsmda_plugin
!  \label{LIS_lsmda_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new LSM
! in a data assimilation mode. The interface mandates that
! a number of routines be implemented and registered for
! each of the LSM that is used in a data assimilation setup.
! Currently two algorithms are supported, Direct Insertion (DI),
! and Ensemble Kalman Filter (EnKF).
!
!  For use with EnKF, the following routines need to be implemented.
!   \begin{description}
!    \item[Return the state prognostic variables]
!      Routine that retrieves the specified state prognostic variable(s)
!      (to be registered using {\tt registerlsmdagetstatevar} and later called
!       through {\tt lsmdagetstatevar} method)
!    \item[Set the state prognostic variables]
!      Routine that sets the specified state prognostic variable(s)
!      (to be registered using {\tt registerlsmdasetstatevar} and later called
!       through {\tt lsmdasetstatevar} method)
!     \item[Retrieve ``Obspred'']
!      Routine that retrieves the model's estimate of observations.
!      (to be registered using {\tt registerlsmdagetobspred} and later
!      called through {\tt lsmdagetobspred} method)
!     \item[QC the LSM state]
!      Routine that QCs the given LSM state for physical consistency.
!      (to be registered using {\tt registerlsmdaqcstate} and later
!      called through {\tt lsmdaqcstate} method)
!     \item[QC the OBS state]
!      Routine that QCs the given OBS state based on land model states
!      (e.g. filter out observations when dense vegetation is present)
!      (to be registered using {\tt registerlsmdaqcobsstate} and later
!      called through {\tt lsmdaqcobsstate} method)
!     \item[Scale LSM state]
!      Scale the LSM state variables in order to change the variables
!      to similar scales so that the matrices are well-conditioned
!      (to be registered using {\tt registerscalelsmstate} and later
!      called through {\tt scalelsmstate} method)
!     \item[Descale LSM state]
!      Descale the LSM state variables from the scaled state.
!      (to be registered using \newline
!      {\tt registerdescalelsmstate} and later
!      called through {\tt descalelsmstate} method)
!     \item[Update LSM state]
!      updates the state variables using the values of the Increment
!      object (note that the set method actually sets the variables to the
!      lsm state, whereas the update method simply changes the LSM\_State
!      object).
!      (to be registered using {\tt registerlsmdaupdatestate} and later
!      called through {\tt lsmdaupdatestate} method)
!    \end{description}
!
!   For use with DI, the following routines are required.
!    \begin{description}
!    \item[Return the state prognostic variables]
!      Routine that retrieves the specified state prognostic variable(s)
!      (to be registered using {\tt registerlsmdagetstatevar} and later called
!       through {\tt lsmdagetstatevar} method)
!    \item[Set the state prognostic variables]
!      Routine that sets the specified state prognostic variable(s)
!      (to be registered using {\tt registerlsmdasetstatevar} and later called
!       through {\tt lsmdasetstatevar} method)
!     \item[Transform Observations]
!      Routine that transforms observations to the LSM state space
!      (could be as simple as a unit conversion).
!      (to be registered using {\tt registerlsmdaobstransform} and later
!      called through {\tt lsmdaobstransform} method)
!     \item[Map observations to LSM space]
!      Routine that maps observations to the LSM state variables.
!      (to be registered using {\tt registerlsmdamapobstolsm} and later
!      called through {\tt lsmdamapobstolsm} method)
!     \end{description}
!
!  The user-defined functions are included in the registry using a
!  user-selected indices. For example, consider the Noah LSM is
!  incorporated in the registry by the following calls. In case of
!  registry functions defined with two indices, the first index
!  refers to Noah LSM and the second index refers to the ``assimilation
!  set'' (in the following example, soil moisture assimilation using
!  synthetic soil moisture). For the definition of plugin indices,
!  please see LIS\_pluginIndices.F90
!
!  \begin{verbatim}
!    call registerlsmdagetstatevar(1,1,noah271_getsoilm)
!    call registerlsmdasetstatevar(1,1,noah271_setsoilm)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call lsmdagetstatevar(1,1) - calls noah271_getsoilm
!    call lsmdasetstatevar(1,1) - calls noah271_setsoilm
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call lsmdagetstatevar(LIS_rc%lsm, LIS_rc%daset)
!    call lsmdasetstatevar(LIS_rc%lsm, LIS_rc%daset)
!    call lsmdaobstransform(LIS_rc%lsm, LIS_rc%daset)
!   \end{verbatim}
!   where $LIS\_rc\%lsm$ and $LIS\_rc\%daset$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
subroutine LIS_lsmda_plugin
!EOP

#if ( ( defined DA_DIRECT_INSERTION ) || \
      ( defined DA_ENKS )             || \
      ( defined DA_ENKF ) )

   use LIS_pluginIndices

#if ( defined SM_NOAH_2_7_1 )
#if ( defined DA_OBS_SNODEP )
   use noah271_dasnow_Mod
#endif
#endif

#if ( defined SM_NOAH_3_3 )
   use noah33_dasoilm_Mod
   use noah33_dasnow_Mod
#endif

#if ( defined SM_NOAH_3_6 )
   use noah36_dasoilm_Mod
   use noah36_dasnow_Mod
#endif

#if ( defined SM_NOAH_3_9 )
   use noah39_dasoilm_Mod
#if ( defined DA_OBS_LDTSI )
   use noah39_dasnow_Mod
#endif
#if ( defined DA_OBS_SNODEP )
   use noah39_dasnow_Mod
#endif
#endif

#if ( defined SM_NOAHMP_3_6 )
   use noahmp36_dasoilm_Mod 
   use noahmp36_dasnow_Mod
   use noahmp36_dasnodep_Mod
   use noahmp36_tws_DAlogMod, only : noahmp36_tws_DAlog
   use noahmp36_datws_Mod
   use noahmp36_daveg_Mod
   use noahmp36_daalbedo_Mod
#endif

#if ( defined SM_NOAHMP_4_0_1 )
   use NoahMP401_dasoilm_Mod
   use noahmp401_dasnodep_Mod
   use noahmp401_daldtsi_Mod
#endif


#if ( defined SM_CLSM_F2_5 )
   use clsmf25_tws_DAlogMod, only : clsmf25_tws_DAlog
   use clsmf25_dasoilm_Mod
   use clsmf25_dasnow_Mod
   use clsmf25_datws_Mod
#endif

#if ( defined SM_RUC_3_7 )
   use ruc37_dasoilm_Mod
   use ruc37_dasnow_Mod
#endif

#if ( defined SM_JULES_4_3 )
   use jules43_dasoilm_Mod
#endif

#if ( defined JULES_5_0_DEV )
! Disable; JULES.5.0 DA not ready yet.
#if ( defined SM_JULES_5_0 )
   use jules50_dasoilm_Mod
   use jules50_dasnodep_Mod
#endif

#if ( defined SM_JULES_5_2 )
   use jules52_dasoilm_Mod
#endif

#if 0
#if ( defined SM_JULES_5_3 )
   use jules53_dasoilm_Mod
#endif
#endif

#endif

#if ( defined SM_NOAH_2_7_1 )
#if ( defined DA_OBS_SNODEP )
! Noah-2.7.1 snow depth
   external noah271_getsnodepvars
   external noah271_transform_snodep
   external noah271_map_snodep
   external noah271_updatesnodep
   external noah271_qcsnodep
   external noah271_setsnodepvars
#endif
#endif

#if ( defined SM_NOAH_3_3 )
! Noah-3.3 soil moisture
   external noah33_getsoilm
   external noah33_setsoilm
   external noah33_getsmpred
   external noah33_getLbandTbPred
   external noah33_qcsoilm
   external noah33_qc_soilmobs
   external noah33_scale_soilm
   external noah33_descale_soilm
   external noah33_updatesoilm

! Noah-3.3 SWE
   external noah33_getsnowvars
   external noah33_setsnowvars
   external noah33_transform_snow
   external noah33_map_snow
   external noah33_getswepred
!    external noah33_getsnowpred_AMSREsnow  !yliu
   external noah33_getsnowpred_PMWsnow  !yliu
   external noah33_getsnwdpred
   external noah33_qcsnow
   external noah33_qc_snowobs
!    external noah33_qc_AMSREsnowobs !yliu
   external noah33_qc_PMWsnowobs !yliu
   external noah33_scale_snow
   external noah33_descale_snow
   external noah33_updatesnowvars

#if ( defined DA_OBS_SNODEP )
! Noah-3.3 snow depth
   external noah33_getsnodepvars
   external noah33_transform_snodep
   external noah33_map_snodep
   external noah33_updatesnodep
   external noah33_qcsnodep
   external noah33_setsnodepvars
   external noah33_getsnodeppred
   external noah33_scale_snodep
   external noah33_descale_snodep
   external noah33_qc_snodepobs
#endif

! Noah-3.3 snow-covered fraction
   external noah33_map_snow_DI          !yliu
   external noah33_updatesnowvars_scfda !yliu
   external noah33_getscfpred           !yliu
   external noah33_qc_scfobs            !yliu
#endif

#if ( defined SM_NOAH_3_6 )
! Noah-3.6 soil moisture
   external noah36_getsoilm
   external noah36_setsoilm
   external noah36_getsmpred
   external noah36_getLbandTbPred
   external noah36_qcsoilm
   external noah36_qc_soilmobs
   external noah36_scale_soilm
   external noah36_descale_soilm
   external noah36_updatesoilm

! Noah-3.6 SWE
   external noah36_getsnowvars
   external noah36_setsnowvars
   external noah36_transform_snow
   external noah36_map_snow
   external noah36_getswepred
!    external noah36_getsnowpred_AMSREsnow  !yliu
   external noah36_getsnowpred_PMWsnow  !yliu
   external noah36_getsnwdpred
   external noah36_qcsnow
   external noah36_qc_snowobs
!    external noah36_qc_AMSREsnowobs !yliu
   external noah36_qc_PMWsnowobs !yliu
   external noah36_scale_snow
   external noah36_descale_snow
   external noah36_updatesnowvars

#if ( defined DA_OBS_SNODEP )
! Noah-3.6 snow depth
   external noah36_getsnodepvars
   external noah36_transform_snodep
   external noah36_map_snodep
   external noah36_updatesnodep
   external noah36_qcsnodep
   external noah36_setsnodepvars
   external noah36_getsnodeppred
   external noah36_scale_snodep
   external noah36_descale_snodep
   external noah36_qc_snodepobs
#endif

! Noah-3.6 snow-covered fraction
   external noah36_map_snow_DI          !yliu
   external noah36_updatesnowvars_scfda !yliu
   external noah36_getscfpred           !yliu
   external noah36_qc_scfobs            !yliu
#endif

#if ( defined SM_NOAH_3_9 )
! Noah-3.6 soil moisture
   external noah39_getsoilm
   external noah39_setsoilm
   external noah39_getsmpred
   external noah39_getLbandTbPred
   external noah39_qcsoilm
   external noah39_qc_soilmobs
   external noah39_scale_soilm
   external noah39_descale_soilm
   external noah39_updatesoilm

#if ( defined DA_OBS_LDTSI )
! Noah-3.9 snow depth
   external noah39_getldtsivars
   external noah39_transform_ldtsi
   external noah39_map_ldtsi
   external noah39_updateldtsi
   external noah39_qcldtsi
   external noah39_setldtsivars
   external noah39_getldtsipred
   external noah39_scale_ldtsi
   external noah39_descale_ldtsi
   external noah39_qc_ldtsiobs
#endif

#if ( defined DA_OBS_SNODEP )
! Noah-3.9 snow depth
   external noah39_getsnodepvars
   external noah39_transform_snodep
   external noah39_map_snodep
   external noah39_updatesnodep
   external noah39_qcsnodep
   external noah39_setsnodepvars
   external noah39_getsnodeppred
   external noah39_scale_snodep
   external noah39_descale_snodep
   external noah39_qc_snodepobs
#endif

#endif

#if ( defined SM_NOAHMP_3_6 )
! MN: Noahmp-3.6 soil moisture
   external noahmp36_getsoilm           
   external noahmp36_setsoilm              
   external noahmp36_getsmpred
!   external noahmp36_getLbandTbPred   !I think we need this for real Lband DA
   external noahmp36_qcsoilm
   external noahmp36_qc_soilmobs
   external noahmp36_scale_soilm
   external noahmp36_descale_soilm
   external noahmp36_updatesoilm

   external noahmp36_getsnowvars
   external noahmp36_setsnowvars
   external noahmp36_getsnwdpred
   external noahmp36_getswepred
   external noahmp36_qcsnow
   external noahmp36_qc_snowobs
   external noahmp36_scale_snow
   external noahmp36_descale_snow
   external noahmp36_updatesnowvars

#if ( defined DA_OBS_SNODEP )
! NoahMP-3.6 SNODEP
   external noahmp36_getsnodepvars
   external noahmp36_transform_snodep
   !external noahmp36_map_snodep
   external noahmp36_updatesnodepvars
   external noahmp36_qcsnodep
   external noahmp36_setsnodepvars
   external noahmp36_getsnodeppred
   external noahmp36_scale_snodep
   external noahmp36_descale_snodep
   external noahmp36_qc_snodepobs
#endif

!NOAHMP3.6 TWS
   external noahmp36_gettws
   external noahmp36_settws
   external noahmp36_qctws
   external noahmp36_gettwspred
   external noahmp36_scale_tws
   external noahmp36_descale_tws
   external noahmp36_updatetws

   external noahmp36_getvegvars          
   external noahmp36_setvegvars  
   external noahmp36_transform_veg
   external noahmp36_map_veg
   external noahmp36_updatevegvars
   external noahmp36_qcveg
   external noahmp36_getLAIpred
   external noahmp36_qc_LAIobs
   external noahmp36_scale_veg
   external noahmp36_descale_veg

   external noahmp36_getalbedovars          
   external noahmp36_setalbedovars  
   external noahmp36_updatealbedovars
   external noahmp36_qcalbedo
   external noahmp36_getalbedopred
   external noahmp36_qc_albedoobs
   external noahmp36_scale_albedo
   external noahmp36_descale_albedo

   external noahmp36_transform_albedo
   external noahmp36_map_albedo
#endif

#if ( defined SM_NOAHMP_4_0_1 ) 
! MN NoahMP.4.0.1 Soil moisture DA
   external NoahMP401_getsoilm           
   external NoahMP401_setsoilm              
   external NoahMP401_getsmpred
!   external NoahMP401_getLbandTbPred   !we need this for Lband DA
   external NoahMP401_qcsoilm
   external NoahMP401_qc_soilmobs
   external NoahMP401_scale_soilm
   external NoahMP401_descale_soilm
   external NoahMP401_updatesoilm

#if ( defined DA_OBS_SNODEP )
! NoahMP-4.0.1 SNODEP
   external noahmp401_getsnodepvars
   external noahmp401_transform_snodep
   external noahmp401_map_snodep
   external noahmp401_updatesnodepvars
   external noahmp401_qcsnodep
   external noahmp401_setsnodepvars
   external noahmp401_getsnodeppred
   external noahmp401_scale_snodep
   external noahmp401_descale_snodep
   external noahmp401_qc_snodepobs
#endif

#if ( defined DA_OBS_LDTSI )
! NoahMP-4.0.1 LDTSI
   external noahmp401_getldtsivars
   external noahmp401_transform_ldtsi
   external noahmp401_map_ldtsi
   external noahmp401_updateldtsivars
   external noahmp401_qcldtsi
   external noahmp401_setldtsivars
   external noahmp401_getldtsipred
   external noahmp401_scale_ldtsi
   external noahmp401_descale_ldtsi
   external noahmp401_qc_ldtsiobs
#endif

#endif

#if ( defined SM_CLSM_F2_5 )
! CLSM-F2.5 soil moisture
   external clsmf25_getsoilm
   external clsmf25_setsoilm
   external clsmf25_getsmpred
   external clsmf25_qcsoilm
   external clsmf25_qc_soilmobs
   external clsmf25_scale_soilm
   external clsmf25_descale_soilm
   external clsmf25_updatesoilm
   external clsmf25_soilm_DAlog

! CLSM-F2.5 SWE
   external clsmf25_getsnowvars
   external clsmf25_setsnowvars
   external clsmf25_getsnwdpred
   external clsmf25_qcsnow
   external clsmf25_qc_snowobs
   external clsmf25_scale_snow
   external clsmf25_descale_snow
   external clsmf25_updatesnowvars
   external clsmf25_snow_DAlog

! CLSM-F2.5 TWS
   external clsmf25_gettws
   external clsmf25_settws
   external clsmf25_qctws
   external clsmf25_gettwspred
   external clsmf25_scale_tws
   external clsmf25_descale_tws
   external clsmf25_update_tws
#endif

#if ( defined SM_RUC_3_7 )
! RUC-3.7 soil moisture
   external ruc37_getsoilm
   external ruc37_setsoilm
   external ruc37_getsmpred
   external ruc37_getLbandTbPred
   external ruc37_qcsoilm
   external ruc37_qc_soilmobs
   external ruc37_scale_soilm
   external ruc37_descale_soilm
   external ruc37_updatesoilm
   external ruc37_getsnowvars
   external ruc37_setsnowvars
   external ruc37_transform_snow
   external ruc37_map_snow
   external ruc37_getswepred
!    external ruc37_getsnowpred_AMSREsnow  !yliu
   external ruc37_getsnowpred_PMWsnow  !yliu
   external ruc37_getsnwdpred
   external ruc37_qcsnow
   external ruc37_qc_snowobs
!    external ruc37_qc_AMSREsnowobs !yliu
   external ruc37_qc_PMWsnowobs !yliu
   external ruc37_scale_snow
   external ruc37_descale_snow
   external ruc37_updatesnowvars

#if 0
#if ( defined DA_OBS_SNODEP )
! RUC-3.7 snow depth
   external ruc37_getsnodepvars
   external ruc37_transform_snodep
   external ruc37_map_snodep
   external ruc37_updatesnodep
   external ruc37_qcsnodep
   external ruc37_setsnodepvars
   external ruc37_getsnodeppred
   external ruc37_scale_snodep
   external ruc37_descale_snodep
   external ruc37_qc_snodepobs
#endif
#endif

! RUC-3.7 snow-covered fraction
   external ruc37_map_snow_DI          !yliu
   external ruc37_updatesnowvars_scfda !yliu
   external ruc37_getscfpred           !yliu
   external ruc37_qc_scfobs            !yliu
#endif

#if ( defined SM_JULES_4_3 )
!MN Jules 43 soil moisture
   external jules43_getsoilm
   external jules43_setsoilm
   external jules43_getsmpred
   !external jules43_getLbandTbPred
   external jules43_qcsoilm
   external jules43_qc_soilmobs
   external jules43_scale_soilm
   external jules43_descale_soilm
   external jules43_updatesoilm
#endif

#if ( defined SM_JULES_5_0 )
!MN Jules 5.0 soil moisture
   external jules50_getsoilm
   external jules50_setsoilm
   external jules50_getsmpred
   !external jules50_getLbandTbPred
   external jules50_qcsoilm
   external jules50_qc_soilmobs
   external jules50_scale_soilm
   external jules50_descale_soilm
   external jules50_updatesoilm

! Yeosang Yoon SNODEP DA
#if ( defined DA_OBS_SNODEP )
   external jules50_getsnodepvars
   external jules50_transform_snodep
   external jules50_map_snodep
   external jules50_updatesnodep
   external jules50_qcsnodep
   external jules50_setsnodepvars
   external jules50_getsnodeppred
   external jules50_scale_snodep
   external jules50_descale_snodep
   external jules50_qc_snodepobs
#endif

#endif

#if ( defined SM_JULES_5_2 )
!MN Jules 5.0 soil moisture
   external jules52_getsoilm
   external jules52_setsoilm
   external jules52_getsmpred
   !external jules52_getLbandTbPred
   external jules52_qcsoilm
   external jules52_qc_soilmobs
   external jules52_scale_soilm
   external jules52_descale_soilm
   external jules52_updatesoilm
#endif

#if 0
#if ( defined SM_JULES_5_3)
!MN Jules 5.0 soil moisture
   external jules53_getsoilm
   external jules53_setsoilm
   external jules53_getsmpred
   !external jules53_getLbandTbPred
   external jules53_qcsoilm
   external jules53_qc_soilmobs
   external jules53_scale_soilm
   external jules53_descale_soilm
   external jules53_updatesoilm
#endif
#endif

#if ( defined SM_NOAH_2_7_1 )
#if ( defined DA_OBS_SNODEP )
! Noah-2.7.1 snow depth
   call registerlsmdainit(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_noah271Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah271_setsnodepvars)
#endif
#endif

#if ( defined SM_NOAH_3_3 )
! Noah-3.3 no observations soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah33_setsoilm)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah33_qcsoilm)

! Noah-3.3 NASA AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 LPRM AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 + SMOS NESDIS soil moisture

   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSNESDISsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 SMOS L2 soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOSL2smobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 ESACCI soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 RT SMOPS soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 RT SMOPS ASCAT soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 RT GCOMW AMSR2 soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3smobsId)//char(0),noah33_updatesoilm)
!MN
! Noah-3.3 NASA SMAP sm obs
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 ASCAT TU Wein soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 synthetic L-band brightness (soil moisture)
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_getLbandTbPred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah33_updatesoilm)

! Noah-3.3 synthetic soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah33_updatesoilm)

!MN: Noah-3.3 PILDAS soil moisture
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noah33_updatesoilm)



! Noah-3.3 ANSA SWE
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_getswepred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah33_qc_snowobs)


! Noah-3.3 PMW SWE and snow depth
! yliu, PMW-based SWE & snow depth assimilation
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_getsnowpred_PMWsnow)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_qc_PMWsnowobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah33_updatesnowvars)
! end ----------------- yliu

! yliu, for AMSR-E/ANSA (and other PMW-based) SWE & snow depth assimilation
!    call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_dasnow_init)
!    call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_getsnowvars)
!    call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_setsnowvars)
!    call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_getsnowpred_AMSREsnow)
!    call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_qcsnow)
!    call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_qc_AMSREsnowobs)
!    call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_scale_snow)
!    call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_descale_snow)
!    call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah33_updatesnowvars)
! end ----------------- yliu

#if ( defined DA_OBS_SNODEP )
! Noah-3.3 snow depth
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_getsnodeppred)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_descale_snodep)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah33_qc_snodepobs)
#endif

! Noah-3.3 ANSA snow depth
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah33_qc_snowobs)

   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_GCOMW_AMSR2L3sndobsId)//char(0),noah33_qc_snowobs)

! Noah-3.3 SSMR snow depth
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah33_qc_snowobs)

   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_qc_snowobs) 
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_synsndId)//char(0),noah33_qc_snowobs)

! Noah-3.3 SSMI snow depth
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah33_qc_snowobs)

! Noah-3.3 ANSA snow-covered fraction
   call registerlsmdainit(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_setsnowvars)
   call registerlsmdaupdatestate(trim(LIS_noah33Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),noah33_updatesnowvars)      !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_updatesnowvars_scfda) !yliu
   call registerlsmdaobstransform(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_transform_snow)
   call registerlsmdaqcstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_qcsnow)
   call registerlsmdamapobstolsm(trim(LIS_noah33Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),noah33_map_snow)           !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_map_snow_DI)         !yliu
!yliu -------------- adding files for SCF assimilation w/ EnKF
   call registerlsmdagetobspred(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_getscfpred)
   call registerlsmdascalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_descale_snow)
   call registerlsmdaqcobsstate(trim(LIS_noah33Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah33_qc_scfobs)
!yliu----------------------------------------------------------
#endif

#if ( defined SM_NOAH_3_6 )
! Noah-3.6 no observations soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah36_setsoilm)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_noobsId)//char(0),noah36_qcsoilm)

! Noah-3.6 NASA AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),noah36_updatesoilm)

! Noah-3.6 LPRM AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),noah36_updatesoilm)

! Noah-3.6 ESACCI soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),noah36_updatesoilm)

! Noah-3.6 RT SMOPS soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noah36_updatesoilm)

! MN 
! Noah-3.6 RT SMOPS ASCAT soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah36_updatesoilm)
!MN
! Noah-3.6 RT SMOPS SMOS soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noah36_updatesoilm)

!MN
! Noah-3.6 RT SMOPS AMSR2 soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noah36_updatesoilm)

!MN
! Noah-3.6 RT SMOPS SMAP soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noah36_updatesoilm)


! Noah-3.6 SMAP(NRT) soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah36_updatesoilm)
!MN
! Noah-3.6 SMAP(NASA) soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah36_updatesoilm)

! Noah-3.6 ASCAT TU Wein soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),noah36_updatesoilm)

! Noah-3.6 synthetic L-band brightness (soil moisture)
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_getLbandTbPred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),noah36_updatesoilm)

! Noah-3.6 synthetic soil moisture
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noah36_updatesoilm)

! Noah-3.6 ANSA SWE
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_getswepred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),noah36_qc_snowobs)

! Noah-3.6 PMW SWE and snow depth
! yliu, PMW-based SWE & snow depth assimilation
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_getsnowpred_PMWsnow)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_qc_PMWsnowobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),noah36_updatesnowvars)
! end ----------------- yliu

! yliu, for AMSR-E/ANSA (and other PMW-based) SWE & snow depth assimilation
!    call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_dasnow_init)
!    call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_getsnowvars)
!    call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_setsnowvars)
!    call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_getsnowpred_AMSREsnow)
!    call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_qcsnow)
!    call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_qc_AMSREsnowobs)
!    call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_scale_snow)
!    call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_descale_snow)
!    call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
!         trim(LIS_AMSREsnowobsId)//char(0),noah36_updatesnowvars)
! end ----------------- yliu

#if ( defined DA_OBS_SNODEP )
! Noah-3.6 snow depth
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_getsnodeppred)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_descale_snodep)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah36_qc_snodepobs)
#endif

! Noah-3.6 ANSA snow depth
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noah36_qc_snowobs)

! Noah-3.6 SSMR snow depth
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),noah36_qc_snowobs)

! Noah-3.6 SSMI snow depth
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),noah36_qc_snowobs)

! Noah-3.6 ANSA snow-covered fraction
   call registerlsmdainit(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_setsnowvars)
   call registerlsmdaupdatestate(trim(LIS_noah36Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),noah36_updatesnowvars)      !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_updatesnowvars_scfda) !yliu
   call registerlsmdaobstransform(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_transform_snow)
   call registerlsmdaqcstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_qcsnow)
   call registerlsmdamapobstolsm(trim(LIS_noah36Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),noah36_map_snow)           !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_map_snow_DI)         !yliu
!yliu -------------- adding files for SCF assimilation w/ EnKF
   call registerlsmdagetobspred(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_getscfpred)
   call registerlsmdascalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_descale_snow)
   call registerlsmdaqcobsstate(trim(LIS_noah36Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),noah36_qc_scfobs)
!yliu----------------------------------------------------------
#endif

#if ( defined SM_NOAH_3_9 )

! Noah-3.9 RT SMOPS ASCAT soil moisture! MN 
   call registerlsmdainit(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noah39_updatesoilm)

! Noah-3.9 SMAP(NRT) soil moisture !MN
   call registerlsmdainit(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noah39_updatesoilm)

! Noah-3.9 SMAP(NASA) soil moisture!MN
   call registerlsmdainit(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noah39_updatesoilm)

#if ( defined DA_OBS_LDTSI )
! Noah-3.9 snow depth
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_getldtsivars)
   call registerlsmdaobstransform(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_transform_ldtsi)
   call registerlsmdamapobstolsm(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_map_ldtsi)
   call registerlsmdaupdatestate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_updateldtsi)
   call registerlsmdaqcstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_qcldtsi)
   call registerlsmdasetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_setldtsivars)
   call registerlsmdagetobspred(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_getldtsipred)
   call registerlsmdascalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_scale_ldtsi)
   call registerlsmdadescalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_descale_ldtsi)
   call registerlsmdaqcobsstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noah39_qc_ldtsiobs)
#endif

#if ( defined DA_OBS_SNODEP )
! Noah-3.6 snow depth
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_getsnodeppred)
   call registerlsmdascalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_descale_snodep)
   call registerlsmdaqcobsstate(trim(LIS_noah39Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noah39_qc_snodepobs)
#endif

#endif


#if ( defined SM_NOAHMP_3_6 )
! Noahmp-3.6 synthetic soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsmId)//char(0),noahmp36_updatesoilm)

!MN: Noahmp-3.6 PILDAS soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),noahmp36_updatesoilm)

   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),noahmp36_updatesoilm)

   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_qc_snowobs) 
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsndId)//char(0),noahmp36_qc_snowobs)


! Yeosang Yoon, SNODEP
#if ( defined DA_OBS_SNODEP )
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_dasnodep_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_getsnodepvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_getsnodeppred)
   call registerlsmdaobstransform(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_transform_snodep)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_qcsnodep)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_qc_snodepobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_descale_snodep)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp36_updatesnodepvars)
#endif

   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_getswepred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_qc_snowobs) 
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_synsweId)//char(0),noahmp36_qc_snowobs)

#if ( defined DA_OBS_ASO_SWE)
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_getswepred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_qc_snowobs) 
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ASOsweobsId)//char(0),noahmp36_qc_snowobs)

#endif

! ANSA snow depth
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),noahmp36_qc_snowobs)

! Noah-MP.3.6 RT SMOPS soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),noahmp36_updatesoilm)

! MN
! Noah-MP.3.6 RT SMOPS_ASCAT soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),noahmp36_updatesoilm)

!MN
! Noah-MP.3.6 RT SMOPS_SMOS soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMOSsmobsId)//char(0),noahmp36_updatesoilm)

!MN
! Noah-MP.3.6 RT SMOPS_AMSR2 soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_AMSR2smobsId)//char(0),noahmp36_updatesoilm)

!MN
! Noah-MP.3.6 RT SMOPS_SMAP soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMOPS_SMAPsmobsId)//char(0),noahmp36_updatesoilm)

!MN
! Noah-MP.3.6 SMAP(NRT) soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),noahmp36_updatesoilm)
!MN
! Noah-MP.3.6 SMAP(NASA) soil moisture
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),noahmp36_updatesoilm)

!  TWS
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_datws_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_gettws)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_settws)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_gettwspred)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_qctws)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_scale_tws)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_descale_tws)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0), noahmp36_updatetws)
   call registerlsmdadiagnosevars(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),noahmp36_tws_DAlog)

!LAI
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_daveg_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_getvegvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_setvegvars)
   call registerlsmdaobstransform(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_transform_veg)
   call registerlsmdamapobstolsm(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_map_veg)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_updatevegvars)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_qcveg)

   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_getLAIpred)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_qc_LAIobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_scale_veg)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSlaiobsId)//char(0),noahmp36_descale_veg)

   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_daveg_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_getvegvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_setvegvars)
   call registerlsmdaobstransform(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_transform_veg)
   call registerlsmdamapobstolsm(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_map_veg)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_updatevegvars)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_qcveg)

   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_getLAIpred)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_qc_LAIobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_scale_veg)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_NASASMAPvodobsId)//char(0),noahmp36_descale_veg)


   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_daveg_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_getvegvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_setvegvars)
   call registerlsmdaobstransform(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_transform_veg)
   call registerlsmdamapobstolsm(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_map_veg)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_updatevegvars)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_qcveg)

   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_getLAIpred)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_qc_LAIobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_scale_veg)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_MODISsportLAIobsId)//char(0),noahmp36_descale_veg)

!albedo
   call registerlsmdainit(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_daalbedo_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_getalbedovars)
   call registerlsmdasetstatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_setalbedovars)
   call registerlsmdaupdatestate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_updatealbedovars)
   call registerlsmdaqcstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_qcalbedo)
   call registerlsmdagetobspred(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_getalbedopred)
   call registerlsmdaqcobsstate(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_qc_albedoobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_scale_albedo)
   call registerlsmdadescalestatevar(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_descale_albedo)

   call registerlsmdaobstransform(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_transform_albedo)
   call registerlsmdamapobstolsm(trim(LIS_noahmp36Id)//"+"//&
        trim(LIS_GLASSalbedoobsId)//char(0),noahmp36_map_albedo)

#endif


#if ( defined SM_NOAHMP_4_0_1 )
! MN
! Noah-MP.4.0.1 RT SMOPS ASCAT soil moisture
   call registerlsmdainit(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),NoahMP401_updatesoilm)

! Noah-MP.4.0.1 SMAP(NRT) soil moisture
   call registerlsmdainit(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),NoahMP401_updatesoilm)
!MN
! Noah-MP.4.0.1 SMAP(NASA) soil moisture
   call registerlsmdainit(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_setsoilm)
   call registerlsmdagetobspred(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_getsmpred)
   call registerlsmdaqcstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_NASASMAPsmobsId )//char(0),NoahMP401_updatesoilm)
! Yeosang Yoon, SNODEP DA
#if ( defined DA_OBS_SNODEP )
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_dasnodep_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_transform_snodep)
   call registerlsmdaupdatestate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_updatesnodepvars)
   call registerlsmdasetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_getsnodeppred)
   call registerlsmdaqcstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_qcsnodep)
   call registerlsmdaqcobsstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_qc_snodepobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),noahmp401_descale_snodep)
#endif

! Yeosang Yoon, LDTSI DA
#if ( defined DA_OBS_LDTSI )
! DA + ldtsi wirings
   call registerlsmdainit(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_daldtsi_init)
   call registerlsmdagetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_getldtsivars)
   call registerlsmdaobstransform(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_transform_ldtsi)
   call registerlsmdaupdatestate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_updateldtsivars)
   call registerlsmdasetstatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_setldtsivars)
   call registerlsmdagetobspred(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_getldtsipred)
   call registerlsmdaqcstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_qcldtsi)
   call registerlsmdaqcobsstate(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_qc_ldtsiobs)
   call registerlsmdascalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_scale_ldtsi)
   call registerlsmdadescalestatevar(trim(LIS_noahmp401Id)//"+"//&
        trim(LIS_ldtsiobsId)//char(0),noahmp401_descale_ldtsi)
#endif
#endif

#if ( defined SM_CLSM_F2_5 )
! CLSM-F2.5 synthetic soil moisture
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_setsoilm)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_getsmpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_synsmId)//char(0),clsmf25_updatesoilm)

! CLSM-F2.5 LPRM AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_setsoilm)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_getsmpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_updatesoilm)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),clsmf25_soilm_DAlog)

! CLSM-F2.5 RT SMOPS soil moisture
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_setsoilm)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_getsmpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_updatesoilm)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),clsmf25_soilm_DAlog)

! CLSM-F2.5 ESACCI soil moisture
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_setsoilm)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_getsmpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_updatesoilm)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ESACCIsmobsId)//char(0),clsmf25_soilm_DAlog)

! CLSM-F2.5 SSMR snow depth
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),clsmf25_snow_DAlog)

! CLSM-F2.5 SSMI snow depth
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),clsmf25_snow_DAlog)

! CLSM-F2.5 ANSA snow depth
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_qc_snowobs)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),clsmf25_snow_DAlog)

! CLSM-F2.5 TWS
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_datws_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_gettws)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_settws)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_gettwspred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_qctws)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_scale_tws)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_descale_tws)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0), clsmf25_update_tws)
!    call registerwritestatevar(trim(LIS_clsmf25Id)//"+"trim(LIS_GRACEtwsobsId)//char(0),clsmf25_write_tws)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_GRACEtwsobsId)//char(0),clsmf25_tws_DAlog)

! CLSM-F2.5 TWS
   call registerlsmdainit(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_datws_init)
   call registerlsmdagetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_gettws)
   call registerlsmdasetstatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_settws)
   call registerlsmdagetobspred(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_gettwspred)
   call registerlsmdaqcstate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_qctws)
   call registerlsmdascalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_scale_tws)
   call registerlsmdadescalestatevar(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_descale_tws)
   call registerlsmdaupdatestate(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0), clsmf25_update_tws)
!    call registerwritestatevar(trim(LIS_clsmf25Id)//"+"trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_write_tws)
   call registerlsmdadiagnosevars(trim(LIS_clsmf25Id)//"+"//&
        trim(LIS_simGRACEJPLobsId)//char(0),clsmf25_tws_DAlog)
#endif

#if ( defined SM_RUC_3_7 )
! RUC-3.7 no observations soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_noobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_noobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_noobsId)//char(0),ruc37_setsoilm)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_noobsId)//char(0),ruc37_qcsoilm)

! RUC-3.7 NASA AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_getsmpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_NASA_AMSREsmobsId)//char(0),ruc37_updatesoilm)

! RUC-3.7 LPRM AMSR-E soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_getsmpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_LPRM_AMSREsmobsId)//char(0),ruc37_updatesoilm)

! RUC-3.7 RT SMOPS soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_getsmpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMOPSsmobsId)//char(0),ruc37_updatesoilm)

#if 0
! RUC-3.7 ASCAT TU Wein soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_getsmpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ASCAT_TUWsmobsId)//char(0),ruc37_updatesoilm)
#endif

#if 0
! RUC-3.7 synthetic L-band brightness (soil moisture)
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_getLbandTbPred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synLbandTbobsId)//char(0),ruc37_updatesoilm)
#endif

#if 0
! RUC-3.7 synthetic soil moisture
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_setsoilm)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_getsmpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_synsmId)//char(0),ruc37_updatesoilm)
#endif

#if 0
! RUC-3.7 ANSA SWE
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_getswepred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASWEsnowobsId)//char(0),ruc37_qc_snowobs)
#endif

! RUC-3.7 PMW SWE and snow depth
! yliu, PMW-based SWE & snow depth assimilation
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_getsnowpred_PMWsnow)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_qc_PMWsnowobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_PMWsnowobsId)//char(0),ruc37_updatesnowvars)

#if 0
#if ( defined DA_OBS_SNODEP )
! RUC-3.7 snow depth
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_getsnodeppred)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_descale_snodep)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),ruc37_qc_snodepobs)
#endif
#endif

#if 0
! RUC-3.7 ANSA snow depth
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASNWDsnowobsId)//char(0),ruc37_qc_snowobs)
#endif

#if 0
! RUC-3.7 SSMR snow depth
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SMMRSNWDsnowobsId)//char(0),ruc37_qc_snowobs)
#endif

#if 0
! RUC-3.7 SSMI snow depth
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_getsnwdpred)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_qc_snowobs)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_updatesnowvars)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_SSMISNWDsnowobsId)//char(0),ruc37_qc_snowobs)
#endif

! RUC-3.7 ANSA snow-covered fraction
   call registerlsmdainit(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_dasnow_init)
   call registerlsmdagetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_getsnowvars)
   call registerlsmdasetstatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_setsnowvars)
   call registerlsmdaupdatestate(trim(LIS_ruc37Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_updatesnowvars)      !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_updatesnowvars_scfda) !yliu
   call registerlsmdaobstransform(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_transform_snow)
   call registerlsmdaqcstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_qcsnow)
   call registerlsmdamapobstolsm(trim(LIS_ruc37Id)//"+"//&
        !trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_map_snow)           !yliu
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_map_snow_DI)         !yliu
!yliu -------------- adding files for SCF assimilation w/ EnKF
   call registerlsmdagetobspred(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_getscfpred)
   call registerlsmdascalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_scale_snow)
   call registerlsmdadescalestatevar(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_descale_snow)
   call registerlsmdaqcobsstate(trim(LIS_ruc37Id)//"+"//&
        trim(LIS_ANSASCFsnowobsId)//char(0),ruc37_qc_scfobs)
!yliu----------------------------------------------------------
#endif

!MN: Jules43 PILDAS soil moisture
#if ( defined SM_JULES_4_3 )
   call registerlsmdainit(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules43Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules43_updatesoilm)
#endif

#if 0
!MN: Jules 5.0 PILDAS soil moisture
#if ( defined SM_JULES_5_0 )
   call registerlsmdainit(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_pildassmobsId)//char(0),jules50_updatesoilm)
#endif
#endif

!MN: Jules 5.0 SMAP(NASA) soil moisture
#if ( defined JULES_5_0_DEV )
!Disable JULES.5.0 DA; not ready yet.
#if ( defined SM_JULES_5_0 )
   call registerlsmdainit(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules50_updatesoilm)
#endif

! MN: Jules 5.0 SMAP(NRT) soil moisture
#if ( defined SM_JULES_5_0 )
   call registerlsmdainit(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules50_updatesoilm)
#endif

! MN 
! Jules 5.0 SMOPS ASCAT soil moisture
#if ( defined SM_JULES_5_0 )
   call registerlsmdainit(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules50_updatesoilm)
#endif

#if ( defined DA_OBS_SNODEP )
! Jules 5.0 snow depth, Yeosang Yoon
! DA + snodep wirings
   call registerlsmdainit(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_dasnodep_init)
   call registerlsmdagetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_getsnodepvars)
   call registerlsmdaobstransform(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_transform_snodep)
   call registerlsmdamapobstolsm(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_map_snodep)
   call registerlsmdaupdatestate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_updatesnodep)
   call registerlsmdaqcstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_qcsnodep)
   call registerlsmdasetstatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_setsnodepvars)
   call registerlsmdagetobspred(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_getsnodeppred)
   call registerlsmdascalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_scale_snodep)
   call registerlsmdadescalestatevar(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_descale_snodep)
   call registerlsmdaqcobsstate(trim(LIS_jules50Id)//"+"//&
        trim(LIS_snodepobsId)//char(0),jules50_qc_snodepobs)
#endif

#endif


!MN: Jules 5.2 SMAP(NASA) soil moisture
#if ( defined JULES_5_0_DEV )
!Disable JULES.5.2 DA; not ready yet.
#if ( defined SM_JULES_5_2)
   call registerlsmdainit(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules52_updatesoilm)
#endif

! MN: Jules 5.2 SMAP(NRT) soil moisture
#if ( defined SM_JULES_5_2 )
   call registerlsmdainit(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules52_updatesoilm)
#endif

! MN 
! Jules 5.2 SMOPS ASCAT soil moisture
#if ( defined SM_JULES_5_2 )
   call registerlsmdainit(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules52Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules52_updatesoilm)
#endif
#endif


#if 0 
!MN: Jules 5.3 SMAP(NASA) soil moisture
#if ( defined SM_JULES_5_3 )
   call registerlsmdainit(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_NASASMAPsmobsId)//char(0),jules53_updatesoilm)
#endif

! MN: Jules 5.3 SMAP(NRT) soil moisture
#if ( defined SM_JULES_5_3 )
   call registerlsmdainit(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMAPNRTsmobsId)//char(0),jules53_updatesoilm)
#endif

! MN 
! Jules 5.3 SMOPS ASCAT soil moisture
#if ( defined SM_JULES_5_3 )
   call registerlsmdainit(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_dasoilm_init)
   call registerlsmdagetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_getsoilm)
   call registerlsmdasetstatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_setsoilm)
   call registerlsmdagetobspred(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_getsmpred)
   call registerlsmdaqcstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_qcsoilm)
   call registerlsmdaqcobsstate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_qc_soilmobs)
   call registerlsmdascalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_scale_soilm)
   call registerlsmdadescalestatevar(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_descale_soilm)
   call registerlsmdaupdatestate(trim(LIS_jules53Id)//"+"//&
        trim(LIS_SMOPS_ASCATsmobsId)//char(0),jules53_updatesoilm)
#endif

#endif ! if 0 for jules53

#endif
end subroutine LIS_lsmda_plugin
end module LIS_lsmda_pluginMod
