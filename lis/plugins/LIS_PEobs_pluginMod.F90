!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_PEobs_pluginMod
!BOP
!
! !MODULE: LIS_PEobs_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions that are
!   used to read observation data that constitute the observation space
!   for parameter estimation. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  16 Jul 2009;   Sujay Kumar  Initial Specification
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_PEobs_plugin    
contains
!BOP
! !ROUTINE: LIS_PEobs_plugin
!  \label{LIS_PEobs_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing routines to handle the 
!  data that constitute the observation space for parameter estimation. 
!
! !INTERFACE:
subroutine LIS_PEobs_plugin
!EOP

#if ( ( defined OPTUE_ALG_ES )        || \
      ( defined OPTUE_ALG_LM )        || \
      ( defined OPTUE_ALG_GA )        || \
      ( defined OPTUE_ALG_SCEUA )     || \
      ( defined OPTUE_ALG_MCSIM )     || \
      ( defined OPTUE_ALG_RWMCMC )    || \
      ( defined OPTUE_ALG_DEMC )      || \
      ( defined OPTUE_ALG_DEMCZ ) )

   use LIS_pluginIndices

#if ( defined PE_OBS_TEMPLATE )
   use templateObs_module,    only : templateObs_setup
#endif

#if ( defined PE_OBS_PESYNSM1 )
   use pesynsm1data_module,   only : pesynsm1data_setup
#endif

#if ( defined PE_OBS_ISCCP_TSKIN )
   use ISCCP_Tskinobs_module, only : ISCCP_Tskinobs_setup
#endif

#if ( defined PE_OBS_WGPBMRSM )
   use wgPBMRsmobs_module,    only : wgPBMRsmdata_setup
#endif

#if ( defined PE_OBS_MACON_LS_DATA )
   use MaconLSDataMod,        only : maconlsobs_setup
#endif

#if ( defined PE_OBS_GLOBAL_LS_DATA )
   use GlobalLSDataMod,       only : globallsobs_setup
#endif

#if ( defined PE_OBS_AMERIFLUX )
   use AmerifluxObsMod,       only : Amerifluxobs_setup
#endif

#if ( defined PE_OBS_CNRS )
   use CNRS_em_obsMod,        only : CNRS_em_obs_setup
#endif

#if ( defined PE_OBS_USDA_ARSSM )
   use USDA_ARSsm_obsMod,     only : USDA_ARSsm_obs_setup
#endif

#if ( defined PE_OBS_ARSSM )
   use ARSsm_obsMod,     only : ARSsm_obs_setup
#endif

#if ( defined PE_OBS_ISMNSM )
   use ISMNsm_obsMod,     only : ISMNsm_obs_setup
#endif

#if ( defined PE_OBS_SMAPSM )
   use SMAPsm_obsMod,     only : SMAPsm_obs_setup
#endif

#if ( defined PE_OBS_AMSRE_SR )
   use AMSRE_SR_em_obsMod,    only : AMSRE_SR_em_obs_setup
#endif

#if ( defined PE_OBS_LPRM_AMSRESM )
   use LPRM_AMSREsm_obsMod,   only : LPRM_AMSREsm_obs_setup
#endif

#if ( defined PE_OBS_ARM )
   use ARMdata_module,        only : ARMdata_setup
#endif

#if ( defined PE_OBS_FLUXNET )
   use FLUXNETdata_module,    only : FLUXNETdata_setup
#endif

#if ( defined PE_OBS_EMPTYOBS )
   use Empty_obsMod
#endif

#if ( defined PE_OBS_TEMPLATE )
   external read_templateObs, write_templateObs, reset_templateObs
#endif

#if ( defined PE_OBS_PESYNSM1 )
   external read_pesynsm1data
#endif

#if ( defined PE_OBS_WGPBMRSM )
   external read_wgPBMRsmdata, write_wgPBMRsmdata, reset_wgPBMRsmdata
#endif

#if ( defined PE_OBS_ISCCP_TSKIN )
   external read_ISCCP_Tskindata
#endif

#if ( defined PE_OBS_MACON_LS_DATA )
   external read_MaconLSObsdata, write_MaconLSObsdata
#endif

#if ( defined PE_OBS_GLOBAL_LS_DATA )
   external read_GlobalLSObsdata, write_GlobalLSObsdata, reset_globallsObsdata
#endif

#if ( defined PE_OBS_AMERIFLUX )
   external read_AmerifluxObs, write_AmerifluxObs
#endif

#if ( defined PE_OBS_CNRS )
   external read_CNRS_em_obsdata, write_CNRS_em_obs_data,reset_CNRS_em_obs_data
#endif

#if ( defined PE_OBS_USDA_ARSSM )
   external read_USDA_ARSsm_obsdata, write_USDA_ARS_sm_obs_data, &
            reset_USDA_ARSsm_obs_data
#endif

#if ( defined PE_OBS_ARSSM )
   external read_ARSsmobs, write_ARSsmobs, reset_ARSsmobs
#endif

#if ( defined PE_OBS_ISMNSM )
   external read_ISMNsmobs, write_ISMNsmobs, reset_ISMNsmobs
#endif

#if ( defined PE_OBS_SMAPSM )
   external read_SMAPsmobs, write_SMAPsmobs, reset_SMAPsmobs
#endif

#if ( defined PE_OBS_AMSRE_SR )
   external read_AMSRE_SR_em_obsdata, write_AMSRE_SR_em_obs_data, &
            reset_AMSRE_SR_em_obs_data
#endif

#if ( defined PE_OBS_CNRS )
   external read_CNRS_MPDI_em_obsdata
#endif

#if ( defined PE_OBS_ARM )
   external read_ARMdata, write_ARMdata, reset_ARMdata
#endif

#if ( defined PE_OBS_LPRM_AMSRESM )
   external read_LPRM_AMSREsm_obsdata, write_LPRM_AMSREsm_obs_data, &
            reset_LPRM_AMSREsm_obs_data
#endif

#if ( defined PE_OBS_FLUXNET )
   external read_FLUXNETdata, write_FLUXNETdata, reset_FLUXNETdata
#endif

#if ( defined PE_OBS_TEMPLATE )
   call registerpeobssetup(trim(LIS_templateObsId)//char(0),templateObs_setup)
   call registergetpeobs(trim(LIS_templateObsId)//char(0),read_templateObs)
   call registerwritepeobs(trim(LIS_templateObsId)//char(0),write_templateObs)
   call registerpeobsreset(trim(LIS_templateObsId)//char(0),reset_templateObs)
#endif

#if ( defined PE_OBS_PESYNSM1 )
   call registerpeobssetup(trim(LIS_pesynsm1Id)//char(0),pesynsm1data_setup)
   call registergetpeobs(trim(LIS_pesynsm1Id)//char(0),read_pesynsm1data)
#endif

#if ( defined PE_OBS_WGPBMRSM )
   call registerpeobssetup(trim(LIS_wgPBMRsmId)//char(0),wgPBMRsmdata_setup)
   call registergetpeobs(trim(LIS_wgPBMRsmId)//char(0),read_wgPBMRsmdata)
   call registerwritepeobs(trim(LIS_wgPBMRsmId)//char(0),write_wgPBMRsmdata)
   call registerpeobsreset(trim(LIS_wgPBMRsmId)//char(0), reset_wgPBMRsmdata)
#endif

#if ( defined PE_OBS_PESYNSM1 )
   call registerpeobssetup(trim(LIS_pesynsm2Id)//char(0),pesynsm1data_setup)
   call registergetpeobs(trim(LIS_pesynsm2Id)//char(0),read_pesynsm1data)
#endif

#if ( defined PE_OBS_ISCCP_TSKIN )
   call registerpeobssetup(trim(LIS_isccpTskinId)//char(0),ISCCP_Tskinobs_setup)
   call registergetpeobs(trim(LIS_isccpTskinId)//char(0),read_ISCCP_Tskindata)
#endif

#if ( defined PE_OBS_MACON_LS_DATA )
   call registerpeobssetup(trim(LIS_maconlandslideObsId)//char(0), &
                           maconlsobs_setup)
   call registergetpeobs(trim(LIS_maconlandslideObsId)//char(0), &
                         read_maconlsObsdata)
   call registerwritepeobs(trim(LIS_maconlandslideObsId)//char(0), &
                           write_maconlsObsdata)
#endif

#if ( defined PE_OBS_GLOBAL_LS_DATA )
   call registerpeobssetup(trim(LIS_globallandslideObsId)//char(0), &
                           globallsobs_setup)
   call registergetpeobs(trim(LIS_globallandslideObsId)//char(0), &
                         read_globallsObsdata)
   call registerwritepeobs(trim(LIS_globallandslideObsId)//char(0), &
                           write_globallsObsdata)
   call registerpeobsreset(trim(LIS_globallandslideObsId)//char(0), &
                           reset_globallsObsdata)
#endif

#if ( defined PE_OBS_AMERIFLUX )
   call registerpeobssetup(trim(LIS_AmerifluxObsId)//char(0),AmerifluxObs_setup)
   call registergetpeobs(trim(LIS_AmerifluxObsId)//char(0),read_AmerifluxObs)
   call registerwritepeobs(trim(LIS_AmerifluxObsId)//char(0),write_AmerifluxObs)
#endif

#if ( defined PE_OBS_CNRS )
   call registerpeobssetup(trim(LIS_CNRSObsId)//char(0),CNRS_em_obs_setup)
   call registergetpeobs(trim(LIS_CNRSObsId)//char(0),read_CNRS_em_obsdata)
   call registerwritepeobs(trim(LIS_CNRSObsId)//char(0),write_CNRS_em_obs_data)
   call registerpeobsreset(trim(LIS_CNRSObsId)//char(0),reset_CNRS_em_obs_data)

   call registerpeobssetup(trim(LIS_CNRS_MPDIObsId)//char(0),CNRS_em_obs_setup)
   call registergetpeobs(trim(LIS_CNRS_MPDIObsId)//char(0), &
                         read_CNRS_MPDI_em_obsdata)
   call registerwritepeobs(trim(LIS_CNRS_MPDIObsId)//char(0), &
                           write_CNRS_em_obs_data)
   call registerpeobsreset(trim(LIS_CNRS_MPDIObsId)//char(0), &
                           reset_CNRS_em_obs_data)
#endif

#if ( defined PE_OBS_AMSRE_SR )
   call registerpeobssetup(trim(LIS_AMSRE_SRObsId)//char(0), &
                           AMSRE_SR_em_obs_setup)
   call registergetpeobs(trim(LIS_AMSRE_SRObsId)//char(0), &
                         read_AMSRE_SR_em_obsdata)
   call registerwritepeobs(trim(LIS_AMSRE_SRObsId)//char(0), &
                           write_AMSRE_SR_em_obs_data)
   call registerpeobsreset(trim(LIS_AMSRE_SRObsId)//char(0), &
                           reset_AMSRE_SR_em_obs_data)
#endif

#if ( defined PE_OBS_LPRM_AMSRESM )
   call registerpeobssetup(trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                           LPRM_AMSREsm_obs_setup)
   call registergetpeobs(trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                         read_LPRM_AMSREsm_obsdata)
   call registerwritepeobs(trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                           write_LPRM_AMSREsm_obs_data)
   call registerpeobsreset(trim(LIS_LPRM_AMSREsmpeObsId)//char(0), &
                           reset_LPRM_AMSREsm_obs_data)
#endif

#if ( defined PE_OBS_FLUXNET )
   call registerpeobssetup(trim(LIS_FLUXNETpeObsId)//char(0),FLUXNETdata_setup)
   call registergetpeobs(trim(LIS_FLUXNETpeObsId)//char(0),read_FLUXNETdata)
   call registerwritepeobs(trim(LIS_FLUXNETpeObsId)//char(0),write_FLUXNETdata)
   call registerpeobsreset(trim(LIS_FLUXNETpeObsId)//char(0),reset_FLUXNETdata)
#endif

#if ( defined PE_OBS_EMPTYOBS )
   call registerpeobssetup(trim(LIS_EmptyObsId)//char(0),Empty_obs_setup)
   call registergetpeobs(trim(LIS_EmptyObsId)//char(0),read_Empty_em_obsdata)
   call registerwritepeobs(trim(LIS_EmptyObsId)//char(0),write_EmptyObsdata)
#endif

#if ( defined PE_OBS_ARM )
   call registerpeobssetup(trim(LIS_ARMObsId)//char(0),ARMdata_setup)
   call registergetpeobs(trim(LIS_ARMObsId)//char(0),read_ARMdata)
   call registerwritepeobs(trim(LIS_ARMObsId)//char(0),write_ARMdata)
   call registerpeobsreset(trim(LIS_ARMObsId)//char(0),reset_ARMdata)
#endif

#if ( defined PE_OBS_USDA_ARSSM )
   call registerpeobssetup(trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                           USDA_ARSsm_obs_setup)
   call registergetpeobs(trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                         read_USDA_ARSsm_obsdata)
   call registerwritepeobs(trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                           write_USDA_ARS_sm_obs_data)
   call registerpeobsreset(trim(LIS_USDA_ARSsmpeObsId)//char(0), &
                           reset_USDA_ARSsm_obs_data)
#endif

#if ( defined PE_OBS_ARSSM )
   call registerpeobssetup(trim(LIS_ARSsmobsId)//char(0), &
                           ARSsm_obs_setup)
   call registergetpeobs(trim(LIS_ARSsmobsId)//char(0), &
                         read_ARSsmobs)
   call registerwritepeobs(trim(LIS_ARSsmobsId)//char(0), &
                           write_ARSsmobs)
   call registerpeobsreset(trim(LIS_ARSsmobsId)//char(0), &
                           reset_ARSsmobs)
#endif

#if ( defined PE_OBS_ISMNSM )
   call registerpeobssetup(trim(LIS_ISMNsmobsId)//char(0), &
                           ISMNsm_obs_setup)
   call registergetpeobs(trim(LIS_ISMNsmobsId)//char(0), &
                         read_ISMNsmobs)
   call registerwritepeobs(trim(LIS_ISMNsmobsId)//char(0), &
                           write_ISMNsmobs)
   call registerpeobsreset(trim(LIS_ISMNsmobsId)//char(0), &
                           reset_ISMNsmobs)
#endif

#if ( defined PE_OBS_ISMNSM )
   call registerpeobssetup(trim(LIS_ISMNsmobsId)//char(0), &
                           ISMNsm_obs_setup)
   call registergetpeobs(trim(LIS_ISMNsmobsId)//char(0), &
                         read_ISMNsmobs)
   call registerwritepeobs(trim(LIS_ISMNsmobsId)//char(0), &
                           write_ISMNsmobs)
   call registerpeobsreset(trim(LIS_ISMNsmobsId)//char(0), &
                           reset_ISMNsmobs)
#endif

#if ( defined PE_OBS_SMAPSM )
   call registerpeobssetup(trim(LIS_SMAPsmobsId)//char(0), &
                           SMAPsm_obs_setup)
   call registergetpeobs(trim(LIS_SMAPsmobsId)//char(0), &
                         read_SMAPsmobs)
   call registerwritepeobs(trim(LIS_SMAPsmobsId)//char(0), &
                           write_SMAPsmobs)
   call registerpeobsreset(trim(LIS_SMAPsmobsId)//char(0), &
                           reset_SMAPsmobs)
#endif

#endif
end subroutine LIS_PEobs_plugin
end module LIS_PEobs_pluginMod
