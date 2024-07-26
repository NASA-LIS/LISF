!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DAobs_pluginMod
!BOP
!
! !MODULE: LDT_DAobs_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LDT-obss. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
!   
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_DAobs_plugin  
contains
!BOP
! !ROUTINE: LDT_DAobs_plugin
!  \label{LDT_DAobs_plugin}
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LDT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LDT-obs. 
!
!
! !INTERFACE:
  subroutine LDT_DAobs_plugin
    use LDT_pluginIndices
!EOP

    use LISlsmSM_obsMod,           only : LISlsmSM_obsInit
    use LISlsmTEFF_obsMod,         only : LISlsmTEFF_obsInit  !Y.Kwon
    use syntheticsm_obsMod,        only : syntheticsm_obsinit
    use NASA_AMSREsm_obsMod,       only : NASA_AMSREsm_obsInit
    use LPRM_AMSREsm_obsMod,       only : LPRM_AMSREsm_obsInit
    use ESACCIsm_obsMod,           only : ESACCIsm_obsInit
    use GRACEtws_obsMod,           only : GRACEtws_obsInit
    use GRACEQLtws_obsMod,         only : GRACEQLtws_obsInit
    use WindSatsm_obsMod,          only : WindSatsm_obsInit
    use SMOPSsm_obsMod,            only : SMOPSsm_obsInit
    use ASCATTUWsm_obsMod,         only : ASCATTUWsm_obsInit
    use SMOSL2sm_obsMod,           only : SMOSL2sm_obsInit
    use SMOSNESDISsm_obsMod,       only : SMOSNESDISsm_obsInit
    use AquariusL2sm_obsMod,       only : AquariusL2sm_obsInit
    use GCOMW_AMSR2L3sm_obsMod,    only : GCOMW_AMSR2L3sm_obsinit
    use simGRACEJPL_obsMod,        only : simGRACEJPL_obsinit
    use SMMRSNWDsnow_obsMod,       only : SMMRSNWDsnow_obsInit
    use SSMISNWDsnow_obsMod,       only : SSMISNWDsnow_obsInit
    use ANSASNWDsnow_obsMod,       only : ANSASNWDsnow_obsInit
    use GCOMW_AMSR2L3snd_obsMod,   only : GCOMW_AMSR2L3snd_obsInit
    use NASASMAPsm_obsMod,         only : NASASMAPsm_obsinit
    use SMOSNRTNNL2sm_obsMod,      only : SMOSNRTNNL2sm_obsinit   !Y.Kwon
    use NASASMAPvod_obsMod,        only : NASASMAPvod_obsinit
    use GLASSlai_obsMod,           only : GLASSlai_obsinit
    use LPRMvod_obsMod,            only : LPRMvod_obsinit
    use MCD15A2Hlai_obsMod,        only : MCD15A2Hlai_obsinit
    use THySM_obsMod,              only : THySM_obsinit
    use LISlsmPrecip_obsMod,       only : LISlsmPrecip_obsInit 
    use VIIRSGVFobsMod,            only : VIIRSGVFobsinit        !Y.Kwon   
    use CDFSGVFobsMod,             only : CDFSGVFobsinit         !Y.Kwon
    use GEOSTEFF_obsMod,           only : GEOSTeffobsinit        !Y.Kwon    
    use SMAPEOPLSMobsMod,          only : SMAPEOPLSMobsinit       !Y.Kwon

    external readLISlsmSMObs
    external readLISlsmTEFFObs    !Y.Kwon
    external readsyntheticsmobs
    external readNASA_AMSREsmObs
    external readLPRM_AMSREsmObs
    external readESACCIsmObs
    external readWindSatsmObs
    external readGRACEtwsObs
    external readGRACEQLtwsObs
    external readSMOPSsmObs
    external readASCATTUWsmObs
    external readSMOSL2smObs
    external readSMOSNESDISsmObs
    external readAquariusL2smObs
    external readGCOMW_AMSR2L3smObs
    external readsimGRACEJPLObs
    external readSMMRSNWDsnowObs
    external readSSMISNWDsnowObs
    external readANSASNWDsnowObs
    external readGCOMW_AMSR2L3sndObs
    external readNASASMAPsmObs
    external readSMOSNRTNNL2smObs      !Y.Kwon
    external readNASASMAPvodObs
    external readGLASSlaiObs
    external readLPRMvodObs
    external readMCD15A2HlaiObs
    external readTHySMobs
    external readLISlsmPrecipObs
    external readVIIRS_GVFObs       !Y.Kwon
    external readCDFS_GVFObs        !Y.Kwon
    external readGEOSTEFFObs        !Y.Kwon
    external readSMAPEOPL_SMObs      !Y.Kwon

    call registerdaobssetup(trim(LDT_LISlsmSMobsId)//char(0), LISlsmSM_obsInit)
    call registerdaobsread(trim(LDT_LISlsmSMobsId)//char(0), readLISlsmSMObs)

    call registerdaobssetup(trim(LDT_LISlsmPrecipobsId)//char(0), LISlsmPrecip_obsInit)
    call registerdaobsread(trim(LDT_LISlsmPrecipobsId)//char(0), readLISlsmPrecipObs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_LISlsmTEFFobsId)//char(0), LISlsmTEFF_obsInit)
    call registerdaobsread(trim(LDT_LISlsmTEFFobsId)//char(0), readLISlsmTEFFObs)
 
    call registerdaobssetup(trim(LDT_syntheticSMobsId)//char(0), &
         syntheticSM_obsinit)
    call registerdaobsread(trim(LDT_syntheticSMobsId)//char(0),&
         readsyntheticSMObs)

    call registerdaobssetup(trim(LDT_NASA_AMSREsmobsId)//char(0), &
         NASA_AMSREsm_obsinit)
    call registerdaobsread(trim(LDT_NASA_AMSREsmobsId)//char(0),&
         readNASA_AMSREsmObs)

    call registerdaobssetup(trim(LDT_LPRM_AMSREsmobsId)//char(0),&
         LPRM_AMSREsm_obsinit)
    call registerdaobsread(trim(LDT_LPRM_AMSREsmobsId)//char(0),&
         readLPRM_AMSREsmObs)

    call registerdaobssetup(trim(LDT_ESACCIsmobsId)//char(0),&
         ESACCIsm_obsinit)
    call registerdaobsread(trim(LDT_ESACCIsmobsId)//char(0),&
         readESACCIsmObs)

    call registerdaobssetup(trim(LDT_GRACEtwsobsId)//char(0),&
         GRACEtws_obsinit)
    call registerdaobsread(trim(LDT_GRACEtwsobsId)//char(0),&
         readGRACEtwsObs)

    call registerdaobssetup(trim(LDT_GRACEQLtwsobsId)//char(0),&
         GRACEQLtws_obsinit)
    call registerdaobsread(trim(LDT_GRACEQLtwsobsId)//char(0),&
         readGRACEQLtwsObs)

    call registerdaobssetup(trim(LDT_WindSatsmobsId)//char(0),&
         WindSatsm_obsinit)
    call registerdaobsread(trim(LDT_WindSatsmobsId)//char(0),&
         readWindSatsmObs)

    call registerdaobssetup(trim(LDT_SMOPSsmobsId)//char(0),&
         SMOPSsm_obsinit)
    call registerdaobsread(trim(LDT_SMOPSsmobsId)//char(0),&
         readSMOPSsmObs)

    call registerdaobssetup(trim(LDT_ASCATTUWsmobsId)//char(0),&
         ASCATTUWsm_obsinit)
    call registerdaobsread(trim(LDT_ASCATTUWsmobsId)//char(0),&
         readASCATTUWsmObs)

    call registerdaobssetup(trim(LDT_SMOSL2smobsId)//char(0),&
         SMOSL2sm_obsinit)
    call registerdaobsread(trim(LDT_SMOSL2smobsId)//char(0),&
         readSMOSL2smObs)

    call registerdaobssetup(trim(LDT_SMOSNESDISsmobsId)//char(0),&
         SMOSNESDISsm_obsinit)
    call registerdaobsread(trim(LDT_SMOSNESDISsmobsId)//char(0),&
         readSMOSNESDISsmObs)

    call registerdaobssetup(trim(LDT_GCOMW_AMSR2L3smobsId)//char(0),&
         GCOMW_AMSR2L3sm_obsinit)
    call registerdaobsread(trim(LDT_GCOMW_AMSR2L3smobsId)//char(0),&
         readGCOMW_AMSR2L3smObs)

    call registerdaobssetup(trim(LDT_AquariusL2smobsId)//char(0),&
         AquariusL2sm_obsinit)
    call registerdaobsread(trim(LDT_AquariusL2smobsId)//char(0),&
         readAquariusL2smObs)

    call registerdaobssetup(trim(LDT_simGRACEJPLobsId)//char(0),&
         simGRACEJPL_obsinit)
    call registerdaobsread(trim(LDT_simGRACEJPLobsId)//char(0),&
         readsimGRACEJPLObs)

    call registerdaobssetup(trim(LDT_SMMRSNWDsnowobsId)//char(0),&
         SMMRSNWDsnow_obsInit)
    call registerdaobsread(trim(LDT_SMMRSNWDsnowobsId)//char(0),&
         readSMMRSNWDsnowObs)

    call registerdaobssetup(trim(LDT_SSMISNWDsnowobsId)//char(0),&
         SSMISNWDsnow_obsInit)
    call registerdaobsread(trim(LDT_SSMISNWDsnowobsId)//char(0),&
         readSSMISNWDsnowObs)

    call registerdaobssetup(trim(LDT_ANSASNWDsnowobsId)//char(0),&
         ANSASNWDsnow_obsInit)
    call registerdaobsread(trim(LDT_ANSASNWDsnowobsId)//char(0),&
         readANSASNWDsnowObs)

    call registerdaobssetup(trim(LDT_GCOMWAMSR2L3sndobsId)//char(0),&
         GCOMW_AMSR2L3snd_obsInit)
    call registerdaobsread(trim(LDT_GCOMWAMSR2L3sndobsId)//char(0),&
         readGCOMW_AMSR2L3sndObs)

    call registerdaobssetup(trim(LDT_NASASMAPsmobsId)//char(0),&
         NASASMAPsm_obsinit)
    call registerdaobsread(trim(LDT_NASASMAPsmobsId)//char(0),&
         readNASASMAPsmObs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_SMOSNRTNNsmobsId)//char(0),&
         SMOSNRTNNL2sm_obsinit)
    call registerdaobsread(trim(LDT_SMOSNRTNNsmobsId)//char(0),&
         readSMOSNRTNNL2smObs)

    call registerdaobssetup(trim(LDT_NASASMAPvodobsId)//char(0),&
         NASASMAPvod_obsinit)
    call registerdaobsread(trim(LDT_NASASMAPvodobsId)//char(0),&
         readNASASMAPvodObs)

    call registerdaobssetup(trim(LDT_GLASSlaiobsId)//char(0),&
         GLASSlai_obsinit)
    call registerdaobsread(trim(LDT_GLASSlaiobsId)//char(0),&
         readGLASSlaiObs)


    call registerdaobssetup(trim(LDT_LPRMvodobsId)//char(0),&
         LPRMvod_obsinit)
    call registerdaobsread(trim(LDT_LPRMvodobsId)//char(0),&
         readLPRMvodObs)

    call registerdaobssetup(trim(LDT_MCD15A2HlaiobsId)//char(0),&
         MCD15A2Hlai_obsinit)
    call registerdaobsread(trim(LDT_MCD15A2HlaiobsId)//char(0),&
         readMCD15A2HlaiObs)

    call registerdaobssetup(trim(LDT_THySMobsId)//char(0),&
         THySM_obsinit)
    call registerdaobsread(trim(LDT_THySMobsId)//char(0),&
         readTHySMobs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_VIIRSgvfobsId)//char(0),&
         VIIRSGVFobsinit)
    call registerdaobsread(trim(LDT_VIIRSgvfobsId)//char(0),&
         readVIIRS_GVFObs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_CDFSgvfobsId)//char(0),&
         CDFSGVFobsinit)
    call registerdaobsread(trim(LDT_CDFSgvfobsId)//char(0),&
         readCDFS_GVFObs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_GEOSTeffobsId)//char(0),&
         GEOSTeffobsinit)
    call registerdaobsread(trim(LDT_GEOSTeffobsId)//char(0),&
         readGEOSTEFFObs)

    !Y.Kwon
    call registerdaobssetup(trim(LDT_SMAPEOPLsmobsId)//char(0),&
         SMAPEOPLSMobsinit)
    call registerdaobsread(trim(LDT_SMAPEOPLsmobsId)//char(0),&
         readSMAPEOPL_SMObs)

  end subroutine LDT_DAobs_plugin
end module LDT_DAobs_pluginMod
