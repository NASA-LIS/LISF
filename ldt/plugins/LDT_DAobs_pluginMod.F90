!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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
    use NASASMAPvod_obsMod,        only : NASASMAPvod_obsinit
    use GLASSlai_obsMod,           only : GLASSlai_obsinit
    use LPRMvod_obsMod,            only : LPRMvod_obsinit

    external readLISlsmSMObs
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
    external readNASASMAPvodObs
    external readGLASSlaiObs
    external readLPRMvodObs

    call registerdaobssetup(trim(LDT_LISlsmSMobsId)//char(0), LISlsmSM_obsInit)
    call registerdaobsread(trim(LDT_LISlsmSMobsId)//char(0), readLISlsmSMObs)

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

  end subroutine LDT_DAobs_plugin
end module LDT_DAobs_pluginMod
