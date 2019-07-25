!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: LVT_datastream_pluginMod
!  \label(LVT_datastream_pluginMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   defining routines that initialize various LVT-obss. 
!   The user defined functions are incorporated into 
!   the appropriate registry to be later invoked through generic calls. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  17 Feb 2004;   Sujay Kumar  Initial Specification
!  17 Oct 2018  Mahdi Navari  Enhanced the LVT reader to read the 
!               Veg. Water Content (VWC) from SMAP SM dataset ! 
!  19 Nov 2018  Mahdi Navari added suport to read SMAP_L3 brightness temperature
!
!EOP
module LVT_datastream_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_datastream_plugin  
contains
!BOP
! 
! !ROUTINE: LVT_datastream_plugin
!  \label{LVT_datastream_plugin}
!
! !INTERFACE:
  subroutine LVT_datastream_plugin
    use LVT_pluginIndices
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
!  This is a plugin point for introducing a new LVT-obs. 
!  The interface mandates that the following interfaces be implemented
!  and registered for each LVT-obs. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    use template_obsMod,        only : template_obsinit
    use LISoutputMod,           only : LISoutputInit
    use LIS6outputMod,          only : LIS6outputInit
    use LISDAdiagOutputMod,     only : LISDAdiagoutputInit
    use SCANGMAO_obsMod,        only : SCANGMAO_obsInit
    use SCAN_obsMod,            only : SCAN_obsInit
    use LISda_obsMod,           only : LISda_obsInit
    use CEOP_obsMod,            only : ceop_obsInit
    use ISCCP_TskinobsMod,      only : ISCCP_TskinObsInit
    use MODIS_LSTobsMod,        only : MODIS_LSTobsInit
    use NASMD_obsMod,           only : NASMD_obsInit
    use SURFRAD_obsMod,         only : SURFRAD_obsInit
    use WGPBMRobsMod,           only : WGPBMRobsInit  
    use SNOTEL_obsMod,          only : SNOTEL_obsInit
    use FMISWE_obsMod,          only : FMISWE_obsInit
    use CMCSNWD_obsMod,         only : CMCSNWD_obsInit
    use SNODAS_obsMod,          only : SNODAS_obsInit
    use NASA_AMSREsm_obsMod,    only : NASA_AMSREsm_obsInit
    use LPRM_AMSREsm_obsMod,    only : LPRM_AMSREsm_obsInit
    use Ameriflux_obsMod,       only : Ameriflux_obsInit
    use ARM_obsMod,             only : ARM_obsInit
    use SMOSREX_obsMod,         only : SMOSREX_obsInit
    use AGRMET_dataMod,         only : AGRMET_dataInit
    use GlobSnow_obsMod,        only : GlobSnow_obsInit
    use SNODEP_metobsMod,       only : SNODEP_metobsInit
    use SNODEP_obsMod,          only : SNODEP_obsInit
    use MOD10A1_obsMod,         only : MOD10A1_obsInit
    use MODSCAG_obsMod,         only : MODSCAG_obsInit
    use MOD10A1V6_obsMod,       only : MOD10A1V6_obsInit
    use ANSASNWD_obsMod,        only : ANSASNWD_obsInit
    use ANSASWE_obsMod,         only : ANSASWE_obsInit
    use CPCPRCP_obsMod,         only : CPCPRCP_obsInit
    use USGSSF_obsMod,          only : USGSSF_obsInit
    use USGSSFgrid_obsMod,      only : USGSSFgrid_obsInit
    use NatSF_obsMod,           only : NatSF_obsInit
    use ISMN_obsMod,            only : ISMN_obsInit
    use FLUXNETmte_obsMod,      only : FLUXNETmte_obsInit
    use MOD16A2_obsMod,         only : MOD16A2_obsInit
    use UWET_obsMod,            only : UWET_obsInit
    use ARSsm_obsMod,           only : ARSsm_obsInit
    use NLDAS2_dataMod,         only : NLDAS2_dataInit
    use GHCN_obsMod,            only : GHCN_obsInit
    use ALEXI_obsMod,           only : ALEXI_obsInit
    use GRACE_obsMod,           only : GRACE_obsInit
    use simGRACE_obsMod,        only : simGRACE_obsInit
    use USGSGWwell_obsMod,      only : USGSGWwell_obsInit
    use PBOH2O_obsMod,          only : PBOH2O_obsInit
    use SMOSL2sm_obsMod,        only : SMOSL2sm_obsInit
    use SMOSNESDIS_smobsMod,    only : SMOSNESDIS_smobsinit
    use SMOSCATDS_smobsMod,     only : SMOSCATDS_smobsinit
    use SMOSL1TB_obsMod,        only : SMOSL1TB_obsInit
    use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3sm_obsinit
    use GCOMW_AMSR2L3snd_obsMod,only : GCOMW_AMSR2L3snd_obsinit
    use SMOPSsm_obsMod,         only : SMOPSsm_obsInit
    use ESACCIsm_obsMod,        only : ESACCIsm_obsInit
    use GIMMSAVHRR_NDVIobsMod,  only : GIMMSAVHRR_NDVIobsinit
    use GIMMSMODIS_NDVIobsMod,  only : GIMMSMODIS_NDVIobsinit
    use GLDAS1obsMod,           only : GLDAS1obsinit
    use GLDAS2obsMod,           only : GLDAS2obsinit
    use MERRA2obsMod,           only : MERRA2obsinit
    use MERRA2asmobsMod,        only : MERRA2asmObsinit
    use MERRAlandObsMod,        only : MERRAlandObsinit
    use ERAinterimLandobsMod,   only : ERAinterimLandobsinit
    use SSEBop_obsMod,          only : SSEBop_obsinit
    use GRDC_obsMod,            only : GRDC_obsinit
    use GOES_LSTobsMod,         only : GOES_LSTobsInit
    use GLERL_dataMod,          only : GLERL_obsinit
    use JULES_obsMod,           only : JULES_obsinit
    use JULES2D_obsMod,         only : JULES2D_obsinit
    use LVTbenchmarkOUT_obsMod, only : LVTbenchmarkOUT_obsInit
    use SMAP_smobsMod,          only : SMAP_smobsinit
    use SMAP_vwcobsMod,         only : SMAP_vwcobsinit !MN
    use SMAP_vodobsMod,         only : SMAP_vodobsinit
    use SMAP_L3TBMod,           only : SMAP_L3TBinit  !MN   
    use SMAP_TBobsMod,          only : SMAP_TBobsinit
    use GOME2_SIFobsMod,        only : GOME2_SIFobsinit
    use Daymet_obsMod,          only : Daymet_obsInit
    use CMORPH_dataMod,         only : CMORPH_dataInit    
    use CHIRPSv2_dataMod,       only : CHIRPSv2_dataInit
    use USCRNsm_obsMod,         only : USCRNsm_obsinit
    use FLUXNET2015_obsMod,     only : FLUXNET2015_obsinit
    use FLUXNET2015NC_obsMod,   only : FLUXNET2015NC_obsinit
    use GLEAM_obsMod,           only : GLEAM_obsinit 
    use USDM_obsMod,            only : USDM_obsinit 
    use LVTpercentile_obsMod,   only : LVTpercentile_obsInit
    use IMDPRCP_obsMod,         only : IMDPRCP_obsInit
    use APHROPRCP_obsMod,       only : APHROPRCP_obsInit
    use TRMM3B42V7_dataMod,     only : TRMM3B42V7_dataInit
    use GLASSlaiobsMod,         only : GLASSlaiobsinit
    use GLASSalbedoobsMod,      only : GLASSalbedoobsinit
    use MODISsportLAIobsMod,    only : MODISsportLAIobsinit
    use FLUXCOM_obsMod,         only : FLUXCOM_obsInit
    use HAR_dataMod,            only : HAR_datainit
    use OCO2_SIFobsMod,         only : OCO2_SIFobsinit
    use ECMWFforc_dataMod,      only : ECMWFforc_datainit    
    use GDASforc_dataMod,       only : GDASforc_datainit    
    use ASOSWE_obsMod,          only : ASOSWE_obsinit
    use IMERG_dataMod,          only : IMERG_datainit
    use UASNOW_obsMod,          only : UASNOW_obsinit
 
    external readtemplateObs
    external readLISoutput
    external readLIS6output
    external readLISDAdiagOutput
    external readSCANGMAOObs
    external readSCANObs
    external readLISdaAsObs  
    external readCEOP
    external readISCCP_TskinObs
    external readMODIS_LSTObs
    external readNASMDObs
    external readSURFRADObs
    external readWGPBMRObs
    external readSNOTELObs
    external readFMISWEobs
    external readCMC_SNWDobs
    external readSNODASobs
    external readNASA_AMSREsmObs
    external readLPRM_AMSREsmObs
    external readAmerifluxObs
    external readARMobs
    external readSMOSREXobs
    external readAGRMETdata
    external readGlobSnowObs
    external readSNODEPmetobs
    external readSNODEPobs
    external readMOD10A1obs
    external readMODSCAGobs
    external readMOD10A1V6obs
    external readANSASNWDobs
    external readANSASWEobs
    external readCPCPRCPobs
    external readUSGSSFobs
    external readUSGSSFgridobs
    external readNatSFobs
    external readISMNobs
    external readFLUXNETmteObs
    external readMOD16A2Obs
    external readUWETObs
    external readARSsmobs
    external readNLDAS2data
    external readGHCNObs
    external readALEXIobs
    external readGRACEObs
    external readsimGRACEObs
    external readUSGSGWwellobs
    external readPBOH2Oobs
    external readSMOSL2smobs
    external readSMOSNESDISsmobs
    external readSMOSCATDSsmobs
    external readSMOSL1TBobs
    external readGCOMW_AMSR2L3smobs
    external readGCOMW_AMSR2L3sndobs
    external readSMOPSsmobs
    external readESACCIsmobs
    external readGLERLobs
    external readGLDAS1obs
    external readGLDAS2obs
    external readMERRA2obs
    external readMERRA2asmObs
    external readMERRAlandObs
    external readERAinterimLandobs
    external readSSEBopObs
    external readGRDCobs
    external readGOES_LSTobs
    external readJULESobs
    external readJULES2Dobs
    external readGIMMSAVHRR_NDVIobs
    external readGIMMSMODIS_NDVIobs
    external readLVTbenchmarkOUTobs
    external readSMAPsmobs
    external readSMAPvodobs
    external readSMAPvwcobs ! MN vegwtation water content
    external readSMAP_L3TB ! MN Tb from SMAP SM data       
    external readSMAPTBobs
    external readGOME2_SIFobs
    external readDaymetObs
    external readCMORPHdata
    external readCHIRPSv2data
    external readUSCRNsmObs
    external readFLUXNET2015Obs
    external readFLUXNET2015NCObs
    external readGLEAMObs
    external readUSDMobs
    external readLVTpercentileAsObs
    external readIMDPRCPobs
    external readAPHROPRCPobs
    external readTRMM3B42V7data
    external readGLASSlaiObs
    external readGLASSalbedoObs
    external readMODISsportLAIObs
    external readFLUXCOMobs
    external readHARdata
    external readOCO2_SIFobs
    external readECMWFforcdata
    external readGDASforcdata
    external readASOSWEObs
    external readIMERGdata
    external readUASNOWObs

    call registerobsread(trim(LVT_LVTbenchmarkobsId)//char(0),&
         readLVTbenchmarkOUTobs)

    call registerobssetup(trim(LVT_templateobsId)//char(0), template_obsInit)
    call registerobsread(trim(LVT_templateobsId)//char(0),readtemplateObs)

    call registerobssetup(trim(LVT_LISoutputId)//char(0), LISoutputInit)
    call registerobsread(trim(LVT_LISoutputId)//char(0),readLISoutput)

    call registerobssetup(trim(LVT_LIS6outputId)//char(0), LIS6outputInit)
    call registerobsread(trim(LVT_LIS6outputId)//char(0),readLIS6output)

    call registerobssetup(trim(LVT_LISDAdiagoutputId)//char(0),&
         LISDAdiagoutputInit)
    call registerobsread(trim(LVT_LISDAdiagoutputId)//char(0),&
         readLISDAdiagoutput)

    call registerobssetup(trim(LVT_SCANGMAOobsId)//char(0), SCANGMAO_obsinit)
    call registerobsread(trim(LVT_SCANGMAOobsId)//char(0),readSCANGMAOObs)

    call registerobssetup(trim(LVT_SCANobsId)//char(0), SCAN_obsinit)
    call registerobsread(trim(LVT_SCANobsId)//char(0),readSCANObs)

    call registerobssetup(trim(LVT_LISdaobsId)//char(0), LISda_obsInit)
    call registerobsread(trim(LVT_LISdaobsId)//char(0),readLISdaAsObs)

    call registerobssetup(trim(LVT_ceopobsId)//char(0),CEOP_obsInit)
    call registerobsread(trim(LVT_ceopobsId)//char(0),readCEOP)

    call registerobssetup(trim(LVT_ISCCP_TskinobsId)//char(0),&
         ISCCP_TskinobsInit)
    call registerobsread(trim(LVT_ISCCP_TskinobsId)//char(0),&
         readISCCP_TskinObs)


    call registerobssetup(trim(LVT_NASMDobsId)//char(0), NASMD_obsinit)
    call registerobsread(trim(LVT_NASMDobsId)//char(0),readNASMDObs)

    call registerobssetup(trim(LVT_SURFRADobsId)//char(0), SURFRAD_obsinit)
    call registerobsread(trim(LVT_SURFRADobsId)//char(0),readSURFRADObs)

    call registerobssetup(trim(LVT_wgPBMRobsId)//char(0), WGPBMRobsinit)
    call registerobsread(trim(LVT_wgPBMRobsId)//char(0),readWGPBMRObs)

    call registerobssetup(trim(LVT_SNOTELobsId)//char(0), SNOTEL_obsinit)
    call registerobsread(trim(LVT_SNOTELobsId)//char(0),readSNOTELObs)

    call registerobssetup(trim(LVT_FMISWEobsId)//char(0), FMISWE_obsinit)
    call registerobsread(trim(LVT_FMISWEobsId)//char(0),readFMISWEobs)

    call registerobssetup(trim(LVT_CMCSNWDobsId)//char(0), CMCSNWD_obsinit)
    call registerobsread(trim(LVT_CMCSNWDobsId)//char(0),readCMC_SNWDobs)

    call registerobssetup(trim(LVT_SNODASobsId)//char(0), SNODAS_obsinit)
    call registerobsread(trim(LVT_SNODASobsId)//char(0),readSNODASobs)

    call registerobssetup(trim(LVT_NASAAMSREsmobsId)//char(0), &
         NASA_AMSREsm_obsinit)
    call registerobsread(trim(LVT_NASAAMSREsmobsId)//char(0),&
         readNASA_AMSREsmObs)

    call registerobssetup(trim(LVT_LPRMAMSREsmobsId)//char(0), &
         LPRM_AMSREsm_obsinit)
    call registerobsread(trim(LVT_LPRMAMSREsmobsId)//char(0),&
         readLPRM_AMSREsmObs)

    call registerobssetup(trim(LVT_AmerifluxobsId)//char(0), Ameriflux_obsinit)
    call registerobsread(trim(LVT_AmerifluxobsId)//char(0),readAmerifluxObs)

    call registerobssetup(trim(LVT_ARMobsId)//char(0), ARM_obsinit)
    call registerobsread(trim(LVT_ARMobsId)//char(0),readARMObs)

    call registerobssetup(trim(LVT_SMOSREXobsId)//char(0), SMOSREX_obsinit)
    call registerobsread(trim(LVT_SMOSREXobsId)//char(0),readSMOSREXObs)

    call registerobssetup(trim(LVT_AGRMETdataId)//char(0), AGRMET_datainit)
    call registerobsread(trim(LVT_AGRMETdataId)//char(0),readAGRMETdata)

    call registerobssetup(trim(LVT_GlobSnowObsId)//char(0), GlobSnow_obsinit)
    call registerobsread(trim(LVT_GlobSnowObsId)//char(0),readGlobSnowObs)

    call registerobssetup(trim(LVT_SNODEPmetObsId)//char(0), SNODEP_metobsinit)
    call registerobsread(trim(LVT_SNODEPmetObsId)//char(0),readSNODEPmetObs)

    call registerobssetup(trim(LVT_SNODEPobsId)//char(0), SNODEP_obsinit)
    call registerobsread(trim(LVT_SNODEPobsId)//char(0),readSNODEPobs)

    call registerobssetup(trim(LVT_MOD10A1obsId)//char(0), MOD10A1_obsinit)
    call registerobsread(trim(LVT_MOD10A1obsId)//char(0),readMOD10A1obs)

    call registerobssetup(trim(LVT_MODSCAGobsId)//char(0), MODSCAG_obsinit)
    call registerobsread(trim(LVT_MODSCAGobsId)//char(0),readMODSCAGobs)

    call registerobssetup(trim(LVT_MOD10A1V6obsId)//char(0), MOD10A1V6_obsinit)
    call registerobsread(trim(LVT_MOD10A1V6obsId)//char(0),readMOD10A1V6obs)

    call registerobssetup(trim(LVT_ANSASNWDobsId)//char(0), ANSASNWD_obsinit)
    call registerobsread(trim(LVT_ANSASNWDobsId)//char(0),readANSASNWDobs)

    call registerobssetup(trim(LVT_ANSASWEobsId)//char(0), ANSASWE_obsinit)
    call registerobsread(trim(LVT_ANSASWEobsId)//char(0),readANSASWEobs)

    call registerobssetup(trim(LVT_CPCPRCPobsId)//char(0), CPCPRCP_obsinit)
    call registerobsread(trim(LVT_CPCPRCPobsId)//char(0),readCPCPRCPobs)

    call registerobssetup(trim(LVT_USGSSFobsId)//char(0), USGSSF_obsinit)
    call registerobsread(trim(LVT_USGSSFobsId)//char(0),readUSGSSFobs)

    call registerobssetup(trim(LVT_USGSSFgridobsId)//char(0), &
         USGSSFgrid_obsinit)
    call registerobsread(trim(LVT_USGSSFgridobsId)//char(0),&
         readUSGSSFgridobs)

    call registerobssetup(trim(LVT_NatSFobsId)//char(0), NatSF_obsinit)
    call registerobsread(trim(LVT_NatSFobsId)//char(0),readNatSFobs)

    call registerobssetup(trim(LVT_ISMNobsId)//char(0), ISMN_obsinit)
    call registerobsread(trim(LVT_ISMNobsId)//char(0),readISMNobs)

    call registerobssetup(trim(LVT_FLUXNETmteobsId)//char(0), &
         FLUXNETmte_obsinit)
    call registerobsread(trim(LVT_FLUXNETmteobsId)//char(0),&
         readFLUXNETmteobs)

    call registerobssetup(trim(LVT_MOD16A2obsId)//char(0), MOD16A2_obsinit)
    call registerobsread(trim(LVT_MOD16A2obsId)//char(0),readMOD16A2obs)

    call registerobssetup(trim(LVT_UWETobsId)//char(0), UWET_obsinit)
    call registerobsread(trim(LVT_UWETobsId)//char(0),readUWETobs)

    call registerobssetup(trim(LVT_ARSsmobsId)//char(0), ARSsm_obsinit)
    call registerobsread(trim(LVT_ARSsmobsId)//char(0),readARSsmobs)

    call registerobssetup(trim(LVT_NLDAS2obsId)//char(0), NLDAS2_datainit)
    call registerobsread(trim(LVT_NLDAS2obsId)//char(0),readNLDAS2data)

    call registerobssetup(trim(LVT_GHCNobsId)//char(0), GHCN_obsinit)
    call registerobsread(trim(LVT_GHCNobsId)//char(0),readGHCNobs)

    call registerobssetup(trim(LVT_ALEXIobsId)//char(0), ALEXI_obsinit)
    call registerobsread(trim(LVT_ALEXIobsId)//char(0),readALEXIobs)

    call registerobssetup(trim(LVT_GRACEobsId)//char(0), GRACE_obsinit)
    call registerobsread(trim(LVT_GRACEobsId)//char(0),readGRACEObs)

    call registerobssetup(trim(LVT_simGRACEobsId)//char(0), simGRACE_obsinit)
    call registerobsread(trim(LVT_simGRACEobsId)//char(0),readsimGRACEObs)

    call registerobssetup(trim(LVT_USGSGWwellobsId)//char(0), USGSGWwell_obsinit)
    call registerobsread(trim(LVT_USGSGWwellobsId)//char(0),readUSGSGWwellobs)

    call registerobssetup(trim(LVT_PBOH2OobsId)//char(0), PBOH2O_obsinit)
    call registerobsread(trim(LVT_PBOH2OobsId)//char(0),readPBOH2Oobs)

    call registerobssetup(trim(LVT_SMOSL2smobsId)//char(0), SMOSL2sm_obsinit)
    call registerobsread(trim(LVT_SMOSL2smobsId)//char(0),readSMOSL2smobs)

    call registerobssetup(trim(LVT_SMOSNESDISsmobsId)//char(0), &
         SMOSNESDIS_smobsinit)
    call registerobsread(trim(LVT_SMOSNESDISsmobsId)//char(0),&
         readSMOSNESDISsmobs)

    call registerobssetup(trim(LVT_SMOSCATDSsmobsId)//char(0), &
         SMOSCATDS_smobsinit)
    call registerobsread(trim(LVT_SMOSCATDSsmobsId)//char(0),&
         readSMOSCATDSsmobs)

    call registerobssetup(trim(LVT_SMOSL1TBobsId)//char(0), SMOSL1TB_obsinit)
    call registerobsread(trim(LVT_SMOSL1TBobsId)//char(0),readSMOSL1TBobs)

    call registerobssetup(trim(LVT_GCOMW_AMSR2L3smobsId)//char(0), &
         GCOMW_AMSR2L3sm_obsinit)
    call registerobsread(trim(LVT_GCOMW_AMSR2L3smobsId)//char(0),&
         readGCOMW_AMSR2L3smobs)

    call registerobssetup(trim(LVT_GCOMW_AMSR2L3sndobsId)//char(0), &
         GCOMW_AMSR2L3snd_obsinit)
    call registerobsread(trim(LVT_GCOMW_AMSR2L3sndobsId)//char(0),&
         readGCOMW_AMSR2L3sndobs)

    call registerobssetup(trim(LVT_SMOPSsmobsId)//char(0), &
         SMOPSsm_obsinit)
    call registerobsread(trim(LVT_SMOPSsmobsId)//char(0),&
         readSMOPSsmobs)

    call registerobssetup(trim(LVT_ESACCIsmobsId)//char(0), &
         ESACCIsm_obsinit)
    call registerobsread(trim(LVT_ESACCIsmobsId)//char(0),&
         readESACCIsmobs)

    call registerobssetup(trim(LVT_GIMMSAVHRR_NDVIobsId)//char(0), &
         GIMMSAVHRR_NDVIobsinit)
    call registerobsread(trim(LVT_GIMMSAVHRR_NDVIobsId)//char(0),&
         readGIMMSAVHRR_NDVIobs)

    call registerobssetup(trim(LVT_GIMMSMODIS_NDVIobsId)//char(0), &
         GIMMSMODIS_NDVIobsinit)
    call registerobsread(trim(LVT_GIMMSMODIS_NDVIobsId)//char(0),&
         readGIMMSMODIS_NDVIobs)

    call registerobssetup(trim(LVT_MODIS_LSTobsId)//char(0), &
         MODIS_LSTobsinit)
    call registerobsread(trim(LVT_MODIS_LSTobsId)//char(0),&
         readMODIS_LSTobs)

    call registerobssetup(trim(LVT_GLDAS2obsId)//char(0), &
         GLDAS2obsinit)
    call registerobsread(trim(LVT_GLDAS2obsId)//char(0),&
         readGLDAS2obs)

    call registerobssetup(trim(LVT_GLDAS1obsId)//char(0), &
         GLDAS1obsinit)
    call registerobsread(trim(LVT_GLDAS1obsId)//char(0),&
         readGLDAS1obs)

    call registerobssetup(trim(LVT_MERRA2obsId)//char(0), &
         MERRA2obsinit)
    call registerobsread(trim(LVT_MERRA2obsId)//char(0),&
         readMERRA2obs)

    call registerobssetup(trim(LVT_MERRA2asmObsId)//char(0), &
         MERRA2asmObsinit)
    call registerobsread(trim(LVT_MERRA2asmObsId)//char(0),&
         readMERRA2asmObs)

    call registerobssetup(trim(LVT_MERRAlandObsId)//char(0), &
         MERRAlandObsinit)
    call registerobsread(trim(LVT_MERRAlandobsId)//char(0),&
         readMERRAlandObs)

    call registerobssetup(trim(LVT_ERAIlandobsId)//char(0), &
         ERAinterimLandobsinit)
    call registerobsread(trim(LVT_ERAIlandobsId)//char(0),&
         readERAinterimLandobs)

    call registerobssetup(trim(LVT_SSEBopobsId)//char(0), &
         SSEBop_obsinit)
    call registerobsread(trim(LVT_SSEBopobsId)//char(0),&
         readSSEBopObs)

    call registerobssetup(trim(LVT_GRDCobsId)//char(0), &
         GRDC_obsinit)
    call registerobsread(trim(LVT_GRDCobsId)//char(0),&
         readGRDCObs)

    call registerobssetup(trim(LVT_GOESLSTobsId)//char(0), &
         GOES_LSTobsinit)
    call registerobsread(trim(LVT_GOESLSTobsId)//char(0),&
         readGOES_LSTobs)

    call registerobssetup(trim(LVT_GLERLobsId)//char(0), &
         GLERL_obsinit)
    call registerobsread(trim(LVT_GLERLobsId)//char(0),&
         readGLERLobs)

    call registerobssetup(trim(LVT_JULESobsId)//char(0), &
         JULES_obsinit)
    call registerobsread(trim(LVT_JULESobsId)//char(0),&
         readJULESObs)

    call registerobssetup(trim(LVT_JULES2DobsId)//char(0), &
         JULES2D_obsinit)
    call registerobsread(trim(LVT_JULES2DobsId)//char(0),&
         readJULES2DObs)

    call registerobssetup(trim(LVT_LVTbenchmarkobsId)//char(0), &
         LVTbenchmarkOUT_obsInit)
    call registerobsread(trim(LVT_LVTbenchmarkobsId)//char(0),&
         readLVTbenchmarkOUTobs)


    call registerobssetup(trim(LVT_SMAPsmobsId)//char(0), &
         SMAP_smobsInit)
    call registerobsread(trim(LVT_SMAPsmobsId)//char(0),&
         readSMAPsmobs)

    call registerobssetup(trim(LVT_SMAPvodobsId)//char(0), &
         SMAP_vodobsInit)
    call registerobsread(trim(LVT_SMAPvodobsId)//char(0),&
         readSMAPvodobs)


! MN: SMAP vegetation water content
!    call registerobssetup(trim(LVT_SMAPvwcobsId)//char(0), &
!         SMAP_vwcobsInit)
!    call registerobsread(trim(LVT_SMAPvwcobsId)//char(0),&
!         readSMAPvwcobs)

! MN: SMAP L3 Tb
    call registerobssetup(trim(LVT_SMAP_L3TbId)//char(0), &
            SMAP_L3TBinit)
    call registerobsread(trim(LVT_SMAP_L3TbId)//char(0),&
            readSMAP_L3TB)

    call registerobssetup(trim(LVT_SMAPTBobsId)//char(0), &
         SMAP_TBobsInit)
    call registerobsread(trim(LVT_SMAPTBobsId)//char(0),&
         readSMAPTBobs)


    call registerobssetup(trim(LVT_GOME2SIFobsId)//char(0), &
         GOME2_SIFobsinit)
    call registerobsread(trim(LVT_GOME2SIFobsId)//char(0),&
         readGOME2_SIFobs)

    call registerobssetup(trim(LVT_DaymetobsId)//char(0), &
         Daymet_obsinit)
    call registerobsread(trim(LVT_DaymetobsId)//char(0),&
         readDaymetObs)

    call registerobssetup(trim(LVT_CMORPHdataId)//char(0), CMORPH_datainit)
    call registerobsread(trim(LVT_CMORPHdataId)//char(0),readCMORPHdata)

    call registerobssetup(trim(LVT_CHIRPSv2dataId)//char(0), CHIRPSv2_datainit)
    call registerobsread(trim(LVT_CHIRPSv2dataId)//char(0),readCHIRPSv2data)
 
    call registerobssetup(trim(LVT_USCRNsmdataId)//char(0), USCRNsm_obsinit)
    call registerobsread(trim(LVT_USCRNsmdataId)//char(0),readUSCRNsmObs)

    call registerobssetup(trim(LVT_FLUXNET2015ObsId)//char(0), &
         FLUXNET2015_obsinit)
    call registerobsread(trim(LVT_FLUXNET2015ObsId)//char(0), &
         readFLUXNET2015Obs)

    call registerobssetup(trim(LVT_FLUXNET2015NCObsId)//char(0), &
         FLUXNET2015NC_obsinit)
    call registerobsread(trim(LVT_FLUXNET2015NCObsId)//char(0), &
         readFLUXNET2015NCObs)

    call registerobssetup(trim(LVT_GLEAMdataId)//char(0), &
         GLEAM_obsinit)
    call registerobsread(trim(LVT_GLEAMdataId)//char(0), &
         readGLEAMobs)

    call registerobssetup(trim(LVT_USDMdataId)//char(0), &
         USDM_obsinit)
    call registerobsread(trim(LVT_USDMdataId)//char(0), &
         readUSDMobs)

    call registerobssetup(trim(LVT_LVTpercentiledataId)//char(0), &
         LVTpercentile_obsInit)
    call registerobsread(trim(LVT_LVTpercentiledataId)//char(0), &
         readLVTpercentileAsObs)

    call registerobssetup(trim(LVT_IMDPRCPdataId)//char(0), IMDPRCP_obsinit)
    call registerobsread(trim(LVT_IMDPRCPdataId)//char(0),readIMDPRCPobs)

    call registerobssetup(trim(LVT_APHROPRCPdataId)//char(0), APHROPRCP_obsinit)
    call registerobsread(trim(LVT_APHROPRCPdataId)//char(0),readAPHROPRCPobs)

    call registerobssetup(trim(LVT_TRMM3B42V7dataId)//char(0), &
         TRMM3B42V7_datainit)
    call registerobsread(trim(LVT_TRMM3B42V7dataId)//char(0),&
         readTRMM3B42V7data)

    call registerobssetup(trim(LVT_GLASSlaiobsId)//char(0), GLASSlaiobsinit)
    call registerobsread(trim(LVT_GLASSlaiobsId)//char(0),readGLASSlaiobs)

    call registerobssetup(trim(LVT_GLASSalbedoobsId)//char(0), GLASSalbedoobsinit)
    call registerobsread(trim(LVT_GLASSalbedoobsId)//char(0),readGLASSalbedoobs)

    call registerobssetup(trim(LVT_MODISsportLAIobsId)//char(0), &
         MODISsportLAIobsinit)
    call registerobsread(trim(LVT_MODISsportLAIobsId)//char(0),&
         readMODISsportLAIobs)

    call registerobssetup(trim(LVT_FLUXCOMobsId)//char(0), &
         FLUXCOM_obsinit)
    call registerobsread(trim(LVT_FLUXCOMobsId)//char(0),&
         readFLUXCOMobs)

    call registerobssetup(trim(LVT_HARdataId)//char(0), &
         HAR_datainit)
    call registerobsread(trim(LVT_HARdataId)//char(0),&
         readHARdata)

    call registerobssetup(trim(LVT_OCO2SIFobsId)//char(0), &
         OCO2_SIFobsinit)
    call registerobsread(trim(LVT_OCO2SIFobsId)//char(0),&
         readOCO2_SIFobs)

    call registerobssetup(trim(LVT_ECMWFdataId)//char(0), &
        ECMWFforc_datainit) 
    call registerobsread(trim(LVT_ECMWFdataId)//char(0),&
         readECMWFforcdata)         

    call registerobssetup(trim(LVT_GDASdataId)//char(0), &
         GDASforc_datainit) 
    call registerobsread(trim(LVT_GDASdataId)//char(0),&
         readGDASforcdata)         

    call registerobssetup(trim(LVT_ASOSWEdataId)//char(0), &
         ASOSWE_obsinit) 
    call registerobsread(trim(LVT_ASOSWEdataId)//char(0),&
         readASOSWEObs)

    call registerobssetup(trim(LVT_IMERGdataId)//char(0), IMERG_datainit)
    call registerobsread(trim(LVT_IMERGdataId)//char(0) , readIMERGdata)

    call registerobssetup(trim(LVT_UASNOWdataId)//char(0), UASNOW_obsinit)
    call registerobsread(trim(LVT_UASNOWdataId)//char(0) , readUASNOWObs)

  end subroutine LVT_datastream_plugin
end module LVT_datastream_pluginMod
