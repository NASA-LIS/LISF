!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: LVT_pluginIndices
!  \label(LVT_pluginIndices)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   The code in this file provides values of indices used to 
!   to register functions in the plugin modules
!
!   The index definitions are simply a convention
!   The user may change these options, and the lis.config 
!   should be changed appropriately to ensure that the correct function
!   is called at run time
! 
!   NOTES: The indices for metrics should be in increasing order, whereas
!   the indices for other plugin sets need not be. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Oct 2008    Sujay Kumar  Initial Specification
! 14 Nov 2017 Jossy Jacob Added MOD10A1V6 (MOD10A1_V006)
!  17 Oct 2018  Mahdi Navari  Enhanced the LVT reader to read the 
!               Veg. Water Content (VWC) from SMAP SM dataset ! 
!  19 Nov 2018  Mahdi Navari added suport to read SMAP_L3 brightness temperature
!  10 Jan 2023  Mahdi Navari added suport for COAMPSout 
!
!EOP
module LVT_pluginIndices

  PRIVATE
   
!BOC
!-------------------------------------------------------------------------
! supported metrics 
!------------------------------------------------------------------------- 
  integer, public,  parameter :: LVT_METRIC_SINDEX = 1
  integer, public,  parameter :: LVT_MEANid        = 1
  integer, public,  parameter :: LVT_Minid        = 2
  integer, public,  parameter :: LVT_Maxid        = 3
  integer, public,  parameter :: LVT_Sumid        = 4
  integer, public,  parameter :: LVT_Stdevid      = 5
  integer, public,  parameter :: LVT_Varianceid    = 6
  integer, public,  parameter :: LVT_Anomalyid     = 7
  integer, public,  parameter :: LVT_RMSEid        = 8
  integer, public,  parameter :: LVT_BIASid        = 9
  integer, public,  parameter :: LVT_MAEid         = 10
  integer, public,  parameter :: LVT_PODYid        = 11
  integer, public,  parameter :: LVT_PODNid        = 12
  integer, public,  parameter :: LVT_FARid         = 13
  integer, public,  parameter :: LVT_POFDid        = 14
  integer, public,  parameter :: LVT_CSIid         = 15
  integer, public,  parameter :: LVT_ACCid         = 16
  integer, public,  parameter :: LVT_FBIASid       = 17
  integer, public,  parameter :: LVT_ETSid         = 18
  integer, public,  parameter :: LVT_Rcorrid       = 19
  integer, public,  parameter :: LVT_Rnkcorrid     = 20
  integer, public,  parameter :: LVT_ARnkcorrid    = 21
  integer, public,  parameter :: LVT_Acorrid       = 22
  integer, public,  parameter :: LVT_ARMSEid       = 23
  integer, public,  parameter :: LVT_NSEid         = 24
  integer, public,  parameter :: LVT_ubRMSEid      = 25
  integer, public,  parameter :: LVT_AREAid        = 26
  integer, public,  parameter :: LVT_waveletStatId = 27
  integer, public,  parameter :: LVT_hnId          = 28
  integer, public,  parameter :: LVT_spiId         = 29
  integer, public,  parameter :: LVT_sriId         = 30
  integer, public,  parameter :: LVT_sswiId        = 31
  integer, public,  parameter :: LVT_sgwiId        = 32
  integer, public,  parameter :: LVT_percentileId  = 33
  integer, public,  parameter :: LVT_RFVid         = 34
  integer, public,  parameter :: LVT_MinTimeid       = 35
  integer, public,  parameter :: LVT_MaxTimeid       = 36
  integer, public,  parameter :: LVT_Tendencyid      = 37
  integer, public,  parameter :: LVT_TendencyCorrid  = 38
  integer, public,  parameter :: LVT_Zscoreid        = 39
  integer, public,  parameter :: LVT_Trendid        = 40
  integer, public,  parameter :: LVT_SdSIId         = 41
  integer, public,  parameter :: LVT_TCId           = 42
  integer, public,  parameter :: LVT_DFRid         = 43   ! EMK
  integer, public,  parameter :: LVT_EFid          = 44    ! EMK
  integer, public,  parameter :: LVT_FFid          = 45    ! EMK
  integer, public,  parameter :: LVT_HSSid          = 46    ! EMK
  integer, public,  parameter :: LVT_PSSid          = 47    ! EMK
  integer, public,  parameter :: LVT_CSSid          = 48    ! EMK
  integer, public,  parameter :: LVT_RELid          = 49
  integer, public,  parameter :: LVT_RESid          = 50
  integer, public,  parameter :: LVT_VULid          = 51
  integer, public,  parameter :: LVT_KMEANSid       = 52

  ! Tian decomposition of mean error...EMK
  integer, public,  parameter :: LVT_THBid        = 53
  integer, public,  parameter :: LVT_TMBid        = 54
  integer, public,  parameter :: LVT_TFBid        = 55
  integer, public,  parameter :: LVT_IEid        = 56
  integer, public,  parameter :: LVT_CEid        = 57

  integer, public,  parameter :: LVT_REid        = 58
  integer, public,  parameter :: LVT_JEid        = 59
  integer, public,  parameter :: LVT_MIid         = 60
  integer, public,  parameter :: LVT_METRIC_EINDEX   = 60

!Information content metrics
!EMK...These are always registered, so they must have unique values

  integer, public,  parameter :: LVT_ICMETRIC_SINDEX = 61
  integer, public,  parameter :: LVT_mentropyid      = 61
  integer, public,  parameter :: LVT_igainid         = 62
  integer, public,  parameter :: LVT_fcomplexityid   = 63
  integer, public,  parameter :: LVT_ecomplexityid   = 64
  integer, public,  parameter :: LVT_ICMETRIC_EINDEX = 65

!ensemble metrics
!EMK...These are currently disabled
  integer, public,  parameter :: LVT_ENSMETRIC_SINDEX = 65
  integer, public,  parameter :: LVT_ENSMETRIC_EINDEX = 65

  integer, public,  parameter :: LVT_NMETRICS        = 65

!-------------------------------------------------------------------------
! Run modes
!------------------------------------------------------------------------- 
   character*50, public,  parameter :: LVT_DataCompId = "Data intercomparison"
   character*50, public,  parameter :: LVT_dastatId = "DA statistics processing"
   character*50, public,  parameter :: LVT_benchMarkId = "Benchmarking"
   character*50, public,  parameter :: LVT_daobsId = "DA observation processing"
   character*50, public,  parameter :: LVT_optUEId  = "OPTUE output processing"
   character*50, public,  parameter :: LVT_rtmrunId  = "RTM output processing"
   character*50, public,  parameter :: LVT_557postId = "557 post"
   character*50, public,  parameter :: LVT_usafsipostId = "USAFSI post"
   character*50, public,  parameter :: LVT_LISpostId = "LIS postprocessing"

!-------------------------------------------------------------------------
! Domains
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LVT_latlonId    = "latlon"
   character*50, public,  parameter :: LVT_mercId      = "mercator"
   character*50, public,  parameter :: LVT_lambertId   = "lambert"
   character*50, public,  parameter :: LVT_gaussId     = "gaussian"
   character*50, public,  parameter :: LVT_polarId     = "polar"
   character*50, public,  parameter :: LVT_utmId       = "UTM"
!-------------------------------------------------------------------------
! Observations
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LVT_templateobsId      = "none" 
   character*50, public,  parameter :: LVT_LISoutputId        = "LIS output"
   character*50, public,  parameter :: LVT_LISdaobsId         = "LIS DAOBS"
   character*50, public,  parameter :: LVT_LIS6outputId       = "LIS6 output"
   character*50, public,  parameter :: LVT_LISDAdiagoutputId  = "LIS DA diagnostics"
   character*50, public,  parameter :: LVT_ceopobsId          = "CEOP"
   character*50, public,  parameter :: LVT_ISCCP_TskinobsId   = "ISCCP LST"
   character*50, public,  parameter :: LVT_MODIS_LSTobsId     = "MODIS LST"
   character*50, public,  parameter :: LVT_SCANGMAOobsId      = "SCAN(GMAO)"
   character*50, public,  parameter :: LVT_SCANobsId          = "SCAN"
   character*50, public,  parameter :: LVT_NASMDobsId         = "NASMD"
   character*50, public,  parameter :: LVT_ISMNobsId          = "ISMN"
   character*50, public,  parameter :: LVT_SURFRADobsId       = "SURFRAD"
   character*50, public,  parameter :: LVT_wgPBMRobsId        = "WG PBMRsm"
   character*50, public,  parameter :: LVT_SNOTELobsId        = "SNOTEL"
   character*50, public,  parameter :: LVT_LSWG_TbobsId       = "LSWG Tb"
   character*50, public,  parameter :: LVT_FMISWEobsId        = "FMI SWE"
   character*50, public,  parameter :: LVT_CMCSNWDobsId       = "CMC"
   character*50, public,  parameter :: LVT_SNODASobsId        = "SNODAS"
   character*50, public,  parameter :: LVT_NASAAMSREsmobsId   = "AMSR-E NASA soil moisture"
   character*50, public,  parameter :: LVT_LPRMAMSREsmobsId   = "AMSR-E LPRM soil moisture"
   character*50, public,  parameter :: LVT_AMMAobsId          = "AMMA"
   character*50, public,  parameter :: LVT_AmerifluxobsId     = "Ameriflux"
   character*50, public,  parameter :: LVT_ARMobsId           = "ARM"
   character*50, public,  parameter :: LVT_SMOSREXobsId       = "SMOSREX"
   character*50, public,  parameter :: LVT_AGRMETdataId       = "AGRMET"
   character*50, public,  parameter :: LVT_GlobSnowObsId      = "Globsnow"
   character*50, public,  parameter :: LVT_SNODEPmetobsId     = "SNODEP metobs"
   character*50, public,  parameter :: LVT_SNODEPobsId        = "SNODEP"
   character*50, public,  parameter :: LVT_MOD10A1obsId       = "MOD10A1"
   character*50, public,  parameter :: LVT_MODSCAGobsId       = "MODSCAG"
   character*50, public,  parameter :: LVT_MOD10A1V6obsId     = "MOD10A1V6"
   character*50, public,  parameter :: LVT_ANSASNWDobsId      = "ANSA snow depth"
   character*50, public,  parameter :: LVT_ANSASWEobsId       = "ANSA SWE"
   character*50, public,  parameter :: LVT_CPCPRCPobsId       = "CPC precipitation"
   character*50, public,  parameter :: LVT_USGSSFobsId        = "USGS streamflow"
   character*50, public,  parameter :: LVT_USGSSFgridobsId    = "USGS streamflow gridded"
   character*50, public,  parameter :: LVT_NatSFobsId         = "Naturalized streamflow"
   character*50, public,  parameter :: LVT_FLUXNETmteObsId    = "FLUXNET MTE"
   character*50, public,  parameter :: LVT_FLUXNET2015ObsId   = "FLUXNET2015 dataset"
   character*50, public,  parameter :: LVT_FLUXNET2015NCObsId = "FLUXNET2015 (NetCDF) dataset"
   character*50, public,  parameter :: LVT_MOD16A2obsId       = "MOD16A2"
   character*50, public,  parameter :: LVT_UWETobsId          = "UW ET"
   character*50, public,  parameter :: LVT_ARSsmobsId         = "USDA ARS soil moisture"  
   character*50, public,  parameter :: LVT_NLDAS2obsId        = "NLDAS2"
   character*50, public,  parameter :: LVT_GHCNobsId          = "GHCN"
   character*50, public,  parameter :: LVT_ALEXIobsId         = "ALEXI"
   character*50, public,  parameter :: LVT_ALEXIesiobsId      = "ALEXI ESI"
   character*50, public,  parameter :: LVT_GRACEobsId         = "GRACE"
   character*50, public,  parameter :: LVT_simGRACEobsId      = "simulated GRACE"
   character*50, public,  parameter :: LVT_USGSGWwellobsId    = "USGS ground water well data"
   character*50, public,  parameter :: LVT_PBOH2OobsId        = "PBO H2O"
   character*50, public,  parameter :: LVT_SMOSL2smobsId      = "SMOS L2 soil moisture"
   character*50, public,  parameter :: LVT_SMOSNESDISsmobsId  = "SMOS (NESDIS) soil moisture"
   character*50, public,  parameter :: LVT_SMOSCATDSsmobsId   = "SMOS (CATDS) soil moisture"
   character*50, public,  parameter :: LVT_SMOSL1TBobsId      = "SMOS L1 TB"
   character*50, public,  parameter :: LVT_GCOMW_AMSR2L3smobsId  = "GCOMW AMSR2 L3 soil moisture"
   character*50, public,  parameter :: LVT_GCOMW_AMSR2L3sndobsId  = "GCOMW AMSR2 L3 snow depth"
   character*50, public,  parameter :: LVT_SMOPSsmobsId  = "SMOPS soil moisture"
   character*50, public,  parameter :: LVT_ESACCIsmobsId = "ESA CCI soil moisture"
   character*50, public,  parameter :: LVT_GIMMSAVHRR_NDVIobsId = "GIMMS AVHRR NDVI"
   character*50, public,  parameter :: LVT_GIMMSMODIS_NDVIobsId = "GIMMS MODIS NDVI"
   character*50, public,  parameter :: LVT_GLDAS1obsId = "GLDAS1"
   character*50, public,  parameter :: LVT_GLDAS2obsId = "GLDAS2"

   character*50, public,  parameter :: LVT_MERRA2obsId    = "MERRA2"
   character*50, public,  parameter :: LVT_MERRA2asmObsId = "MERRA2asm"
   character*50, public,  parameter :: LVT_MERRAlandobsId = "MERRA-Land"
   character*50, public,  parameter :: LVT_ERAIlandobsId  = "ERA interim land"
   character*50, public,  parameter :: LVT_SSEBopobsId    = "SSEB"
   character*50, public,  parameter :: LVT_GRDCobsId      = "GRDC"
   character*50, public,  parameter :: LVT_GOESLSTobsId   = "GOES LST"
   character*50, public,  parameter :: LVT_GLERLobsId     = "GLERL hydro data"
   character*50, public,  parameter :: LVT_JULESobsId     = "JULES data"
   character*50, public,  parameter :: LVT_JULES2dobsId     = "JULES 2d data"
   character*50, public,  parameter :: LVT_SMAPsmobsId    = "SMAP soil moisture"
   character*50, public,  parameter :: LVT_SMAPvodobsId    = "SMAP vegetation optical depth" 
   character*50, public,  parameter :: LVT_LPRMvodobsId  = "LPRM vegetation optical depth"
   character*50, public,  parameter :: LVT_SMAPvwcobsId    = "SMAP vegetation water content" ! MN
   character*50, public,  parameter :: LVT_SMAP_L3TbId    = "SMAP L3 Tb" ! MN
   character*50, public,  parameter :: LVT_SMAPTBobsId    = "SMAP TB"
   character*50, public,  parameter :: LVT_GOME2SIFobsId  = "GOME2 SIF"
   character*50, public,  parameter :: LVT_DaymetobsId    = "Daymet"
   character*50, public,  parameter :: LVT_LVTbenchmarkobsId = "LVT benchmark"
   character*50, public,  parameter :: LVT_CMORPHdataId   = "CMORPH"
   character*50, public,  parameter :: LVT_TRMM3B42V7dataId   = "3B42V7"
   character*50, public,  parameter :: LVT_CHIRPSv2dataId = "CHIRPSv2"
   character*50, public,  parameter :: LVT_USCRNsmdataId  = "USCRN soil moisture"
   character*50, public,  parameter :: LVT_GLEAMdataId    = "GLEAM"
   character*50, public,  parameter :: LVT_USDMdataId     = "USDM"
   character*50, public,  parameter :: LVT_IMDprcpdataId     = "IMD gridded precipitation"
   character*50, public,  parameter :: LVT_APHROprcpdataId     = "APHRODITE precipitation"
   character*50, public,  parameter :: LVT_LVTpercentiledataId = "LVT percentile"
   character*50, public,  parameter :: LVT_GLASSlaiobsId = "GLASS LAI"
   character*50, public,  parameter :: LVT_GLASSalbedoobsId = "GLASS ALBEDO"
   character*50, public,  parameter :: LVT_MODISsportLAIobsId = "MODIS SPORT LAI"
   character*50, public,  parameter :: LVT_FLUXCOMobsId = "FLUXCOM"
   character*50, public,  parameter :: LVT_HARdataId = "HAR"
   character*50, public,  parameter :: LVT_OCO2SIFobsId  = "OCO2 SIF"
   character*50, public,  parameter :: LVT_ECMWFdataId = "ECMWF"
   character*50, public,  parameter :: LVT_GDASdataId = "GDAS"
   character*50, public,  parameter :: LVT_ASOSWEdataId = "ASO SWE"
   character*50, public,  parameter :: LVT_IMERGdataId = "GPM IMERG"
   character*50, public,  parameter :: LVT_IMERGMonthlydataId = &
        "GPM IMERG Monthly"
   character*50, public,  parameter :: LVT_UASNOWdataId = "UA SNOW"
   character*50, public,  parameter :: LVT_ozFluxdataId = "OzFlux"
   character*50, public,  parameter :: LVT_JASMINsmobsId = "JASMIN soil moisture"
   character*50, public,  parameter :: LVT_MCD15A2HobsId = "MCD15A2H LAI"
   character*50, public,  parameter :: LVT_ERA5obsId      = "ERA5"
   character*50, public,  parameter :: LVT_FluxSatobsId = "FluxSAT GPP"
   character*50, public,  parameter :: LVT_THySMobsId = "THySM"
   character*50, public,  parameter :: LVT_UASMAPobsId = "UA SMAP"
   character*50, public,  parameter :: LVT_GRUNobsId = "GRUN runoff"
   character*50, public,  parameter :: LVT_COAMPSoutId = "COAMPSout"
   character*50, public,  parameter :: LVT_SMAP_E_OPLId = "OPL E SMAP soil moisture retrieval"  
!-------------------------------------------------------------------------
! Training algorithms
!------------------------------------------------------------------------- 
   character*50, public,  parameter :: LVT_LinearRegressionId = "Linear regression"

!EOC
 end module LVT_pluginIndices
