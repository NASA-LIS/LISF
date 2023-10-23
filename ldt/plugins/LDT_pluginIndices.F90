!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_pluginIndices
!BOP

! !MODULE: LDT_pluginIndices
!
! !DESCRIPTION:
!  The code in this file provides values of indices
!  used to register functions in the plugin modules.
!
!  The index definitions are simply a convention.
!  The user may change these options, and the lis.config
!  should be changed appropriately to ensure that the
!  correct function is called at run time.
!
! !REVISION HISTORY:
!  23 Oct 2008: Sujay Kumar  -- Initial Specification
!  17 Jul 2012: KR Arsenault -- Updated entries with capitalization rules
!  01 Mar 2020: Yeosang Yoon -- Added MERIT DEM
!  28 Jun 2022: Eric Kemp -- Added NAFPA background precipitation
!
!EOP
  PRIVATE

!BOC
!-------------------------------------------------------------------------
! Run modes
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_LSMparamprocId   = "LSM parameter processing"
   character*50, public,  parameter :: LDT_DApreprocId      = "DA preprocessing"
   character*50, public,  parameter :: LDT_EnsRstpreprocId  = "Ensemble restart processing"
   character*50, public,  parameter :: LDT_climoRstProcId   = "Climatological restart processing"
   character*50, public,  parameter :: LDT_rstTransformProcId = "Restart transformation processing"
   character*50, public,  parameter :: LDT_NUWRFpreprocId   = "NUWRF preprocessing for real"
   character*50, public,  parameter :: LDT_ANNprocId        = "ANN processing"
   character*50, public,  parameter :: LDT_MetForcprocId    = "Metforce processing"
   character*50, public,  parameter :: LDT_MetTDscaleprocId = "Metforce temporal downscaling"
   character*50, public,  parameter :: LDT_StatDscaleMetforcprocId = "Statistical downscaling of met forcing"
   character*50, public,  parameter :: LDT_usafsiId = "USAFSI analysis"
   character*50, public,  parameter :: LDT_OPTUEparamprocId   = "OPTUE parameter processing"
   character*50, public,  parameter :: LDT_obsSimprocId   = "Observation simulator"
   character*50, public,  parameter :: LDT_LISHydropreprocId  = "LISHydro preprocessing for WRFHydro"
   character*50, public,  parameter :: LDT_SMAP_E_OPLId       = "OPL E SMAP soil moisture retrieval"  !Y.Kwon

!-------------------------------------------------------------------------
! Domains
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_latlonId    = "latlon"     ! 0
   character*50, public,  parameter :: LDT_mercId      = "mercator"   ! 1
   character*50, public,  parameter :: LDT_lambertId   = "lambert"    ! 3
   character*50, public,  parameter :: LDT_gaussId     = "gaussian"   ! 4
   character*50, public,  parameter :: LDT_polarId     = "polar"      ! 5
   character*50, public,  parameter :: LDT_easev2Id    = "ease V2"    !
   character*50, public,  parameter :: LDT_utmId       = "UTM"        ! ?
   character*50, public,  parameter :: LDT_hrapId      = "hrap"       ! 10
   character*50, public,  parameter :: LDT_catdomainId = "catchment"

!-------------------------------------------------------------------------
! DA Observations
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_LISlsmSMobsId              &
        = "LIS LSM soil moisture"
   character*50, public,  parameter :: LDT_LISlsmTEFFobsId            &
        = "LIS LSM effective soil temperature"                               !Y.Kwon
   character*50, public,  parameter :: LDT_syntheticSMobsId           &
        = "Synthetic soil moisture"
   character*50, public,  parameter :: LDT_NASA_AMSREsmobsId          &
        = "AMSR-E(NASA) soil moisture"
   character*50, public,  parameter :: LDT_LPRM_AMSREsmobsId          &
        = "AMSR-E(LPRM) soil moisture"
   character*50, public,  parameter :: LDT_ESACCIsmobsId              &
        = "ESA CCI soil moisture"
   character*50, public,  parameter :: LDT_WindSatsmobsId             &
        = "WindSat soil moisture"
   character*50, public,  parameter :: LDT_SMOPSsmobsId               &
        = "SMOPS soil moisture"
   character*50, public,  parameter :: LDT_ASCATTUWsmobsId            &
        = "ASCAT TUW soil moisture"
   character*50, public,  parameter :: LDT_GRACEtwsobsId              &
        = "GRACE TWS"
   character*50, public,  parameter :: LDT_GRACEQLtwsobsId            &
        = "GRACE QL TWS"
   character*50, public,  parameter :: LDT_SMOSL2smobsId              &
        = "SMOS L2 soil moisture"
   character*50, public,  parameter :: LDT_SMOSNESDISsmobsId          &
        = "SMOS NESDIS soil moisture"
   character*50, public,  parameter :: LDT_GCOMW_AMSR2L3smobsId       &
        = "GCOMW AMSR2 L3 soil moisture"
   character*50, public,  parameter :: LDT_AquariusL2smobsId          &
        = "Aquarius L2 soil moisture"
   character*50, public,  parameter :: LDT_simGRACEJPLobsId           &
        = "Simulated GRACE (JPL)"
   character*50, public,  parameter :: LDT_SMMRSNWDsnowobsId          &
        = "SMMR snow depth"
   character*50, public,  parameter :: LDT_SSMISNWDsnowobsId          &
        = "SSMI snow depth"
   character*50, public,  parameter :: LDT_ANSASNWDsnowobsId          &
        = "ANSA snow depth"
   character*50, public,  parameter :: LDT_GCOMWAMSR2L3sndobsId       &
        = "GCOMW AMSR2 L3 snow depth"
   character*50, public,  parameter :: LDT_NASASMAPsmobsId            &
        = "NASA SMAP soil moisture"
   character*50, public,  parameter :: LDT_SMAPEOPLsmobsId            &
        = "SMAP_E_OPL soil moisture"                                        !Y.Kwon
   character*50, public,  parameter :: LDT_THySMobsId            &
        = "THySM soil moisture"
   character*50, public,  parameter :: LDT_SMOSNRTNNsmobsId            &
        = "SMOS NRT NN soil moisture"                                        !Y.Kwon
   character*50, public,  parameter :: LDT_NASASMAPvodobsId            &
        = "NASA SMAP vegetation optical depth"
   character*50, public,  parameter :: LDT_GLASSlaiobsId            &
        = "GLASS LAI"
   character*50, public,  parameter :: LDT_LPRMvodobsId            &
        = "LPRM vegetation optical depth"
   character*50, public,  parameter :: LDT_MCD15A2HlaiobsId            &
        = "MCD15A2H LAI"
   character*50, public,  parameter :: LDT_LISlsmPrecipobsId          &
        = "LIS LSM total precipitation"
   character*50, public,  parameter :: LDT_VIIRSgvfobsId            &
        = "VIIRS GVF"                                                    !Y.Kwon
   character*50, public,  parameter :: LDT_CDFSgvfobsId            &
        = "CDFS GVF"                                                     !Y.Kwon
   character*50, public,  parameter :: LDT_GEOSTeffobsId            &
        = "GEOS effective soil temperature"                              !Y.Kwon
!-------------------------------------------------------------------------
! Meteorological forcings
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_metForcTemplateId = "none"
   character*50, public,  parameter :: LDT_agrmetId       = "AGRMET"
   character*50, public,  parameter :: LDT_agrradpsId     = "AGRMET radiation (polar stereographic)"
   character*50, public,  parameter :: LDT_agrradId       = "AGRMET radiation (latlon)"
   character*50, public,  parameter :: LDT_gdasId         = "GDAS"
   character*50, public,  parameter :: LDT_geos5fcstId    = "GEOS5 forecast"
   character*50, public,  parameter :: LDT_ecmwfId        = "ECMWF"
   character*50, public,  parameter :: LDT_princetonId    = "PRINCETON"
   character*50, public,  parameter :: LDT_merra2Id       = "MERRA2"
   character*50, public,  parameter :: LDT_era5Id         = "ERA5"
   character*50, public,  parameter :: LDT_gswp1Id        = "GSWP1"
   character*50, public,  parameter :: LDT_gswp2Id        = "GSWP2"
   character*50, public,  parameter :: LDT_nldas2Id       = "NLDAS2"
   character*50, public,  parameter :: LDT_gldasId        = "GLDAS"
   character*50, public,  parameter :: LDT_gdasLSWGId     = "GDAS(LSWG)"
   character*50, public,  parameter :: LDT_gfsId          = "GFS"
   character*50, public,  parameter :: LDT_narrId         = "NARR"
   character*50, public,  parameter :: LDT_nam242Id       = "NAM242"
   character*50, public,  parameter :: LDT_wrfoutv2Id     = "WRFoutv2"
   character*50, public,  parameter :: LDT_WRFakId        = "WRF AK"
   character*50, public,  parameter :: LDT_cmapId         = "CMAP"
!   character*50, public,  parameter :: LDT_TRMM3B42RTId   = "TRMM 3B42RT"
   character*50, public,  parameter :: LDT_TRMM3B42V6Id   = "TRMM 3B42V6"
   character*50, public,  parameter :: LDT_TRMM3B42V7Id   = "TRMM 3B42V7"
   character*50, public,  parameter :: LDT_TRMM3B42RTV7Id = "TRMM 3B42RTV7"
   character*50, public,  parameter :: LDT_cmorphId       = "CPC CMORPH"
   character*50, public,  parameter :: LDT_stg2Id         = "CPC STAGEII"
   character*50, public,  parameter :: LDT_stg4Id         = "CPC STAGEIV"
   character*50, public,  parameter :: LDT_RFE2DailyId    = "RFE2(daily)"
   character*50, public,  parameter :: LDT_RFE2gdasId     = "RFE2(GDAS bias-corrected)"
   character*50, public,  parameter :: LDT_chirps2Id      = "CHIRPS2"
   character*50, public,  parameter :: LDT_ceopId         = "CEOP"
   character*50, public,  parameter :: LDT_ALMIPIIId      = "ALMIPII"
   character*50, public,  parameter :: LDT_scanId         = "SCAN"
   character*50, public,  parameter :: LDT_armsId         = "ARMS"
   character*50, public,  parameter :: LDT_d2pcpcarId     = "D2PCPCAR"
   character*50, public,  parameter :: LDT_d2pcpoklId     = "D2PCPOKL"
   character*50, public,  parameter :: LDT_Noah31BondId   = "Noah Bondville"
   character*50, public,  parameter :: LDT_FASSTsingleId  = "FASST test"
   character*50, public,  parameter :: LDT_TRIGRSseattleId ="TRIGRS test"
   character*50, public,  parameter :: LDT_snotelId       = "SNOTEL"
   character*50, public,  parameter :: LDT_coopId         = "COOP"
   character*50, public,  parameter :: LDT_rhoneAGGId     = "Rhone AGG"
   character*50, public,  parameter :: LDT_vicforcingId   = "VIC processed forcing"

!-------------------------------------------------------------------------
! LSMS
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_templateLSMId = "none"
   character*50, public,  parameter :: LDT_noah271Id     = "Noah.2.7.1"
   character*50, public,  parameter :: LDT_noah32Id      = "Noah.3.2"
   character*50, public,  parameter :: LDT_noah33Id      = "Noah.3.3"
   character*50, public,  parameter :: LDT_noah36Id      = "Noah.3.6"
   character*50, public,  parameter :: LDT_noah39Id      = "Noah.3.9"
   character*50, public,  parameter :: LDT_noahmp36Id    = "Noah-MP.3.6"
   character*50, public,  parameter :: LDT_noahmp401Id   = "Noah-MP.4.0.1"
   character*50, public,  parameter :: LDT_clm2Id        = "CLM.2"
   character*50, public,  parameter :: LDT_clm45Id       = "CLM.4.5"
   character*50, public,  parameter :: LDT_vic411Id      = "VIC.4.1.1"
   character*50, public,  parameter :: LDT_vic412Id      = "VIC.4.1.2"
   character*50, public,  parameter :: LDT_mosaicId      = "Mosaic"
   character*50, public,  parameter :: LDT_hyssibId      = "HySSIB"
   character*50, public,  parameter :: LDT_sib2Id        = "SiB2"
   character*50, public,  parameter :: LDT_tessId        = "HTESSEL"
   character*50, public,  parameter :: LDT_jules50Id     = "JULES.5.0"
   character*50, public,  parameter :: LDT_cableId       = "CABLE"
   character*50, public,  parameter :: LDT_fasstId       = "FASST"
   character*50, public,  parameter :: LDT_sheelsId      = "SHEELS"
   character*50, public,  parameter :: LDT_clsmf25Id     = "CLSMF2.5"
   character*50, public,  parameter :: LDT_GeoWRSI2Id    = "GeoWRSI.2"
   character*50, public,  parameter :: LDT_rdhm356Id     = "RDHM.3.5.6"
   character*50, public,  parameter :: LDT_sachtet356Id  = "SACHTET.3.5.6"
   character*50, public,  parameter :: LDT_snow17Id      = "SNOW17"
   character*50, public,  parameter :: LDT_ruc37Id       = "RUC.3.7"

!-------------------------------------------------------------------------
! Lake models and data
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_flakeId    = "FLake"
   character*50, public,  parameter :: LDT_gldbv1Id   = "GLDBv1"
   character*50, public,  parameter :: LDT_gldbv2Id   = "GLDBv2"
   character*50, public,  parameter :: LDT_glwdId     = "GLWD"

!-------------------------------------------------------------------------
! Snow models and data
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_Crocus81Id    = "Crocus8.1"  ! this is SURFEX version comes from surf_version.F90
   character*50, public,  parameter :: LDT_snowmodelId   = "SnowModel"

!-------------------------------------------------------------------------
! Landcover sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_avhrrlcLISId   = "AVHRR"
   character*50, public,  parameter :: LDT_avhrrlcGFSId   = "AVHRR_GFS"
   character*50, public,  parameter :: LDT_usgslcLISId    = "USGS_LIS"
   character*50, public,  parameter :: LDT_usgslcNATId    = "USGS_Native"
   character*50, public,  parameter :: LDT_modislcLISId   = "MODIS_LIS"
   character*50, public,  parameter :: LDT_modislcNATId   = "MODIS_Native"
   character*50, public,  parameter :: LDT_mcd12q1Id      = "MCD12Q1"
   character*50, public,  parameter :: LDT_modislcPFTId   = "MODIS_Native_PFT"
   character*50, public,  parameter :: LDT_ukmoigbpPFTId  = "UKMO_IGBP_Native_PFT"
   character*50, public,  parameter :: LDT_UM_ancillaryId = "UM_Native_Ancillary"
   character*50, public,  parameter :: LDT_ALMIPIIlcId    = "ALMIPII"
   character*50, public,  parameter :: LDT_isalcId        = "ISA"
   character*50, public,  parameter :: LDT_clsmf25lcId    = "CLSMF2.5"
   character*50, public,  parameter :: LDT_vic411lcId     = "VIC411"
   character*50, public,  parameter :: LDT_vic412lcId     = "VIC412"
   character*50, public,  parameter :: LDT_clm45lcId      = "CLM45"
   character*50, public,  parameter :: LDT_nalcmsSMlcId   = "NALCMS_SM"
   character*50, public,  parameter :: LDT_nalcmsSMIGBPlcId  = "NALCMS_SM_IGBPNCEP"
   character*50, public,  parameter :: LDT_constId        = "CONSTANT"

!-------------------------------------------------------------------------
! Croptype sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_umdcropmapId   = "UMDCROPMAP"
   character*50, public,  parameter :: LDT_monfredacropId = "Monfreda08"

!------------------------------------------------------------------------
! Mask sources
!-------------------------------------------------------------------------
! Global/Continent -
   character*50, public,  parameter :: LDT_modis44WmaskId = "MOD44W"
   character*50, public,  parameter :: LDT_vic411maskId   = "VIC411"
! Regional -
   character*50, public,  parameter :: LDT_gismaskId      = "ESRI"
   character*50, public,  parameter :: LDT_wrsimaskId     = "WRSI"
   character*50, public,  parameter :: LDT_maskmaskId     = "file"

!-------------------------------------------------------------------------
! Topography sources
!-------------------------------------------------------------------------
   character*50, public, parameter :: LDT_gtopoLISId = "GTOPO30_LIS"
   character*50, public, parameter :: LDT_gtopoGFSId = "GTOPO30_GFS"
   character*50, public, parameter :: LDT_gtopoNATId = "GTOPO30_Native"
   character*50, public, parameter :: LDT_srtmLISId  = "SRTM_LIS"
   character*50, public, parameter :: LDT_srtmNATId  = "SRTM_Native"
   character*50, public, parameter :: LDT_merit1KId  = "MERIT_1K"
   character*50, public, parameter :: LDT_nedSMId    = "NED_SM"

!-------------------------------------------------------------------------
! Soils sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_faoSoilId         = "FAO"
   character*50, public,  parameter :: LDT_statsgo1FAOLISId  = "STATSGOFAO_LIS"  ! STATSGO+FAO
   character*50, public,  parameter :: LDT_zoblerGFSId       = "ZOBLER_GFS"
   character*50, public,  parameter :: LDT_statsgo1FAONATId  = "STATSGOFAO_Native"
   character*50, public,  parameter :: LDT_statsgov1LISId    = "STATSGO_LIS"     ! STATSGO v1 only
   character*50, public,  parameter :: LDT_statsgov1NATId    = "STATSGOv1"       ! STATSGO v1 only
   character*50, public,  parameter :: LDT_statsgov2Id       = "STATSGOv2"       ! STATSGO v2 only
   character*50, public,  parameter :: LDT_ALMIPIIsoilId     = "ALMIPII"
   character*50, public,  parameter :: LDT_specialsoilId     = "Special"
   character*50, public,  parameter :: LDT_ISRICsoilId       = "ISRIC"
   character*50, public,  parameter :: LDT_ukmofracId        = "UKMOFRAC"        ! UKMO

!-------------------------------------------------------------------------
! LAI/SAI sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_avhrrlaiId    = "AVHRR"
   character*50, public,  parameter :: LDT_modislaiId    = "MODIS"
   character*50, public,  parameter :: LDT_clsmf25laiId  = "CLSMF2.5"
!-------------------------------------------------------------------------
! Greenness data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_gfracClimLISId = "NCEP_LIS"
   character*50, public,  parameter :: LDT_gfracClimNATId = "NCEP_Native"
   character*50, public,  parameter :: LDT_gfracClsmf25Id = "CLSMF2.5"
   character*50, public,  parameter :: LDT_gfracSACHTETId = "SACHTET.3.5.6"
!-------------------------------------------------------------------------
! Maximum greenness data sources
!-------------------------------------------------------------------------
  character*50, public,  parameter :: LDT_ncepshdmaxLISId = "NCEP_LIS"
!-------------------------------------------------------------------------
! Minimum greenness data sources
!-------------------------------------------------------------------------
  character*50, public,  parameter :: LDT_ncepshdminLISId = "NCEP_LIS"
!-------------------------------------------------------------------------
! PET data sources
!-------------------------------------------------------------------------
  character*50, public,  parameter :: LDT_petSACHTETId  = "SACHTET.3.5.6"

!-------------------------------------------------------------------------
! Irrigation sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_modisirrigId  = "MODIS"
   character*50, public,  parameter :: LDT_modOGirrigId  = "MODIS_OG"
   character*50, public,  parameter :: LDT_gripcirrigId  = "GRIPC"
   character*50, public,  parameter :: LDT_irriggwratioId  = "USGS_Native"

   character*50, public,  parameter :: LDT_userinputirrigId = "UserDerived"

!-------------------------------------------------------------------------
! Albedo data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_albedoClimLISId    = "NCEP_LIS"
   character*50, public,  parameter :: LDT_albedoClimNATId    = "NCEP_Native"
   character*50, public,  parameter :: LDT_albedoClimNATQtrId = "NCEP_NativeQtr"
   character*50, public,  parameter :: LDT_mxsnalbMODISId     = "MODIS"
   character*50, public,  parameter :: LDT_mxsnalbGFSId       = "NCEP_GFS"
   character*50, public,  parameter :: LDT_mxsnalbBarlageId   = "Barlage_Native"
   character*50, public,  parameter :: LDT_mxsnalbSACHTETId   = "SACHTET.3.5.6"

   character*50, public,  parameter :: LDT_albnirclsmf25Id = "CLSMF2.5"
   character*50, public,  parameter :: LDT_albvisclsmf25Id = "CLSMF2.5"

!-------------------------------------------------------------------------
! Bottom temperature (Tbot) data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_nceptbotId    = "NCEP_LIS"
   character*50, public,  parameter :: LDT_ncepgfstbotId = "NCEP_GFS"
   character*50, public,  parameter :: LDT_islscp1tbotId = "ISLSCP1"
   character*50, public,  parameter :: LDT_tbotSACHTETId = "SACHTET.3.5.6"

!-------------------------------------------------------------------------
! Slope type data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_slopetypeLISId   = "NCEP_LIS"
   character*50, public,  parameter :: LDT_slopetypeGFSId   = "NCEP_GFS"
   character*50, public,  parameter :: LDT_slopetypeNATId   = "NCEP_Native"
   character*50, public,  parameter :: LDT_slopetypeConstId = "CONSTANT"

!-------------------------------------------------------------------------
! Climate downscaling data sources
!-------------------------------------------------------------------------
!- PPT:
   character*50, public,  parameter :: LDT_prismpptId     = "PRISM"
   character*50, public,  parameter :: LDT_worldclimpptId = "WORLDCLIM"
   character*50, public,  parameter :: LDT_nafpabackgfspptId = "NAFPA_BACK_GFS"
   character*50, public,  parameter :: LDT_nafpabackgalwempptId = &
        "NAFPA_BACK_GALWEM"

!- TMIN:
   character*50, public,  parameter :: LDT_prismtminId     = "PRISM"
   character*50, public,  parameter :: LDT_worldclimtminId = "WORLDCLIM"
!- TMAX:
   character*50, public,  parameter :: LDT_prismtmaxId     = "PRISM"
   character*50, public,  parameter :: LDT_worldclimtmaxId = "WORLDCLIM"

!-------------------------------------------------------------------------
! Slope type data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_GLIMSId = "GLIMS"

!-------------------------------------------------------------------------
! Routing model data sources
!-------------------------------------------------------------------------
   character*50, public, parameter  :: LDT_HYMAPId  = "HYMAP"
   character*50, public, parameter  :: LDT_HYMAP2Id = "HYMAP2"

!-------------------------------------------------------------------------
! ANN data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_ANNLISlsmSMId     &
        = "LIS LSM soil moisture"
   character*50, public,  parameter :: LDT_ANNsynSMId     &
        = "Synthetic soil moisture"
   character*50, public,  parameter :: LDT_ANNLPRMAMSREsmobsId &
        = "AMSR-E(LPRM) soil moisture"
   character*50, public,  parameter :: LDT_ANNMODISlstobsId &
        = "MODIS LST"
   character*50, public,  parameter :: LDT_ANNGHCNsnwdobsId &
        = "GHCN snow depth"
   character*50, public,  parameter :: LDT_ANNMOD10A1obsId &
        = "MOD10A1 snow cover"
   character*50, public,  parameter :: LDT_ANNGCOMWAMSR2TbobsId &
        = "GCOMW AMSR2 Tb"

!-------------------------------------------------------------------------
! Met forcing scaling (up/down) options
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_simplewgtId = "Simple weighting"

!-------------------------------------------------------------------------
!  Statistical downscaling options
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_forcingClimoId  = "Climatology"
   character*50, public,  parameter :: LDT_bayesianMergeId = "Bayesian merging"

!-------------------------------------------------------------------------
!  obs simulator nature run source
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_LISoutNatureRunDataId = "LIS output"

!-------------------------------------------------------------------------
!  obs simulator OSSE mask
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LDT_LISoutOSSEmaskDataId = "LIS output"
   character*50, public,  parameter :: LDT_AMSR2OSSEmaskDataId = "AMSR2"
   character*50, public,  parameter :: LDT_MODISOSSEmaskDataId = "MODIS"
   character*50, public,  parameter :: LDT_Sentinel1AOSSEmaskDataId = "Sentinel1A"
   character*50, public,  parameter :: LDT_TSMMOSSEmaskDataId = "TSMM"
   

!EOC
 end module LDT_pluginIndices
