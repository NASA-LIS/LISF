!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_pluginIndices
!BOP
!
!  !MODULE: LIS_pluginIndices
! 
!  !DESCRIPTION: 
!   The code in this file provides values of indices used to 
!   to register functions in the plugin modules
!  
!  The index definitions are simply a convention
!  The user may change these options, and the lis.config 
!  should be changed appropriately to ensure that the correct function
!  is called at run time
! 
!  !REVISION HISTORY: 
!  23 Oct 2006: Sujay Kumar  Initial Specification
!  17 Jan 2011: David Mocko, added max/min greenness & slope type
!  08 Feb 2011: Yudong Tian, added cmem3 rtm. 
!  01 Jun 2012: Sujay Kumar, changed Ids to character strings
!  22 May 2013: Shugong Wang, added FLAKE.1.0  
!  27 Jan 2014: Shugong Wang, added HRAP projection
!   4 Nov 2014: Jonathan Case, added support for daily NESDIS/VIIRS GVF for Noah
!  16 Aug 2016: Mahdi Navari, added PILDAS  
!
!EOP
  PRIVATE
   
!BOC
!-------------------------------------------------------------------------
! Run modes
!------------------------------------------------------------------------- 
   character*50, public,  parameter :: LIS_retroId     = "retrospective"
   character*50, public,  parameter :: LIS_RTMforwardId  = "RTM forward"
   character*50, public,  parameter :: LIS_agrmetrunId = "AGRMET ops"
   character*50, public,  parameter :: LIS_wrfcplId    = "WRF coupling"
   character*50, public,  parameter :: LIS_gcecplId    = "GCE coupling"
   character*50, public,  parameter :: LIS_gfscplId    = "GFS coupling"
   character*50, public,  parameter :: LIS_paramEstimRunId = "parameter estimation"
   character*50, public,  parameter :: LIS_smootherDAId    = "ensemble smoother"
   character*50, public,  parameter :: LIS_forecastrunId      = "forecast"
!-------------------------------------------------------------------------
! Domains
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_latlonId    = "latlon"
   character*50, public,  parameter :: LIS_mercId      = "mercator"
   character*50, public,  parameter :: LIS_lambertId   = "lambert"
   character*50, public,  parameter :: LIS_gaussId     = "gaussian"
   character*50, public,  parameter :: LIS_polarId     = "polar"
   character*50, public,  parameter :: LIS_utmId       = "UTM"
   character*50, public,  parameter :: LIS_catdomainId = "catchment"
   character*50, public,  parameter :: LIS_hrapId      = "hrap"
!-------------------------------------------------------------------------
! LSMS
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_templateLSMId = "none"
   character*50, public,  parameter :: LIS_noah271Id   = "Noah.2.7.1"
   character*50, public,  parameter :: LIS_noah32Id    = "Noah.3.2"
   character*50, public,  parameter :: LIS_noah33Id    = "Noah.3.3"
   character*50, public,  parameter :: LIS_noah36Id    = "Noah.3.6"
   character*50, public,  parameter :: LIS_noah39Id    = "Noah.3.9"
   character*50, public,  parameter :: LIS_noahmp36Id  = "NoahMP.3.6"
   character*50, public,  parameter :: LIS_noahmp401Id = "Noah-MP.4.0.1"
   character*50, public,  parameter :: LIS_ruc37Id     = "RUC.3.7"
   character*50, public,  parameter :: LIS_clm2Id      = "CLM.2"
   character*50, public,  parameter :: LIS_vic411Id    = "VIC.4.1.1"
   character*50, public,  parameter :: LIS_vic412Id    = "VIC.4.1.2"
   character*50, public,  parameter :: LIS_mosaicId    = "Mosaic"
   character*50, public,  parameter :: LIS_hyssibId    = "HySSIB"
   !character*50, public,  parameter :: LIS_sib2Id      = "SiB2"
   !character*50, public,  parameter :: LIS_tessId      = "HTESSEL"
   character*50, public,  parameter :: LIS_jules43Id     = "JULES.4.3"
   character*50, public,  parameter :: LIS_jules50Id     = "JULES.5.0"
   character*50, public,  parameter :: LIS_jules52Id     = "JULES.5.2"
   character*50, public,  parameter :: LIS_jules53Id     = "JULES.5.3"
   character*50, public,  parameter :: LIS_cableId     = "CABLE"
   character*50, public,  parameter :: LIS_fasstId     = "FASST"
   !character*50, public,  parameter :: LIS_sheelsId    = "SHEELS"
   character*50, public,  parameter :: LIS_clsmf25Id   = "CLSM F2.5"
   character*50, public,  parameter :: LIS_geowrsi2Id  = "GeoWRSI.2"
   character*50, public,  parameter :: LIS_rdhm356lsmId = "RDHM.3.5.6"
   character*50, public,  parameter :: LIS_summa1Id     = "SUMMA.1.0"

!-------------------------------------------------------------------------
! Lake models
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_FLAKE1Id   = "FLAKE.1.0"

!-------------------------------------------------------------------------
! Glacier models
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_templateGLId = "template glacier"
   character*50, public,  parameter :: LIS_noahmpglacier3911Id  = "NoahMP-GL.3.9.1.1"

!-------------------------------------------------------------------------
! Open water models
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_templateOpenWaterId    = "template open water"

!-------------------------------------------------------------------------
! Met forcings 
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_metForcTemplateId = "none"
   character*50, public,  parameter :: LIS_generatedForcId   = "LDT-generated"
   character*50, public,  parameter :: LIS_climstdId         = "CLIM-Standard"
   character*50, public,  parameter :: LIS_genensfcstId      = "GenEnsFcst"
   character*50, public,  parameter :: LIS_pptensfcstId      = "PPTEnsFcst"

   character*50, public,  parameter :: LIS_gdasId            = "GDAS"
   character*50, public,  parameter :: LIS_gdasT1534Id       = "GDAS T1534"
   character*50, public,  parameter :: LIS_geosId            = "GEOS"
   character*50, public,  parameter :: LIS_geos5fcstId       = "GEOS5 forecast"
   character*50, public,  parameter :: LIS_ecmwfId           = "ECMWF"
   character*50, public,  parameter :: LIS_gswp1Id           = "GSWP1"
   character*50, public,  parameter :: LIS_gswp2Id           = "GSWP2"
   character*50, public,  parameter :: LIS_ecmwfreanalId     = "ECMWF reanalysis"
   character*50, public,  parameter :: LIS_agrmetId          = "AGRMET"
   character*50, public,  parameter :: LIS_princetonId       = "PRINCETON"
   character*50, public,  parameter :: LIS_nldas1Id          = "NLDAS1"
   character*50, public,  parameter :: LIS_nldas2Id          = "NLDAS2"

   character*50, public,  parameter :: LIS_gldasId           = "GLDAS"
   character*50, public,  parameter :: LIS_gfsId             = "GFS"
   character*50, public,  parameter :: LIS_merralandId       = "MERRA-Land"
   character*50, public,  parameter :: LIS_merra2Id          = "MERRA2"

   character*50, public,  parameter :: LIS_cmapId            = "CMAP"
   character*50, public,  parameter :: LIS_chirps2Id         = "CHIRPS2"

   character*50, public,  parameter :: LIS_TRMM3B42RTId      = "TRMM 3B42RT"
   character*50, public,  parameter :: LIS_TRMM3B42RTV7Id    = "TRMM 3B42RTV7"
   character*50, public,  parameter :: LIS_TRMM3B42V6Id      = "TRMM 3B42V6"
   character*50, public,  parameter :: LIS_TRMM3B42V7Id      = "TRMM 3B42V7" ! SY
   character*50, public,  parameter :: LIS_cmorphId          = "CPC CMORPH"
   character*50, public,  parameter :: LIS_imergId           = "GPM IMERG"
   character*50, public,  parameter :: LIS_stg2Id            = "CPC STAGEII"
   character*50, public,  parameter :: LIS_stg4Id            = "CPC STAGEIV"

   character*50, public,  parameter :: LIS_narrId            = "NARR"
   character*50, public,  parameter :: LIS_ALMIPIIId         = "ALMIPII"
   character*50, public,  parameter :: LIS_RFE2DailyId       = "RFE2(daily)"
   character*50, public,  parameter :: LIS_ceopId            = "CEOP"
   character*50, public,  parameter :: LIS_scanId            = "SCAN"
   character*50, public,  parameter :: LIS_armsId            = "ARMS"
   character*50, public,  parameter :: LIS_gdasLSWGId        = "GDAS(LSWG)"
   !character*50, public,  parameter :: LIS_d2pcpcarId        = "D2PCPCAR"
   !character*50, public,  parameter :: LIS_d2pcpoklId        = "D2PCPOKL"
   character*50, public,  parameter :: LIS_rdhm356Id         = "RDHM.3.5.6"
   character*50, public,  parameter :: LIS_gdas3dId          = "GDAS(3d)"
   character*50, public,  parameter :: LIS_agrradpsId        = "AGRMET radiation (polar stereographic)"
   character*50, public,  parameter :: LIS_agrradId          = "AGRMET radiation (latlon)"
   character*50, public,  parameter :: LIS_BondvilleId       = "Bondville"
   character*50, public,  parameter :: LIS_LoobosId          = "Loobos"
   character*50, public,  parameter :: LIS_FASSTsingleId     = "FASST test"
   character*50, public,  parameter :: LIS_TRIGRSseattleId   = "TRIGRS test"
   character*50, public,  parameter :: LIS_snotelId          = "SNOTEL"
   character*50, public,  parameter :: LIS_coopId            = "COOP"
   character*50, public,  parameter :: LIS_rhoneAGGId        = "Rhone AGG"
   character*50, public,  parameter :: LIS_RFE2gdasId        = "RFE2(GDAS bias-corrected)"
   character*50, public,  parameter :: LIS_vicforcingId      = "VIC processed forcing"
   character*50, public,  parameter :: LIS_PALSmetforcId     = "PALS station forcing"
   character*50, public,  parameter :: LIS_PILDASmetforcId   = "PILDAS"
   character*50, public,  parameter :: LIS_USGSPETforcId     = "PET USGS"
   character*50, public,  parameter :: LIS_capaId            = "CaPA"
   character*50, public,  parameter :: LIS_nam242Id          = "NAM242"
   character*50, public,  parameter :: LIS_WRFoutId          = "WRFout"
   character*50, public,  parameter :: LIS_AWAPforcId        = "AWAP"
   character*50, public,  parameter :: LIS_HiMATGMUforcId    = "HiMAT GMU"
!-------------------------------------------------------------------------
! land surface parameters
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! LAI/SAI sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_modis_RT_laiId  = "MODIS real-time"
   character*50, public,  parameter :: LIS_ALMIPIIlaiId    = "ALMIPII"
!-------------------------------------------------------------------------
! greenness data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_NESDISgfracId ="NESDIS weekly"
   character*50, public,  parameter :: LIS_SPORTgfracId = "SPORT"
   character*50, public,  parameter :: LIS_VIIRSgfracId = "VIIRS"
   character*50, public,  parameter :: LIS_ALMIPIIgfracId = "ALMIPII"
!-------------------------------------------------------------------------
! roughness data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_ALMIPIIroughnessId = "ALMIPII"
!-------------------------------------------------------------------------
! albedo data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_ALMIPIIalbedoId = "ALMIPII"
!-------------------------------------------------------------------------
! emissivity data sources
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_ALMIPIIemissId = "ALMIPII"
!-------------------------------------------------------------------------
! data assimilation algorithms
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_diId         = "Direct insertion"
   character*50, public,  parameter :: LIS_ekfId        = "EKF"
   character*50, public,  parameter :: LIS_enkfId       = "EnKF"
   character*50, public,  parameter :: LIS_ensrfId      = "EnSRF"
   character*50, public,  parameter :: LIS_enksId       = "EnKS"
   character*50, public,  parameter :: LIS_pfId         = "PF"
!-------------------------------------------------------------------------
! perturbation algorithms
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_uniformpertId = "uniform"
   character*50, public,  parameter :: LIS_gmaopertId = "GMAO scheme"
!-------------------------------------------------------------------------
! Assimilation set
! DA variable being updated with observations
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_noobsId             = "none"
   character*50, public,  parameter :: LIS_synsmId             = "Synthetic SM"
   character*50, public,  parameter :: LIS_synsweId            = "Synthetic SWE"
   character*50, public,  parameter :: LIS_synlstId            = "Synthetic LST"
   character*50, public,  parameter :: LIS_synsndId            = "Synthetic SND"
   character*50, public,  parameter :: LIS_synSnowTbID         = "Synthetic Snow TB"              

   character*50, public,  parameter :: LIS_multisynsmobsId     = "Synthetic(Multilayer) sm"
   character*50, public,  parameter :: LIS_isccpTskinId        = "ISCCP LST"
   character*50, public,  parameter :: LIS_NASA_AMSREsmobsId   = "AMSR-E(NASA) soil moisture"
   character*50, public,  parameter :: LIS_LPRM_AMSREsmobsId   = "AMSR-E(LPRM) soil moisture"
   character*50, public,  parameter :: LIS_ESACCIsmobsId       = "ESA CCI soil moisture"
   character*50, public,  parameter :: LIS_WindSatsmobsId      = "Windsat"
   character*50, public,  parameter :: LIS_WindSatCsmobsId     = "Windsat C-band"
   character*50, public,  parameter :: LIS_snodepobsId         = "SNODEP"
   character*50, public,  parameter :: LIS_ldtsiobsId          = "LDTSI"
   character*50, public,  parameter :: LIS_ANSASWEsnowobsId    = "ANSA SWE"
   character*50, public,  parameter :: LIS_ANSASCFsnowobsId    = "ANSA SCF"
   character*50, public,  parameter :: LIS_ANSASNWDsnowobsId   = "ANSA snow depth"
   character*50, public,  parameter :: LIS_SMMRSNWDsnowobsId   = "SMMR snow depth"
   character*50, public,  parameter :: LIS_SSMISNWDsnowobsId   = "SSMI snow depth"
   character*50, public,  parameter :: LIS_AMSREsweobsId       = "AMSR-E SWE"
!   character*50, public,  parameter :: LIS_AMSREsnowobsId      = "AMSR-E snow" !yliu
   character*50, public,  parameter :: LIS_PMWsnowobsId        = "PMW snow" !yliu
   character*50, public,  parameter :: LIS_modisscfId          = "MODIS SCF"
   character*50, public,  parameter :: LIS_GRACEtwsobsId       = "GRACE TWS"
   character*50, public,  parameter :: LIS_simGRACEJPLobsId    = "Simulated GRACE (JPL)"
   character*50, public,  parameter :: LIS_synLbandTbobsId     = "Synthetic L-band Tb"
   character*50, public,  parameter :: LIS_SMOPSsmobsId        = "SMOPS soil moisture"
   character*50, public,  parameter :: LIS_SMOPS_ASCATsmobsId  = "SMOPS-ASCAT soil moisture" ! MN
   character*50, public,  parameter :: LIS_SMOPS_SMOSsmobsId   = "SMOPS-SMOS soil moisture"  ! MN
   character*50, public,  parameter :: LIS_SMOPS_AMSR2smobsId  = "SMOPS-AMSR2 soil moisture" ! MN
   character*50, public,  parameter :: LIS_SMOPS_SMAPsmobsId   = "SMOPS-SMAP soil moisture"  ! MN
   character*50, public,  parameter :: LIS_ASCAT_TUWsmobsId    = "ASCAT (TUW) soil moisture"
   character*50, public,  parameter :: LIS_IMSscaobsId         = "IMS snow cover"
   character*50, public,  parameter :: LIS_GCOMW_AMSR2L3smobsId = &
        "GCOMW AMSR2 L3 soil moisture"
   character*50, public,  parameter :: LIS_GCOMW_AMSR2L3sndobsId = &
        "GCOMW AMSR2 L3 snow depth"
   character*50, public,  parameter :: LIS_SMOSL2smobsId         = &
        "SMOS L2 soil moisture"
   character*50, public,  parameter :: LIS_pildassmobsId         = &
        "PILDAS SM"
   character*50, public,  parameter :: LIS_SMOSNESDISsmobsId     = &
        "SMOS(NESDIS) soil moisture"
   character*50, public,  parameter :: LIS_NASASMAPsmobsId       = &
        "SMAP(NASA) soil moisture"
   character*50, public,  parameter :: LIS_NASASMAPvodobsId      = &
        "SMAP(NASA) vegetation optical depth"
   character*50, public,  parameter :: LIS_GLASSlaiobsId         = &
        "GLASS LAI"
   character*50, public,  parameter :: LIS_MODISsportLAIobsId    = &
        "MODIS SPoRT LAI"
   character*50, public,  parameter :: LIS_GLASSalbedoobsId      = &
        "GLASS Albedo"
   character*50, public,  parameter :: LIS_SMAPNRTsmobsId        = &
        "SMAP(NRT) soil moisture"
   character*50, public,  parameter :: LIS_ASOsweobsId           = &
        "ASO SWE"
!-------------------------------------------------------------------------
! Bias Estimation Algorithms
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_gmaobiasId = "Adaptive bias correction"
!-------------------------------------------------------------------------
! Optimization / Uncertainty estimation Algorithms
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_ESOptId ="Enumerated search"
   character*50, public,  parameter :: LIS_LMOptId ="Levenberg marquardt"
   character*50, public,  parameter :: LIS_GAOptId = "Genetic algorithm"
   character*50, public,  parameter :: LIS_SCEUAOptId = "Shuffled complex evolution"
   character*50, public,  parameter :: LIS_MCSIMId = "Monte carlo sampling"
   character*50, public,  parameter :: LIS_RWMCMCId = "Random walk markov chain monte carlo"
   character*50, public,  parameter :: LIS_DEMCId = "Differential evolution markov chain"
   character*50, public,  parameter :: LIS_DEMCzId = "Differential evolution markov chain z"
!-------------------------------------------------------------------------
! Optimization Set (parameter estimation calibration dataset) definitions
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_EmptyObsId = "No obs"
   character*50, public,  parameter :: LIS_templateObsId = "NONE"
   character*50, public,  parameter :: LIS_wgPBMRsmId = "WG PBMR sm"
   character*50, public,  parameter :: LIS_pesynsm1Id = "Synthetic sm1"
   character*50, public,  parameter :: LIS_pesynsm2Id = "Synthetic sm2"
   character*50, public,  parameter :: LIS_AmerifluxObsId = "Ameriflux obs"
   character*50, public,  parameter :: LIS_ARMObsId = "ARM obs"
   character*50, public,  parameter :: LIS_maconlandslideObsId = "Macon landslide obs"
   character*50, public,  parameter :: LIS_globallandslideObsId = "Global landslide obs"
   character*50, public,  parameter :: LIS_CNRSObsId = "CNRS"
   character*50, public,  parameter :: LIS_CNRS_MPDIObsId = "CNRS MPDI"
   character*50, public,  parameter :: LIS_AMSRE_SRObsId = "AMSRE SR"
   character*50, public,  parameter :: LIS_LPRM_AMSREsmpeObsId = "AMSR-E(LPRM) pe soil moisture"
   character*50, public,  parameter :: LIS_FLUXNETpeObsId = "Gridded FLUXNET"
   character*50, public,  parameter :: LIS_USDA_ARSsmpeObsId = "USDA ARSsm"
   character*50, public,  parameter :: LIS_ARSsmobsId = "ARS sm" ! SY
   character*50, public,  parameter :: LIS_ISMNsmobsId = "ISMN sm" 
   character*50, public,  parameter :: LIS_SMAPsmobsId = "SMAP sm"
!-------------------------------------------------------------------------
! Objective Function Evaluation Criteria
!-------------------------------------------------------------------------
   character*50, public,  parameter :: LIS_LMestimateId = "LM"
   character*50, public,  parameter :: LIS_LSestimateId = "Least squares"
   character*50, public,  parameter :: LIS_LLestimateId = "Likelihood"
   character*50, public,  parameter :: LIS_PestimateId  = "Probability"
!-------------------------------------------------------------------------
! Radiative Transfer Models
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_crtmId = "CRTM"
   character*50, public, parameter :: LIS_crtm2Id = "CRTM2"
   character*50, public, parameter :: LIS_crtm2EMId = "CRTM2EM"
   character*50, public, parameter :: LIS_cmem3Id = "CMEM"
   character*50, public, parameter :: LIS_tauomegaRTMId = "Tau Omega"
!-------------------------------------------------------------------------
! Land Slide Models
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_GLSId = "GLS"
   character*50, public, parameter :: LIS_TRIGRSId = "TRIGRS"
!-------------------------------------------------------------------------
! Routing Models
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_NLDASrouterId = "NLDAS router"
   character*50, public, parameter :: LIS_HYMAProuterId = "HYMAP router"
   character*50, public, parameter :: LIS_HYMAP2routerId = "HYMAP2 router"
!-------------------------------------------------------------------------
! Runoff data support
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_LISrunoffdataId   = "LIS runoff output"
   character*50, public, parameter :: LIS_GLDAS1runoffdataId = "GLDAS1 runoff data"
   character*50, public, parameter :: LIS_GLDAS2runoffdataId = "GLDAS2 runoff data"
   character*50, public, parameter :: LIS_NLDAS2runoffdataId = "NLDAS2 runoff data"
   character*50, public, parameter :: LIS_MERRA2runoffdataId = "MERRA2 runoff data"
   character*50, public, parameter :: LIS_ERAIlandrunoffdataId = "ERA interim land runoff data"
   character*50, public, parameter :: LIS_GWBMIPrunoffdataId = "GWB MIP runoff data"
!-------------------------------------------------------------------------
!  Irrigation models
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_sprinklerIrrigationId = "Sprinkler"
   character*50, public, parameter :: LIS_floodIrrigationId = "Flood"
   character*50, public, parameter :: LIS_dripIrrigationId  = "Drip"
!-------------------------------------------------------------------------
!  Forecasting algorithms
!-------------------------------------------------------------------------
   character*50, public, parameter :: LIS_ESPbootId = "ESP boot"
   character*50, public, parameter :: LIS_ESPconvId = "ESP conventional"

!EOC
 end module LIS_pluginIndices
