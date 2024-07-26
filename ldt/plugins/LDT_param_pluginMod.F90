!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_param_pluginMod
!BOP
!
! !MODULE: LDT_param_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions used for
!   defining routines to read various sources of parameters maps.
!   The user defined functions are incorporated into the
!   appropriate registry to be later invoked through generic calls.
!
!
! !REVISION HISTORY:
!  11 Dec 2003:  Sujay Kumar  - Initial Specification
!  11 Feb 2013:  KR Arsenault - Updated to accommodate new parameter types and options
!  01 Mar 2020:  Yeosang Yoon - Added MERIT DEM
!  29 Jun 2020:  Mahdi Navari - Glacier fraction added 
!  12 Apr 2021:  Wanshu Nie   - groundwater irrigation ratio added
!  28 Jun 2022:  Eric Kemp    - Added NAFPA background precipitation
!EOP

  use LDT_pluginIndices

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------

  PUBLIC :: LDT_landcover_plugin
  PUBLIC :: LDT_soils_plugin
  PUBLIC :: LDT_topo_plugin
  PUBLIC :: LDT_laisai_plugin
  PUBLIC :: LDT_irrigation_plugin
  PUBLIC :: LDT_gfrac_plugin
  PUBLIC :: LDT_alb_plugin

  PUBLIC :: LDT_LSMparam_plugin
  PUBLIC :: LDT_routingparam_plugin
  PUBLIC :: LDT_lakeparam_plugin
  PUBLIC :: LDT_climate_plugin
  PUBLIC :: LDT_forcingparams_plugin

  PUBLIC :: LDT_glacier_plugin

contains

!BOP
! !ROUTINE: LDT_LSMparam_plugin
!  \label{LDT_LSMparam_plugin}
!
! !DESCRIPTION:
!
! !INTERFACE:
  subroutine LDT_LSMparam_plugin
!EOP

    use Noah_parmsMod
    use CLSMF25_parmsMod
    use RDHM_parmsMod
    use SACHTET_parmsMod
    use Snow17_parmsMod
    use GeoWRSI_parmsMod
    use VIC_parmsMod
    use SiB2_parmsMod
    use CLM2_parmsMod
    use CLM45_parmsMod
    use Mosaic_parmsMod
    use RUC_parmsMod
    use JULES50_parmsMod
    use Crocus_parmsMod
    use SnowModel_parmsMod

    external :: registerlsmparamprocinit
    external :: registerlsmparamprocwriteheader
    external :: registerlsmparamprocwritedata
    
  ! Noah 2.7.1 LSM:
    call registerlsmparamprocinit(trim(LDT_noah271Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noah271Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noah271Id)//char(0),&
         NoahParms_writeData)

  ! Noah 3.2 LSM:
    call registerlsmparamprocinit(trim(LDT_noah32Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noah32Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noah32Id)//char(0),&
         NoahParms_writeData)

  ! Noah 3.3 LSM:
    call registerlsmparamprocinit(trim(LDT_noah33Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noah33Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noah33Id)//char(0),&
         NoahParms_writeData)

  ! Noah 3.6 LSM:
    call registerlsmparamprocinit(trim(LDT_noah36Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noah36Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noah36Id)//char(0),&
         NoahParms_writeData)

  ! Noah 3.9 LSM:
    call registerlsmparamprocinit(trim(LDT_noah39Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noah39Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noah39Id)//char(0),&
         NoahParms_writeData)

  ! Noah-MP (v3.6) LSM:
    call registerlsmparamprocinit(trim(LDT_noahmp36Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noahmp36Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noahmp36Id)//char(0),&
         NoahParms_writeData)

  ! Noah-MP (v4.0.1) LSM:
    call registerlsmparamprocinit(trim(LDT_noahmp401Id)//char(0),&
         NoahParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_noahmp401Id)//char(0),&
         NoahParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_noahmp401Id)//char(0),&
         NoahParms_writeData)

  ! CLSM F2.5 LSM:
    call registerlsmparamprocinit(trim(LDT_clsmf25Id)//char(0),&
         catchmentParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_clsmf25Id)//char(0),&
         catchmentParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_clsmf25Id)//char(0),&
         catchmentParms_writeData)

  ! RDHM:

    ! ONLY SAC-HTET 3.5.6 LSM:
    call registerlsmparamprocinit(trim(LDT_sachtet356Id)//char(0),&
         SACHTETParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_sachtet356Id)//char(0),&
         SACHTETParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_sachtet356Id)//char(0),&
         SACHTETParms_writeData)

    ! ONLY Snow-17:
    call registerlsmparamprocinit(trim(LDT_snow17Id)//char(0),&
         Snow17Parms_init)
    call registerlsmparamprocwriteheader(trim(LDT_snow17Id)//char(0),&
         Snow17Parms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_snow17Id)//char(0),&
         Snow17Parms_writeData)

    ! RDHM: SAC-HTET 3.5.6 and Snow-17 combined:
    call registerlsmparamprocinit(trim(LDT_rdhm356Id)//char(0),&
         RDHMParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_rdhm356Id)//char(0),&
         RDHMParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_rdhm356Id)//char(0),&
         RDHMParms_writeData)

  ! GeoWRSI 2.0:
    call registerlsmparamprocinit(trim(LDT_GeoWRSI2Id)//char(0),&
         GeoWRSIParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_GeoWRSI2Id)//char(0),&
         GeoWRSIParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_GeoWRSI2Id)//char(0),&
         GeoWRSIParms_writeData)

  ! VIC LSM versions:
    call registerlsmparamprocinit(trim(LDT_vic411Id)//char(0),&
         VICParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_vic411Id)//char(0),&
         VICParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_vic411Id)//char(0),&
         VICParms_writeData)

    call registerlsmparamprocinit(trim(LDT_vic412Id)//char(0),&
         VICParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_vic412Id)//char(0),&
         VICParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_vic412Id)//char(0),&
         VICParms_writeData)

  ! SiB2:
    call registerlsmparamprocinit(trim(LDT_SiB2Id)//char(0),&
         SiB2Parms_init)
    call registerlsmparamprocwriteheader(trim(LDT_SiB2Id)//char(0),&
         SiB2Parms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_SiB2Id)//char(0),&
         SiB2Parms_writeData)

  ! CLM v2 LSM:
    call registerlsmparamprocinit(trim(LDT_clm2Id)//char(0),&
         CLM2Parms_init)
    call registerlsmparamprocwriteheader(trim(LDT_clm2Id)//char(0),&
         CLM2Parms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_clm2Id)//char(0),&
         CLM2Parms_writeData)

  ! CLM v4.5 LSM:
    call registerlsmparamprocinit(trim(LDT_clm45Id)//char(0),&
         CLM45Parms_init)
    call registerlsmparamprocwriteheader(trim(LDT_clm45Id)//char(0),&
         CLM45Parms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_clm45Id)//char(0),&
         CLM45Parms_writeData)

  ! Mosaic LSM:
    call registerlsmparamprocinit(trim(LDT_mosaicId)//char(0),&
         MosaicParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_mosaicId)//char(0),&
         MosaicParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_mosaicId)//char(0),&
         MosaicParms_writeData)

  ! RUC 3.7 LSM:
    call registerlsmparamprocinit(trim(LDT_ruc37Id)//char(0),&
         RUCParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_ruc37Id)//char(0),&
         RUCParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_ruc37Id)//char(0),&
         RUCParms_writeData)

  ! JULES 5.0 LSM:
    call registerlsmparamprocinit(trim(LDT_jules50Id)//char(0),&
         JULES50Parms_init)
    call registerlsmparamprocwriteheader(trim(LDT_jules50Id)//char(0),&
         JULES50Parms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_jules50Id)//char(0),&
         JULES50Parms_writeData)

  ! Crocus 8.1 :
    call registerlsmparamprocinit(trim(LDT_Crocus81Id)//char(0),&
        CrocusParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_Crocus81Id)//char(0),&
         CrocusParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_Crocus81Id)//char(0),&
         CrocusParms_writeData)

  ! SnowModel LSM:
    call registerlsmparamprocinit(trim(LDT_snowmodelId)//char(0),&
         SnowModelParms_init)
    call registerlsmparamprocwriteheader(trim(LDT_snowmodelId)//char(0),&
         SnowModelParms_writeHeader)
    call registerlsmparamprocwritedata(trim(LDT_snowmodelId)//char(0),&
         SnowModelParms_writeData)

  end subroutine LDT_LSMparam_plugin

!BOP
! !ROUTINE: LDT_routingparam_plugin
!  \label{LDT_Routingparam_plugin}
!
! !DESCRIPTION:
!
! !INTERFACE:
  subroutine LDT_routingparam_plugin
!EOP

    use HYMAP_parmsMod   ! Set for both HYMAP 1 and 2

    ! HYMAP - version 1
    call registerroutingparamprocinit(trim(LDT_HYMAPId)//char(0),&
         HYMAPParms_init)
    call registerroutingparamprocwriteheader(trim(LDT_HYMAPId)//char(0),&
         HYMAPParms_writeHeader)
    call registerroutingparamprocwritedata(trim(LDT_HYMAPId)//char(0),&
         HYMAPParms_writeData)

    ! HYMAP - version 2 (No difference to version 1 at this time ...)
    call registerroutingparamprocinit(trim(LDT_HYMAP2Id)//char(0),&
         HYMAPParms_init)
    call registerroutingparamprocwriteheader(trim(LDT_HYMAP2Id)//char(0),&
         HYMAPParms_writeHeader)
    call registerroutingparamprocwritedata(trim(LDT_HYMAP2Id)//char(0),&
         HYMAPParms_writeData)

  end subroutine LDT_routingparam_plugin

!BOP
! !ROUTINE: LDT_lakeparam_plugin
!  \label{LDT_Lakeparam_plugin}
!
! !DESCRIPTION:
!
! !INTERFACE:
  subroutine LDT_lakeparam_plugin
!EOP
    use FLAKE_parmsMod

    call registerlakeparamprocinit(trim(LDT_flakeId)//char(0),&
         FLAKEparms_init)
    call registerlakeparamprocwriteheader(trim(LDT_flakeId)//char(0),&
         FLAKEparms_writeHeader)
    call registerlakeparamprocwritedata(trim(LDT_flakeId)//char(0),&
         FLAKEparms_writeData)

  end subroutine LDT_lakeparam_plugin

!BOP
! !ROUTINE: LDT_landcover_plugin
!  \label{LDT_landcover_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new landcover datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the landmask data]
!      Routines to retrieve the landmask data
!      (to be registered using {\tt registerreadmask} and later called
!       using the generic {\tt readmask} method)
!  \item[read the regional mask data]
!      Routines to retrieve the regional mask data
!      (to be registered using {\tt registerreadregmask} and later called
!       using the generic {\tt readregmask} method)
!  \item[read the landcover data]
!      Routines to retrieve the landcover data
!      (to be registered using {\tt registerreadlc} and later called
!       using the generic {\tt readlandcover} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the landcover
!  datasets from the vegetation source with LIS using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LIS uses and index of 1
!
!  \begin{verbatim}
!    call registerreadmask(1,1,read_latlon_avhrrmask)
!    call registerreadlc(1,1,read_latlon_avhrrlc)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readmask(1,1) - calls read_latlon_avhrrmask
!    call readlc(1,1)   - calls read_latlon_avhrrlc
!  \end{verbatim}
!
!   In the LIS code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readmask(ldt%source,ldt%nest,lis%datavalues)
!     call readregmask(ldt%source,ldt%nest,lis%datavalues)
!     call readlandcover(ldt%source,ldt%nest,lis%datavalues)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%vegsrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_landcover_plugin
!EOP

    use Monfredaetal08_crops_module, only : read_Monfredaetal08_croptype

    external set_AVHRR_lc_attribs
    external read_avhrr_lc

    external read_avhrr_gfs_lc

    external read_USGS_lc
    external read_USGSNative_lc

    external set_MODISNative_lc_attribs
    external read_MODISNative_lc
    external read_MCD12Q1_lc
    external read_MODISNative_PFT
    external read_UKMO_IGBP_PFT
    external read_UM_ancillary
    external read_MODIS_lc

    external read_ALMIPII_lc
    external read_CLSMF25_lc
    external read_VIC411_lc
    external read_VIC412_lc
    external read_ISA_lc
    external read_CONSTANT_lc
    external read_SACHTET356_lc
    external read_CLM45_lc
    external read_NALCMS_SM_lc
    external read_NALCMS_SM_IGBPNCEP_lc

    external read_regmask_gis
    external read_regmask_wrsi
    external read_regmask_mask
!    external read_VIC411_maskfile

    external read_UMDCROPMAP_croptype
!    external read_Monfredaetal08_croptype

    external read_ALMIPII_droot
!    external read_UMDCROPMAP_rootdepth
! _______________________________________

! - Landcover sources:
  ! AVHRR-UMD (LIS-based):
    call registerreadlc(trim(LDT_avhrrlcLISId)//char(0),read_AVHRR_lc)

  ! AVHRR-UMD (GFS-domain):
    call registerreadlc(trim(LDT_avhrrlcGFSId)//char(0),read_AVHRR_GFS_lc)

  ! MODIS-IGBP/NCEP (Native):
    call registerreadlc(trim(LDT_modislcNATId)//char(0), read_MODISNative_lc)

  ! MODIS-IGBP/NCEP (Native) to JULES PFT:
    call registerreadlc(trim(LDT_modislcPFTId)//char(0), read_MODISNative_PFT)

  ! UKMO IGBP to JULES PFT:
    call registerreadlc(trim(LDT_ukmoigbpPFTId)//char(0), read_UKMO_IGBP_PFT)

  ! UM ancillary
    call registerreadlc(trim(LDT_UM_ancillaryId)//char(0), read_UM_ancillary)

  ! MODIS-IGBP/NCEP (LIS-based):
    call registerreadlc(trim(LDT_modislcLISId)//char(0), read_MODIS_lc)
    call registerreadlc(trim(LDT_mcd12q1Id)//char(0), read_MCD12Q1_lc)

  ! USGS (Native):
    call registerreadlc(trim(LDT_usgslcNATId)//char(0), read_USGSNative_lc)
  ! USGS (LIS-based):
    call registerreadlc(trim(LDT_usgslcLISId)//char(0), read_USGS_lc)

  ! ALMIPII (ECONOMAP):
    call registerreadlc(trim(LDT_ALMIPIIlcId)//char(0), read_ALMIPII_lc)
  ! Catchment LSM - F2.5:
    call registerreadlc(trim(LDT_clsmf25lcId)//char(0), read_CLSMF25_lc)
  ! VIC-4.1.1:
    call registerreadlc(trim(LDT_vic411lcId)//char(0), read_VIC411_lc)
  ! VIC-4.1.2:
    call registerreadlc(trim(LDT_vic412lcId)//char(0), read_VIC412_lc)
  ! ISA:
    call registerreadlc(trim(LDT_isalcId)//char(0), read_ISA_lc)
  ! SACHTET Landcover:
    call registerreadlc(trim(LDT_sachtet356Id)//char(0), read_SACHTET356_lc)
  ! CLM-4.5 Landcover:
    call registerreadlc(trim(LDT_clm45lcId)//char(0), read_CLM45_lc)

  ! SnowModel-based NALCMS Landcover:
    call registerreadlc(trim(LDT_nalcmsSMlcId)//char(0), read_NALCMS_SM_lc)
  ! NALCMS/SnowModel-mapped to IGBP/NCEP Landcover:
    call registerreadlc(trim(LDT_nalcmsSMIGBPlcId)//char(0), read_NALCMS_SM_IGBPNCEP_lc)

  ! Constant Landcover:
    call registerreadlc(trim(LDT_constId)//char(0), read_CONSTANT_lc)

! - Mask sources:
    call registerreadregmask(trim(LDT_gismaskId)//char(0), read_regmask_gis)
    call registerreadregmask(trim(LDT_wrsimaskId)//char(0), read_regmask_wrsi)
    call registerreadregmask(trim(LDT_maskmaskId)//char(0), read_regmask_mask)
!    call registerreadregmask(trim(LDT_vic411maskId)//char(0), read_VIC411_maskfile)

! - Crop type sources:
    call registerreadcroptype(trim(LDT_umdcropmapId)//char(0), read_UMDCROPMAP_croptype)
    call registerreadcroptype(trim(LDT_monfredacropId)//char(0), read_Monfredaetal08_croptype)

! - Root depth:
    call registerreadrootdepth(trim(LDT_ALMIPIIlcId)//char(0),read_ALMIPII_droot)
!    call registerreadrootdepth(trim(LDT_umdcropmapId)//char(0),read_UMDCROPMAP_droot)

  end subroutine LDT_landcover_plugin

!BOP
! !ROUTINE: LDT_topo_plugin
!  \label{LDT_topo_plugin}
!
! !DESCRIPTION:
!
! This is a plugin point for introducing new topography datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the elevation data]
!      Routines to retrieve the elevation data
!      (to be registered using {\tt registerreadelev} and later called
!       using the generic {\tt readelev} method)
!  \item[read the slope  data]
!      Routines to retrieve the slope data
!      (to be registered using {\tt registerreadslope} and later called
!       using the generic {\tt readslope} method)
!  \item[read the aspect data]
!      Routines to retrieve the aspect data
!      (to be registered using {\tt registerreadaspect} and later called
!       using the generic {\tt readaspect} method)
!  \item[read the curvature data]
!      Routines to retrieve the curvature data
!      (to be registered using {\tt registerreadcurv} and later called
!       using the generic {\tt readcurv} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the topography
!  datasets from the GTOPO30 source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows,
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadelev(1,1,read_elev_gtopo30)
!    call registerreadslope(1,1,read_slope_gtopo30)
!    call registerreadaspect(1,1,read_aspect_gtopo30)
!    call registerreadcurv(1,1,read_curv_gtopo30)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readelev(1,1)   - calls read_elev_gtopo30
!    call readslope(1,1)  - calls read_slope_gtopo30
!    call readaspect(1,1)  - calls read_aspect_gtopo30
!    call readcurv(1,1)  - calls read_curv_gtopo30
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readelev(ldt%domain, ldt%toposrc)
!    call readslope(ldt%domain, ldt%toposrc)
!    call readaspect(ldt%domain, ldt%toposrc)
!    call readcurv(ldt%domain, ldt%toposrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%toposrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_topo_plugin
!EOP

    external read_GTOPO30_elev
    external read_GTOPO30_GFS_elev
    external read_GTOPO30Native_elev
    external read_GTOPO30_slope
    external read_GTOPO30_aspect
    external read_GTOPO30_curv

    external read_SRTM_elev
    external read_SRTM_slope
    external read_SRTM_aspect
    external read_SRTM_Native_elev
    external read_SRTM_Native_slope
    external read_SRTM_Native_aspect

    external read_CONSTANT_elev
    external read_CONSTANT_slope
    external read_CONSTANT_aspect

    external read_MERIT1K_elev
    external read_MERIT1K_slope
    external read_MERIT1K_aspect

    external read_NED_SM_elev
    external read_NED_SM_slope
    external read_NED_SM_aspect
    external read_NED_SM_curvature

 !- GTOPO30:
    call registerreadelev(trim(LDT_gtopoLISId)//char(0),read_GTOPO30_elev)
    call registerreadelev(trim(LDT_gtopoGFSId)//char(0),read_GTOPO30_GFS_elev)
    call registerreadelev(trim(LDT_gtopoNATId)//char(0),read_GTOPO30Native_elev)
    call registerreadslope(trim(LDT_gtopoLISId)//char(0),read_GTOPO30_slope)
    call registerreadaspect(trim(LDT_gtopoLISId)//char(0),read_GTOPO30_aspect)
    call registerreadcurv(trim(LDT_gtopoLISId)//char(0),read_GTOPO30_curv)

 !- SRTM:
    call registerreadelev(trim(LDT_srtmLISId)//char(0),read_SRTM_elev)
    call registerreadslope(trim(LDT_srtmLISId)//char(0),read_SRTM_slope)
    call registerreadaspect(trim(LDT_srtmLISId)//char(0),read_SRTM_aspect)
    call registerreadelev(trim(LDT_srtmNATId)//char(0),read_SRTM_Native_elev)
    call registerreadslope(trim(LDT_srtmNATId)//char(0),read_SRTM_Native_slope)
    call registerreadaspect(trim(LDT_srtmNATId)//char(0),read_SRTM_Native_aspect)

 !- Constant value:
    call registerreadelev(trim(LDT_constId)//char(0),read_CONSTANT_elev)
    call registerreadslope(trim(LDT_constId)//char(0),read_CONSTANT_slope)
    call registerreadaspect(trim(LDT_constId)//char(0),read_CONSTANT_aspect)

!- MERIT:
    call registerreadelev(trim(LDT_merit1KId)//char(0),read_MERIT1K_elev)
    call registerreadslope(trim(LDT_merit1KId)//char(0),read_MERIT1K_slope)
    call registerreadaspect(trim(LDT_merit1KId)//char(0),read_MERIT1K_aspect)

!- NED (SnowModel file version) elevation:
    call registerreadelev(trim(LDT_nedSMId)//char(0),read_NED_SM_elev)
    call registerreadslope(trim(LDT_nedSMId)//char(0),read_NED_SM_slope)
    call registerreadaspect(trim(LDT_nedSMId)//char(0),read_NED_SM_aspect)
    call registerreadcurv(trim(LDT_nedSMId)//char(0),read_NED_SM_curvature)

  end subroutine LDT_topo_plugin


!BOP
! !ROUTINE: LDT_soils_plugin
!  \label{LDT_soils_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new soil parameter datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the texture data]
!      Routines to retrieve the texture data
!      (to be registered using {\tt registerreadsoiltexture} and later called
!       using the generic {\tt readtexture} method)
!  \item[read the sand fraction data]
!      Routines to retrieve the sand fraction data
!      (to be registered using {\tt registerreadsand} and later called
!       using the generic {\tt readsand} method)
!  \item[read the clay fraction data]
!      Routines to retrieve the clay fraction data
!      (to be registered using {\tt registerreadclay} and later called
!       using the generic {\tt readclay} method)
!  \item[read the silt fraction data]
!      Routines to retrieve the silt fraction data
!      (to be registered using {\tt registerreadsilt} and later called
!       using the generic {\tt readsilt} method)
!  \item[read the soil color data]
!      Routines to retrieve the soil color data
!      (to be registered using {\tt registerreadcolor} and later called
!       using the generic {\tt readcolor} method)
!  \item[read the soil porosity data]
!      Routines to retrieve the soil porosity data
!      (to be registered using {\tt registerreadporosity} and later called
!       using the generic {\tt readporosity} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the soil
!  datasets from the FAO source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows (only a
!  subset of the parameters is shown below, others can be defined in
!  the same mannter),
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadsand(1,1,read_FAOsand)
!    call registerreadclay(1,1,read_FAOclay)
!    call registerreadsilt(1,1,read_FAOsilt)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readsand(1,1)  - calls read_FAOsand
!    call readclay(1,1)  - calls read_FAOclay
!    call readsilt(1,1)  - calls read_FAOsilt
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readsand(ldt%domain, ldt%soilsrc)
!    call readclay(ldt%domain, ldt%soilsrc)
!    call readsilt(ldt%domain, ldt%soilsrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%soilsrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_soils_plugin
!EOP

  ! FAO-only
    external read_FAO_soilfractions, read_FAO_texture, read_FAO_color, &
             read_FAO_porosity

  ! STATSGO+FAO blended global product (NCAR/NCEP)
    external set_STATSGOFAO_LIS_texture_attribs
    external read_STATSGOFAO_texture

    external set_STATSGOFAO_Native_texture_attribs
    external set_ZOBLER_GFS_texture_attribs
    external read_STATSGOFAO_Native_texture
    external read_ZOBLER_GFS_texture

  ! Original LIS-STATSGO v1 (USDA/PSU)
    external set_STATSGOv1_texture_attribs
    external set_STATSGOv1_hsg_attribs
    external read_STATSGOv1_soilfractions, read_STATSGOv1_texture
    external read_STATSGOv1_porosity, read_STATSGOv1_bedrockdepth
    external read_STATSGOv1_hydrosoilgroup, ed_STATSGOv1_bulkdensity
    external read_STATSGOv1_domrockfrag, read_STATSGOv1_rockvolfrag
    external read_STATSGOv1_availwatcap

  ! Native STATSGO v2 (USDA)
!    external read_STATSGOv2_texture, read_STATSGOv2_soilfractions

  ! ALMIPII
    external read_ALMIPII_soilfractions, read_ALMIPII_dsoil !, read_ALMIPII_drootlayer

  ! CLSM Fortuna 2.5 readers
    external read_CLSMF25_porosity

  ! Constant cases
    external set_CONSTANT_texture_attribs
    external read_CONSTANT_texture, read_CONSTANT_soilfractions, read_CONSTANT_porosity

  ! Special cases
    external set_Special_texture_attribs
    external read_Special_texture
    external read_Special_soilfractions

    external set_ISRIC_texture_attribs
    external read_ISRIC_texture
    external read_ISRIC_soilfractions

!== FAO (global; original LIS-processed):
    call registerreadsoilfrac(trim(LDT_faoSoilId)//char(0),read_FAO_soilfractions)
    call registerreadsoiltexture(trim(LDT_faoSoilId)//char(0), read_FAO_texture)
    call registerreadcolor(trim(LDT_faoSoilId)//char(0),read_FAO_color)
    call registerreadporosity(trim(LDT_faoSoilId)//char(0),read_FAO_porosity)

!== STATSGO.v1 and FAO (merged STATSGOv1+FAO global maps): ==
  ! Original LIS domain:
    call registersettextureattribs(trim(LDT_statsgo1FAOLISId)//char(0),&
                                 set_STATSGOFAO_LIS_texture_attribs)
    call registerreadsoiltexture(trim(LDT_statsgo1FAOLISId)//char(0),&
                                 read_STATSGOFAO_texture)

  ! Native domain:
    call registersettextureattribs(trim(LDT_statsgo1FAONATId)//char(0),&
         set_STATSGOFAO_Native_texture_attribs)
    call registerreadsoiltexture(trim(LDT_statsgo1FAONATId)//char(0),&
         read_STATSGOFAO_Native_texture)
  ! old GFS files
    call registersettextureattribs(trim(LDT_zoblerGFSId)//char(0),&
         set_ZOBLER_GFS_texture_attribs)
    call registerreadsoiltexture(trim(LDT_zoblerGFSId)//char(0),&
                                 read_ZOBLER_GFS_texture)

!== STATSGO v1 (CONUS data only): ==

  ! Original LIS domain:
    call registerreadsoilfrac(trim(LDT_statsgov1LISId)//char(0),&
                              read_STATSGOv1_soilfractions)

    call registerreadcolor(trim(LDT_statsgov1LISId)//char(0),read_FAO_color)
  ! for now using the same color dataset

  ! Native domain:
    call registersettextureattribs(trim(LDT_statsgov1NATId)//char(0),&
         set_STATSGOv1_texture_attribs)

    call registerreadsoiltexture(trim(LDT_statsgov1NATId)//char(0),&
                                 read_STATSGOv1_texture)

!    call registerreadsoilfrac(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_soilfractions)

!    call registerreadporosity(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_porosity)

!    call registerreadbdrckdpth(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_bedrockdepth)

    call registersethsgattribs(trim(LDT_statsgov1NATId)//char(0),&
         set_STATSGOv1_hsg_attribs)

    call registerreadhsg(trim(LDT_statsgov1NATId)//char(0),&
                              read_STATSGOv1_hydrosoilgroup)

!    call registerreadbulkdensity(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_bulkdensity)

!    call registerreaddomrockfrag(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_domrockfrag)

!    call registerreadrockvolfrag(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_rockvolfrag)

!    call registerreadavailwatcap(trim(LDT_statsgov1NATId)//char(0),&
!                              read_STATSGOv1_availwatcap)

!== STATSGO v2 (CONUS data only for now): ==

!    call registerreadsoiltexture(trim(LDT_statsgov2SoilId)//char(0),&
!                                 read_STATSGOv2_texture)

!== Constant value case:
    call registersettextureattribs(trim(LDT_constId)//char(0),&
                                 set_CONSTANT_texture_attribs)
    call registerreadsoiltexture(trim(LDT_constId)//char(0),read_CONSTANT_texture)
    call registerreadsoilfrac(trim(LDT_constId)//char(0),read_CONSTANT_soilfractions)
    call registerreadporosity(trim(LDT_constId)//char(0),read_CONSTANT_porosity)

 !- ALMIPII:
    call registerreadsoilfrac(trim(LDT_ALMIPIIsoilId)//char(0),read_ALMIPII_soilfractions)
    call registerreadsoildepth(trim(LDT_ALMIPIIsoilId)//char(0),read_ALMIPII_dsoil)
!    call registerreadrootdepthlayer(trim(LDT_ALMIPIIsoilId)//char(0),read_ALMIPII_droot)

!== Catchment-Fortuna 2.5: ==
    call registerreadporosity(trim(LDT_clsmf25Id)//char(0),read_CLSMF25_porosity)

!== Special case study:
    call registersettextureattribs(trim(LDT_specialSoilId)//char(0),set_Special_texture_attribs)
    call registerreadsoiltexture(trim(LDT_specialSoilId)//char(0),read_Special_texture)
    call registerreadsoilfrac(trim(LDT_specialSoilId)//char(0),read_Special_soilfractions)
!=== ISRIC soils

    call registersettextureattribs(trim(LDT_ISRICsoilId)//char(0),&
                                 set_ISRIC_texture_attribs)
    call registerreadsoiltexture(trim(LDT_ISRICSoilId)//char(0), &
         read_ISRIC_texture)
    call registerreadsoilfrac(trim(LDT_ISRICsoilId)//char(0),&
         read_ISRIC_soilfractions)

  end subroutine LDT_soils_plugin


!BOP
! !ROUTINE: LDT_laisai_plugin
!  \label{LDT_laisai_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new LAI/SAI datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the LAI data]
!      Routines to retrieve the LAI data
!      (to be registered using {\tt registerreadlai} and later called
!       using the generic {\tt readlai} method)
!  \item[read the SAI data]
!      Routines to retrieve the SAI data
!      (to be registered using {\tt registerreadsai} and later called
!       using the generic {\tt readsai} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the LAI/SAI
!  datasets from the AVHRR source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadlai(1,1,read_ll_avhrrlai)
!    call registerreadsai(1,1,read_ll_avhrrsai)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readlai(1,1)  - calls read_ll_avhrrlai
!    call readsai(1,1)  - calls read_ll_avhrrsai
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readlai(ldt%domain,ldtlaisrc)
!    call readsai(ldt%domain,ldtlaisrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%laisrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_laisai_plugin
!EOP

    external set_AVHRR_lai_attribs
    external read_AVHRR_lai,read_AVHRR_sai

    external set_CLSMF25_lai_attribs
    external read_CLSMF25_lai

    external read_CLSMF25_laimax,read_CLSMF25_laimin
    external read_CONSTANT_lai, read_CONSTANT_sai

  ! AVHRR:
    call registersetlaiattribs(trim(LDT_avhrrlaiId)//char(0),&
         set_AVHRR_lai_attribs)
    call registerreadlai(trim(LDT_avhrrlaiId)//char(0),read_AVHRR_lai)
    call registerreadsai(trim(LDT_avhrrlaiId)//char(0),read_AVHRR_sai)

  ! CLSM F2.5:
    call registersetlaiattribs(trim(LDT_clsmf25laiId)//char(0),&
         set_CLSMF25_lai_attribs)
    call registerreadlai(trim(LDT_clsmf25laiId)//char(0),read_CLSMF25_lai)
    call registerreadlaimax(trim(LDT_clsmf25laiId)//char(0),read_CLSMF25_laimax)
    call registerreadlaimin(trim(LDT_clsmf25laiId)//char(0),read_CLSMF25_laimin)

  ! Constant values:
    call registerreadlai(trim(LDT_constId)//char(0),read_CONSTANT_lai)
    call registerreadsai(trim(LDT_constId)//char(0),read_CONSTANT_sai)

  end subroutine LDT_laisai_plugin


!
! !ROUTINE: LDT_irrigation_plugin
!  \label{LDT_irrigation_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new irrigation datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the irrigation type data]
!      Routines to retrieve the irrigation type data
!      (to be registered using {\tt registerreadirrigtype} and later called
!       using the generic {\tt readirrigtype} method)
!  \item[read the irrigation fraction data]
!      Routines to retrieve the irrigation fraction data
!      (to be registered using {\tt registerreadirrigfrac} and later called
!       using the generic {\tt readirrigfrac} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the irrigation
!  datasets from the MODIS source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadirrigtype(1,1,read_ll_modisirrigtype)
!    call registerreadirrigfrac(1,1,read_ll_modisirrigtype)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readirrigtype(1,1)  - calls read_ll_modisirrigtype
!    call readirrigfrac(1,1)  - calls read_ll_modisirrigfrac
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readirrigtype(ldt%domain,ldtirrigtypesrc)
!    call readirrigfrac(ldt%domain,ldtirrigfracsrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%irrig[var]src$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_irrigation_plugin
!EOP

    external read_OzdoganGutman_irrigfrac

    external read_GRIPC_irrigtype
    external read_GRIPC_irrigfrac
    external read_UserDerived_irrigfrac

    external read_USGSNative_irriggwratio

    call registerreadirrigfrac(trim(LDT_modOGirrigId)//char(0),&
         read_OzdoganGutman_irrigfrac)

    call registerreadirrigtype(trim(LDT_gripcirrigId)//char(0),read_GRIPC_irrigtype)
    call registerreadirrigfrac(trim(LDT_gripcirrigId)//char(0),read_GRIPC_irrigfrac)

    ! Added user-derived irrigation fraction input option:
    call registerreadirrigfrac(trim(LDT_userinputirrigId)//char(0),read_UserDerived_irrigfrac)
    ! Added irrigation groundwater ratio input option
    call registerreadirriggwratio(trim(LDT_irriggwratioId)//char(0),read_USGSNative_irriggwratio)

  end subroutine LDT_irrigation_plugin



!BOP
! !ROUTINE: LDT_gfrac_plugin
!  \label{LDT_gfrac_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new greenness datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the greenness data]
!      Routines to retrieve the greenness data
!      (to be registered using {\tt registerreadgfrac} and later called
!       using the generic {\tt readgfrac} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the greenness
!  datasets from the NCEP source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadgfrac(1,1,read_llgfrac)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readgfrac(1,1)  - calls read_llgfrac
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readgfrac(ldt%domain,ldt%gfracsrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%gfracsrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_gfrac_plugin
!EOP
    external read_AVHRR_gfrac, read_NCEP_shdmax, read_NCEP_shdmin

    external set_AVHRRNative_gfrac_attribs
    external read_AVHRRNative_gfrac

    external read_NCEPNative_shdmax, read_NCEPNative_shdmin

    external set_CLSMF25_gfrac_attribs
    external read_CLSMF25_gfrac

    external read_CLSMF25_gfracmax,read_CLSMF25_gfracmin
    external read_SACHTET356_gfrac
    external read_CONSTANT_gfrac, read_CONSTANT_gfracmax, read_CONSTANT_gfracmin

! !USES:
 !- Noah LSM:
    call registerreadgfrac(trim(LDT_gfracClimLISId)//char(0),read_AVHRR_gfrac)
    call registerreadgfrac(trim(LDT_gfracClimNATId)//char(0),read_AVHRRNative_gfrac)

    call registerreadshdmax(trim(LDT_gfracClimLISId)//char(0),read_NCEP_shdmax)
    call registerreadshdmax(trim(LDT_gfracClimNATId)//char(0),read_NCEPNative_shdmax)
    call registerreadshdmin(trim(LDT_gfracClimLISId)//char(0),read_NCEP_shdmin)
    call registerreadshdmin(trim(LDT_gfracClimNATId)//char(0),read_NCEPNative_shdmin)

 !- Catchment LSM:
    call registerreadgfrac(trim(LDT_gfracClsmf25Id)//char(0),read_CLSMF25_gfrac)
    call registerreadshdmax(trim(LDT_gfracClsmf25Id)//char(0),read_CLSMF25_gfracmax)
    call registerreadshdmin(trim(LDT_gfracClsmf25Id)//char(0),read_CLSMF25_gfracmin)

 !- SAC-HTET v3.5.6 LSM:
    call registerreadgfrac(trim(LDT_gfracSACHTETId)//char(0),read_SACHTET356_gfrac)

 !- Constant value:
    call registerreadgfrac(trim(LDT_constId)//char(0),read_CONSTANT_gfrac)
    call registerreadshdmax(trim(LDT_constId)//char(0),read_CONSTANT_gfracmax)
    call registerreadshdmin(trim(LDT_constId)//char(0),read_CONSTANT_gfracmin)

  end subroutine LDT_gfrac_plugin

!BOP
! !ROUTINE: LDT_alb_plugin
!  \label{LDT_alb_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new albedo datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the albedo climatology data]
!      Routines to retrieve the albedo climatology data
!      (to be registered using {\tt registerreadalbedo} and later called
!       using the generic {\tt readalbedo} method)
!  \item[read the max snow albedo data]
!      Routines to retrieve the max snow albedo data
!      (to be registered using {\tt registerreadmxsnoalb} and later called
!       using the generic {\tt readmxsnoalb} method)
!  \item[read the albedo NIR factor data]
!      Routines to retrieve the albedo NIR factor data
!      (to be registered using {\tt registerreadalbnir} and later called
!       using the generic {\tt readalbnir} method)
!  \item[read the albedo VIS factor data]
!      Routines to retrieve the albedo VIS factor data
!      (to be registered using {\tt registerreadalbvis} and later called
!       using the generic {\tt readalbvis} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of the albedo
!  datasets from the NCEP source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows
!  if the index of the source is defined to 1 and latlon projection in
!  LDT uses and index of 1
!
!  \begin{verbatim}
!    call registerreadalbedo(1,1,read_llalbedo)
!    call registerreadmxsnoalb(1,1,read_llmxsnoalb)
!  \end{verbatim}
!
!   The functions registered above are invoked using generic calls as
!   follows:
!
!  \begin{verbatim}
!    call readalbedo(1,1)    - calls read_llalbedo
!    call readmxsnoalb(1,1)  - calls read_llmxsnoalb
!  \end{verbatim}
!
!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!     call readalbedo(ldt%domain,ldt%albedoarc)
!     call readmxsnoalb(ldt%domain,ldt%albedosrc)
!   \end{verbatim}
!   where $ldt\%domain$ and $ldt\%albedosrc$ are set through the configuration
!   utility, enabling the user make a selection at runtime.
!
! !INTERFACE:
  subroutine LDT_alb_plugin
!EOP

    external read_Briegleb_albedo

    external read_BrieglebNative_monalbedo
    external read_BrieglebNative_qtralbedo

    external read_RobinsonKukla_mxsnoalb
    external read_GFS_mxsnoalb
    external read_RobinsonKuklaNative_mxsnoalb
    external read_BarlageNative_mxsnoalb
    external read_SACHTET356_mxsnoalb

    external read_CONSTANT_albedo
    external read_CONSTANT_mxsnoalb

 !- Monthly albedo climatologies:
    call registerreadalbedo(trim(LDT_albedoClimLISId)//char(0),&
         read_Briegleb_albedo)

    call registerreadalbedo(trim(LDT_albedoClimNATId)//char(0),&
         read_BrieglebNative_monalbedo)

    call registerreadalbedo(trim(LDT_albedoClimNATQtrId)//char(0),&
         read_BrieglebNative_qtralbedo)

 !- Max snow albedo (Noah LSM)
    call registerreadmxsnoalb(trim(LDT_albedoClimLISId)//char(0),&
         read_RobinsonKukla_mxsnoalb)
    call registerreadmxsnoalb(trim(LDT_albedoClimNATId)//char(0),&
         read_RobinsonKuklaNative_mxsnoalb)
    call registerreadmxsnoalb(trim(LDT_mxsnalbGFSId)//char(0),&
         read_GFS_mxsnoalb)
    call registerreadmxsnoalb(trim(LDT_mxsnalbBarlageId)//char(0),&
         read_BarlageNative_mxsnoalb)

 !- SAC-HTET (3.5.6)
    call registerreadmxsnoalb(trim(LDT_mxsnalbSACHTETId)//char(0),&
         read_SACHTET356_mxsnoalb)

 !- Constant value:
    call registerreadalbedo(trim(LDT_constId)//char(0),read_CONSTANT_albedo)
    call registerreadmxsnoalb(trim(LDT_constId)//char(0),read_CONSTANT_mxsnoalb)

 end subroutine LDT_alb_plugin


!BOP
! !ROUTINE: LDT_climate_plugin
!  \label{LDT_climate_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing new climate downscaling datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the precipitation climate downscaling data]
!      Routines to retrieve the climate downscaling data
!      (to be registered using {\tt registerreadclimppt} and later called
!       using the generic {\tt readclimppt} method)
!  \item[read the min temp climate downscaling data]
!      Routines to retrieve the min temp climatology data
!      (to be registered using {\tt registerreadclimtmin} and later called
!       using the generic {\tt readclimtmin} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of higher res PPT
!  climtology datasets from the PRISM source with LDT using a latlon projection.
!  The methods should be defined in the registry as follows:
!
!    call readclimppt(ldt%source,ldt%nest,lis%datavalues)
!    call readclimtmin(ldt%source,ldt%nest,lis%datavalues)
!    call readclimtmax(ldt%source,ldt%nest,lis%datavalues)
!
!  \end{verbatim}

! !INTERFACE:
  subroutine LDT_climate_plugin
!EOP
    use LDT_NAFPA_back_climpptMod, only: LDT_read_NAFPA_back_gfs_climppt, &
         LDT_read_NAFPA_back_galwem_climppt
    external read_PRISM_climppt
    external read_WorldClim_climppt
    external read_NLDAS_climppt

    external :: registerreadclimppt

! !USES:
!- Precipitation downscaling:
    call registerreadclimppt(trim(LDT_prismpptId)//char(0),&
         read_PRISM_climppt)

    call registerreadclimppt(trim(LDT_worldclimpptId)//char(0),&
         read_WorldClim_climppt)

    call registerreadclimppt(trim(LDT_nafpabackgfspptId)//char(0),&
         LDT_read_NAFPA_back_gfs_climppt)

    call registerreadclimppt(trim(LDT_nafpabackgalwempptId)//char(0),&
         LDT_read_NAFPA_back_galwem_climppt)

!- Temperature downscaling:

  end subroutine LDT_climate_plugin

!BOP
! !ROUTINE: LDT_forcingparams_plugin
!  \label{LDT_forcingparams_plugin}
!
! !DESCRIPTION:
! This is a plugin point for introducing forcing parameters datasets.
! The interface mandates that the following routines be implemented
! and registered for each parameter data source.
!
!  \begin{description}
!  \item[read the forcing parameter elevation data]
!      Routines to retrieve the forcing parameter elevation data
!      (to be registered using {\tt registerreadforcelem} and later called
!       using the generic {\tt readclimppt} method)
!  \end{description}
!
!  The user-defined functions are included in the registry using
!  two indices. For example, consider the incorporation of meteorological
!  terrain height data from the GDAS source using a latlon projection.
!  The methods should be defined in the registry as follows:
!
!    call readforcelev(ldt%source,ldt%nest,ldt%datavalues)
!
!  \end{verbatim}
!
  subroutine LDT_forcingparams_plugin
!
!EOP

    external read_gdas_elev
    external read_nldas2_elev
    external read_nam242_elev
    external read_princeton_elev
    external read_ecmwf_elev
    external read_merra2_elev
    external read_era5_elev
    external read_wrfoutv2_elev
    external read_wrfak_elev
!    external read_geos5_elev

! !USES:
! - Read forcing parameter: Elevation/terrain height

!- GDAS forcings:
    call registerreadforcelev(trim(LDT_gdasId)//char(0),&
         read_gdas_elev)

!- CONUS-only forcings:
    call registerreadforcelev(trim(LDT_nldas2Id)//char(0),&
         read_nldas2_elev)

!- NAM242 forcing:
    call registerreadforcelev(trim(LDT_nam242Id)//char(0),&
         read_nam242_elev)

!- Princeton forcing:
    call registerreadforcelev(trim(LDT_princetonId)//char(0),&
         read_princeton_elev)

!- ECMWF forcing:
    call registerreadforcelev(trim(LDT_ecmwfId)//char(0),&
         read_ecmwf_elev)

!- MERRA2 forcing:
    call registerreadforcelev(trim(LDT_merra2Id)//char(0),&
         read_merra2_elev)

!- ERA5 forcing:
    call registerreadforcelev(trim(LDT_era5Id)//char(0),&
         read_era5_elev)

!- WRFoutv2 forcing:
    call registerreadforcelev(trim(LDT_wrfoutv2Id)//char(0),&
         read_WRFoutv2_elev)

!- WRF-Alaska forcing:
    call registerreadforcelev(trim(LDT_WRFakId)//char(0),&
         read_WRFAK_elev)

!- GEOS5 forcing:
!    call registerreadforcelev(trim(LDT_geos5Id)//char(0),&
!         read_geos5_elev)

  end subroutine LDT_forcingparams_plugin

!BOP
! !ROUTINE: LDT_glacier_plugin
!  \label{LDT_Glacier_plugin}
!
! !DESCRIPTION:
!
! !INTERFACE:
  subroutine LDT_glacier_plugin
!EOP
    external read_GLIMS_glaciermask
    external read_GLIMS_glacierfraction


!   In the LDT code, the above calls are typically invoked in the
!   following manner.
!   \begin{verbatim}
!    call readglacierfrac(ldt%domain,ldtglacierfracsrc)

    call registerreadglaciermask(trim(LDT_GLIMSId)//char(0),&
         read_GLIMS_glaciermask)
    call registerreadglacierfrac(trim(LDT_GLIMSId)//char(0),&
         read_GLIMS_glacierfraction)

  end subroutine LDT_glacier_plugin


end module LDT_param_pluginMod
