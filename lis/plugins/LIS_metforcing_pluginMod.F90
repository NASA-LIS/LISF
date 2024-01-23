!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_metforcing_pluginMod
!BOP
!
! !MODULE: LIS_metforcing_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   incorporating a new met forcing scheme. The user defined functions
!   are incorporated into the appropriate registry to be later invoked
!   through generic calls. 
!   
! !REVISION HISTORY: 
!  11 Dec 2003   Sujay Kumar:  Initial Specification
!  11 Oct 2014   K. Arsenault: Reorganzed for syncing with LDT
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_metforcing_plugin
  
contains
!BOP
! !ROUTINE: LIS_metforcing_plugin
!  \label{LIS_metforcing_plugin}
! 
! !INTERFACE:
subroutine LIS_metforcing_plugin

! !DESCRIPTION:
!
! This is a plugin point for introducing a new met forcing scheme.
! The interface mandates that the following interfaces be implemented
! and registered for each met forcing scheme. 
!
!  \begin{description}
!  \item[retrieval of forcing data]      
!      Routines to retrieve forcing data and to interpolate them. 
!      (to be registered using {\tt registerget} and later invoked through
!       {\tt retrieveforcing} method)
!  \item[definition of native domain]
!      Routines to define the native domain 
!      (to be registered using registerinitmetforc and later
!       invoked through {\tt initmetforc} method)
!  \item[temporal interpolation] 
!      Interpolate forcing data temporally. 
!      (to be registered using {\tt registertimeinterp} and later 
!       invoked through {\tt timeinterp} method)
!  \end{description}
!
!  The index used in the register calls should be used 
!  to select the appropriate met forcing scheme. For example, 
!  assume that the GDAS forcing scheme is incorporated in the registry
!  by the following calls. 
!   \begin{verbatim}
!     call registerinitmetforc(1,init_GDAS)
!     call registerget(1,get_gdas)             
!     call registertimeinterp(1,timeinterp_gdas)
!     call registerforcingfinal(1,finalize_gdas) 
!    \end{verbatim}
!   The index used here to register these functions is 1. To invoke
!   these methods, the corresponding calls will be: 
!   \begin{verbatim}
!     call initmetforc(1)  - calls init_gdas
!     call retrieveforcing(1)  - calls get_gdas
!     call timeinterp(1)  - calls timeinterp_gdas
!     call forcingfinalize(1) - calls finalize_gdas
!   \end{verbatim}
!   
!   In the LIS code, the above calls are typically invoked in the 
!   following manner. 
!   \begin{verbatim}
!     call initmetforc(lis%force)  
!     call retrieveforcing(lis%force,n)
!     call timeinterp(lis%force,n)  
!     call forcingfinalize(lis%force)
!   \end{verbatim}
!   where $lis\%force$ is set through the configuration utility, enabling
!   the user to select any of the met forcing schemes, at runtime. 
!EOP
   use LIS_pluginIndices

#if ( defined MF_MET_TEMPLATE )
   use metForcTemplate_forcingMod
#endif

#if ( defined MF_LDT_GENERATED )
   use metForcGenerated_forcingMod
#endif

#if ( defined MF_CLIMATOLOGY )
   use climatology_forcingMod
#endif

#if ( defined MF_GENENSFCST )
   use genEnsFcst_forcingMod
#endif

#if ( defined MF_PPTENSFCST )
   use pptEnsFcst_forcingMod
#endif

#if ( defined MF_GDAS )
   use gdas_forcingMod
#endif

#if ( defined MF_ECMWF )
   use ecmwf_forcingMod
#endif

#if ( defined MF_PRINCETON )
   use princeton_forcingMod
#endif

#if ( defined MF_GLDAS )
   use gldas_forcingMod
#endif

#if ( defined MF_GFS )
   use gfs_forcingMod
#endif

#if ( defined MF_MERRA2 )
   use merra2_forcingMod
#endif

#if ( defined MF_GEOS_IT )
   use geosit_forcingMod
#endif

#if ( defined MF_ERA5 )
   use era5_forcingMod
#endif

#if ( defined MF_GSWP1 )
   use gswp1_forcingMod
#endif

#if ( defined MF_GSWP2 )
   use gswp2_forcingMod
#endif

#if ( defined MF_AGRMET )
   use AGRMET_forcingMod
#endif

#if ( defined MF_AGRMET_RADIATION_POLAR_STEREOGRAPHIC )
   use agrradps_forcingMod
#endif

#if ( defined MF_GEOS5_FORECAST )
   use geos5fcst_forcingMod
#endif

#if ( defined MF_GEFS_FORECAST )
   use gefs_forcingMod
#endif

#if ( defined MF_TRMM_3B42RT )
   use TRMM3B42RT_forcingMod
#endif

#if ( defined MF_TRMM_3B42RTV7 )
   use TRMM3B42RTV7_forcingMod
#endif

#if ( defined MF_TRMM_3B42V6 )
   use TRMM3B42V6_forcingMod
#endif

#if ( defined MF_TRMM_3B42V7 )
   use TRMM3B42V7_forcingMod 
#endif

#if ( defined MF_CPC_CMORPH )
   use cmorph_forcingMod
#endif

#if ( defined MF_GPM_IMERG )
   use imerg_forcingMod
#endif

#if ( defined MF_CPC_STAGEII )
   use stg2_forcingMod
#endif

#if ( defined MF_CPC_STAGEIV )
   use stg4_forcingMod
#endif

#if ( defined MF_CMAP )
   use cmap_forcingMod
#endif

#if ( defined MF_NLDAS2 )
   use nldas2_forcingMod
#endif

#if ( defined MF_NARR )
   use narr_forcingMod
#endif

#if ( defined MF_NAM242 )
   use nam242_forcingMod
#endif

#if ( defined MF_VIC_PROCESSED_FORCING )
   use vic_forcingMod
#endif

#if ( defined MF_RFE2_DAILY )
   use RFE2Daily_forcingMod
#endif

#if ( defined MF_RFE2_GDAS_BIAS_CORRECTED )
   use RFE2gdas_forcingMod
#endif

#if ( defined MF_PET_USGS )
   use petusgs_forcingMod
#endif

#if ( defined MF_CHIRPS2 )
   use chirps2_forcingMod
#endif

#if ( defined MF_SCAN )
   use scan_forcingMod
#endif

#if ( defined MF_BONDVILLE )
   use Bondville_forcingMod
#endif

#if ( defined MF_PLUMBER2)
   use plumber2_forcingMod
#endif

#if ( defined MF_LOOBOS )
   use Loobos_forcingMod
#endif

#if ( defined MF_SNOTEL )
   use snotel_forcingMod
#endif

#if ( defined MF_PALS_STATION_FORCING )
   use PALSmetdata_forcingMod
#endif

#if ( defined MF_WRFOUT )
   use WRFout_forcingMod
#endif

#if ( defined MF_COAMPSOUT )
   use COAMPSout_forcingMod
#endif

#if ( defined MF_WRFOUTV2 )
   use WRFoutv2_forcingMod
#endif

#if ( defined MF_WRF_AK )
   use WRF_AKdom_forcingMod
#endif

#if ( defined MF_GDAS_T1534 )
   use gdasT1534_forcingMod
#endif

#if ( defined MF_AWAP)
   use AWAP_forcingMod
#endif

#if ( defined MF_AWRAL)
   use AWRAL_forcingMod
#endif

#if ( defined MF_HIMAT_GMU)
   use HiMATGMU_forcingMod
#endif

#if ( defined MF_MRMS )
   use mrms_grib_forcingMod
#endif

#if ( defined MF_GDDP )
   use gddp_forcingMod
#endif   

#if ( defined MF_GALWEM_FORECAST )
   use galwem_forcingMod
#endif

#if ( defined MF_GALWEM_GE_FORECAST )
   use galwemge_forcingMod
#endif

#if ( defined MF_MOGREPS_G_FORECAST )
   use mogrepsg_forcingMod
#endif

#if ( defined MF_MET_TEMPLATE )
   external get_metForcTemplate
   external timeinterp_metForcTemplate
   external finalize_metForcTemplate
   external reset_metForcTemplate
#endif

#if ( defined MF_LDT_GENERATED )
   external get_metForcGenerated
   external timeinterp_metForcGenerated
   external finalize_metForcGenerated
   external reset_metForcGenerated
#endif

#if ( defined MF_CLIMATOLOGY )
   external get_climatology
   external timeinterp_climatology
   external finalize_climatology
   external reset_climatology
#endif

#if ( defined MF_GENENSFCST )
   external get_genEnsFcst
   external timeinterp_genEnsFcst
   external finalize_genEnsFcst
   external reset_genEnsFcst
#endif

#if ( defined MF_PPTENSFCST )
   external get_pptEnsFcst
   external timeinterp_pptEnsFcst
   external finalize_pptEnsFcst
   external reset_pptEnsFcst
#endif


#if ( defined MF_GDAS )
   external get_gdas
   external timeinterp_gdas
   external reset_gdas
   external finalize_gdas
#endif

#if ( defined MF_ECMWF )
   external get_ecmwf
   external timeinterp_ecmwf
   external finalize_ecmwf
   external reset_ecmwf
#endif

#if ( defined MF_PRINCETON )
   external get_princeton
   external timeinterp_princeton
   external finalize_princeton
   external reset_princeton
#endif

#if ( defined MF_GLDAS )
   external get_gldas
   external timeinterp_gldas
   external finalize_gldas
#endif

#if ( defined MF_GFS )
   external get_gfs
   external timeinterp_gfs
   external finalize_gfs
#endif

#if ( defined MF_MERRA2 )
   external get_merra2
   external timeinterp_merra2
   external finalize_merra2
   external reset_merra2
#endif

#if ( defined MF_GEOS_IT )
   external get_geosit
   external timeinterp_geosit
   external finalize_geosit
   external reset_geosit
#endif

#if ( defined MF_ERA5 )
   external get_era5
   external timeinterp_era5
   external finalize_era5
   external reset_era5
#endif

#if ( defined MF_GSWP1 )
   external get_gswp1
   external timeinterp_gswp1
   external finalize_gswp1
#endif

#if ( defined MF_GSWP2 )
   external get_gswp2
   external timeinterp_gswp2
   external finalize_gswp2
#endif

#if ( defined MF_AGRMET )
   external get_agrmet
   external timeinterp_agrmet
   external reset_agrmet
   external finalize_agrmet
#endif

#if ( defined MF_AGRMET_RADIATION_POLAR_STEREOGRAPHIC )
   external get_agrradps
   external timeinterp_agrradps
   external finalize_agrradps
#endif

#if ( defined MF_GEOS5_FORECAST )
   external get_geos5fcst
   external timeinterp_geos5fcst
   external finalize_geos5fcst
   external reset_geos5fcst
#endif

#if ( defined MF_GEFS_FORECAST )
   external get_gefs
   external timeinterp_gefs
   external finalize_gefs
   external reset_gefs
#endif

#if ( defined MF_TRMM_3B42RT )
   external get_TRMM3B42RT
   external timeinterp_TRMM3B42RT
   external finalize_TRMM3B42RT
#endif

#if ( defined MF_TRMM_3B42RTV7 )
   external get_TRMM3B42RTV7
   external timeinterp_TRMM3B42RTV7
   external finalize_TRMM3B42RTV7
#endif

#if ( defined MF_TRMM_3B42V6 )
   external get_TRMM3B42V6
   external timeinterp_TRMM3B42V6
   external finalize_TRMM3B42V6
#endif

#if ( defined MF_TRMM_3B42V7 )
   external get_TRMM3B42V7
   external timeinterp_TRMM3B42V7
   external finalize_TRMM3B42V7
#endif

#if ( defined MF_CPC_CMORPH )
   external get_cmorph
   external timeinterp_cmorph
   external finalize_cmorph
#endif

#if ( defined MF_GPM_IMERG )
   external get_imerg
   external timeinterp_imerg
   external finalize_imerg
   external reset_imerg
#endif

#if ( defined MF_CPC_STAGEII )
   external get_stg2
   external timeinterp_stg2
   external finalize_stg2
#endif

#if ( defined MF_CPC_STAGEIV )
   external get_stg4
   external timeinterp_stg4
   external finalize_stg4
   external reset_stg4
#endif

#if ( defined MF_CMAP )
   external get_cmap
   external timeinterp_cmap
   external finalize_cmap
#endif

#if ( defined MF_NLDAS2 )
   external get_nldas2
   external timeinterp_nldas2
   external finalize_nldas2
   external reset_nldas2
#endif

#if ( defined MF_NARR )
   external get_narr
   external timeinterp_narr
   external finalize_narr
#endif

#if ( defined MF_NAM242 )
   external get_nam242
   external timeinterp_nam242
   external finalize_nam242
#endif

#if ( defined MF_VIC_PROCESSED_FORCING )
   external getvicforcing
   external time_interp_vicforcing
   external vicforcing_finalize
#endif

#if ( defined MF_RFE2_DAILY )
   external get_RFE2Daily
   external timeinterp_RFE2Daily
   external finalize_RFE2Daily
#endif

#if ( defined MF_RFE2_GDAS_BIAS_CORRECTED )
   external get_RFE2gdas
   external timeinterp_RFE2gdas
   external finalize_RFE2gdas
   external reset_RFE2gdas
#endif

#if ( defined MF_PET_USGS )
   external get_petusgs
   external timeinterp_petusgs
   external finalize_petusgs
#endif

#if ( defined MF_CHIRPS2 )
   external get_chirps2
   external timeinterp_chirps2
   external finalize_chirps2
   external reset_chirps2
#endif

#if ( defined MF_SCAN )
   external get_scan
   external timeinterp_scan
   external finalize_scan
#endif

#if ( defined MF_BONDVILLE )
   external get_Bondville
   external timeinterp_Bondville
   external finalize_Bondville
#endif

#if ( defined MF_PLUMBER2)
   external get_plumber2
   external timeinterp_plumber2
   external finalize_plumber2
   external reset_plumber2
#endif


#if ( defined MF_LOOBOS )
   external get_Loobos
   external timeinterp_Loobos
   external finalize_Loobos
#endif

#if ( defined MF_SNOTEL )
   external get_snotel
   external timeinterp_snotel
   external finalize_snotel
#endif

#if ( defined MF_PALS_STATION_FORCING )
   external get_PALSmetdata
   external timeinterp_PALSmetdata
   external finalize_PALSmetdata
   external reset_PALSmetdata
#endif

#if ( defined MF_WRFOUT )
   external get_WRFout
   external timeinterp_WRFout
   external reset_WRFout
   external finalize_WRFout
#endif

#if ( defined MF_COAMPSOUT )
   external get_COAMPSout
   external timeinterp_COAMPSout
   external reset_COAMPSout
   external finalize_COAMPSout
#endif

#if ( defined MF_WRFOUTV2 )
   external get_WRFoutv2
   external timeinterp_WRFoutv2
   external reset_WRFoutv2
   external finalize_WRFoutv2
#endif

#if ( defined MF_WRF_AK )
   external get_WRF_AKdom
   external timeinterp_WRF_AKdom
   external reset_WRF_AKdom
   external finalize_WRF_AKdom
#endif

#if ( defined MF_GDAS_T1534 )
   external get_gdasT1534
   external timeinterp_gdasT1534
   external reset_gdasT1534
   external finalize_gdasT1534
#endif

#if ( defined MF_AWAP )
   external get_AWAP
   external timeinterp_AWAP
   external reset_AWAP
   external finalize_AWAP
#endif

#if ( defined MF_AWRAL )
   external get_AWRAL
   external timeinterp_AWRAL
   external reset_AWRAL
   external finalize_AWRAL
#endif

#if ( defined MF_HIMAT_GMU )
   external get_HiMATGMU
   external timeinterp_HiMATGMU
   external reset_HiMATGMU
   external finalize_HiMATGMU
#endif

#if ( defined MF_MRMS )
   external get_mrms_grib
   external timeinterp_mrms_grib
   external finalize_mrms_grib
   external reset_mrms_grib
#endif

#if ( defined MF_GDDP )
   external get_gddp
   external timeinterp_gddp
   external finalize_gddp
   external reset_gddp
#endif

#if ( defined MF_GALWEM_FORECAST )
   external get_galwem
   external timeinterp_galwem
   external finalize_galwem
   external reset_galwem
#endif
 
#if ( defined MF_GALWEM_GE_FORECAST )
   external get_galwemge
   external timeinterp_galwemge
   external finalize_galwemge
   external reset_galwemge
#endif

#if ( defined MF_MOGREPS_G_FORECAST )
   external get_mogrepsg
   external timeinterp_mogrepsg
   external finalize_mogrepsg
   external reset_mogrepsg
#endif

   external :: registerinitmetforc
   external :: registerretrievemetforc
   external :: registertimeinterpmetforc
   external :: registerfinalmetforc
   external :: registerresetmetforc

#if ( defined MF_MET_TEMPLATE )
! - Meteorological Forcing Template:
   call registerinitmetforc(trim(LIS_metForcTemplateId)//char(0), &
                                 init_metForctemplate)
   call registerretrievemetforc(trim(LIS_metForcTemplateId)//char(0), &
                                get_MetForctemplate)
   call registertimeinterpmetforc(trim(LIS_metForcTemplateId)//char(0), &
                                  timeinterp_MetForctemplate)
   call registerfinalmetforc(trim(LIS_metForcTemplateId)//char(0), &
                             finalize_metForctemplate)
   call registerresetmetforc(trim(LIS_metForcTemplateId)//char(0), &
                             reset_metForctemplate)
#endif

#if ( defined MF_LDT_GENERATED )
! - Generated Meteorological Forcing (e.g., in LDT):
   call registerinitmetforc(trim(LIS_generatedForcId)//char(0), &
                            init_metForcGenerated)
   call registerretrievemetforc(trim(LIS_generatedForcId)//char(0), &
                                get_metForcGenerated)
   call registertimeinterpmetforc(trim(LIS_generatedForcId)//char(0), &
                                  timeinterp_metForcGenerated)
   call registerresetmetforc(trim(LIS_generatedForcId)//char(0), &
                             reset_metForcGenerated)
   call registerfinalmetforc(trim(LIS_generatedForcId)//char(0), &
                             finalize_metForcGenerated)
#endif

#if ( defined MF_CLIMATOLOGY )
! - Climatological Meteorological Forcing (e.g.,generated in LDT):
   call registerinitmetforc(trim(LIS_climstdId)//char(0),init_climatology)
   call registerretrievemetforc(trim(LIS_climstdId)//char(0),get_climatology)
   call registertimeinterpmetforc(trim(LIS_climstdId)//char(0),&
        timeinterp_climatology)
   call registerresetmetforc(trim(LIS_climstdId)//char(0),&
        reset_climatology)
   call registerfinalmetforc(trim(LIS_climstdId)//char(0),&
        finalize_climatology)
#endif

#if ( defined MF_GENENSFCST )
! - Generic Ensemble Forecast Reader (input files generated by user):
   call registerinitmetforc(trim(LIS_genensfcstId)//char(0),init_genEnsFcst)
   call registerretrievemetforc(trim(LIS_genensfcstId)//char(0),get_genEnsFcst)
   call registertimeinterpmetforc(trim(LIS_genensfcstId)//char(0),&
        timeinterp_genEnsFcst)
   call registerresetmetforc(trim(LIS_genensfcstId)//char(0),&
        reset_genEnsFcst)
   call registerfinalmetforc(trim(LIS_genensfcstId)//char(0),&
        finalize_genEnsFcst)
#endif

#if ( defined MF_PPTENSFCST )
! - Generic Ensemble Forecast Reader (input files generated by user):
   call registerinitmetforc(trim(LIS_pptensfcstId)//char(0),init_pptEnsFcst)
   call registerretrievemetforc(trim(LIS_pptensfcstId)//char(0),get_pptEnsFcst)
   call registertimeinterpmetforc(trim(LIS_pptensfcstId)//char(0),&
        timeinterp_pptEnsFcst)
   call registerresetmetforc(trim(LIS_pptensfcstId)//char(0),&
        reset_pptEnsFcst)
   call registerfinalmetforc(trim(LIS_pptensfcstId)//char(0),&
        finalize_pptEnsFcst)
#endif


#if ( defined MF_GDAS )
! - GDAS (GSFC) Forcing:
   call registerinitmetforc(trim(LIS_gdasId)//char(0),init_gdas)
   call registerretrievemetforc(trim(LIS_gdasId)//char(0),get_gdas)
   call registertimeinterpmetforc(trim(LIS_gdasId)//char(0),timeinterp_gdas)
   call registerresetmetforc(trim(LIS_gdasId)//char(0),reset_gdas)
   call registerfinalmetforc(trim(LIS_gdasId)//char(0),finalize_gdas)
#endif

#if ( defined MF_ECMWF )
! - ECMWF Forcing:
   call registerinitmetforc(trim(LIS_ecmwfId)//char(0),init_ECMWF)
   call registerretrievemetforc(trim(LIS_ecmwfId)//char(0),get_ecmwf)
   call registertimeinterpmetforc(trim(LIS_ecmwfId)//char(0),timeinterp_ecmwf)
   call registerfinalmetforc(trim(LIS_ecmwfId)//char(0),finalize_ecmwf)
   call registerresetmetforc(trim(LIS_ecmwfId)//char(0),reset_ecmwf)
#endif

#if ( defined MF_PRINCETON )
! - PRINCETON Reanalysis Forcing:
   call registerinitmetforc(trim(LIS_princetonId)//char(0),init_PRINCETON)
   call registerretrievemetforc(trim(LIS_princetonId)//char(0),get_princeton)
   call registertimeinterpmetforc(trim(LIS_princetonId)//char(0), &
                                  timeinterp_princeton)
   call registerfinalmetforc(trim(LIS_princetonId)//char(0),finalize_princeton)
   call registerresetmetforc(trim(LIS_princetonId)//char(0),reset_princeton)
#endif

#if ( defined MF_GLDAS )
! - GLDAS Reanalysis Forcing:
   call registerinitmetforc(trim(LIS_gldasId)//char(0),init_GLDAS)
   call registerretrievemetforc(trim(LIS_gldasId)//char(0),get_gldas)
   call registertimeinterpmetforc(trim(LIS_gldasId)//char(0), &
                                  timeinterp_gldas)
   call registerfinalmetforc(trim(LIS_gldasId)//char(0),finalize_gldas)
#endif

#if ( defined MF_GFS )
! - GFS Forecast Forcing:
   call registerinitmetforc(trim(LIS_gfsId)//char(0),init_GFS)
   call registerretrievemetforc(trim(LIS_gfsId)//char(0),get_gfs)
   call registertimeinterpmetforc(trim(LIS_gfsId)//char(0),timeinterp_gfs)
   call registerfinalmetforc(trim(LIS_gfsId)//char(0),finalize_gfs)
#endif

#if ( defined MF_MERRA2 )
! - MERRA2 Reanalysis Forcing:
   call registerinitmetforc(trim(LIS_merra2Id)//char(0),init_MERRA2)
   call registerretrievemetforc(trim(LIS_merra2Id)//char(0),get_merra2)
   call registertimeinterpmetforc(trim(LIS_merra2Id)//char(0), &
                                  timeinterp_merra2)
   call registerresetmetforc(trim(LIS_merra2Id)//char(0),reset_merra2)
   call registerfinalmetforc(trim(LIS_merra2Id)//char(0),finalize_merra2)
#endif

#if ( defined MF_GEOS_IT )
! - GEOS-IT Forcing:
   call registerinitmetforc(trim(LIS_geositId)//char(0),init_geosit)
   call registerretrievemetforc(trim(LIS_geositId)//char(0),get_geosit)
   call registertimeinterpmetforc(trim(LIS_geositId)//char(0), &
                                  timeinterp_geosit)
   call registerresetmetforc(trim(LIS_geositId)//char(0),reset_geosit)
   call registerfinalmetforc(trim(LIS_geositId)//char(0),finalize_geosit)
#endif

#if ( defined MF_ERA5)
! - ERA5 Reanalysis Forcing:
   call registerinitmetforc(trim(LIS_era5Id)//char(0),init_ERA5)
   call registerretrievemetforc(trim(LIS_era5Id)//char(0),get_era5)
   call registertimeinterpmetforc(trim(LIS_era5Id)//char(0), &
                                  timeinterp_era5)
   call registerresetmetforc(trim(LIS_era5Id)//char(0),reset_era5)
   call registerfinalmetforc(trim(LIS_era5Id)//char(0),finalize_era5)
#endif

#if ( defined MF_GSWP1 )
! - GWSP1 Forcing:
   call registerinitmetforc(trim(LIS_gswp1Id)//char(0),init_GSWP1)
   call registerretrievemetforc(trim(LIS_gswp1Id)//char(0),get_gswp1)
   call registertimeinterpmetforc(trim(LIS_gswp1Id)//char(0), &
                                  timeinterp_gswp1)
   call registerfinalmetforc(trim(LIS_gswp1Id)//char(0),finalize_gswp1)
#endif

#if ( defined MF_GSWP2 )
! - GSWP2 Forcing:
   call registerinitmetforc(trim(LIS_gswp2Id)//char(0),init_GSWP2)
   call registerretrievemetforc(trim(LIS_gswp2Id)//char(0),get_gswp2)
   call registertimeinterpmetforc(trim(LIS_gswp2Id)//char(0), &
                                  timeinterp_gswp2)
   call registerfinalmetforc(trim(LIS_gswp2Id)//char(0),finalize_gswp2)
#endif

#if ( defined MF_AGRMET )
! - AGRMET Forcing:
   call registerinitmetforc(trim(LIS_agrmetId)//char(0),init_AGRMET)
   call registerretrievemetforc(trim(LIS_agrmetId)//char(0),get_agrmet)
   call registertimeinterpmetforc(trim(LIS_agrmetId)//char(0), &
                                  timeinterp_agrmet)
   call registerresetmetforc(trim(LIS_agrmetId)//char(0),reset_agrmet)
   call registerfinalmetforc(trim(LIS_agrmetId)//char(0),finalize_agrmet)
#endif

#if ( defined MF_AGRMET_RADIATION_POLAR_STEREOGRAPHIC )
! - AGRMET Polar Stereographic Radiation:
   call registerinitmetforc(trim(LIS_agrradpsId)//char(0),init_AGRRADPS)
   call registerretrievemetforc(trim(LIS_agrradpsId)//char(0),get_agrradps)
   call registertimeinterpmetforc(trim(LIS_agrradpsId)//char(0), &
                                  timeinterp_agrradps)
   call registerfinalmetforc(trim(LIS_agrradpsId)//char(0),finalize_agrradps)
#endif

#if ( defined MF_GEOS5_FORECAST )
! - GEOSv5 forecast forcing:
   call registerinitmetforc(trim(LIS_geos5fcstId)//char(0),init_GEOS5FCST)
   call registerretrievemetforc(trim(LIS_geos5fcstId)//char(0),get_geos5fcst)
   call registertimeinterpmetforc(trim(LIS_geos5fcstId)//char(0), &
                                  timeinterp_geos5fcst)
   call registerfinalmetforc(trim(LIS_geos5fcstId)//char(0),finalize_geos5fcst)
   call registerresetmetforc(trim(LIS_geos5fcstId)//char(0),reset_geos5fcst)
#endif

#if ( defined MF_GEFS_FORECAST )
! - GEFS forecast forcing:
   call registerinitmetforc(trim(LIS_gefsId)//char(0),init_GEFS)
   call registerretrievemetforc(trim(LIS_gefsId)//char(0),get_gefs)
   call registertimeinterpmetforc(trim(LIS_gefsID)//char(0), &
                                  timeinterp_gefs)
   call registerfinalmetforc(trim(LIS_gefsId)//char(0),finalize_gefs)
   call registerresetmetforc(trim(LIS_gefsId)//char(0),reset_gefs)
#endif

#if ( defined MF_TRMM_3B42RT )
! - 3B42RT TRMM (Near-) Real-time (added by Yudong)
   call registerinitmetforc(trim(LIS_TRMM3B42RTId)//char(0),init_TRMM3B42RT)
   call registerretrievemetforc(trim(LIS_TRMM3B42RTId)//char(0),get_TRMM3B42RT)
   call registertimeinterpmetforc(trim(LIS_TRMM3B42RTId)//char(0), &
                                  timeinterp_TRMM3B42RT)
   call registerfinalmetforc(trim(LIS_TRMM3B42RTId)//char(0), &
                             finalize_TRMM3B42RT)
#endif

#if ( defined MF_TRMM_3B42RTV7 )
! - 3B42RT-V7 TRMM (Near-) Real-time (K. Arsenault)
   call registerinitmetforc(trim(LIS_TRMM3B42RTV7Id)//char(0),init_TRMM3B42RTV7)
   call registerretrievemetforc(trim(LIS_TRMM3B42RTV7Id)//char(0), &
                                get_TRMM3B42RTV7)
   call registertimeinterpmetforc(trim(LIS_TRMM3B42RTV7Id)//char(0), &
                                  timeinterp_TRMM3B42RTV7)
   call registerfinalmetforc(trim(LIS_TRMM3B42RTV7Id)//char(0), &
                             finalize_TRMM3B42RTV7)
#endif

#if ( defined MF_TRMM_3B42V6 )
! - 3B42V6 TRMM version 6 (added by Yudong)
   call registerinitmetforc(trim(LIS_TRMM3B42V6Id)//char(0),init_TRMM3B42V6)
   call registerretrievemetforc(trim(LIS_TRMM3B42V6Id)//char(0),get_TRMM3B42V6)
   call registertimeinterpmetforc(trim(LIS_TRMM3B42V6Id)//char(0), &
                                  timeinterp_TRMM3B42V6)
   call registerfinalmetforc(trim(LIS_TRMM3B42V6Id)//char(0), &
                             finalize_TRMM3B42V6)
#endif

#if ( defined MF_TRMM_3B42V7 )
! - 3B42V7 TRMM version 7 (added by Soni)
   call registerinitmetforc(trim(LIS_TRMM3B42V7Id)//char(0),init_TRMM3B42V7)
   call registerretrievemetforc(trim(LIS_TRMM3B42V7Id)//char(0),get_TRMM3B42V7)
   call registertimeinterpmetforc(trim(LIS_TRMM3B42V7Id)//char(0), &
                                  timeinterp_TRMM3B42V7)
   call registerfinalmetforc(trim(LIS_TRMM3B42V7Id)//char(0), &
                             finalize_TRMM3B42V7)
#endif

#if ( defined MF_CPC_CMORPH )
! - CMORPH 30min, 8KM, added by Yudong
   call registerinitmetforc(trim(LIS_cmorphId)//char(0),init_CMORPH)
   call registerretrievemetforc(trim(LIS_cmorphId)//char(0),get_cmorph)
   call registertimeinterpmetforc(trim(LIS_cmorphId)//char(0), &
                                  timeinterp_cmorph)
   call registerfinalmetforc(trim(LIS_cmorphId)//char(0),finalize_cmorph)
#endif

#if ( defined MF_GPM_IMERG )
! - IMERG 30min, 0.1-deg, added by J.Case (3/5/2015)
   call registerinitmetforc(trim(LIS_imergId)//char(0),init_IMERG)
   call registerretrievemetforc(trim(LIS_imergId)//char(0),get_imerg)
   call registertimeinterpmetforc(trim(LIS_imergId)//char(0), &
                                  timeinterp_imerg)
   call registerfinalmetforc(trim(LIS_imergId)//char(0),finalize_imerg)
   call registerresetmetforc(trim(LIS_imergId)//char(0),reset_imerg)
#endif

#if ( defined MF_CPC_STAGEII )
! - STAGE 2 HRAP, 4KM, added by K. Arsenault
   call registerinitmetforc(trim(LIS_stg2Id)//char(0),init_STG2)
   call registerretrievemetforc(trim(LIS_stg2Id)//char(0),get_stg2)
   call registertimeinterpmetforc(trim(LIS_stg2Id)//char(0), &
                                  timeinterp_stg2)
   call registerfinalmetforc(trim(LIS_stg2Id)//char(0),finalize_stg2)
#endif

#if ( defined MF_CPC_STAGEIV )
! - STAGE 4 HRAP, 4KM, added by K. Arsenault
   call registerinitmetforc(trim(LIS_stg4Id)//char(0),init_STG4)
   call registerretrievemetforc(trim(LIS_stg4Id)//char(0),get_stg4)
   call registertimeinterpmetforc(trim(LIS_stg4Id)//char(0), &
                                  timeinterp_stg4)
   call registerfinalmetforc(trim(LIS_stg4Id)//char(0),finalize_stg4)
   call registerresetmetforc(trim(LIS_stg4Id)//char(0),reset_stg4)
#endif

#if ( defined MF_CMAP )
! - CMAP/GDAS (analysis) Forcing:
   call registerinitmetforc(trim(LIS_cmapId)//char(0),init_CMAP)
   call registerretrievemetforc(trim(LIS_cmapId)//char(0),get_cmap)
   call registertimeinterpmetforc(trim(LIS_cmapId)//char(0), &
                                  timeinterp_cmap)
   call registerfinalmetforc(trim(LIS_cmapId)//char(0),finalize_cmap)
#endif

#if ( defined MF_NLDAS2 )
! - NLDAS2 Forcing:
   call registerinitmetforc(trim(LIS_nldas2Id)//char(0),init_NLDAS2)
   call registerretrievemetforc(trim(LIS_nldas2Id)//char(0),get_nldas2)
   call registertimeinterpmetforc(trim(LIS_nldas2Id)//char(0), &
                                  timeinterp_nldas2)
   call registerfinalmetforc(trim(LIS_nldas2Id)//char(0),finalize_nldas2)
   call registerresetmetforc(trim(LIS_nldas2Id)//char(0),reset_nldas2)
#endif

#if ( defined MF_NARR )
! - NARR profile data for CRTM 
   call registerinitmetforc(trim(LIS_narrId)//char(0),init_NARR)
   call registerretrievemetforc(trim(LIS_narrId)//char(0),get_narr)
   call registertimeinterpmetforc(trim(LIS_narrId)//char(0), &
                                  timeinterp_narr)
   call registerfinalmetforc(trim(LIS_narrId)//char(0),finalize_narr)
#endif

#if ( defined MF_NAM242 )
! - NAM242 forcing
   call registerinitmetforc(trim(LIS_nam242Id)//char(0),init_nam242)
   call registerretrievemetforc(trim(LIS_nam242Id)//char(0),get_nam242)
   call registertimeinterpmetforc(trim(LIS_nam242Id)//char(0), &
                                  timeinterp_nam242)
   call registerfinalmetforc(trim(LIS_nam242Id)//char(0),finalize_nam242)
#endif

#if ( defined MF_VIC_PROCESSED_FORCING )
! - VIC processed data
   call registerinitmetforc(trim(LIS_vicforcingId)//char(0), &
                            defineNativevicforcing)
   call registerretrievemetforc(trim(LIS_vicforcingId)//char(0),getvicforcing)
   call registertimeinterpmetforc(trim(LIS_vicforcingId)//char(0), &
                                  time_interp_vicforcing)
   call registerfinalmetforc(trim(LIS_vicforcingId)//char(0), &
                             vicforcing_finalize)
#endif

#if ( defined MF_RFE2_DAILY )
! - FEWSNET CPC RFE2.0 Data
   call registerinitmetforc(trim(LIS_RFE2DailyId)//char(0),init_RFE2Daily)
   call registerretrievemetforc(trim(LIS_RFE2DailyId)//char(0),get_RFE2Daily)
   call registertimeinterpmetforc(trim(LIS_RFE2DailyId)//char(0), &
                                  timeinterp_RFE2Daily)
   call registerfinalmetforc(trim(LIS_RFE2DailyId)//char(0),finalize_RFE2Daily)
#endif

#if ( defined MF_RFE2_GDAS_BIAS_CORRECTED )
! - FEWSNET RFE2.0-GDAS (Temporally downscaled) Data
   call registerinitmetforc(trim(LIS_RFE2gdasId)//char(0),init_RFE2gdas)
   call registerretrievemetforc(trim(LIS_RFE2gdasId)//char(0),get_RFE2gdas)
   call registertimeinterpmetforc(trim(LIS_RFE2gdasId)//char(0), &
                                  timeinterp_RFE2gdas)
   call registerfinalmetforc(trim(LIS_RFE2gdasId)//char(0),finalize_RFE2gdas)
   call registerresetmetforc(trim(LIS_RFE2gdasId)//char(0),reset_RFE2gdas)
#endif

#if ( defined MF_PET_USGS )
! - FEWSNET USGS PET Data
   call registerinitmetforc(trim(LIS_USGSPETforcId)//char(0),init_petusgs)
   call registerretrievemetforc(trim(LIS_USGSPETforcId)//char(0),get_petusgs)
   call registertimeinterpmetforc(trim(LIS_USGSPETforcId)//char(0), &
                                  timeinterp_petusgs)
   call registerfinalmetforc(trim(LIS_USGSPETforcId)//char(0),finalize_petusgs)
#endif

#if ( defined MF_CHIRPS2 )
! - CHIRPS 2.0 Precipitation Data 
   call registerinitmetforc(trim(LIS_chirps2Id)//char(0),init_chirps2)
   call registerretrievemetforc(trim(LIS_chirps2Id)//char(0),get_chirps2)
   call registertimeinterpmetforc(trim(LIS_chirps2Id)//char(0), &
                                  timeinterp_chirps2)
   call registerfinalmetforc(trim(LIS_chirps2Id)//char(0),finalize_chirps2)
   call registerresetmetforc(trim(LIS_chirps2Id)//char(0),reset_chirps2)
#endif

#if ( defined MF_SCAN )
! - SCAN Forcing:
   call registerinitmetforc(trim(LIS_scanId)//char(0),init_SCAN)
   call registerretrievemetforc(trim(LIS_scanId)//char(0), get_scan)
   call registertimeinterpmetforc(trim(LIS_scanId)//char(0), &
                                  timeinterp_scan)
   call registerfinalmetforc(trim(LIS_scanId)//char(0),finalize_scan)
#endif

#if (defined MF_PLUMBER2)
   call registerinitmetforc(trim(LIS_plumber2Id)//char(0),init_plumber2)
   call registerretrievemetforc(trim(LIS_plumber2Id)//char(0),get_plumber2)
   call registertimeinterpmetforc(trim(LIS_plumber2Id)//char(0),timeinterp_plumber2)
   call registerresetmetforc(trim(LIS_plumber2Id)//char(0),reset_plumber2)
   call registerfinalmetforc(trim(LIS_plumber2Id)//char(0),finalize_plumber2)
#endif


#if ( defined MF_BONDVILLE )
! - Noah3.1 Bondville test case
   call registerinitmetforc(trim(LIS_BondvilleId)//char(0),init_Bondville)
   call registerretrievemetforc(trim(LIS_BondvilleId)//char(0),get_Bondville)
   call registertimeinterpmetforc(trim(LIS_BondvilleId)//char(0), &
                                  timeinterp_Bondville)
   call registerfinalmetforc(trim(LIS_BondvilleId)//char(0),finalize_Bondville)
#endif

#if ( defined MF_LOOBOS )
! - JULES Loobos test case
   call registerinitmetforc(trim(LIS_LoobosId)//char(0),init_Loobos)
   call registerretrievemetforc(trim(LIS_LoobosId)//char(0),get_Loobos)
   call registertimeinterpmetforc(trim(LIS_LoobosId)//char(0), &
                                  timeinterp_Loobos)
   call registerfinalmetforc(trim(LIS_LoobosId)//char(0),finalize_Loobos)
#endif

#if ( defined MF_SNOTEL )
! - SNOTEL Precipitation Forcing (added by Yuqiong Liu):
   call registerinitmetforc(trim(LIS_snotelId)//char(0),init_SNOTEL)
   call registerretrievemetforc(trim(LIS_snotelId)//char(0),get_snotel)
   call registertimeinterpmetforc(trim(LIS_snotelId)//char(0), &
                                  timeinterp_snotel)
   call registerfinalmetforc(trim(LIS_snotelId)//char(0),finalize_snotel)
#endif

#if ( defined MF_PALS_STATION_FORCING )
! - PALS station data
   call registerinitmetforc(trim(LIS_PALSmetforcId)//char(0),init_PALSMETDATA)
   call registerretrievemetforc(trim(LIS_PALSmetforcId)//char(0), &
                                get_PALSmetdata)
   call registertimeinterpmetforc(trim(LIS_PALSmetforcId)//char(0), &
                                  timeinterp_PALSmetdata)
   call registerfinalmetforc(trim(LIS_PALSmetforcId)//char(0), &
                             finalize_PALSmetdata)
   call registerresetmetforc(trim(LIS_PALSmetforcId)//char(0), &
                             reset_PALSmetdata)
#endif

#if ( defined MF_WRFOUT )
! - WRFout forcing
   call registerinitmetforc(trim(LIS_WRFoutId)//char(0),init_WRFout)
   call registerretrievemetforc(trim(LIS_WRFoutId)//char(0),get_WRFout)
   call registertimeinterpmetforc(trim(LIS_WRFoutId)//char(0), &
                                  timeinterp_WRFout)
   call registerfinalmetforc(trim(LIS_WRFoutId)//char(0),finalize_WRFout)
   call registerresetmetforc(trim(LIS_WRFoutId)//char(0),reset_WRFout)
#endif

#if ( defined MF_COAMPSOUT )
! - COAMPSout forcing
   call registerinitmetforc(trim(LIS_COAMPSoutId)//char(0),init_COAMPSout)
   call registerretrievemetforc(trim(LIS_COAMPSoutId)//char(0),get_COAMPSout)
   call registertimeinterpmetforc(trim(LIS_COAMPSoutId)//char(0), &
                                  timeinterp_COAMPSout)
   call registerfinalmetforc(trim(LIS_COAMPSoutId)//char(0),finalize_COAMPSout)
   call registerresetmetforc(trim(LIS_COAMPSoutId)//char(0),reset_COAMPSout)
#endif

#if ( defined MF_WRFOUTV2 )
! - WRFout-4KM HiRes forcing
   call registerinitmetforc(trim(LIS_WRFoutv2Id)//char(0),init_WRFoutv2)
   call registerretrievemetforc(trim(LIS_WRFoutv2Id)//char(0),get_WRFoutv2)
   call registertimeinterpmetforc(trim(LIS_WRFoutv2Id)//char(0), &
                                  timeinterp_WRFoutv2)
   call registerfinalmetforc(trim(LIS_WRFoutv2Id)//char(0),finalize_WRFoutv2)
   call registerresetmetforc(trim(LIS_WRFoutv2Id)//char(0),reset_WRFoutv2)
#endif

#if ( defined MF_WRF_AK )
! - WRFout-4KM Alaska (AK) domain forcing
   call registerinitmetforc(trim(LIS_WRFakId)//char(0),init_WRF_AKdom)
   call registerretrievemetforc(trim(LIS_WRFakId)//char(0),get_WRF_AKdom)
   call registertimeinterpmetforc(trim(LIS_WRFakId)//char(0), &
                                  timeinterp_WRF_AKdom)
   call registerfinalmetforc(trim(LIS_WRFakId)//char(0),finalize_WRF_AKdom)
   call registerresetmetforc(trim(LIS_WRFakId)//char(0),reset_WRF_AKdom)
#endif

#if ( defined MF_GDAS_T1534 )
! - NCEP Gaussian T1534 forcing
   call registerinitmetforc(trim(LIS_gdasT1534Id)//char(0),init_gdasT1534)
   call registerretrievemetforc(trim(LIS_gdasT1534Id)//char(0),get_gdasT1534)
   call registertimeinterpmetforc(trim(LIS_gdasT1534Id)//char(0), &
                                  timeinterp_gdasT1534)
   call registerresetmetforc(trim(LIS_gdasT1534Id)//char(0),reset_gdasT1534)
   call registerfinalmetforc(trim(LIS_gdasT1534Id)//char(0),finalize_gdasT1534)
#endif

#if ( defined MF_AWAP)
   call registerinitmetforc(trim(LIS_AWAPforcId)//char(0),init_AWAP)
   call registerretrievemetforc(trim(LIS_AWAPforcId)//char(0),get_AWAP)
   call registertimeinterpmetforc(trim(LIS_AWAPforcId)//char(0), &
                                  timeinterp_AWAP)
   call registerresetmetforc(trim(LIS_AWAPforcId)//char(0),reset_AWAP)
   call registerfinalmetforc(trim(LIS_AWAPforcId)//char(0),finalize_AWAP)
#endif

#if ( defined MF_AWRAL)
   call registerinitmetforc(trim(LIS_AWRALforcId)//char(0),init_AWRAL)
   call registerretrievemetforc(trim(LIS_AWRALforcId)//char(0),get_AWRAL)
   call registertimeinterpmetforc(trim(LIS_AWRALforcId)//char(0), &
                                  timeinterp_AWRAL)
   call registerresetmetforc(trim(LIS_AWRALforcId)//char(0),reset_AWRAL)
   call registerfinalmetforc(trim(LIS_AWRALforcId)//char(0),finalize_AWRAL)
#endif

#if ( defined MF_HIMAT_GMU)
   call registerinitmetforc(trim(LIS_HiMATGMUforcId)//char(0),init_HiMATGMU)
   call registerretrievemetforc(trim(LIS_HiMATGMUforcId)//char(0),get_HiMATGMU)
   call registertimeinterpmetforc(trim(LIS_HiMATGMUforcId)//char(0), &
                                  timeinterp_HiMATGMU)
   call registerresetmetforc(trim(LIS_HiMATGMUforcId)//char(0),reset_HiMATGMU)
   call registerfinalmetforc(trim(LIS_HiMATGMUforcId)//char(0),finalize_HiMATGMU)
#endif

#if ( defined MF_MRMS )
! - MRMS operational forcing, 0.1-deg, added by J. Erlingis
   call registerinitmetforc(trim(LIS_mrmsId)//char(0),init_MRMS_grib)
   call registerretrievemetforc(trim(LIS_mrmsId)//char(0),get_mrms_grib)
   call registertimeinterpmetforc(trim(LIS_mrmsId)//char(0), &
                                  timeinterp_mrms_grib)
   call registerfinalmetforc(trim(LIS_mrmsId)//char(0),finalize_mrms_grib)
   call registerresetmetforc(trim(LIS_mrmsId)//char(0),reset_mrms_grib)
#endif
#if ( defined MF_GDDP)
   call registerinitmetforc(trim(LIS_gddpId)//char(0),init_gddp)
   call registerretrievemetforc(trim(LIS_gddpId)//char(0),get_gddp)
   call registertimeinterpmetforc(trim(LIS_gddpId)//char(0), &
                                  timeinterp_gddp)
   call registerfinalmetforc(trim(LIS_gddpId)//char(0),finalize_gddp)
   call registerresetmetforc(trim(LIS_gddpId)//char(0),reset_gddp)
#endif

#if ( defined MF_GALWEM_FORECAST)
   call registerinitmetforc(trim(LIS_galwemId)//char(0),init_galwem)
   call registerretrievemetforc(trim(LIS_galwemId)//char(0),get_galwem)
   call registertimeinterpmetforc(trim(LIS_galwemId)//char(0), &
                                  timeinterp_galwem)
   call registerfinalmetforc(trim(LIS_galwemId)//char(0),finalize_galwem)
   call registerresetmetforc(trim(LIS_galwemId)//char(0),reset_galwem)
#endif

#if ( defined MF_GALWEM_GE_FORECAST)
   call registerinitmetforc(trim(LIS_galwemgeId)//char(0),init_galwemge)
   call registerretrievemetforc(trim(LIS_galwemgeId)//char(0),get_galwemge)
   call registertimeinterpmetforc(trim(LIS_galwemgeId)//char(0), &
                                  timeinterp_galwemge)
   call registerfinalmetforc(trim(LIS_galwemgeId)//char(0),finalize_galwemge)
   call registerresetmetforc(trim(LIS_galwemgeId)//char(0),reset_galwemge)
#endif

#if ( defined MF_MOGREPS_G_FORECAST)
   call registerinitmetforc(trim(LIS_mogrepsgId)//char(0),init_mogrepsg)
   call registerretrievemetforc(trim(LIS_mogrepsgId)//char(0),get_mogrepsg)
   call registertimeinterpmetforc(trim(LIS_mogrepsgId)//char(0), &
                                  timeinterp_mogrepsg)
   call registerfinalmetforc(trim(LIS_mogrepsgId)//char(0),finalize_mogrepsg)
   call registerresetmetforc(trim(LIS_mogrepsgId)//char(0),reset_mogrepsg)
#endif

end subroutine LIS_metforcing_plugin

end module LIS_metforcing_pluginMod
