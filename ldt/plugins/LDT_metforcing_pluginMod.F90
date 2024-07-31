!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_metforcing_pluginMod
!BOP
!
! !MODULE: LDT_metforcing_pluginMod
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   incorporating a new met forcing scheme. The user defined functions
!   are incorporated into the appropriate registry to be later invoked
!   through generic calls. 
!   
! !REVISION HISTORY: 
!  11 Dec 2003   Sujay Kumar:  Initial Specification
!  11 Oct 2014   K. Arsenault: Reorganzing for syncing with LDT  
! 
!EOP  
  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LDT_metforcing_plugin
  
contains
!BOP
! !ROUTINE: LDT_metforcing_plugin
!  \label{LDT_metforcing_plugin}
! 
! !INTERFACE:
  subroutine LDT_metforcing_plugin

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
!   where $ldt\%force$ is set through the congiguration utility, enabling
!   the user to select any of the met forcing schemes, at runtime. 
!EOP

    use LDT_pluginIndices

!    use metForcTemplate_forcingMod
    use gdas_forcingMod
    use ecmwf_forcingMod
    use princeton_forcingMod
!    use rhoneAGG_forcingMod
    use gldas_forcingMod
    use gfs_forcingMod
    use merra2_forcingMod
    use era5_forcingMod
    use gswp1_forcingMod
    use gswp2_forcingMod
#if ( defined MF_AGRMET )
    use AGRMET_forcingMod
#endif
    use agrradps_forcingMod
!    use agrrad_forcingMod
    use geos5fcst_forcingMod
!    use gdasLSWG_forcingMod
    use WRFoutv2_forcingMod
    use WRF_AKdom_forcingMod

    use TRMM3B42RTV7_forcingMod
    use TRMM3B42V6_forcingMod
    use TRMM3B42V7_forcingMod 
    use cmorph_forcingMod
    use stg2_forcingMod
    use stg4_forcingMod
    use cmap_forcingMod
    use nldas2_forcingMod
    use narr_forcingMod
    use nam242_forcingMod
!    use vic_forcingMod
    use RFE2Daily_forcingMod
    use RFE2gdas_forcingMod
    use chirps2_forcingMod
!    use petusgs_forcingMod

!    use scan_forcingMod
!    use ceop_forcingMod
!    use arms_forcingMod
!    use ALMIPII_forcingMod
!    use Bondville_forcingMod
!    use FASSTsingle_forcingMod
!    use snotel_forcingMod
!    use coop_forcingMod
!    use PALSmetdata_forcingMod
!    use pildas_forcingMod
!    use capa_forcingMod
!    use WRFout_forcingMod

!    external get_metForcTemplate
!    external timeinterp_metForcTemplate
!    external finalize_metForcTemplate
!    external reset_metForcTemplate

    external get_gdas
    external timeinterp_gdas
    external reset_gdas
    external finalize_gdas

    external get_ecmwf
    external timeinterp_ecmwf
    external reset_ecmwf
    external finalize_ecmwf

    external get_agrmet
    external timeinterp_agrmet
    external reset_agrmet
    external finalize_agrmet

    external get_princeton
    external timeinterp_princeton
    external finalize_princeton
    external reset_princeton

!    external get_rhoneAGG
!    external timeinterp_rhoneAGG
!    external finalize_rhoneAGG

    external get_gswp2
    external timeinterp_gswp2
    external finalize_gswp2

    external get_gswp1
    external timeinterp_gswp1
    external finalize_gswp1

    external get_gldas
    external timeinterp_gldas
    external finalize_gldas

    external get_gfs
    external timeinterp_gfs
    external finalize_gfs

    external get_merra2
    external timeinterp_merra2
    external finalize_merra2
    external reset_merra2

    external get_era5
    external timeinterp_era5
    external finalize_era5
    external reset_era5

    external get_agrradps
    external timeinterp_agrradps
    external finalize_agrradps

!    external get_agrrad
!    external timeinterp_agrrad
!    external finalize_agrrad

    external get_cmap
    external timeinterp_cmap
    external finalize_cmap
    external reset_cmap

    external get_TRMM3B42RTV7
    external timeinterp_TRMM3B42RTV7
    external finalize_TRMM3B42RTV7

    external get_TRMM3B42V6
    external timeinterp_TRMM3B42V6
    external finalize_TRMM3B42V6

    external get_TRMM3B42V7
    external timeinterp_TRMM3B42V7
    external finalize_TRMM3B42V7

    external get_cmorph
    external timeinterp_cmorph
    external finalize_cmorph

    external get_stg2
    external timeinterp_stg2
    external finalize_stg2

    external get_stg4
    external timeinterp_stg4
    external finalize_stg4
    external reset_stg4

    external get_nldas2
    external timeinterp_nldas2
    external finalize_nldas2
    external reset_nldas2

!    external getvicforcing
!    external time_interp_vicforcing
!    external vicforcing_finalize

    external get_geos5fcst
    external timeinterp_geos5fcst
    external finalize_geos5fcst
    external reset_geos5fcst

!    external get_gdasLSWG
!    external timeinterp_gdasLSWG
!    external finalize_gdasLSWG

    external get_narr
    external timeinterp_narr
    external finalize_narr

    external get_nam242
    external timeinterp_nam242
    external finalize_nam242

    external get_RFE2gdas
    external timeinterp_RFE2gdas
    external finalize_RFE2gdas
    external reset_RFE2gdas

    external get_RFE2Daily
    external timeinterp_RFE2Daily
    external finalize_RFE2Daily
    external reset_RFE2Daily

    external get_chirps2
    external timeinterp_chirps2
    external finalize_chirps2
    external reset_chirps2

!    external get_petusgs
!    external timeinterp_petusgs
!    external finalize_petusgs

!    external get_scan
!    external timeinterp_scan
!    external finalize_scan

!    external get_arms
!    external timeinterp_arms
!    external finalize_arms
!    external reset_arms

!    external get_ceop
!    external timeinterp_ceop
!    external finalize_ceop

!    external get_ALMIPII
!    external timeinterp_ALMIPII
!    external finalize_ALMIPII

!    external get_Bondville
!    external timeinterp_Bondville
!    external finalize_Bondville

!    external get_FASSTsingle
!    external timeinterp_FASSTsingle
!    external finalize_FASSTsingle

!    external get_snotel
!    external timeinterp_snotel
!    external finalize_snotel

!    external get_coop
!    external timeinterp_coop
!    external finalize_coop

!    external get_PALSmetdata
!    external timeinterp_PALSmetdata
!    external finalize_PALSmetdata
!    external reset_PALSmetdata

!    external get_pildas
!    external timeinterp_pildas
!    external finalize_pildas

!    external get_capa
!    external timeinterp_capa
!    external finalize_capa

!    external get_WRFout
!    external timeinterp_WRFout
!    external reset_WRFout
!    external finalize_WRFout


#if 0
! - Meteorological Forcing Template:
    call registerinitmetforc(trim(LDT_metForcTemplateId)//char(0),init_metForctemplate)
    call registerretrievemetforc(trim(LDT_metForcTemplateId)//char(0),get_MetForctemplate)
    call registertimeinterpmetforc(trim(LDT_metForcTemplateId)//char(0),&
         timeinterp_MetForctemplate)
    call registerfinalmetforc(trim(LDT_metForcTemplateId)//char(0),&
         finalize_metForctemplate)
    call registerresetmetforc(trim(LDT_metForcTemplateId)//char(0),&
         reset_metForctemplate)
#endif

! - GDAS (GSFC) Forcing:
    call registerinitmetforc(trim(LDT_gdasId)//char(0),init_gdas)
    call registerretrievemetforc(trim(LDT_gdasId)//char(0),get_gdas)
    call registertimeinterpmetforc(trim(LDT_gdasId)//char(0),timeinterp_gdas)
    call registerresetmetforc(trim(LDT_gdasId)//char(0),reset_gdas)
    call registerfinalmetforc(trim(LDT_gdasId)//char(0),finalize_gdas)

! - ECMWF Forcing:
    call registerinitmetforc(trim(LDT_ecmwfId)//char(0),init_ECMWF)
    call registerretrievemetforc(trim(LDT_ecmwfId)//char(0),get_ecmwf)
    call registertimeinterpmetforc(trim(LDT_ecmwfId)//char(0),timeinterp_ecmwf)
    call registerresetmetforc(trim(LDT_ecmwfId)//char(0),reset_ecmwf)
    call registerfinalmetforc(trim(LDT_ecmwfId)//char(0),finalize_ecmwf)

! - PRINCETON Reanalysis Forcing:
    call registerinitmetforc(trim(LDT_princetonId)//char(0),init_PRINCETON)
    call registerretrievemetforc(trim(LDT_princetonId)//char(0),get_princeton)
    call registertimeinterpmetforc(trim(LDT_princetonId)//char(0),timeinterp_princeton)
    call registerfinalmetforc(trim(LDT_princetonId)//char(0),finalize_princeton)
    call registerresetmetforc(trim(LDT_princetonId)//char(0),reset_princeton)

#if 0
    call registerinitmetforc(trim(LDT_rhoneAGGId)//char(0),init_RHONEAGG)
    call registerretrievemetforc(trim(LDT_rhoneAGGId)//char(0),get_rhoneAGG)
    call registertimeinterpmetforc(trim(LDT_rhoneAGGId)//char(0),timeinterp_rhoneAGG)
    call registerfinalmetforc(trim(LDT_rhoneAGGId)//char(0),finalize_rhoneAGG)
#endif

! - GLDAS Reanalysis Forcing:
    call registerinitmetforc(trim(LDT_gldasId)//char(0),init_GLDAS)
!    call registerretrievemetforc(trim(LDT_gldasId)//char(0),get_gldas)
!    call registertimeinterpmetforc(trim(LDT_gldasId)//char(0),timeinterp_gldas)
!    call registerfinalmetforc(trim(LDT_gldasId)//char(0),finalize_gldas)

! - GFS Forecast Forcing:
    call registerinitmetforc(trim(LDT_gfsId)//char(0),init_GFS)
!    call registerretrievemetforc(trim(LDT_gfsId)//char(0),get_gfs)
!    call registertimeinterpmetforc(trim(LDT_gfsId)//char(0),timeinterp_gfs)
!    call registerfinalmetforc(trim(LDT_gfsId)//char(0),finalize_gfs)

! - MERRA-2 Reanalysis Forcing:
    call registerinitmetforc(trim(LDT_merra2Id)//char(0),init_MERRA2)
    call registerretrievemetforc(trim(LDT_merra2Id)//char(0),get_merra2)
    call registertimeinterpmetforc(trim(LDT_merra2Id)//char(0),timeinterp_merra2)
    call registerresetmetforc(trim(LDT_merra2Id)//char(0),reset_merra2)
    call registerfinalmetforc(trim(LDT_merra2Id)//char(0),finalize_merra2)

! - ERA5 Reanalysis Forcing:
    call registerinitmetforc(trim(LDT_ERA5Id)//char(0),init_ERA5)
    call registerretrievemetforc(trim(LDT_ERA5Id)//char(0),get_ERA5)
    call registertimeinterpmetforc(trim(LDT_ERA5Id)//char(0),timeinterp_ERA5)
    call registerresetmetforc(trim(LDT_ERA5Id)//char(0),reset_ERA5)
    call registerfinalmetforc(trim(LDT_ERA5Id)//char(0),finalize_ERA5)

! - WRFv2 Analysis Forcing:
    call registerinitmetforc(trim(LDT_wrfoutv2Id)//char(0),init_WRFoutv2)

! - WRF Alaska Forcing:
    call registerinitmetforc(trim(LDT_wrfakId)//char(0),init_WRF_AKdom)

! - GSWP2 Forcing:
    call registerinitmetforc(trim(LDT_gswp2Id)//char(0),init_GSWP2)
!    call registerretrievemetforc(trim(LDT_gswp2Id)//char(0),get_gswp2)
!    call registertimeinterpmetforc(trim(LDT_gswp2Id)//char(0),timeinterp_gswp2)
!    call registerfinalmetforc(trim(LDT_gswp2Id)//char(0),finalize_gswp2)
! - GWSP1 Forcing:
    call registerinitmetforc(trim(LDT_gswp1Id)//char(0),init_GSWP1)
!    call registerretrievemetforc(trim(LDT_gswp1Id)//char(0),get_gswp1)
!    call registertimeinterpmetforc(trim(LDT_gswp1Id)//char(0),timeinterp_gswp1)
!    call registerfinalmetforc(trim(LDT_gswp1Id)//char(0),finalize_gswp1)

#if ( defined MF_AGRMET )
! - AGRMET Forcing:
    call registerinitmetforc(trim(LDT_agrmetId)//char(0),init_AGRMET)
    call registerretrievemetforc(trim(LDT_agrmetId)//char(0),get_agrmet)
    call registertimeinterpmetforc(trim(LDT_agrmetId)//char(0),timeinterp_agrmet)
    call registerresetmetforc(trim(LDT_agrmetId)//char(0),reset_agrmet)
    call registerfinalmetforc(trim(LDT_agrmetId)//char(0),finalize_agrmet)
#endif

! - AGRMET Polar Stereographic Radiation:
    call registerinitmetforc(trim(LDT_agrradpsId)//char(0),init_AGRRADPS)
!    call registerretrievemetforc(trim(LDT_agrradpsId)//char(0),get_agrradps)
!    call registertimeinterpmetforc(trim(LDT_agrradpsId)//char(0),timeinterp_agrradps)
!    call registerfinalmetforc(trim(LDT_agrradpsId)//char(0),finalize_agrradps)

#if 0
! - AGRMET Radiation:
    call registerinitmetforc(trim(LDT_agrradId)//char(0),init_AGRRAD)
    call registerretrievemetforc(trim(LDT_agrradId)//char(0),get_agrrad)
    call registertimeinterpmetforc(trim(LDT_agrradId)//char(0),timeinterp_agrrad)
    call registerfinalmetforc(trim(LDT_agrradId)//char(0),finalize_agrrad)
#endif

! - NLDAS2 Forcing:
    call registerinitmetforc(trim(LDT_nldas2Id)//char(0),init_NLDAS2)
    call registerretrievemetforc(trim(LDT_nldas2Id)//char(0),get_nldas2)
    call registertimeinterpmetforc(trim(LDT_nldas2Id)//char(0),timeinterp_nldas2)
    call registerfinalmetforc(trim(LDT_nldas2Id)//char(0),finalize_nldas2)
    call registerresetmetforc(trim(LDT_nldas2Id)//char(0),reset_nldas2)

! - NAM242 forcing
    call registerinitmetforc(trim(LDT_nam242Id)//char(0),init_nam242)
    call registerretrievemetforc(trim(LDT_nam242Id)//char(0),get_nam242)
    call registertimeinterpmetforc(trim(LDT_nam242Id)//char(0),timeinterp_nam242)
    call registerfinalmetforc(trim(LDT_nam242Id)//char(0),finalize_nam242)

!! - PRECIPITATION-ONLY DATASETS - !!

! - CMAP/GDAS (analysis) Forcing:
    call registerinitmetforc(trim(LDT_cmapId)//char(0),init_CMAP)
    call registerretrievemetforc(trim(LDT_cmapId)//char(0),get_cmap)
    call registertimeinterpmetforc(trim(LDT_cmapId)//char(0),timeinterp_cmap)
    call registerfinalmetforc(trim(LDT_cmapId)//char(0), finalize_cmap)
    call registerresetmetforc(trim(LDT_cmapId)//char(0),reset_cmap)

! - 3B42RTV7 TRMM (Near-) Real-time (added by K. Arsenault; original by Yudong)
    call registerinitmetforc(trim(LDT_TRMM3B42RTV7Id)//char(0),init_TRMM3B42RTV7)
    call registerretrievemetforc(trim(LDT_TRMM3B42RTV7Id)//char(0),get_TRMM3B42RTV7)
    call registertimeinterpmetforc(trim(LDT_TRMM3B42RTV7Id)//char(0),timeinterp_TRMM3B42RTV7)
    call registerfinalmetforc(trim(LDT_TRMM3B42RTV7Id)//char(0),finalize_TRMM3B42RTV7)

! - 3B42V6 TRMM version 6 (added by Yudong)
    call registerinitmetforc(trim(LDT_TRMM3B42V6Id)//char(0),init_TRMM3B42V6)
    call registerretrievemetforc(trim(LDT_TRMM3B42V6Id)//char(0),get_TRMM3B42V6)
    call registertimeinterpmetforc(trim(LDT_TRMM3B42V6Id)//char(0),timeinterp_TRMM3B42V6)
    call registerfinalmetforc(trim(LDT_TRMM3B42V6Id)//char(0),finalize_TRMM3B42V6)

! - 3B42V7 TRMM version 7 (added by Soni)
    call registerinitmetforc(trim(LDT_TRMM3B42V7Id)//char(0),init_TRMM3B42V7)
    call registerretrievemetforc(trim(LDT_TRMM3B42V7Id)//char(0),get_TRMM3B42V7)
    call registertimeinterpmetforc(trim(LDT_TRMM3B42V7Id)//char(0),timeinterp_TRMM3B42V7)
    call registerfinalmetforc(trim(LDT_TRMM3B42V7Id)//char(0),finalize_TRMM3B42V7)

! - CMORPH 30min, 8KM, added by Yudong
    call registerinitmetforc(trim(LDT_cmorphId)//char(0),init_CMORPH)
    call registerretrievemetforc(trim(LDT_cmorphId)//char(0),get_cmorph)
    call registertimeinterpmetforc(trim(LDT_cmorphId)//char(0),timeinterp_cmorph)
    call registerfinalmetforc(trim(LDT_cmorphId)//char(0),finalize_cmorph)

! - STAGE 2 HRAP, 4KM, added by K. Arsenault
    call registerinitmetforc(trim(LDT_stg2Id)//char(0),init_STG2)
    call registerretrievemetforc(trim(LDT_stg2Id)//char(0),get_stg2)
    call registertimeinterpmetforc(trim(LDT_stg2Id)//char(0),timeinterp_stg2)
    call registerfinalmetforc(trim(LDT_stg2Id)//char(0),finalize_stg2)
! - STAGE 4 HRAP, 4KM, added by K. Arsenault
    call registerinitmetforc(trim(LDT_stg4Id)//char(0),init_STG4)
    call registerretrievemetforc(trim(LDT_stg4Id)//char(0),get_stg4)
    call registertimeinterpmetforc(trim(LDT_stg4Id)//char(0),timeinterp_stg4)
    call registerfinalmetforc(trim(LDT_stg4Id)//char(0),finalize_stg4)
    call registerresetmetforc(trim(LDT_stg4Id)//char(0),reset_stg4)

! - FEWSNET CPC RFE2.0 Data
    call registerinitmetforc(trim(LDT_RFE2DailyId)//char(0),init_RFE2Daily)
    call registerretrievemetforc(trim(LDT_RFE2DailyId)//char(0),get_RFE2Daily)
    call registertimeinterpmetforc(trim(LDT_RFE2DailyId)//char(0),timeinterp_RFE2Daily)
    call registerfinalmetforc(trim(LDT_RFE2DailyId)//char(0),finalize_RFE2Daily)
    call registerresetmetforc(trim(LDT_RFE2DailyId)//char(0),reset_RFE2Daily)

! - FEWSNET RFE2.0-GDAS (Temporally downscaled) Data
    call registerinitmetforc(trim(LDT_RFE2gdasId)//char(0),init_RFE2gdas)
    call registerretrievemetforc(trim(LDT_RFE2gdasId)//char(0),get_RFE2gdas)
    call registertimeinterpmetforc(trim(LDT_RFE2gdasId)//char(0),timeinterp_RFE2gdas)
    call registerfinalmetforc(trim(LDT_RFE2gdasId)//char(0),finalize_RFE2gdas)
    call registerresetmetforc(trim(LDT_RFE2gdasId)//char(0),reset_RFE2gdas)

! - CHIRPS 2.0 Precipitation Data 
    call registerinitmetforc(trim(LDT_chirps2Id)//char(0),init_chirps2)
    call registerretrievemetforc(trim(LDT_chirps2Id)//char(0),get_chirps2)
    call registertimeinterpmetforc(trim(LDT_chirps2Id)//char(0),timeinterp_chirps2)
    call registerfinalmetforc(trim(LDT_chirps2Id)//char(0),finalize_chirps2)
    call registerresetmetforc(trim(LDT_chirps2Id)//char(0),reset_chirps2)

!! - END OF PRECIPITATION-ONLY DATASETS - !!

#if 0
! - FEWSNET USGS PET Data
    call registerinitmetforc(trim(LDT_USGSPETforcId)//char(0),init_petusgs)
    call registerretrievemetforc(trim(LDT_USGSPETforcId)//char(0),get_petusgs)
    call registertimeinterpmetforc(trim(LDT_USGSPETforcId)//char(0),timeinterp_petusgs)
    call registerfinalmetforc(trim(LDT_USGSPETforcId)//char(0),finalize_petusgs)
#endif

! - NARR profile data for CRTM 
    call registerinitmetforc(trim(LDT_narrId)//char(0),init_NARR)
!    call registerretrievemetforc(trim(LDT_narrId)//char(0),get_narr)
!    call registertimeinterpmetforc(trim(LDT_narrId)//char(0),timeinterp_narr)
!    call registerfinalmetforc(trim(LDT_narrId)//char(0),finalize_narr)

#if 0
! - GDAS LSWG profiles
    call registerinitmetforc(trim(LDT_gdasLSWGId)//char(0),init_gdasLSWG)
    call registerretrievemetforc(trim(LDT_gdasLSWGId)//char(0),get_gdasLSWG)
    call registertimeinterpmetforc(trim(LDT_gdasLSWGId)//char(0),timeinterp_gdasLSWG)
    call registerfinalmetforc(trim(LDT_gdasLSWGId)//char(0),finalize_gdasLSWG)
#endif

! - GEOSv5 forecast forcing:
    call registerinitmetforc(trim(LDT_geos5fcstId)//char(0),init_GEOS5FCST)
    call registerretrievemetforc(trim(LDT_geos5fcstId)//char(0),get_geos5fcst)
    call registertimeinterpmetforc(trim(LDT_geos5fcstId)//char(0),timeinterp_geos5fcst)
    call registerfinalmetforc(trim(LDT_geos5fcstId)//char(0),finalize_geos5fcst)
    call registerresetmetforc(trim(LDT_geos5fcstId)//char(0),reset_geos5fcst)

#if 0
! - VIC processed data
    call registerinitmetforc(trim(LDT_vicforcingId)//char(0),defineNativevicforcing)
    call registerretrievemetforc(trim(LDT_vicforcingId)//char(0),getvicforcing)
    call registertimeinterpmetforc(trim(LDT_vicforcingId)//char(0),time_interp_vicforcing)
    call registerfinalmetforc(trim(LDT_vicforcingId)//char(0),vicforcing_finalize)
#endif

#if 0
! - CEOP
    call registerinitmetforc(trim(LDT_ceopId)//char(0),init_CEOP)
    call registerretrievemetforc(trim(LDT_ceopId)//char(0),get_ceop)
    call registertimeinterpmetforc(trim(LDT_ceopId)//char(0),timeinterp_ceop)
    call registerfinalmetforc(trim(LDT_ceopId)//char(0), finalize_ceop)

! - SCAN Forcing:
    call registerinitmetforc(trim(LDT_scanId)//char(0),init_SCAN)
    call registerretrievemetforc(trim(LDT_scanId)//char(0), get_scan)
    call registertimeinterpmetforc(trim(LDT_scanId)//char(0),timeinterp_scan)
    call registerfinalmetforc(trim(LDT_scanId)//char(0), finalize_scan)

! - SNOTEL Precipitation Forcing (added by Yuqiong Liu):
    call registerinitmetforc(trim(LDT_snotelId)//char(0),init_SNOTEL)
    call registerretrievemetforc(trim(LDT_snotelId)//char(0),get_snotel)
    call registertimeinterpmetforc(trim(LDT_snotelId)//char(0),timeinterp_snotel)
    call registerfinalmetforc(trim(LDT_snotelId)//char(0), finalize_snotel)

! - COOP Precipitation Forcing (added by Yuqiong Liu):
    call registerinitmetforc(trim(LDT_coopId)//char(0),init_COOP)
    call registerretrievemetforc(trim(LDT_coopId)//char(0),get_coop)
    call registertimeinterpmetforc(trim(LDT_coopId)//char(0),timeinterp_coop)
    call registerfinalmetforc(trim(LDT_coopId)//char(0), finalize_coop)

! - ARMS station data
    call registerinitmetforc(trim(LDT_armsId)//char(0),init_ARMS)
    call registerretrievemetforc(trim(LDT_armsId)//char(0),get_arms)
    call registertimeinterpmetforc(trim(LDT_armsId)//char(0),timeinterp_arms)
    call registerfinalmetforc(trim(LDT_armsId)//char(0),finalize_arms)
    call registerresetmetforc(trim(LDT_armsId)//char(0),reset_arms)

! - ALMIPII
    call registerinitmetforc(trim(LDT_ALMIPIIId)//char(0),init_ALMIPII)
    call registerretrievemetforc(trim(LDT_ALMIPIIId)//char(0),get_ALMIPII)
    call registertimeinterpmetforc(trim(LDT_ALMIPIIId)//char(0),timeinterp_ALMIPII)
    call registerfinalmetforc(trim(LDT_ALMIPIIId)//char(0),finalize_ALMIPII)

! - Noah3.1 Bondville test case
    call registerinitmetforc(trim(LDT_BondvilleId)//char(0),init_Bondville)
    call registerretrievemetforc(trim(LDT_BondvilleId)//char(0), get_Bondville)
    call registertimeinterpmetforc(trim(LDT_BondvilleId)//char(0),timeinterp_Bondville)
    call registerfinalmetforc(trim(LDT_BondvilleId)//char(0), finalize_Bondville)

! - FASST single point test case
    call registerinitmetforc(trim(LDT_FASSTsingleId)//char(0),init_FASSTsingle)
    call registerretrievemetforc(trim(LDT_FASSTsingleId)//char(0), get_FASSTsingle)
    call registertimeinterpmetforc(trim(LDT_FASSTsingleId)//char(0),timeinterp_FASSTsingle)
    call registerfinalmetforc(trim(LDT_FASSTsingleId)//char(0), finalize_FASSTsingle)

! - PALS station data
    call registerinitmetforc(trim(LDT_PALSmetforcId)//char(0),init_PALSMETDATA)
    call registerretrievemetforc(trim(LDT_PALSmetforcId)//char(0),get_PALSmetdata)
    call registertimeinterpmetforc(trim(LDT_PALSmetforcId)//char(0),timeinterp_PALSmetdata)
    call registerfinalmetforc(trim(LDT_PALSmetforcId)//char(0),finalize_PALSmetdata)
    call registerresetmetforc(trim(LDT_PALSmetforcId)//char(0),reset_PALSmetdata)

! - PILDAS station data
    call registerinitmetforc(trim(LDT_pildasmetforcId)//char(0),init_pildas)
    call registerretrievemetforc(trim(LDT_pildasmetforcId)//char(0),get_pildas)
    call registertimeinterpmetforc(trim(LDT_pildasmetforcId)//char(0),timeinterp_pildas)
    call registerfinalmetforc(trim(LDT_pildasmetforcId)//char(0),finalize_pildas)

! - CAPA precipitation
    call registerinitmetforc(trim(LDT_capaId)//char(0),init_capa)
    call registerretrievemetforc(trim(LDT_capaId)//char(0),get_capa)
    call registertimeinterpmetforc(trim(LDT_capaId)//char(0),timeinterp_capa)
    call registerfinalmetforc(trim(LDT_capaId)//char(0),finalize_capa)

! - WRFout forcing
    call registerinitmetforc(trim(LDT_WRFoutId)//char(0),init_WRFout)
    call registerretrievemetforc(trim(LDT_WRFoutId)//char(0),get_WRFout)
    call registertimeinterpmetforc(trim(LDT_WRFoutId)//char(0),timeinterp_WRFout)
    call registerfinalmetforc(trim(LDT_WRFoutId)//char(0),finalize_WRFout)
    call registerresetmetforc(trim(LDT_WRFoutId)//char(0),reset_WRFout)
#endif

  end subroutine LDT_metforcing_plugin

end module LDT_metforcing_pluginMod
