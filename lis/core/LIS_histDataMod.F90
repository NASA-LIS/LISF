!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_histDataMod
!BOP
!
!  !MODULE: LIS_histDataMod
! 
!  !DESCRIPTION: 
!   This module is used by the user to define the metadata associated with 
!   model output. The module lists a superset of the land surface model
!   variables. The user can choose a subset from this list through the 
!   lis configuration file for model output. Currently this list includes
!   the variable definitions from ALMA specification. 
!
!   \textsl{http://www.lmd.jussieu.fr/ALMA/}
!   
!  !REVISION HISTORY: 
!  21 Oct 2005    Sujay Kumar  Initial Specification
!  19 Jan 2007    Chuck Alonge Added Snow Depth Option (Future use in WPS)
!   4 Jul 2008    Sujay Kumar Redesigned the I/O to enable generic I/O for  
!                     all LSMs. 
!  11 May 2011    Sujay Kumar, Updated to be generic across all LIS I/O and
!                     not just land surface model output
!  26 May 2011    Soni Yatheendradas: Added Potential ET forcings for FEWSNET
!  11 Jan 2012    Jiarui Dong, Added the Sac-HT/Snow17 implementation
!  28 Jan 2014    David Mocko; Updates for AFWA GRIB-1 files using GRIBAPI library
!  21 Feb 2014    David Mocko; More updates for matching AFWA GRIB-1 files
!                     to the LIS-6 GRIB output
!  14 Mar 2014    David Mocko: Added CAPE, CRAINFFORC, CMFORC, and CHFORC outputs
!  28 Mar 2014    David Mocko: More refinements/fixes to AFWA GRIB-1 CONFIGS
!   6 Apr 2015    Hiroko Beaudoing: Added GRIB2 specific surface type indexes
!                 106 and 114 in place for 112 (some variables left alone) 
!  01 Jun 2017    Augusto Getirana: Add/update HyMAP2 outputs [SWS, differetial and 
!                      potential evaporation and deep water infiltration (DWI)]
!  
!  
!EOP
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_mpiMod

  implicit none

  PRIVATE 

  public :: LIS_histDataInit
  public :: LIS_diagnoseSurfaceOutputVar
  public :: LIS_diagnoseRoutingOutputVar
  public :: LIS_diagnoseRTMOutputVar
  public :: LIS_diagnoseIrrigationOutputVar
  public :: LIS_resetOutputVars
  public :: LIS_rescaleCount

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
! MOC - Model Output Convention
!-----------------------------------------------------------------------------
  public :: LIS_histData
  ! LSM
  public :: LIS_MOC_SWNET   
  public :: LIS_MOC_LWNET   
  public :: LIS_MOC_QLE     
  public :: LIS_MOC_QH      
  public :: LIS_MOC_QG      
  public :: LIS_MOC_QF      
  public :: LIS_MOC_QV      
  public :: LIS_MOC_QTAU    
  public :: LIS_MOC_QA         
  public :: LIS_MOC_DELSURFHEAT
  public :: LIS_MOC_DELCOLDCONT
  public :: LIS_MOC_BR
  public :: LIS_MOC_EF
  public :: LIS_MOC_SNOWF     
  public :: LIS_MOC_RAINF     
  public :: LIS_MOC_EVAP      
  public :: LIS_MOC_QS        
  public :: LIS_MOC_QREC      
  public :: LIS_MOC_QSB       
  public :: LIS_MOC_QSM       
  public :: LIS_MOC_QFZ       
  public :: LIS_MOC_QST       
  public :: LIS_MOC_DELSOILMOIST
  public :: LIS_MOC_DELSWE    
  public :: LIS_MOC_DELSURFSTOR
  public :: LIS_MOC_DELINTERCEPT
  public :: LIS_MOC_SNOWT     
  public :: LIS_MOC_VEGT      
  public :: LIS_MOC_BARESOILT 
  public :: LIS_MOC_AVGSURFT  
  public :: LIS_MOC_RADT      
  public :: LIS_MOC_ALBEDO    
  public :: LIS_MOC_ALBEDO_DIR_V    
  public :: LIS_MOC_ALBEDO_DIF_V    
  public :: LIS_MOC_ALBEDO_DIR_N    
  public :: LIS_MOC_ALBEDO_DIF_N    
  public :: LIS_MOC_SWE       
  public :: LIS_MOC_SNOWDENSITY
  public :: LIS_MOC_SNOWGRAIN
  public :: LIS_MOC_SWEVEG    
  public :: LIS_MOC_SURFSTOR  
  public :: LIS_MOC_SOILMOIST 
  public :: LIS_MOC_SOILTEMP  
  public :: LIS_MOC_SMLIQFRAC
  public :: LIS_MOC_SMFROZFRAC
  public :: LIS_MOC_SOILWET   
  public :: LIS_MOC_MATRICPOTENTIAL
  public :: LIS_MOC_POTEVAP   
  public :: LIS_MOC_ECANOP    
  public :: LIS_MOC_TVEG      
  public :: LIS_MOC_ESOIL     
  public :: LIS_MOC_EWATER    
  public :: LIS_MOC_ROOTMOIST 
  public :: LIS_MOC_ROOTTEMP
  public :: LIS_MOC_CANOPINT  
  public :: LIS_MOC_EVAPSNOW     
  public :: LIS_MOC_SUBSNOW   
  public :: LIS_MOC_SUBSURF   
  public :: LIS_MOC_ACOND   
  public :: LIS_MOC_WATERTABLED
  public :: LIS_MOC_TWS
  public :: LIS_MOC_GWS
  public :: LIS_MOC_SNOWCOVER
  public :: LIS_MOC_SALBEDO
  public :: LIS_MOC_SNOWTPROF
  public :: LIS_MOC_SNOWDEPTH
  public :: LIS_MOC_SNOWTHICK
  public :: LIS_MOC_SLIQFRAC
  public :: LIS_MOC_LAYERSNOWDEPTH
  public :: LIS_MOC_LAYERSNOWDENSITY
  public :: LIS_MOC_LWUP
  public :: LIS_MOC_GPP
  public :: LIS_MOC_NPP
  public :: LIS_MOC_NEE
  public :: LIS_MOC_AUTORESP
  public :: LIS_MOC_HETERORESP
  public :: LIS_MOC_LEAFRESP
  public :: LIS_MOC_TOTSOILCARB
  public :: LIS_MOC_TOTLIVBIOM
  
  public :: LIS_MOC_WINDFORC  
  public :: LIS_MOC_RAINFFORC 
  public :: LIS_MOC_SNOWFFORC 
  public :: LIS_MOC_CRAINFFORC 
  public :: LIS_MOC_LSRAINFFORC 
  public :: LIS_MOC_CSNOWFFORC 
  public :: LIS_MOC_LSSNOWFFORC 
  public :: LIS_MOC_TAIRFORC  
  public :: LIS_MOC_QAIRFORC  
  public :: LIS_MOC_PSURFFORC 
  public :: LIS_MOC_SWDOWNFORC
  public :: LIS_MOC_LWDOWNFORC
  public :: LIS_MOC_DIRECTSWFORC
  public :: LIS_MOC_DIFFUSESWFORC
  public :: LIS_MOC_NWINDFORC  
  public :: LIS_MOC_EWINDFORC  
  public :: LIS_MOC_FHEIGHTFORC  
  public :: LIS_MOC_CHFORC  
  public :: LIS_MOC_CMFORC  
  public :: LIS_MOC_EMISSFORC  
  public :: LIS_MOC_MIXRATIOFORC  
  public :: LIS_MOC_COSZENFORC  
  public :: LIS_MOC_ALBEDOFORC  
  public :: LIS_MOC_PARDRFORC  
  public :: LIS_MOC_PARDFFORC  
!<for vic>
  public :: LIS_MOC_SNOWFLAGFORC
  public :: LIS_MOC_DENSITYFORC
  public :: LIS_MOC_VAPORPRESSFORC
  public :: LIS_MOC_VAPORPRESSDEFICITFORC
  public :: LIS_MOC_ARESIST
!</for vic>

  public :: LIS_MOC_LANDMASK  
  public :: LIS_MOC_LANDCOVER 
  public :: LIS_MOC_SOILTYPE  
  public :: LIS_MOC_SANDFRAC
  public :: LIS_MOC_CLAYFRAC
  public :: LIS_MOC_SILTFRAC
  public :: LIS_MOC_POROSITY
  public :: LIS_MOC_SOILCOLOR 
  public :: LIS_MOC_ELEVATION
  public :: LIS_MOC_SLOPE
  public :: LIS_MOC_LAI       
  public :: LIS_MOC_SAI       
  public :: LIS_MOC_SNFRALBEDO
  public :: LIS_MOC_MXSNALBEDO
  public :: LIS_MOC_GREENNESS 
  public :: LIS_MOC_TEMPBOT  

  public :: LIS_MOC_CCOND
  public :: LIS_MOC_RELSMC
  public :: LIS_MOC_RHMIN
  public :: LIS_MOC_TOTALPRECIP     
  public :: LIS_MOC_CRAINF
  PUBLIC :: LIS_MOC_LSRAINF
  PUBLIC :: LIS_MOC_LSSNOWF
  PUBLIC :: LIS_MOC_CSNOWF
  public :: LIS_MOC_SOILET
  public :: LIS_MOC_Z0BRD
  public :: LIS_MOC_ROUGHNESS
  public :: LIS_MOC_THERMAL_ROUGHNESS
  public :: LIS_MOC_T2DIAG
  public :: LIS_MOC_Q2DIAG
  public :: LIS_MOC_RNET   
  public :: LIS_MOC_CH
  public :: LIS_MOC_CM
  public :: LIS_MOC_MIXRATIO

  public :: LIS_MOC_SACUZTWC
  public :: LIS_MOC_SACUZFWC
  public :: LIS_MOC_SACLZTWC
  public :: LIS_MOC_SACLZFSC
  public :: LIS_MOC_SACLZFPC
  public :: LIS_MOC_SACADIMPC
  public :: LIS_MOC_SNOW17SWE
  public :: LIS_MOC_SNOW17LIQW
  public :: LIS_MOC_SNOW17NEGHS
  public :: LIS_MOC_SNOW17ACCMAX
  public :: LIS_MOC_SNOW17AEADJ
  public :: LIS_MOC_SNOW17RMLT

!</for vic>
  public :: LIS_MOC_VIC_PET_SATSOIL
  public :: LIS_MOC_VIC_PET_H2OSURF
  public :: LIS_MOC_VIC_PET_SHORT
  public :: LIS_MOC_VIC_PET_TALL
  public :: LIS_MOC_VIC_PET_NATVEG
  public :: LIS_MOC_VIC_PET_VEGNOCR
!</for vic>

!FLDAS
  public :: LIS_MOC_PETFORC
  public :: LIS_MOC_REFETFORC
  ! FLDAS-WRSI OUTPUTS LIST
  public :: LIS_MOC_SOS
  public :: LIS_MOC_WRSI
  public :: LIS_MOC_KF2
  public :: LIS_MOC_SumWR
  public :: LIS_MOC_SumET
  public :: LIS_MOC_SWI
  public :: LIS_MOC_SOSa
  public :: LIS_MOC_TotalSurplusWater
  public :: LIS_MOC_MaxSurplusWater
  public :: LIS_MOC_TotalWaterDeficit
  public :: LIS_MOC_MaxWaterDeficit
  public :: LIS_MOC_TotalAETInitial
  public :: LIS_MOC_TotalWRInitial
  public :: LIS_MOC_TotalSurplusWaterInitial
  public :: LIS_MOC_TotalWaterDeficitInitial
  public :: LIS_MOC_TotalAETVeg
  public :: LIS_MOC_TotalWRVeg
  public :: LIS_MOC_TotalSurplusWaterVeg
  public :: LIS_MOC_TotalWaterDeficitVeg
  public :: LIS_MOC_TotalAETFlower
  public :: LIS_MOC_TotalWRFlower
  public :: LIS_MOC_TotalSurplusWaterFlower
  public :: LIS_MOC_TotalWaterDeficitFlower
  public :: LIS_MOC_TotalAETRipe
  public :: LIS_MOC_TotalWRRipe
  public :: LIS_MOC_TotalSurplusWaterRipe
  public :: LIS_MOC_TotalWaterDeficitRipe
  public :: LIS_MOC_PermWiltDate
  public :: LIS_MOC_Wilting1
  public :: LIS_MOC_Wilting2
  public :: LIS_MOC_WRSIa
  public :: LIS_MOC_growing_season
  public :: LIS_MOC_WHC
  public :: LIS_MOC_LGP
  public :: LIS_MOC_WR_TimeStep ! SY
  public :: LIS_MOC_AET_TimeStep ! SY
  public :: LIS_MOC_WRSI_TimeStep ! SY
  public :: LIS_MOC_SurplusWater_TimeStep ! SY

! NLDAS
  public :: LIS_MOC_CAPEFORC

  ! Routing

  public :: LIS_MOC_RIVSTO
  public :: LIS_MOC_RIVDPH
  public :: LIS_MOC_RIVVEL
  public :: LIS_MOC_STREAMFLOW
  public :: LIS_MOC_FLDOUT
  public :: LIS_MOC_FLDEVAP
  public :: LIS_MOC_FLDSTO
  public :: LIS_MOC_FLDDPH
  public :: LIS_MOC_FLDVEL
  public :: LIS_MOC_FLDFRC
  public :: LIS_MOC_FLDARE  
  public :: LIS_MOC_SFCELV
  public :: LIS_MOC_RNFSTO
  public :: LIS_MOC_BSFSTO
  public :: LIS_MOC_RNFDWI
  public :: LIS_MOC_BSFDWI
  public :: LIS_MOC_SURFWS
  public :: LIS_MOC_EWAT
  public :: LIS_MOC_EDIF

  ! RTM
  public :: LIS_MOC_RTM_EMISSIVITY
  public :: LIS_MOC_RTM_TB
  public :: LIS_MOC_RTM_SM
  ! Irrigation
  public :: LIS_MOC_IRRIGATEDWATER  

  public :: LIS_MOC_LSM_COUNT
  public :: LIS_MOC_ROUTING_COUNT
  public :: LIS_MOC_RTM_COUNT
  public :: LIS_MOC_IRRIG_COUNT

  ! FLAKE 2003, Added by Shugong Wang 05/20/2013 
  public ::   LIS_MOC_LAKE_T_SNOW
  public ::   LIS_MOC_LAKE_T_ICE
  public ::   LIS_MOC_LAKE_T_MNW
  public ::   LIS_MOC_LAKE_T_WML
  public ::   LIS_MOC_LAKE_T_BOT
  public ::   LIS_MOC_LAKE_T_B1
  public ::   LIS_MOC_LAKE_C_T
  public ::   LIS_MOC_LAKE_H_SNOW
  public ::   LIS_MOC_LAKE_H_ICE
  public ::   LIS_MOC_LAKE_H_ML
  public ::   LIS_MOC_LAKE_H_B1
  public ::   LIS_MOC_LAKE_T_SFC
  public ::   LIS_MOC_LAKE_ALBEDO_WATER
  public ::   LIS_MOC_LAKE_ALBEDO_ICE
  public ::   LIS_MOC_LAKE_ALBEDO_SNOW
  public ::   LIS_MOC_LAKE_UFR_A
  public ::   LIS_MOC_LAKE_UFR_W
  public ::   LIS_MOC_LAKE_WCONV
  public ::   LIS_MOC_LAKE_Q_SE
  public ::   LIS_MOC_LAKE_Q_LA
  public ::   LIS_MOC_LAKE_I_W
  public ::   LIS_MOC_LAKE_Q_LWA
  public ::   LIS_MOC_LAKE_Q_LWW
  public ::   LIS_MOC_LAKE_Q_BOT 

  ! SAC-HTET and Snow-17

  public ::   LIS_MOC_SACUZTWH
  public ::   LIS_MOC_SACUZFWH
  public ::   LIS_MOC_SACLZTWH
  public ::   LIS_MOC_SACLZFSH
  public ::   LIS_MOC_SACLZFPH
  
  public ::   LIS_MOC_SACSWINT
  public ::   LIS_MOC_SACTSINT
  public ::   LIS_MOC_SACSWHINT
  public ::   LIS_MOC_SACFROST
  
  ! NoahMP
  public ::   LIS_MOC_CANOPY_TEMP 
  public ::   LIS_MOC_CANOPY_VP   
  public ::   LIS_MOC_CANOPY_WF   
  public ::   LIS_MOC_CANOPY_INTL   
  public ::   LIS_MOC_GROUNDAVGT   
  public ::   LIS_MOC_GROUNDVEGT   
  public ::   LIS_MOC_SOWN_NLAYER   
  public ::   LIS_MOC_SNOW_LBDFSS   
  public ::   LIS_MOC_SOIL_LBDFSS   
  public ::   LIS_MOC_SNOWICE   
  public ::   LIS_MOC_SNOWLIQ   
  public ::   LIS_MOC_WT_AQUI_SATSOIL
  public ::   LIS_MOC_LAKEWATER   
  public ::   LIS_MOC_LEAFMASS    
  public ::   LIS_MOC_ROOTMASS    
  public ::   LIS_MOC_STEMMASS    
  public ::   LIS_MOC_WOODMASS    
  public ::   LIS_MOC_CARBON_DEEPSOIL   
  public ::   LIS_MOC_CARBON_SHALLOWSOIL    
  public ::   LIS_MOC_SNOWAGE   
  public ::   LIS_MOC_BETWEENWATER    
  public ::   LIS_MOC_QRECTOGW    
  public ::   LIS_MOC_QRECFROMGW    
  public ::   LIS_MOC_FSR   
  public ::   LIS_MOC_FCEV    
  public ::   LIS_MOC_FGEV    
  public ::   LIS_MOC_FCTR    
  public ::   LIS_MOC_VEGE2MT   
  public ::   LIS_MOC_BARE2MT   
  public ::   LIS_MOC_VEGE2MQ2    
  public ::   LIS_MOC_BARE2MQ2    
  public ::   LIS_MOC_APAR    
  public ::   LIS_MOC_PSCO2   
  public ::   LIS_MOC_SAV   
  public ::   LIS_MOC_SAG   
  public ::   LIS_MOC_PONDING    
  public ::   LIS_MOC_PONDING1   
  public ::   LIS_MOC_PONDING2   
  public ::   LIS_MOC_RSSUN   
  public ::   LIS_MOC_RSSHA   
  public ::   LIS_MOC_BGAP    
  public ::   LIS_MOC_WGAP    
  public ::   LIS_MOC_CHV   
  public ::   LIS_MOC_CHB   
  public ::   LIS_MOC_SHG   
  public ::   LIS_MOC_SHC   
  public ::   LIS_MOC_SHB   
  public ::   LIS_MOC_EVG   
  public ::   LIS_MOC_EVB   
  public ::   LIS_MOC_GHV   
  public ::   LIS_MOC_GHB   
  public ::   LIS_MOC_IRV   
  public ::   LIS_MOC_IRC   
  public ::   LIS_MOC_IRB   
  public ::   LIS_MOC_HTR    
  public ::   LIS_MOC_HEVC   
  public ::   LIS_MOC_CHLEAF    
  public ::   LIS_MOC_CHUC    
  public ::   LIS_MOC_CHV2    
  public ::   LIS_MOC_CHB2    
  public ::   LIS_MOC_FPICE   
  ! end Noahmp
 
  ! RUC 
  public :: LIS_MOC_QVG
  public :: LIS_MOC_QCG
  public :: LIS_MOC_QSG
  public :: LIS_MOC_SNNOT75CM
  public :: LIS_MOC_FRZPREC_DEN
  public :: LIS_MOC_FRZPREC
  public :: LIS_MOC_ACC_SNOWF
  public :: LIS_MOC_ACC_FRZPREC
  public :: LIS_MOC_ACC_EVAP
  public :: LIS_MOC_QSFC
  public :: LIS_MOC_DEW_FROST
  public :: LIS_MOC_DRIP
  public :: LIS_QH_SNOW
  public :: LIS_MOC_SNOWTHRESH
  ! end ruc 

  ! JULES
  public :: LIS_MOC_GS
  public :: LIS_MOC_GC 
  PUBLIC :: LIS_MOC_JULES_STHZW
  PUBLIC :: LIS_MOC_JULES_STHU
  PUBLIC :: LIS_MOC_JULES_STHU_MIN
  PUBLIC :: LIS_MOC_JULES_STHF
  PUBLIC :: LIS_MOC_JULES_SMVCCL
  PUBLIC :: LIS_MOC_JULES_SMVCST
  PUBLIC :: LIS_MOC_JULES_SMVCWT
  PUBLIC :: LIS_MOC_JULES_FSAT
  PUBLIC :: LIS_MOC_JULES_FWETL 
  public :: LIS_MOC_JULES_ESOIL     
  
  integer :: LIS_MOC_JULES_STHZW = -9999
  integer :: LIS_MOC_JULES_STHU = -9999
  integer :: LIS_MOC_JULES_STHU_MIN = -9999
  integer :: LIS_MOC_JULES_STHF = -9999
  integer :: LIS_MOC_JULES_SMVCCL = -9999
  integer :: LIS_MOC_JULES_SMVCST = -9999
  integer :: LIS_MOC_JULES_SMVCWT = -9999
  integer :: LIS_MOC_JULES_FSAT  = -9999
  integer :: LIS_MOC_JULES_FWETL = -9999 
  integer :: LIS_MOC_JULES_ESOIL = -9999 

   ! ALMA ENERGY BALANCE COMPONENTS
  integer :: LIS_MOC_SWNET      = -9999
  integer :: LIS_MOC_LWNET      = -9999
  integer :: LIS_MOC_QLE        = -9999
  integer :: LIS_MOC_QH         = -9999
  integer :: LIS_MOC_QG         = -9999
  integer :: LIS_MOC_QF         = -9999
  integer :: LIS_MOC_QV         = -9999
  integer :: LIS_MOC_QTAU       = -9999
  integer :: LIS_MOC_QA         = -9999
  integer :: LIS_MOC_DELSURFHEAT = -9999
  integer :: LIS_MOC_DELCOLDCONT = -9999

   ! ALMA WATER BALANCE COMPONENTS
  integer :: LIS_MOC_SNOWF      = -9999
  integer :: LIS_MOC_RAINF      = -9999
  integer :: LIS_MOC_EVAP       = -9999
  integer :: LIS_MOC_QS         = -9999
  integer :: LIS_MOC_QREC       = -9999
  integer :: LIS_MOC_QSB        = -9999
  integer :: LIS_MOC_QSM        = -9999
  integer :: LIS_MOC_QFZ        = -9999
  integer :: LIS_MOC_QST        = -9999
  integer :: LIS_MOC_DELSOILMOIST = -9999
  integer :: LIS_MOC_DELSWE     = -9999
  integer :: LIS_MOC_DELSURFSTOR  = -9999
  integer :: LIS_MOC_DELINTERCEPT = -9999

   ! ALMA SURFACE STATE VARIABLES
  integer :: LIS_MOC_SNOWT      = -9999
  integer :: LIS_MOC_VEGT       = -9999
  integer :: LIS_MOC_BARESOILT  = -9999
  integer :: LIS_MOC_AVGSURFT   = -9999
  integer :: LIS_MOC_RADT       = -9999
  integer :: LIS_MOC_ALBEDO     = -9999
  integer :: LIS_MOC_ALBEDO_DIR_V     = -9999
  integer :: LIS_MOC_ALBEDO_DIF_V     = -9999
  integer :: LIS_MOC_ALBEDO_DIR_N     = -9999
  integer :: LIS_MOC_ALBEDO_DIF_N     = -9999
  integer :: LIS_MOC_SWE        = -9999
  integer :: LIS_MOC_SNOWDENSITY = -9999
  integer :: LIS_MOC_SNOWGRAIN  = -9999
  integer :: LIS_MOC_SWEVEG     = -9999 
  integer :: LIS_MOC_SURFSTOR   = -9999
   
   ! ALMA SUBSURFACE STATE VARIABLES
   integer :: LIS_MOC_SOILMOIST  = -9999
   integer :: LIS_MOC_SOILTEMP   = -9999
   integer :: LIS_MOC_SMLIQFRAC  = -9999
   integer :: LIS_MOC_SMFROZFRAC = -9999
   integer :: LIS_MOC_SOILWET    = -9999
   integer :: LIS_MOC_MATRICPOTENTIAL    = -9999

   ! ALMA EVAPORATION COMPONENTS
   integer :: LIS_MOC_POTEVAP    = -9999
   integer :: LIS_MOC_ECANOP     = -9999
   integer :: LIS_MOC_TVEG       = -9999
   integer :: LIS_MOC_ESOIL      = -9999
   integer :: LIS_MOC_EWATER     = -9999
   integer :: LIS_MOC_ROOTMOIST  = -9999
   integer :: LIS_MOC_CANOPINT   = -9999
   integer :: LIS_MOC_EVAPSNOW   = -9999
   integer :: LIS_MOC_SUBSNOW    = -9999
   integer :: LIS_MOC_SUBSURF    = -9999
   integer :: LIS_MOC_ACOND      = -9999

   ! ALMA OTHER HYDROLOGIC VARIABLES
  integer :: LIS_MOC_WATERTABLED= -9999
  integer :: LIS_MOC_TWS        = -9999
  integer :: LIS_MOC_GWS        = -9999

   ! ALMA COLD SEASON PROCESSES
  integer :: LIS_MOC_SNOWCOVER  = -9999
  integer :: LIS_MOC_SALBEDO    = -9999
  integer :: LIS_MOC_SNOWTPROF  = -9999
  integer :: LIS_MOC_SNOWDEPTH  = -9999
  integer :: LIS_MOC_SNOWTHICK  = -9999
  integer :: LIS_MOC_SLIQFRAC   = -9999
  integer :: LIS_MOC_SNOWTHRESH = -9999
  integer :: LIS_MOC_LAYERSNOWDEPTH = -9999
  integer :: LIS_MOC_LAYERSNOWDENSITY = -9999

   ! ALMA VARIABLES TO BE COMPARED WITH REMOTE SENSED DATA
   integer :: LIS_MOC_LWUP       = -9999

   ! ALMA CARBON VARIABLES
   integer :: LIS_MOC_GPP        = -9999
   integer :: LIS_MOC_NPP        = -9999
   integer :: LIS_MOC_NEE        = -9999
   integer :: LIS_MOC_AUTORESP   = -9999
   integer :: LIS_MOC_HETERORESP = -9999
   integer :: LIS_MOC_LEAFRESP   = -9999
   integer :: LIS_MOC_TOTSOILCARB= -9999
   integer :: LIS_MOC_TOTLIVBIOM = -9999

   ! ALMA FORCING VARIABLES
   integer :: LIS_MOC_WINDFORC   = -9999
   integer :: LIS_MOC_RAINFFORC  = -9999
   integer :: LIS_MOC_SNOWFFORC  = -9999
   integer :: LIS_MOC_CRAINFFORC  = -9999
   integer :: LIS_MOC_LSRAINFFORC = -9999
   integer :: LIS_MOC_CSNOWFFORC  = -9999
   integer :: LIS_MOC_LSSNOWFFORC = -9999
   integer :: LIS_MOC_TAIRFORC   = -9999
   integer :: LIS_MOC_QAIRFORC   = -9999
   integer :: LIS_MOC_PSURFFORC  = -9999
   integer :: LIS_MOC_SWDOWNFORC = -9999
   integer :: LIS_MOC_LWDOWNFORC = -9999

   ! CLSM FORCING VARIABLES
   integer :: LIS_MOC_PARDRFORC  = -9999
   integer :: LIS_MOC_PARDFFORC  = -9999

   ! PARAMETER OUTPUT - EXPERIMENTAL (USE W/WRF-WPS)
   integer :: LIS_MOC_LANDMASK   = -9999
   integer :: LIS_MOC_LANDCOVER  = -9999
   integer :: LIS_MOC_SOILTYPE   = -9999
   integer :: LIS_MOC_SANDFRAC   = -9999
   integer :: LIS_MOC_CLAYFRAC   = -9999
   integer :: LIS_MOC_SILTFRAC   = -9999
   integer :: LIS_MOC_POROSITY   = -9999
   integer :: LIS_MOC_SOILCOLOR  = -9999
   integer :: LIS_MOC_ELEVATION  = -9999
   integer :: LIS_MOC_SLOPE      = -9999
   integer :: LIS_MOC_LAI        = -9999
   integer :: LIS_MOC_SAI        = -9999
   integer :: LIS_MOC_SNFRALBEDO = -9999
   integer :: LIS_MOC_MXSNALBEDO = -9999
   integer :: LIS_MOC_GREENNESS  = -9999
   integer :: LIS_MOC_TEMPBOT   = -9999

   ! NLDAS OUTPUT
   integer :: LIS_MOC_CCOND    = -9999

   ! ADDITIONAL AFWA VARIABLES
   integer :: LIS_MOC_RELSMC       = -9999
   integer :: LIS_MOC_RHMIN        = -9999
   integer :: LIS_MOC_ROOTTEMP  = -9999
   integer :: LIS_MOC_TOTALPRECIP = -9999
   integer :: LIS_MOC_CRAINF  = -9999
   integer :: LIS_MOC_LSRAINF = -9999
   integer :: LIS_MOC_CSNOWF  = -9999
   integer :: LIS_MOC_LSSNOWF = -9999

   ! multivariate diagnostics (Bowen Ratio, Evaporative fraction)
   integer :: LIS_MOC_BR = -9999
   integer :: LIS_MOC_EF = -9999

   ! ADDITIONAL COUPLING FORCING VARIABLES
   integer :: LIS_MOC_DIRECTSWFORC  = -9999
   integer :: LIS_MOC_DIFFUSESWFORC = -9999
   integer :: LIS_MOC_NWINDFORC     = -9999
   integer :: LIS_MOC_EWINDFORC     = -9999
   integer :: LIS_MOC_FHEIGHTFORC   = -9999
   integer :: LIS_MOC_CHFORC        = -9999
   integer :: LIS_MOC_CMFORC        = -9999
   integer :: LIS_MOC_EMISSFORC     = -9999
   integer :: LIS_MOC_MIXRATIOFORC  = -9999
   integer :: LIS_MOC_COSZENFORC    = -9999
   integer :: LIS_MOC_ALBEDOFORC    = -9999

   ! ADDITIONAL Noah3.x variables
   integer :: LIS_MOC_SOILET  = -9999
   integer :: LIS_MOC_Z0BRD   = -9999
   integer :: LIS_MOC_ROUGHNESS   = -9999
   integer :: LIS_MOC_THERMAL_ROUGHNESS   = -9999

   !t2,q2 diagnostics
   integer :: LIS_MOC_T2DIAG = -9999
   integer :: LIS_MOC_Q2DIAG = -9999
   integer :: LIS_MOC_RNET = -9999
   integer :: LIS_MOC_CH     = -9999
   integer :: LIS_MOC_CM     = -9999
   integer :: LIS_MOC_MIXRATIO = -9999

!<for vic>
   !Additional VIC forcing variables
   integer :: LIS_MOC_SNOWFLAGFORC          = -9999
   integer :: LIS_MOC_DENSITYFORC           = -9999
   integer :: LIS_MOC_VAPORPRESSFORC        = -9999
   integer :: LIS_MOC_VAPORPRESSDEFICITFORC = -9999
   integer :: LIS_MOC_ARESIST               = -9999
!</for vic>

   integer :: LIS_MOC_SACUZTWC    = -9999
   integer :: LIS_MOC_SACUZFWC    = -9999
   integer :: LIS_MOC_SACLZTWC    = -9999
   integer :: LIS_MOC_SACLZFSC    = -9999
   integer :: LIS_MOC_SACLZFPC    = -9999
   integer :: LIS_MOC_SACADIMPC    = -9999
   integer :: LIS_MOC_SNOW17SWE    = -9999
   integer :: LIS_MOC_SNOW17LIQW    = -9999
   integer :: LIS_MOC_SNOW17NEGHS    = -9999
   integer :: LIS_MOC_SNOW17ACCMAX    = -9999
   integer :: LIS_MOC_SNOW17AEADJ    = -9999
   integer :: LIS_MOC_SNOW17RMLT    = -9999

!<for vic>
   integer :: LIS_MOC_VIC_PET_SATSOIL   = -9999
   integer :: LIS_MOC_VIC_PET_H2OSURF   = -9999
   integer :: LIS_MOC_VIC_PET_SHORT     = -9999
   integer :: LIS_MOC_VIC_PET_TALL      = -9999
   integer :: LIS_MOC_VIC_PET_NATVEG    = -9999
   integer :: LIS_MOC_VIC_PET_VEGNOCR   = -9999
!</for vic>

   !FLDAS
   integer :: LIS_MOC_PETFORC         = -9999
   integer :: LIS_MOC_REFETFORC       = -9999

   ! FLDAS-WRSI OUTPUTS LIST
   integer :: LIS_MOC_SOS = -9999
   integer :: LIS_MOC_WRSI = -9999
   integer :: LIS_MOC_KF2 = -9999
   integer :: LIS_MOC_SumWR = -9999
   integer :: LIS_MOC_SumET = -9999
   integer :: LIS_MOC_SWI = -9999
   integer :: LIS_MOC_SOSa = -9999
   integer :: LIS_MOC_TotalSurplusWater = -9999
   integer :: LIS_MOC_MaxSurplusWater = -9999
   integer :: LIS_MOC_TotalWaterDeficit = -9999
   integer :: LIS_MOC_MaxWaterDeficit = -9999
   integer :: LIS_MOC_TotalAETInitial = -9999
   integer :: LIS_MOC_TotalWRInitial = -9999
   integer :: LIS_MOC_TotalSurplusWaterInitial = -9999
   integer :: LIS_MOC_TotalWaterDeficitInitial = -9999
   integer :: LIS_MOC_TotalAETVeg = -9999
   integer :: LIS_MOC_TotalWRVeg = -9999
   integer :: LIS_MOC_TotalSurplusWaterVeg = -9999
   integer :: LIS_MOC_TotalWaterDeficitVeg = -9999
   integer :: LIS_MOC_TotalAETFlower = -9999
   integer :: LIS_MOC_TotalWRFlower = -9999
   integer :: LIS_MOC_TotalSurplusWaterFlower = -9999
   integer :: LIS_MOC_TotalWaterDeficitFlower = -9999
   integer :: LIS_MOC_TotalAETRipe = -9999
   integer :: LIS_MOC_TotalWRRipe = -9999
   integer :: LIS_MOC_TotalSurplusWaterRipe = -9999
   integer :: LIS_MOC_TotalWaterDeficitRipe = -9999
   integer :: LIS_MOC_PermWiltDate = -9999
   integer :: LIS_MOC_Wilting1 = -9999
   integer :: LIS_MOC_Wilting2 = -9999
   integer :: LIS_MOC_WRSIa = -9999
   integer :: LIS_MOC_growing_season = -9999
   integer :: LIS_MOC_WHC                      = -9999
   integer :: LIS_MOC_LGP                      = -9999
   integer :: LIS_MOC_WR_TimeStep              = -9999 ! SY
   integer :: LIS_MOC_AET_TimeStep             = -9999 ! SY
   integer :: LIS_MOC_WRSI_TimeStep            = -9999 ! SY
   integer :: LIS_MOC_SurplusWater_TimeStep    = -9999 ! SY

   !NLDAS
   integer :: LIS_MOC_CAPEFORC         = -9999

   !Variables related to streamflow routing
   integer :: LIS_MOC_STREAMFLOW = -9999
   integer :: LIS_MOC_RIVSTO = -9999
   integer :: LIS_MOC_RIVDPH = -9999
   integer :: LIS_MOC_RIVVEL = -9999
   integer :: LIS_MOC_FLDOUT = -9999
   integer :: LIS_MOC_FLDEVAP = -9999
   integer :: LIS_MOC_FLDSTO = -9999
   integer :: LIS_MOC_FLDDPH = -9999
   integer :: LIS_MOC_FLDVEL = -9999
   integer :: LIS_MOC_FLDFRC = -9999
   integer :: LIS_MOC_FLDARE = -9999
   integer :: LIS_MOC_SFCELV = -9999
   integer :: LIS_MOC_RNFSTO = -9999
   integer :: LIS_MOC_BSFSTO = -9999
   integer :: LIS_MOC_RNFDWI = -9999
   integer :: LIS_MOC_BSFDWI = -9999
   integer :: LIS_MOC_SURFWS = -9999

   integer :: LIS_MOC_EWAT = -9999
   integer :: LIS_MOC_EDIF = -9999
   !Variables related to RTMs
   integer :: LIS_MOC_RTM_EMISSIVITY = -9999
   integer :: LIS_MOC_RTM_TB = -9999
   integer :: LIS_MOC_RTM_SM = -9999

   integer :: LIS_MOC_IRRIGATEDWATER

   integer :: LIS_MOC_LSM_COUNT
   integer :: LIS_MOC_ROUTING_COUNT
   integer :: LIS_MOC_RTM_COUNT
   integer :: LIS_MOC_IRRIG_COUNT
   
   ! <- for FLAKE 2013->
   integer :: LIS_MOC_LAKE_T_SNOW    =   -9999
   integer :: LIS_MOC_LAKE_T_ICE =   -9999
   integer :: LIS_MOC_LAKE_T_MNW =   -9999
   integer :: LIS_MOC_LAKE_T_WML =   -9999
   integer :: LIS_MOC_LAKE_T_BOT =   -9999
   integer :: LIS_MOC_LAKE_T_B1  =   -9999
   integer :: LIS_MOC_LAKE_C_T   =   -9999
   integer :: LIS_MOC_LAKE_H_SNOW    =   -9999
   integer :: LIS_MOC_LAKE_H_ICE =   -9999
   integer :: LIS_MOC_LAKE_H_ML  =   -9999
   integer :: LIS_MOC_LAKE_H_B1  =   -9999
   integer :: LIS_MOC_LAKE_T_SFC =   -9999
   integer :: LIS_MOC_LAKE_ALBEDO_WATER  =   -9999
   integer :: LIS_MOC_LAKE_ALBEDO_ICE    =   -9999
   integer :: LIS_MOC_LAKE_ALBEDO_SNOW   =   -9999
   integer :: LIS_MOC_LAKE_UFR_A =   -9999
   integer :: LIS_MOC_LAKE_UFR_W =   -9999
   integer :: LIS_MOC_LAKE_WCONV =   -9999
   integer :: LIS_MOC_LAKE_Q_SE  =   -9999
   integer :: LIS_MOC_LAKE_Q_LA  =   -9999
   integer :: LIS_MOC_LAKE_I_W   =   -9999
   integer :: LIS_MOC_LAKE_Q_LWA =   -9999
   integer :: LIS_MOC_LAKE_Q_LWW =   -9999
   integer :: LIS_MOC_LAKE_Q_BOT =   -9999
   ! <- end for FLAKE 2013 ->


   ! <- SACHTET ->

   integer :: LIS_MOC_SACUZTWH = -9999
   integer :: LIS_MOC_SACUZFWH = -9999
   integer :: LIS_MOC_SACLZTWH = -9999
   integer :: LIS_MOC_SACLZFSH = -9999
   integer :: LIS_MOC_SACLZFPH = -9999
   integer :: LIS_MOC_SACSWINT = -9999
   integer :: LIS_MOC_SACTSINT = -9999
   integer :: LIS_MOC_SACSWHINT = -9999
   integer :: LIS_MOC_SACFROST = -9999

   ! <- NoahMP ->
    integer ::  LIS_MOC_CANOPY_TEMP = -9999
    integer ::  LIS_MOC_CANOPY_VP   = -9999
    integer ::  LIS_MOC_CANOPY_WF   = -9999
    integer ::  LIS_MOC_CANOPY_INTL   = -9999
    integer ::  LIS_MOC_GROUNDAVGT   = -9999
    integer ::  LIS_MOC_GROUNDVEGT   = -9999
    integer ::  LIS_MOC_SOWN_NLAYER   = -9999
    integer ::  LIS_MOC_SNOW_LBDFSS   = -9999
    integer ::  LIS_MOC_SOIL_LBDFSS   = -9999
    integer ::  LIS_MOC_SNOWICE   = -9999
    integer ::  LIS_MOC_SNOWLIQ   = -9999
    integer ::  LIS_MOC_WT_AQUI_SATSOIL    = -9999
    integer ::  LIS_MOC_LAKEWATER   = -9999
    integer ::  LIS_MOC_LEAFMASS    = -9999
    integer ::  LIS_MOC_ROOTMASS    = -9999
    integer ::  LIS_MOC_STEMMASS    = -9999
    integer ::  LIS_MOC_WOODMASS    = -9999
    integer ::  LIS_MOC_CARBON_DEEPSOIL   = -9999
    integer ::  LIS_MOC_CARBON_SHALLOWSOIL    = -9999
    integer ::  LIS_MOC_SNOWAGE   = -9999
    integer ::  LIS_MOC_BETWEENWATER    = -9999
    integer ::  LIS_MOC_QRECTOGW    = -9999
    integer ::  LIS_MOC_QRECFROMGW    = -9999
    integer ::  LIS_MOC_FSR   = -9999
    !integer ::  LIS_MOC_LWUP    = -9999
    integer ::  LIS_MOC_FCEV    = -9999
    integer ::  LIS_MOC_FGEV    = -9999
    integer ::  LIS_MOC_FCTR    = -9999
    integer ::  LIS_MOC_VEGE2MT   = -9999
    integer ::  LIS_MOC_BARE2MT   = -9999
    integer ::  LIS_MOC_VEGE2MQ2    = -9999
    integer ::  LIS_MOC_BARE2MQ2    = -9999
    integer ::  LIS_MOC_APAR    = -9999
    integer ::  LIS_MOC_PSCO2   = -9999
    integer ::  LIS_MOC_SAV   = -9999
    integer ::  LIS_MOC_SAG   = -9999
    integer ::  LIS_MOC_PONDING    = -9999
    integer ::  LIS_MOC_PONDING1   = -9999
    integer ::  LIS_MOC_PONDING2   = -9999
    integer ::  LIS_MOC_RSSUN   = -9999
    integer ::  LIS_MOC_RSSHA   = -9999
    integer ::  LIS_MOC_BGAP    = -9999
    integer ::  LIS_MOC_WGAP    = -9999
    integer ::  LIS_MOC_CHV   = -9999
    integer ::  LIS_MOC_CHB   = -9999
    integer ::  LIS_MOC_SHG   = -9999
    integer ::  LIS_MOC_SHC   = -9999
    integer ::  LIS_MOC_SHB   = -9999
    integer ::  LIS_MOC_EVG   = -9999
    integer ::  LIS_MOC_EVB   = -9999
    integer ::  LIS_MOC_GHV   = -9999
    integer ::  LIS_MOC_GHB   = -9999
    integer ::  LIS_MOC_IRV   = -9999
    integer ::  LIS_MOC_IRC   = -9999
    integer ::  LIS_MOC_IRB   = -9999
    integer ::  LIS_MOC_HTR    = -9999
    integer ::  LIS_MOC_HEVC   = -9999
    integer ::  LIS_MOC_CHLEAF    = -9999
    integer ::  LIS_MOC_CHUC    = -9999
    integer ::  LIS_MOC_CHV2    = -9999
    integer ::  LIS_MOC_CHB2    = -9999
    integer ::  LIS_MOC_FPICE   = -9999
!  <- end Noah MP  ->

!   <- RUC -> 
   integer :: LIS_MOC_QVG = -9999
   integer :: LIS_MOC_QCG = -9999
   integer :: LIS_MOC_QSG = -9999
   integer :: LIS_MOC_SNNOT75CM = -9999
   integer :: LIS_MOC_FRZPREC_DEN = -9999
   integer :: LIS_MOC_FRZPREC = -9999
   integer :: LIS_MOC_ACC_SNOWF = -9999
   integer :: LIS_MOC_ACC_FRZPREC = -9999
   integer :: LIS_MOC_ACC_EVAP = -9999
   integer :: LIS_MOC_QSFC = -9999
   integer :: LIS_MOC_DEW_FROST = -9999
   integer :: LIS_MOC_DRIP = -9999
   integer :: LIS_QH_SNOW = -9999 
!   <- end RUC ->  

    ! <- JULES -> 
    integer :: LIS_MOC_GS = -9999
    integer :: LIS_MOC_GC = -9999 

#if 0
   ! SPECIAL CASE INDICES
   ! These are required because Min/Max support cannot be generically
   ! handled for GRIB-1 output.  The routine writeSingleGrib1Var maps
   ! these two entries to LIS_MOC_TAIRFORC.
   ! They should not be counted in the LIS_MOC_COUNT total count.
   integer, parameter :: LIS_MOC_TAIRFORC_MIN = 112
   integer, parameter :: LIS_MOC_TAIRFORC_MAX = 113

   ! READ ABOVE NOTE ABOUT SPECIAL CASE INDICES
   integer, parameter :: LIS_MOC_COUNT      = 145
   ! Add the special cases.  LIS_MOC_GRIB_COUNT should be used only in
   ! LIS_gribMod.F90.
   integer, parameter :: LIS_MOC_GRIB_COUNT = 145

   TODO: implement support for the special case of LIS_MOC_TAIRFORC_MIN
         and LIS_MOC_TAIRFORC_MAX.
#endif

   real, parameter :: LIS_MOC_MAX_NUM =  999999.0
   real, parameter :: LIS_MOC_MIN_NUM = -999999.0
  
  type, public :: LIS_metadataEntry
     character(len=100) :: long_name
     character(len=100) :: standard_name
     character(len=50)  :: short_name
     character(len=20)  :: units
     integer       :: nunits
     character(len=20)  :: dir
     integer       :: ndirs
     character(len=1)   :: format          ! (scientific - E, else - F)
     integer       :: form
     integer       :: index           ! LIS_MOC_INDEX
     integer       :: vlevels
     integer       :: varId_def         
     integer       :: varId_opt1
     integer       :: varId_opt2
     integer :: varId_max ! EMK
     integer :: varId_min ! EMK
     integer       :: gribSF          ! GRIB scale factor
     integer       :: gribSfc         ! GRIB surface
     integer       :: gribLvl         ! GRIB level
     integer       :: gribDis         ! GRIB2 discipline
     integer       :: gribCat         ! GRIB2 category
     !integer       :: gribminId      ! GRIB id for minimum field
     !integer       :: gribmaxId      ! GRIB id for maximum field
     integer       :: timeAvgOpt
     integer       :: selectOpt
     integer       :: minMaxOpt
     integer       :: stdOpt
     real          :: valid_min
     real          :: valid_max
     character(len=20), allocatable :: unittypes(:)
     character(len=20), allocatable :: dirtypes(:)
     integer, allocatable :: count(:,:)
     integer              :: diagFlag
     real, allocatable :: minimum(:,:) ! ntiles, vlevels
     real, allocatable :: maximum(:,:) ! ntiles, vlevels
     real, allocatable :: modelOutput(:,:,:) !timeavg, ntiles, vlevels

     type(LIS_metadataEntry), pointer :: next
  end type LIS_metadataEntry

  ! To create an array of pointers, you must create a derived
  ! datatype containing the pointer.  Then you create
  ! an array of this datatype.
  type, public :: dep
     type(LIS_metadataEntry), pointer :: dataEntryPtr
  end type dep

  type, public :: output_meta
!BOC
     real*8  :: time         ! time to begin writing output
     integer :: syear        ! year to begin writing output
     integer :: smonth       ! month to begin writing output
     integer :: sday         ! day to begin writing output
     integer :: shour        ! hour to begin writing output
     integer :: smin         ! minutes to begin writing output
     integer :: ssec         ! seconds to begin writing output

     integer :: month       ! specific output writing time (month)
     integer :: day         ! specific output writing time (day)
     integer :: hour        ! specific output writing time (hour)
     integer :: min         ! specific output writing time (min)
     integer :: sec         ! specific output writing time (sec)

     type(LIS_metadataEntry), pointer :: head_lsm_list
     type(LIS_metadataEntry), pointer :: head_routing_list
     type(LIS_metadataEntry), pointer :: head_rtm_list
     type(LIS_metadataEntry), pointer :: head_irrig_list

     type(dep), allocatable, dimension(:) :: ptr_into_lsm_list
     type(dep), allocatable, dimension(:) :: ptr_into_routing_list
     type(dep), allocatable, dimension(:) :: ptr_into_rtm_list
     type(dep), allocatable, dimension(:) :: ptr_into_irrig_list
!EOC
  end type output_meta

  type(output_meta),     allocatable :: LIS_histData(:)

contains
 
!BOP
!  !ROUTINE: LIS_histDataInit
! \label{LIS_histDataInit}
! 
! !INTERFACE: 
  subroutine LIS_histDataInit(n, ntiles)
! !USES: 

    implicit none

! !ARGUMENTS: 
    integer,  intent(IN)   :: n 
    integer,  intent(IN)   :: ntiles
! 
! !DESCRIPTION: 
!  This routine initializes the required linked lists to hold the selected
!  list of LSM variables
!
!   The arguments are: 
!   \begin{description}
!    \item[n]  index of the nest \newline
!    \item[ntiles]  size of the tilespace \newline
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[get\_moc\_attributes] (\ref{get_moc_attributes}) \newline
!      adds history output objects to the linked lists and sets the
!      user-definable elements
!    \item[register\_dataEntry] (\ref{register_dataEntry}) \newline
!      completes the setting of the elements in the history output objects
!      and allocates the data structures related to the variable being output
!    \item[LIS\_resetOutputVars] (\ref{LIS_resetOutputVars}) \newline
!      resets the arrays storing the variable values. 
!   \end{description}
!EOP
    type(ESMF_Config) :: modelSpecConfig
    logical           :: file_exists
    integer           :: rc

    integer           :: grib_depthlvl
    integer           :: grib_snowlvl

    !hkb--GRIB2 specific Depth below land surface = 106
    !hkb--GRIB2 specific Snow level = 114
    if (LIS_rc%wout == "grib2") then
       grib_depthlvl = 106
       grib_snowlvl  = 114
    else
       grib_depthlvl = 112
       grib_snowlvl  = 112
    endif

    LIS_histData(n)%head_lsm_list     => null()
    LIS_histData(n)%head_routing_list => null()
    LIS_histData(n)%head_rtm_list     => null()
    LIS_histData(n)%head_irrig_list   => null()

    LIS_MOC_LSM_COUNT     = 0
    LIS_MOC_ROUTING_COUNT = 0
    LIS_MOC_RTM_COUNT     = 0
    LIS_MOC_IRRIG_COUNT   = 0 

    call ESMF_ConfigFindLabel(LIS_config,"Model output attributes file:", &
                              rc=rc)
    call LIS_verify(rc,'Model output attributes file: not specified')

    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%outputSpecFile(n),rc=rc)
    write(LIS_logunit,*) '[INFO] Opening Model Output Attributes File', &
                         LIS_rc%outputSpecFile(n)

    inquire(file=LIS_rc%outputSpecFile(n),exist=file_exists)
    if(.not.file_exists) then 
       write(LIS_logunit,*) '[ERR] Model output attributes file does not exist...'
       write(LIS_logunit,*) '[ERR] Program stopping... '
       call LIS_endrun()
    endif

    modelSpecConfig = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(modelSpecConfig,trim(LIS_rc%outputSpecFile(n)), &
         rc=rc)     
!-------------------------------------------------------------------------
! read the meta data attributes for each variable
!-------------------------------------------------------------------------
  !!! JULES

    call ESMF_ConfigFindLabel(modelSpecConfig,"sthu:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sthu",&
         "sthu",&
         "sthu",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_STHU,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    ! soil moist fraction in deep (water table) layer.
    call ESMF_ConfigFindLabel(modelSpecConfig,"sthzw:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sthzw",&
         "sthzw",&
         "sthzw",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_STHZW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    ! 
    call ESMF_ConfigFindLabel(modelSpecConfig,"fsat:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "fsat",&
         "fsat",&
         "fsat",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_FSAT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    ! 
    call ESMF_ConfigFindLabel(modelSpecConfig,"fwetl:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "fwetl",&
         "fwetl",&
         "fwetl",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_FWETL,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"JESoil:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "JESoil",&
         "evaporation_from_soil_water_storage_JULES_only",&
         "evaporation from soil water storage JULES only",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_ESOIL,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"sthf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sthf",&
         "sthf",&
         "sthf",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_STHF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif


    call ESMF_ConfigFindLabel(modelSpecConfig,"sthu_min:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sthu_min",&
         "sthu_min",&
         "sthu_min",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_STHU_MIN,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"smvccl:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "smvccl",&
         "smvccl",&
         "smvccl",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_SMVCCL,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"smvcst:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "smvcst",&
         "smvcst",&
         "smvcst",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_SMVCST,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"smvcwt:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "smvcwt",&
         "smvcwt",&
         "smvcwt",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_JULES_SMVCWT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    ! Arguments to register_dataEntry:
    ! LIS_MOC_INDEX, tail dataEntry, nest index, number of units, 
    ! number of tiles, units, stats form, 
    ! GRIB surface type, GRIB levels type

    ! AIX requires the strings in these array contrustors 
    ! to have the same length.  So I have padded the shorter string.
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"Swnet:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Swnet",&
         "surface_net_downward_shortwave_flux",&
         "net downward shortwave radiation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SWNET,&
            LIS_histData(n)%head_lsm_list,&
            n,nunits=1,ntiles=ntiles,unittypes=(/"W/m2"/),&
            ndirs=2,dirtypes=(/"UP","DN"/),&
            form=1,gribSFC=1,gribLvl=1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lwnet:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Lwnet",&
         "surface_net_downward_longwave_flux",&
         "net downward longwave radiation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LWNET,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Rnet:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Rnet",&
         "net_radiation_flux",&
         "total net radiation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RNET,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qle:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qle",&
         "surface_upward_latent_heat_flux",&
         "latent heat flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QLE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qh",&
         "surface_upward_sensible_heat_flux",&
         "sensible heat flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qg",&
         "downward_heat_flux_in_soil",&
         "soil heat flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qf",&
         "energy_of_fusion",&
         "energy of fusion",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"S2L","L2S"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qv:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qv",&
         "surface_snow_sublimation_heat_flux",&
         "energy of sublimation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QV,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"S2V","V2S"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qtau:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qtau",&
         "momentum_flux",&
         "momentum flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QTAU,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"N/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qa:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qa",&
         "advective_energy",&
         "advective energy",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QA,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelSurfHeat:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DelSurfHeat",&
         "change_in_heat_storage",&
         "change in heat storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELSURFHEAT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"J/m2"/),2,(/"INC","DEC"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelColdCont:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DelColdCont",&
         "change_in_cold_content",&
         "change in cold content",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELCOLDCONT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"J/m2"/),2,(/"INC","DEC"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"BR:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "BR",&
         "bowen_ratio",&
         "bowen ratio",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_BR,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EF:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EF",&
         "evaporative_fraction",&
         "evaporative fraction",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Snowf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Snowf",&
         "snowfall_rate",&
         "snowfall rate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Rainf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Rainf",&
         "precipitation_rate",&
         "precipitation rate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RAINF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CRainf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CRainf",&
         "convective_rainfall_rate",&
         "convective rainfall",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CRAINF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Evap:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Evap",&
         "total_evapotranspiration",&
         "total evapotranspiration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EVAP,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s","mm/hr ","W/m2  "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qs:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qs",&
         "surface_runoff_amount",&
         "surface runoff",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QS,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"IN ","OUT"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qsb:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qsb",&
         "subsurface_runoff_amount",&
         "subsurface runoff amount",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QSB,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
#if (defined AFWA_GRIB_CONFIGS)
! Hard-code baseflow surface and level to AFWA's specifications
! to make the LIS-7 output match the LIS-6 style. - dmm
            2,(/"IN ","OUT"/),2,112,200,&
#else
            2,(/"IN ","OUT"/),2,1,1,&
#endif
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qrec:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qrec",&
         "recharge",&
         "recharge",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QREC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),2,(/"IN ","OUT"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qsm:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qsm",&
         "snowmelt",&
         "snowmelt",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QSM,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"S2L","L2S"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qfz:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qfz",&
         "refreezing_of_water",&
         "refreezing water",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QFZ,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),2,(/"L2S","S2L"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qst:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qst",&
         "snow_throughfall",&
         "snow throughfall",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QST,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),1,("-"),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelSoilMoist:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DelSoilMoist",&
         "change_in_soil_moisture",&
         "change in soil moisture",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELSOILMOIST,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),2,(/"INC","DEC"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelSWE:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DelSWE",&
         "change_in_swe",&
         "change in snow water equivalent",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELSWE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),2,(/"INC","DEC"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelSurfStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "DelSurfStor",&
         "DelSurfStor",&
         "change in surface water storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELSURFSTOR,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2 ","kg/m2s"/),2,(/"INC","DEC"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DelIntercept:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "DelIntercept",&
         "change_in_interception_storage",&
         "change in interception storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DELINTERCEPT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),2,(/"INC","DEC"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowT",&
         "temperature_in_surface_snow",&
         "surface snow temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"VegT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "VegT",&
         "canopy_temperature",&
         "canopy temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VEGT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"BareSoilT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "BareSoilT",&
         "bare_soil_temperature",&
         "bare soil temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_BARESOILT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"AvgSurfT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AvgSurfT",&
         "surface_temperature",&
         "surface temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_AVGSURFT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"AvgGrndT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AvgGrndT",&
         "average_ground_surface_temperature",&
         "average ground surface temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GROUNDAVGT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"VegGrndT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "VegGrndT",&
         "vegetated_ground_surface_temperature",&
         "vegetated ground surface temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GROUNDVEGT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RadT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RadT",&
         "surface_radiative_temperature",&
         "surface radiative temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RADT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Albedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Albedo",&
         "surface_albedo",&
         "surface albedo",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDO,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"AlbDirVis:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AlbDirVis",&
         "surface_albedo_for_direct_beam_visible_light",&
         "surface albedo for direct beam visible light",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDO_DIR_V,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"AlbDifVis:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AlbDifVis",&
         "surface_albedo_for_diffuse_visible_light",&
         "surface albedo for diffuse visible light",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDO_DIF_V,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"AlbDirNIR:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AlbDirNIR",&
         "surface_albedo_for_direct_beam_NIR_light",&
         "surface albedo for direct beam NIR light",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDO_DIR_N,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"AlbDifNIR:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AlbDifNIR",&
         "surface_albedo_for_diffuse_NIR_light",&
         "surface albedo for diffuse NIR light",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDO_DIF_N,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    

    call ESMF_ConfigFindLabel(modelSpecConfig,"SWE:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SWE",&
         "liquid_water_content_of_surface_snow",&
         "snow water equivalent",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SWE,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2","m    "/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowDensity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowDensity",&
         "snow_density_for_each_layer",&
         "snow density for each layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWDENSITY,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m3"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
   

    call ESMF_ConfigFindLabel(modelSpecConfig,"LayerSnowDensity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LayerSnowDensity",&
         "snow_density_for_each_layer",&
         "snow_density_for_each_layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAYERSNOWDENSITY,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m3"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
 
    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowGrain:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowGrain",&
         "snow_grain_size_for_each_layer",&
         "snow grain size for each layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWGRAIN,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"micron"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowDepth:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowDepth",&
         "snow_depth",&
         "snow depth",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWDEPTH,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"m ", "cm", "mm"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
       ! cm is added for VIC, Shugong Wang 02/20/2012
    endif
    
    ! added by Shugong Wang 05/02/2018 for JULES 
    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowThick:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowThick",&
         "snow_thick",&
         "thickness of snow layers",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWTHICK,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LayerSnowDepth:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LayerSnowDepth",&
         "snow_depth_for_each_layer",&
         "snow_depth_for_each_layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAYERSNOWDEPTH,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"m ", "cm", "mm"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowIce:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowIce",&
         "snow_ice",&
         "snow ice",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWICE,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2","mm   "/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SWEVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SWEVeg",&
         "swe_intercepted_by_vegetation",&
         "swe intercepted by vegetation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SWEVEG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowAge:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowAge",&
         "Snow_Age",&
         "Snow Age",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SNOWAGE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    ! Added by Zhuo Wang on 11/11/2018
    call ESMF_ConfigFindLabel(modelSpecConfig,"Smcwtd:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
          "SMCWTD",&
          "smcwtd",&
          "SM content in the layer to the water table when deep",rc)
    if ( rc == 1 ) then
         call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_BETWEENWATER,&
              LIS_histData(n)%head_lsm_list,&
              n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,1,1,&
              model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SurfStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SurfStor",&
         "surface_water_storage",&
         "surface water storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SURFSTOR,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SoilMoist:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SoilMoist",&
         "soil_moisture_content",&
         "soil moisture content",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILMOIST,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2","m3/m3"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SoilTemp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SoilTemp",&
         "soil_temperature",&
         "soil temperature",rc)
    if ( rc == 1 ) then
!hkb-- GRIB2 specific Depth below land surface = 106
       if (LIS_rc%wout .eq. "grib2") then
        call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILTEMP,&
             LIS_histData(n)%head_lsm_list,&
             n,1,ntiles,(/"K"/),1,(/"-"/),1,grib_depthlvl,0,&
             model_patch=.true.)
       else
        call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILTEMP,&
             LIS_histData(n)%head_lsm_list,&
             n,1,ntiles,(/"K"/),1,(/"-"/),1,112,0,&
             model_patch=.true.)
       endif
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SmLiqFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SmLiqFrac",&
         "liquid_fraction_of_soil_moisture",&
         "average layer fraction of liquid moisture",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SMLIQFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-    ","m3/m3"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SmFrozFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SmFrozFrac",&
         "frozen_fraction_of_soil_moisture",&
         "average layer fraction of frozen moisture",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SMFROZFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-    ","m3/m3"/),1,(/"-"/),1,grib_depthlvl,0)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SoilWet:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SoilWet",&
         "total_soil_wetness",&
         "total soil wetness",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILWET,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"- ","% ","mm"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"MatricPotential:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "MatricPotential",&
         "soil_matric_potential",&
         "soil matric potential",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_MATRICPOTENTIAL,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"m ","mm"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"PotEvap:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "PotEvap",&
         "potential_evapotranspiration",&
         "potential evapotranspiration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_POTEVAP,&
            LIS_histData(n)%head_lsm_list,&
            n,4,ntiles,(/"kg/m2s","mm/hr ","W/m2  ","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"ECanop:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ECanop",&
         "interception_evaporation",&
         "interception evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ECANOP,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s","mm/hr ","W/m2  "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "TVeg",&
         "vegetation_transpiration",&
         "vegetation transpiration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TVEG,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s","mm/hr ","W/m2  "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"ESoil:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ESoil",&
         "bare_soil_evaporation",&
         "bare soil evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ESOIL,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s","mm/hr ","W/m2  "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EWater:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EWater",&
         "open_water_evaporation",&
         "open water evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EWATER,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"m3    ", "kg/m2s"/),2,&
            (/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RootMoist:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RootMoist",&
         "root_zone_soil_moisture",&
         "root zone soil moisture",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ROOTMOIST,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"m3/m3","kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CanopInt:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CanopInt",&
         "total_canopy_water_storage",&
         "total canopy water storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CANOPINT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EvapSnow:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EvapSnow",&
         "snow_evaporation",&
         "snow evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EVAPSNOW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SubSnow:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SubSnow",&
         "snow_sublimation",&
         "snow sublimation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SUBSNOW,&
            LIS_histData(n)%head_lsm_list,n,&
            4,ntiles,(/"kg/m2s","mm/hr ","W/m2  ","mm    "/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SubSurf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SubSurf",&
         "sublimation_of_the_snow_free_area",&
         "sublimation of the snow free area",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SUBSURF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"ACond:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ACond",&
         "aerodynamic_conductance",&
         "aerodynamic conductance",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ACOND,&
            LIS_histData(n)%head_lsm_list,n,&
            1,ntiles,(/"m/s"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CCond:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CCond",&
         "canopy_conductance",&
         "canopy conductance",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CCOND,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WaterTableD:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WaterTableD",&
         "water_table_depth",&
         "water table depth",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WATERTABLED,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TWS:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "TWS",&
         "terrestrial_water_storage",&
         "terrestrial water storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TWS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"GWS:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "GWS",&
         "ground_water_storage",&
         "ground water storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GWS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,grib_depthlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Snowcover:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Snowcover",&
         "surface_snow_area_fraction",&
         "snow cover",rc)
    if ( rc == 1 ) then
       ! EMK...Added percentage
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWCOVER,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SAlbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SAlbedo",&
         "snow_albedo",&
         "snow albedo",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SALBEDO,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowTProf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowTProf",&
         "snow_temperature_profile",&
         "snow temperature profile",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWTPROF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,grib_snowlvl,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SLiqFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SLiqFrac",&
         "snow_liquid_fraction_on_ground",&
         "snow liquid fraction on ground",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SLIQFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LWup:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LWup",&
         "longwave_radiation_up",&
         "longwave radiation up",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LWUP,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"UP"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"GPP:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "GPP",&
         "gross_primary_production",&
         "gross primary production",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GPP,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s2 ","umol/m2s", "g/m2s   "/),&
            2,(/"IN ","OUT"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"NPP:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "NPP",&
         "net_primary_productivity",&
         "net primary productivity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_NPP,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s2 ","umol/m2s", "g/m2s   "/),&
            2,(/"IN ","OUT"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"NEE:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "NEE",&
         "net_ecosystem_exchange",&
         "net ecosystem exchange",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_NEE,&
            LIS_histData(n)%head_lsm_list,&
            n,3,ntiles,(/"kg/m2s2 ","umol/m2s", "g/m2s   "/),&
            2,(/"IN ","OUT"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"AutoResp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "AutoResp",&
         "autoresp",&
         "",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_AUTORESP,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s2 ","umol/m2s"/),&
            2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"HeteroResp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "HeteroResp",&
         "heteroresp",&
         "",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_HETERORESP,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s2 ","umol/m2s"/),&
            2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LeafResp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LeafResp",&
         "leafresp",&
         "",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LEAFRESP,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s2 ","umol/m2s"/),&
            2,(/"UP","DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotSoilCarb:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "total_soil_and_litter_carbon_content",&
         "total soil and litter carbon content",&
         "",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TOTSOILCARB,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),&
            1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotLivBiom:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "total_living_biomass_carbon_content",&
         "total living biomass carbon content",&
         "",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TOTLIVBIOM,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),&
            1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SoilET:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SoilET",&
         "soil_evaporation",&
         "soil evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILET,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s", "W/m2  "/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Z0brd:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "z0brd",&
         "z0brd",&
         "z0brd",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_Z0BRD,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"Gs:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Gs",&
         "Gridbox_surface_conductance_to_evaporation",&
         "Gridbox surface conductance to evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Gc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Gc",&
         "tile_surface_conductance_to_evaporation",&
         "tile surface conductance to evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Ch:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Ch",&
         "heat_exchange_coefficient",&
         "heat exchange coefficient",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Cm:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Cm",&
         "momentum_exchange_coefficient",&
         "momentum exchange coefficient",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CM,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"T2diag:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "T2diag",&
         "diagnostic_t2",&
         "diagnostic t2",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_T2DIAG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Q2diag:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Q2diag",&
         "diagnostic_q2",&
         "diagnostic q2",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_Q2DIAG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RootTemp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RootTemp",&
         "root_zone_temperature",&
         "root zone temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ROOTTEMP,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Wind_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Wind_f",&
         "wind_speed",&
         "wind speed",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WINDFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
#if (defined AFWA_GRIB_CONFIGS)
! Hard-code wind surface and level to AFWA's specifications
! to make the LIS-7 output match the LIS-6 style. - dmm
            ntiles,(/"m/s   ","km/day"/),1,(/"-"/),1,105,10,&
#else
            ntiles,(/"m/s   ","km/day"/),1,(/"-"/),1,1,1,&
#endif
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Rainf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Rainf_f",&
         "rainfall_flux",&
         "rainfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RAINFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Snowf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Snowf_f",&
         "snowfall_flux",&
         "snowfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CRainf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CRainf_f",&
         "convective_rainfall_flux",&
         "convective rainfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CRAINFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"LSRainf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LSRainf_f",&
         "large_scale_rainfall_flux",&
         "large scale rainfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LSRAINFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"CSnowf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CSnowf_f",&
         "convective_snowfall_flux",&
         "convective snowfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CSNOWFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"LSSnowf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LSSnowf_f",&
         "large_scale_snowfall_flux",&
         "large scale snowfall flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LSSNOWFFORC,&
            LIS_histData(n)%head_lsm_list,n,2,&
            ntiles,(/"kg/m2s","kg/m2 "/),2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Tair_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Tair_f",&
         "air_temperature",&
         "air temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TAIRFORC,&
            LIS_histData(n)%head_lsm_list,&
#if (defined AFWA_GRIB_CONFIGS)
! Hard-code air temperature surface and level to AFWA's specifications
! to make the LIS-7 output match the LIS-6 style. - dmm
            n,1,ntiles,(/"K"/),1,(/"-"/),1,105,2,&
#else
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
#endif
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Qair_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qair_f",&
         "specific_humidity",&
         "specific humidity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QAIRFORC,&
            LIS_histData(n)%head_lsm_list,&
#if (defined AFWA_GRIB_CONFIGS)
! Hard-code specific humidity surface and level to AFWA's specifications
! to make the LIS-7 output match the LIS-6 style. - dmm
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),2,105,2,&
#else
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),2,1,1,&
#endif
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Psurf_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Psurf_f",&
         "surface_air_pressure",&
         "surface pressure",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_PSURFFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"Pa"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SWdown_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SWdown_f",&
         "surface_downwelling_shortwave_flux_in_air",&
         "surface downward shortwave radiation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SWDOWNFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LWdown_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "LWdown_f",&
         "surface_downwelling_longwave_flux_in_air",&
         "surface downward longwave radiation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LWDOWNFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),2,(/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DirectSW_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DirectSW_f",&
         "surface_direct_downwelling_shortwave_flux_in_air",&
         "direct shortwave flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DIRECTSWFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"DiffuseSW_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DiffuseSW_f",&
         "surface_diffuse_downwelling_shortwave_flux_in_air",&
         "diffuse shortwave flux", rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DIFFUSESWFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"NWind_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "NWind_f",&
         "northward_wind",&
         "northward wind",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_NWINDFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),2,(/"E", "N"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EWind_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EWind_f",&
         "eastward_wind",&
         "eastward wind",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EWINDFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),2,(/"E", "N"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FHeight_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "FHeight_f",&
         "forcing_height",&
         "forcing height",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_FHEIGHTFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Ch_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Ch_f",&
         "heat_exchange_coefficient",&
         "heat exchange coefficient",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CHFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Cm_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Cm_f",&
         "momentum_exchange_coefficient",&
         "momentum exchange coefficient",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CMFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"MixRatio_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MixRatio_f",&
         "mixing_ratio",&
         "mixing ratio",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_MIXRATIOFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CosZenith_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CosZenith_f",&
         "cosine_solar_zenith_angle",&
         "cosine of solar zenith angle",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_COSZENFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Albedo_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Albedo_f",&
         "albedo",&
         "albedo",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ALBEDOFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"PARDR_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "PARDR_f",&
         "surface_downward_PAR_direct",&
         "surface downward PAR direct",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_PARDRFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"PARDF_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "PARDF_f",&
         "surface_downward_PAR_diffuse",&
         "surface downward PAR diffuse",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_PARDFFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Landmask:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Landmask",&
         "landmask",&
         "landmask",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LANDMASK,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Landcover:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Landcover",&
         "landcover",&
         "landcover",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LANDCOVER,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Soiltype:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Soiltype",&
         "soiltype",&
         "soiltype",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILTYPE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SandFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SandFrac",&
         "sandfrac",&
         "fraction of sand",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SANDFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"ClayFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ClayFrac",&
         "clayfrac",&
         "fraction of clay",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CLAYFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SiltFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SiltFrac",&
         "siltfrac",&
         "fraction of silt",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SILTFRAC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Porosity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Porosity",&
         "porosity",&
         "porosity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_POROSITY,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Soilcolor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SoilColor",&
         "soilcolor",&
         "soil color",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOILCOLOR,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Elevation:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Elevation",&
         "elevation",&
         "elevation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ELEVATION,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Slope:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Slope",&
         "slope",&
         "slope",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SLOPE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LAI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LAI",&
         "leaf_area_index",&
         "leaf area index",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAI,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SAI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SAI",&
         "stem_area_index",&
         "stem area index",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SAI,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Snfralbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Snfralbedo",&
         "snow_free_albedo",&
         "snow free albedo",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNFRALBEDO,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Mxsnalbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Mxsnalbedo",&
         "maximum_snow_free_albedo",&
         "maximum snow free albedo",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_MXSNALBEDO,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Greenness:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Greenness",&
         "green_vegetation_fraction",&
         "green vegetation fraction",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_GREENNESS,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Tempbot:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Tempbot",&
         "bottom_temperature",&
         "bottom temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TEMPBOT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"PET_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "PET_f",&
         "potential_evaporation",&
         "potential evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_PETFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RefET_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RefET_f",&
         "reference_ET",&
         "reference ET",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_REFETFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"CAPE_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CAPE_f",&
         "convective_available_potential_energy",&
         "Convective Available Potential Energy",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CAPEFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"J/kg"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SOS:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SOS",&
         "start_of_season",&
         "start of season",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WRSI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WRSI",&
         "water_requirements_satisfaction_index",&
         "water requirements satisfaction index",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WRSI,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"KF2:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "KF2",&
         "KF2",&
         "KF2",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_KF2,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SumWR:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SumWR",&
         "SumWR",&
         "SumWR",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SumWR,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SumET:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SumET",&
         "SumET",&
         "SumET",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SumET,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SWI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SWI",&
         "soil_water_index",&
         "soil water index",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SWI,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"%"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SOSa:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SOSa",&
         "SOSa",&
         "SOSa",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SOSa,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWater:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TotalSurplusWater",&
         "total_surplus_water",&
         "total surplus water", rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalSurplusWater,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"MaxSurplusWater:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MaxSurplusWater",&
         "max_surplus_water",&
         "max surplus water", rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_MaxSurplusWater,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficit:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TotalWaterDeficit",&
         "total_water_deficit",&
         "total water deficit", rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWaterDeficit,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"MaxWaterDeficit:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MaxWaterDeficit",&
         "max_water_deficit",&
         "max water deficit",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_MaxWaterDeficit,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETInitial:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TotalAETInitial",&
         "total_AET_initial",&
         "total AET initial",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalAETInitial,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRInitial:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TotalWRInitial",&
         "total_WR_initial",&
         "total WR initial",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWRInitial,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterInitial:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TotalSurplusWaterInitial",&
         "total_surplus_water_initial",&
         "total surplus water initial", rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalSurplusWaterInitial,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitInitial:",&
         rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list, &
         "TotalWaterDeficitInitial",&
         "total_water_deficit_initial",&
         "total water deficit initial",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWaterDeficitInitial,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "TotalAETVeg",&
         "total_AET_veg",&
         "total AET veg",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalAETVeg,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "TotalWRVeg",&
         "total_WR_veg",&
         "total WR veg",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWRVeg,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalSurplusWaterVeg",&
         "total_surplus_water_veg",&
         "total surplus water veg",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalSurplusWaterVeg,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitVeg:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalWaterDeficitVeg",&
         "total_water_deficit_veg",&
         "total water deficit veg",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWaterDeficitVeg,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETFlower:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalAETFlower",&
         "totalAETFlower",&
         "totalAETFlower",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalAETFlower,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRFlower:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalWRFlower",&
         "TotalWRFlower",&
         "TotalWRFlower",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWRFlower,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterFlower:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list, &
         "TotalSurplusWaterFlower",&
         "TotalSurplusWaterFlower",&
         "TotalSurplusWaterFlower",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalSurplusWaterFlower,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitFlower:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list, &
         "TotalWaterDeficitFlower",&
         "TotalWaterDeficitFlower",&
         "TotalWaterDeficitFlower",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWaterDeficitFlower,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalAETRipe:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalAETRipe",&
         "TotalAETRipe",&
         "TotalAETRipe",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalAETRipe,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWRRipe:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalWRRipe",&
         "TotalWRRipe",&
         "TotalWRRipe",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWRRipe,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalSurplusWaterRipe:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "TotalSurplusWaterRipe",&
         "TotalSurplusWaterRipe",&
         "TotalSurplusWaterRipe",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalSurplusWaterRipe,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalWaterDeficitRipe:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list, &
         "TotalWaterDeficitRipe",&
         "TotalWaterDeficitRipe",&
         "TotalWaterDeficitRipe",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TotalWaterDeficitRipe,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"PermWiltDate:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "PermWiltDate",&
         "PermWiltDate",&
         "PermWiltDate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_PermWiltDate,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Wilting1:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Wilting1",&
         "Wilting1",&
         "Wilting1",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_Wilting1,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Wilting2:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Wilting2",&
         "Wilting2",&
         "Wilting2",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_Wilting2,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WRSIa:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WRSIa",&
         "WRSIa",&
         "WRSIa",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WRSIa,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"growing_season:",rc=rc)
    call get_moc_attributes(modelSpecConfig,LIS_histData(n)%head_lsm_list,&
         "growing_season",&
         "growing_season",&
         "growing season",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_growing_season,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WHC:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "WHC",&
         "WHC",&
         "WHC",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WHC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"LGP:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "LGP",&
         "LGP",&
         "LGP",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LGP,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WR_TimeStep:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "WR_TimeStep",&
         "WR_TimeStep",&
         "WR_TimeStep",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WR_TimeStep,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"AET_TimeStep:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "AET_TimeStep",&
         "AET_TimeStep",&
         "AET_TimeStep",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_AET_TimeStep,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"WRSI_TimeStep:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "WRSI_TimeStep",&
         "WRSI_TimeStep",&
         "WRSI_TimeStep",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_WRSI_TimeStep,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SurplusWater_TimeStep:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "SurplusWater_TimeStep",&
         "SurplusWater_TimeStep",&
         "SurplusWater_TimeStep",rc)
    if (rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SurplusWater_TimeStep,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif

    !<for vic>
    call ESMF_ConfigFindLabel(modelSpecConfig,"Snowflag_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Snowflag",&
         "Snowflag",&
         "Snowflag",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOWFLAGFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Density_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Atmospheric_Density",&
         "Atmospheric_Density",&
         "Atmospheric Density",rc)
    if ( rc == 1 ) then
       ! the unit was '-' for air density. Chaged to 'kg/m3' by Shugong Wang on 02/17/2012     
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DENSITYFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m3"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"VaporPress_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Vapor_Pressure",&
         "Vapor_Pressure",&
         "Vapor Pressure",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VAPORPRESSFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"VaporPressDeficit_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Vapor_Pressure_Deficit",&
         "Vapor_Pressure_Deficit",&
         "Vapor Pressure Deficit",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VAPORPRESSDEFICITFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"AResist:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Aerodynamic_Resistance",&
         "Aerodynamic_Resistance",&
         "Aerodynamic Resistance",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ARESIST,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"s/m"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_tsint:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sac_tsint",&
         "sac_tsint",&
         "sac soil temperature of intended layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACTSINT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_swint:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sac_swint",&
         "sac_swint",&
         "sac total volumetric soil moisture content of intented layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACSWINT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_swhint:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "sac_swhint",&
         "sac_swhint",&
         "sac liquid volumetric soil moisture content of intented layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACSWHINT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_frost:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_frost",&
         "sac_frost",&
         "sac_frost",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACFROST,&
         LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),2,1,1)
    endif

    
    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uztwc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_uztwc",&
         "sac_uztwc",&
         "sac_uztwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACUZTWC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uzfwc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_uzfwc",&
         "sac_uzfwc",&
         "sac_uzfwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACUZFWC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lztwc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lztwc",&
         "sac_lztwc",&
         "sac_lztwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZTWC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfsc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lzfsc",&
         "sac_lzfsc",&
         "sac_lzfsc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZFSC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfpc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lzfpc",&
         "sac_lzfpc",&
         "sac_lzfpc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZFPC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_adimpc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_adimpc",&
         "sac_adimpc",&
         "sac_adimpc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACADIMPC,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"mm", "- "/),1,(/"-"/),2,1,1)
    endif


    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uztwh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_uztwc",&
         "sac_uztwc",&
         "sac_uztwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACUZTWH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_uzfwh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_uzfwc",&
         "sac_uzfwc",&
         "sac_uzfwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACUZFWH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lztwh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lztwc",&
         "sac_lztwc",&
         "sac_lztwc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZTWH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfsh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lzfsc",&
         "sac_lzfsc",&
         "sac_lzfsc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZFSH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"sac_lzfph:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "sac_lzfpc",&
         "sac_lzfpc",&
         "sac_lzfpc",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SACLZFPH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_swe:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_swe",&
         "snow17_swe",&
         "snow17_swe",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17SWE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_aeadj:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_aeadj",&
         "snow17_aeadj",&
         "snow17_aeadj",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17AEADJ,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_neghs:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_neghs",&
         "snow17_neghs",&
         "snow17_neghs",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17NEGHS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_liqw:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_liqw",&
         "snow17_liqw",&
         "snow17_liqw",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17LIQW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_accmax:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_accmax",&
         "snow17_accmax",&
         "snow17_accmax",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17ACCMAX,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),2,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"snow17_rmlt:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "snow17_rmlt",&
         "snow17_rmlt",&
         "snow17_rmlt",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNOW17RMLT,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),2,1,1)
    endif
!<for vic>

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_satsoil:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_satsoil",&
         "vic_pet_satsoil",&
         "vic_pet_satsoil",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_SATSOIL,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_h2osurf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_h2osurf",&
         "vic_pet_h2osurf",&
         "vic_pet_h2osurf",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_H2OSURF,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_short:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_short",&
         "vic_pet_short",&
         "vic_pet_short",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_SHORT,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_tall:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_tall",&
         "vic_pet_tall",&
         "vic_pet_tall",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_TALL,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_natveg:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_natveg",&
         "vic_pet_natveg",&
         "vic_pet_natveg",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_NATVEG,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"vic_pet_vegnocr:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "vic_pet_vegnocr",&
         "vic_pet_vegnocr",&
         "vic_pet_vegnocr",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_VIC_PET_VEGNOCR,&
                               LIS_histData(n)%head_lsm_list,&
                               n,2,ntiles,(/"kg/m2 ","kg/m2s"/),1,(/"-"/),&
                               2,1,1,model_patch=.true.)
    endif
!</for vic>

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tsnow:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTsnow",&
         "Lake_temperature_at_air_snow_interface",&
         "Lake temperature at the air snow interface",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_SNOW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tice:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTice",&
         "Lake_temperature_at_snow_ice_interface",&
         "Lake temperature at snow ice interface",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_ICE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tmnw:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTmnw",&
         "Mean_temperature_of_the_water_column",&
         "Mean temperature of the water column",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_MNW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Twml:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTwml",&
         "Lake_temperature_of_the_mixed_layer",&
         "Lake temperature of the mixed layer",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_WML,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tbot:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTbot",&
         "Lake_temperature_at_the_water_bottom",&
         "Lake temperature at the water bottom",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_BOT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Tb1:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeTb1",&
         "Temperature_at_the_bottom_of_upper_layer_of_sediments",&
         "Temperature at the bottom of upper layer of sediments",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_T_B1,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_CT:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeCT",&
         "Thermocline_shape_factor_of_lake",&
         "Thermocline shape factor of lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_C_T,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hice:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeHice",&
         "Ice_thickness_above_lake",&
         "Ice thickness above lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_H_ICE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hml:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeHml",&
         "Thickness_of_mixed_layer_of_lake",&
         "Thickness of mixed layer of lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_H_ML,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Hb1:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeHb1",&
         "Thickness_of_upper_layer_of_bottom_sediments",&
         "Thickness of upper layer of bottom sediments",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_H_B1,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Walbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeWalbedo",&
         "Water_surface_albedo_over_lake",&
         "Water surface albedo over lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_ALBEDO_WATER,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_IceAlbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeIceAlbedo",&
         "Ice_surface_albedo_over_lake",&
         "Ice surface albedo over lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_ALBEDO_ICE,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_SnowAlbedo:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeSnowAlbedo",&
         "Snow_surface_albedo_over_lake",&
         "Snow surface albedo over lake",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_ALBEDO_SNOW,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif
    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_UFRa:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeUFRa",&
         "Lake_friction_velocity_in_air",&
         "Lake friction velocity in air",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_UFR_A,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_UFRw:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeUFRw",&
         "Lake_friction_velocity_in_surface_water",&
         "Lake friction velocity in surface water",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_UFR_W,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_WConv:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeWConv",&
         "Lake_convective_velocity_scale",&
         "Lake convective velocity scale",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_WCONV,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m/s"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_IW:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeIW",&
         "Lake_radiation_flux_at_the_interface",&
         "Lake radiation flux at the interface",rc)

    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_I_W,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Lake_Qbot:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LakeQbot",&
         "Lake_heat_flux_across_water_sediment_boundary",&
         "Lake heat flux across water sediment boundary",rc)

    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LAKE_Q_BOT,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RelSMC:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RelSMC",&
         "relative_soil_moisture",&
         "relative soil moisture",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RELSMC,&
            LIS_histData(n)%head_lsm_list,n,&
            2,ntiles,(/"-","%"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RHMin:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RHMin",&
         "min_relative_humidity",&
         "minimum relative humidity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_RHMIN,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"-","%"/),1,(/"-"/),1,105,2,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"TotalPrecip:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "TotalPrecip",&
         "total_precipitation_amount",&
         "total precipitation amount",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_TOTALPRECIP,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Emiss_f:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Emiss_f",&
         "emissivity",&
         "emissivity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_EMISSFORC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Roughness:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Roughness",&
         "surface_roughness_length",&
         "surface roughness",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_ROUGHNESS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    ! JULES 
    call ESMF_ConfigFindLabel(modelSpecConfig,"ThermalRoughness:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ThermalRoughness",&
         "surface_thermal_roughness_length",&
         "surface thermal roughness",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_THERMAL_ROUGHNESS,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,112,0,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"LSRainf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LSRainf",&
         "large_scale_rainfall_rate",&
         "large scale rainfall rate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LSRAINF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"LSSnowf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LSSnowf",&
         "large_scale_snowfall_rate",&
         "large scale snowfall rate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LSSNOWF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif
    
    call ESMF_ConfigFindLabel(modelSpecConfig,"CSnowf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CSnowf",&
         "convective_snowfall_rate",&
         "convective snowfall rate",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_CSNOWF,&
            LIS_histData(n)%head_lsm_list,&
            n,2,ntiles,(/"kg/m2s","kg/m2 "/),&
            2,(/"UP","DN"/),2,1,1,&
            model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Streamflow:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "Streamflow",&
         "streamflow",&
         "streamflow",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_STREAMFLOW,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m3/s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RiverStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "RiverStor",&
         "River_Water_Storage",&
         "River Water Storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_RIVSTO,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m3"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RiverDepth:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "RiverDepth",&
         "River_Depth",&
         "River Depth",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_RIVDPH,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RiverVelocity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "RiverFlowVelocity",&
         "River_Flow_Velocity",&
         "River Flow Velocity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_RIVVEL,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodQ:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodQ",&
         "Floodplain_Water_Discharge",&
         "Floodplain Water Discharge",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldout,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m3/s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodEvap:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodEvap",&
         "Floodplain_evaporation",&
         "Floodplain evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldevap,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodStor",&
         "Floodplain_Water_Storage",&
         "Floodplain Water Storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldsto,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m3"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodDepth:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodDepth",&
         "Floodplain_Depth",&
         "Floodplain Depth",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_flddph,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodVelocity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodVelocity",&
         "Floodplain Flow Velocity",&
         "Floodplain_Flow_Velocity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldvel,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m/s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodedFrac:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodedFrac",&
         "Flooded Fraction",&
         "Flooded_Fraction",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldfrc,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"FloodedArea:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "FloodedArea",&
         "Flooded Area",&
         "Flooded_Area",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_fldare,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m2"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SurfElev:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "SurfElev",&
         "Surface Water Elevation",&
         "Surface_Water_Elevation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_sfcelv,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RunoffStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "RunoffStor",&
         "Runoff Reservoir Storage",&
         "Runoff_Reservoir_Storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_rnfsto,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"BaseflowStor:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "BaseflowStor",&
         "Baseflow Reservoir Storage",&
         "Baseflow_Reservoir_Storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_bsfsto,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif


    call ESMF_ConfigFindLabel(modelSpecConfig,"RunoffDWI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "RunoffDWI",&
         "Runoff Deep Water Infiltration",&
         "Runoff_Deep_Water_Infiltration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_RNFDWI,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"BaseflowDWI:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "BaseflowDWI",&
         "Baseflow Deep Water Infiltration",&
         "Baseflow_Deep_Water_Infiltration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_BSFDWI,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"SWS:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "SWS",&
         "Surface Water Storage",&
         "Surface_Water_Storage",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_SURFWS,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"mm"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EvapWater:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "EvapWater",&
         "Evaporation_open_water",&
         "Evaporation_open_water",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_ewat,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"EvapDif:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_routing_list,&
         "EvapDif",&
         "Differential_evaporation",&
         "Differential_evaporation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_ROUTING_COUNT,LIS_MOC_edif,&
            LIS_histData(n)%head_routing_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),1,1,1,model_patch=.true.)
    endif


    call ESMF_ConfigFindLabel(modelSpecConfig,"RTM emissivity:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_rtm_list,&
         "RTM_emissivity",&
         "rtm_emissivity",&
         "rtm emissivity",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_RTM_COUNT,LIS_MOC_RTM_EMISSIVITY, &
            LIS_histData(n)%head_rtm_list,&
            n,1,ntiles,(/"-"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RTM Tb:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_rtm_list,&
         "RTM_Tb",&
         "rtm_brightness_temperature",&
         "rtm brightness temperature",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_RTM_COUNT,LIS_MOC_RTM_TB,&
            LIS_histData(n)%head_rtm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"RTM SoilMoist:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_rtm_list,&
         "RTM_SoilMoist",&
         "RTM_SoilMoist",&
         "RTM SoilMoist",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_RTM_COUNT,LIS_MOC_RTM_SM,&
            LIS_histData(n)%head_rtm_list,&
            n,1,ntiles,(/"m3/m3"/),1,(/"-"/),1,1,1)
    endif

    call ESMF_ConfigFindLabel(modelSpecConfig,"Irrigated water:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_irrig_list, &
         "IrrigatedWater",&
         "total_irrigated_water_amount",&
         "irrigated water amount",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_IRRIG_COUNT,LIS_MOC_IRRIGATEDWATER,&
            LIS_histData(n)%head_irrig_list,&
            n,1,ntiles,(/"kg/m2s"/),&
            1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif

    !<- NoahMP ->
!    call ESMF_ConfigFindLabel(modelSpecConfig,"LwUP:",rc=rc)
!    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
!         "LwUP",&
!         "total_net_longwave_radiation_to_atmosphere",&
!         "total net longwave radiation to atmosphere",rc)
!    if ( rc == 1 ) then
!       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_LWUP,&
!            LIS_histData(n)%head_lsm_list,&
!            n,1,ntiles,(/"W/m2"/),2,(/"UP","DN"/),1,1,1,&
!            model_patch=.true.)
!    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "VegCanopT:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "VegCanopT", &
         "canopy_air_temperature",  &
         "canopy air temperature",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CANOPY_TEMP, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"K"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "CanopVP:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CanopVP", &
         "canopy_air_vapor_pressure",  &
         "canopy air vapor pressure",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CANOPY_VP, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"Pa"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "CanopIntLiq:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CanopIntLiq", &
         "intercepted_liquid_water",  &
         "intercepted liquid water",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CANOPY_INTL, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"mm"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "CanopWetFrac:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "CanopWetFrac", &
         "canopy_wet_fraction",  &
         "canopy wet fraction",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CANOPY_WF, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"-"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    ! Added by Zhuo Wang on 11/11/2018
    Call ESMF_ConfigFindLabel(modelSpecConfig, "Wslake:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Wslake",  &
         "Wslake",  &
         "lake water storage",rc)
    if ( rc == 1 ) then
         call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_LAKEWATER, &
              LIS_histData(n)%head_lsm_list,&
              n, 1, ntiles,(/"mm"/), 1, (/"-"/),1,1,1,&
              model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "ActSnowNL:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ActSnowNL", &
         "actual_number_of_snow_layers",  &
         "actual number of snow layers",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SOWN_NLAYER, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"-"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "z_snow:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "z_snow", &
         "snow_layer-bottom_depth_from_snow_surface",  &
         "snow layer-bottom depth from snow surface",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SNOW_LBDFSS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "z_soil:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "z_soil", &
         "soil_layer-bottom_depth_from_snow_surface",  &
         "soil layer-bottom depth from snow surface",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SOIL_LBDFSS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "SnowLiq:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowLiq", &
         "snow-layer_liquid_water",  &
         "snow-layer liquid water",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SNOWLIQ, &
            LIS_histData(n)%head_lsm_list,&
            n, 2, ntiles,(/"kg/m2","mm   "/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "WT:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WT", &
         "water_in_aquifer_and_saturated_soil",  &
         "water in aquifer and saturated soil",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_WT_AQUI_SATSOIL, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"mm"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "LeafMass:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "LeafMass", &
         "leaf_mass",  &
         "leaf mass",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_LEAFMASS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "RootMass:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RootMass", &
         "mass_of_fine_roots",  &
         "mass of fine roots",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_ROOTMASS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    Call ESMF_ConfigFindLabel(modelSpecConfig, "StemMass:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "StemMass", &
         "mass_of_wood_stem",  &
         "mass of wood stem", rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_STEMMASS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "WoodMass:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WoodMass", &
         "mass_of_wood_including_woody_roots",  &
         "mass of wood including woody roots",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_WOODMASS, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "DeepSoilCarbon:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "DeepSoilCarbon", &
         "stable_carbon_in_deep_soil",  &
         "stable carbon in deep soil",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CARBON_DEEPSOIL, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "ShallowSoilCarbon:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ShallowSoilCarbon", &
         "short-lived_carbon_in_shallow_soil",  &
         "short-lived carbon in shallow soil",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CARBON_SHALLOWSOIL, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"g/m2"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "RechToGW:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RechToGW", &
         "recharge_to_the_water_table_when_groundwater_is_deep",  &
         "recharge to the water table when groundwater is deep",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_QRECTOGW, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "RechFromGW:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RechFromGW", &
         "recharge_from_the_water_table_when_groundwater_is_shallow",  &
         "recharge from the water table when groundwater is shallow",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_QRECFROMGW, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "SwReflect:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SwReflect", &
         "total_reflected_solar_radiation",   &
         "total reflected solar radiation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_FSR, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "fcev:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "fcev", &
         "canopy_evaporative_heat_to_atmosphere",   &
         "canopy evaporative heat to atmosphere",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_FCEV, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "fctr:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "fctr", &
         "transpiration_heat_to_atmosphere",   &
         "transpiration heat to atmosphere",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_FCTR, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "VegT2m:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "VegT2m", &
         "2-m_air_temperature_over_vegetated_part",  &
         "2-m air temperature over vegetated part",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_VEGE2MT, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"K"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "QairT2m:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "QairT2m", &
         "2-m_specific_humidity_over_vegetation",  &
         "2-m specific humidity over vegetation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_VEGE2MQ2, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"kg/kg"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "APAR:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "APAR", &
         "absorbed_photosynthesis_active_energy_by_canopy",   &
         "absorbed photosynthesis active radiation energy by canopy",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_APAR, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"IN ", "OUT"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "PSCO2:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "PSCO2", &
         "total_photosynthesis_of_CO2",   &
         "total photosynthesis of CO2",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_PSCO2, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"umol/m2s"/), 2, (/"IN ", "OUT"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "SAV:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SAV", &
         "solar_radiation_absorbed_by_vegetation",   &
         "solar radiation absorbed by vegetation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SAV, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"IN ", "OUT"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "SAG:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SAG", &
         "solar_radiation_absorbed_by_ground",   &
         "solar radiation absorbed by ground",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SAG, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"IN ", "OUT"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "ponding:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ponding", &
         "surface_ponding",  &
         "surface ponding",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_PONDING, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"mm"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "ponding2:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ponding2", &
         "surface_ponding2",  &
         "surface ponding2",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_PONDING2, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"mm"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "RsShaded:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RsShaded", &
         "shaded_stomatal_resistance",  &
         "shaded stomatal resistance",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_RSSHA, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"s/m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "RsSunlit:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RsSunlit", &
         "sunlit_stomatal_resistance",  &
         "sunlit stomatal resistance",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_RSSUN, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"s/m"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    Call ESMF_ConfigFindLabel(modelSpecConfig, "WCanoGapFrac:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "WCanoGapFrac", &
         "within-canopy_gap_fraction_for_beam",  &
         "within-canopy gap fraction for beam",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_WGAP, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"-"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "BCanoGapFrac:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "BCanoGapFrac", &
         "between-canopy_gap_fraction_for_beam",  &
         "between-canopy gap fraction for beam",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_BGAP, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"-"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    call esmf_configfindlabel(modelspecconfig, "ChVeg:", rc = rc)
    call get_moc_attributes(modelspecconfig, lis_histdata(n)%head_lsm_list, &
         "ChVeg", &
         "sensible_heat_exchange_coefficient_over_vegetated_fraction",  &
         "sensible heat exchange coefficient over vegetated fraction",rc)
    if ( rc == 1 ) then
        call register_dataentry(lis_moc_lsm_count, LIS_MOC_CHV, &
            lis_histdata(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    call esmf_configfindlabel(modelspecconfig, "ChBare:", rc = rc)
    call get_moc_attributes(modelspecconfig, lis_histdata(n)%head_lsm_list, &
         "ChBare", &
         "sensible_heat_exchange_coefficient_over_bare-ground_fraction",  &
         "sensible heat exchange coefficient over bare-ground fraction",rc)
    if ( rc == 1 ) then
        call register_dataentry(lis_moc_lsm_count, LIS_MOC_CHB, &
            lis_histdata(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "QhBare:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "QhBare", &
         "bare_ground_sensible_heat",   &
         "bare ground sensible heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SHB, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "QhGrnd:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "QhGrnd", &
         "ground_sensible_heat",   &
         "ground sensible heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SHG, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "QhCano:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "QhCano", &
         "canopy_sensible_heat",   &
         "canopy sensible heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SHC, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "EvapHBare:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EvapBare", &
         "bare_ground_evaporation_heat",   &
         "bare ground evaporation heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_EVB, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "EvapHGrnd:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "EvapHGrnd", &
         "ground_evaporation_heat",   &
         "ground evaporation heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_EVG, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "GrndHBare:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "GrndHBare", &
         "bare_ground_heat_flux",   &
         "bare ground heat flux",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_GHB, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "GrndHVeg:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "GrndHVeg", &
         "vegetated_ground_heat_flux",   &
         "vegetated ground heat flux",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_GHV, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "IRC:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "IRC", &
         "canopy_net_long_wave_radiation",   &
         "canopy net long wave radiation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_IRC, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "IRV:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "IRV", &
         "vegetated_ground_net_long_wave_radiation",   &
         "vegetated ground net long wave radiation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_IRV, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "IRB:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "IRB", &
         "bare_ground_net_long_wave_radiation",   &
         "bare ground net long wave radiation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_IRB, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "HeatTR:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "HeatTR", &
         "transpiration_heat",   &
         "transpiration heat",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_HTR, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "HeatEVC:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "HeatEVC", &
         "canopy_evaporation_heat",   &
         "canopy evaporation heat", rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_HEVC, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"W/m2"/), 2, (/"UP", "DN"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "ChLeaf:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ChLeaf", &
         "leaf_exchange_coefficient",  &
         "leaf exchange coefficient",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CHLEAF, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    Call ESMF_ConfigFindLabel(modelSpecConfig, "ChUC:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ChUC", &
         "under_canopy_exchange_coefficient",  &
         "under canopy exchange coefficient",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CHUC, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "ChV2:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ChV2", &
         "sensible_heat_exchange_coefficient_over_vegetated_fraction",  &
         "sensible heat exchange coefficient over vegetated fraction",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CHV2, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "ChB2:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "ChB2", &
         "sensible_heat_exchange_coefficient_over_bare_ground",  &
         "sensible heat exchange coefficient over bare ground",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_CHB2, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"m/s"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif


    Call ESMF_ConfigFindLabel(modelSpecConfig, "fpice:", rc = rc)
    Call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "fpice", &
         "snow_fraction_in_precipitation",  &
         "snow fraction in precipitation",rc)
    if ( rc == 1 ) then
        call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_FPICE, &
            LIS_histData(n)%head_lsm_list,&
            n, 1, ntiles,(/"-"/), 1, (/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    !<- end NoahMP ->

    !<- RUC addition ->
    ! RUC  density of frozen precipitation (kg m{-3})
    call ESMF_ConfigFindLabel(modelSpecConfig,"Density_FrzRain:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Density_FrzRain",&
         "frozen_rain_density",&
         "frozen rain density",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_FRZPREC_DEN,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m3"/),1,(/"-"/),1,1,1)
    endif

    ! RUC   time-step frozen precipitation (kg m{-2})
    call ESMF_ConfigFindLabel(modelSpecConfig,"TimeStep_FrzRain:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "TimeStep_FrzRain",&
         "Time_Step_Frozen_Rain",&
         "time step frozen rain",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_FRZPREC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1)
    endif


    !  RUC: effective cloud water mixing ratio at the surface 
    call ESMF_ConfigFindLabel(modelSpecConfig,"MixRatio_QCG:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MixRatio_QCG",&
         "effective_cloud_water_mixing_ratio_at_the_surface",&
         "effective cloud water mixing ratio at the surface",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QCG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    !  RUC:  effective mixing ratio at the surface ( kg kg{-1} )
    call ESMF_ConfigFindLabel(modelSpecConfig,"MixRatio_QVG:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MixRatio_QVG",&
         "effective_mixing_ratio_at_the_surface",&
         "effective mixing ratio at the surface",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QVG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    ! RUC:  surface water vapor mixing ratio at satration (kg/kg) 
    call ESMF_ConfigFindLabel(modelSpecConfig,"MixRatio_QSG:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "MixRatio_QSG",&
         "surface_water_vapor_mixing_ratio_at_satration",&
         "surface water vapor mixing ratio at satration",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QCG,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

   ! RUC  specific humidity at the surface [kg/kg]
    call ESMF_ConfigFindLabel(modelSpecConfig,"Qair_sfc:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "Qair_sfc",&
         "specific_humidity_at_the_surface",&
         "specific humidity at the surface",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_QSFC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/kg"/),1,(/"-"/),2,1,1,&
            model_patch=.true.)
    endif


    ! RUC  snow temperature at 7.5 cm depth (k)
    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowT7.5cm:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowT7.5cm",&
         "snow_temperature_at_7.5_cm_depth",&
         "snow temperature at 7.5 cm depth",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_SNNOT75CM,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"K"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    ! RUC  dewfall (or frostfall for t<273.15) ( m )
    call ESMF_ConfigFindLabel(modelSpecConfig,"Dew_Frost:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Dew_Frost",&
         "dew_or_frost",&
         "dew or frost",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_MOC_DEW_FROST,&
         LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),2,1,1)
    endif

    ! RUC  throughfall of precipitation from canopy (kg m{-2} s{-1})
    call ESMF_ConfigFindLabel(modelSpecConfig,"Throughfall_Prcp:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Throughfall_Prcp",&
         "throughfall_of_precipitation_from_canopy",&
         "throughfall of precipitation from canopy",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_DRIP,&
         LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2s"/),1,(/"-"/),2,1,1)
    endif

    ! RUC  snow heat flux (w/m^2: negative, if downward from surface)
    call ESMF_ConfigFindLabel(modelSpecConfig,"Qh_Snow:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list,&
         "Qh_Snow",&
         "snow_heat_flux",&
         "snow heat flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT,LIS_QH_SNOW,&
         LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"W/m2"/),1,(/"-"/),2,1,1)
    endif

    ! RUC  run total snowfall accumulation (kg m{-2})
    call ESMF_ConfigFindLabel(modelSpecConfig,"RunTotal_Snowf:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RunTotal_Snowf",&
         "Run_Total_Snow_Fall",&
         "run total snow fall",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_ACC_SNOWF,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    ! RUC  run total frozen precipitation accumulation (kg m{-2})
    call ESMF_ConfigFindLabel(modelSpecConfig,"RunTotal_FrzRain:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RunTotal_FrzRain",&
         "Run_Total_Frozen_Precipitation",&
         "run total frozen precipitation",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_ACC_FRZPREC,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    ! RUC  run total evaporation flux accumulation (kg m{-2})
    call ESMF_ConfigFindLabel(modelSpecConfig,"RunTotal_Evap:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "RunTotal_Evap",&
         "Run_Total_Evaporation_Flux",&
         "run total evaporation flux",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_ACC_EVAP,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"kg/m2"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif

    ! RUC  snow depth threshold ( m )
    call ESMF_ConfigFindLabel(modelSpecConfig,"SnowThresh:",rc=rc)
    call get_moc_attributes(modelSpecConfig, LIS_histData(n)%head_lsm_list, &
         "SnowThresh",&
         "Snow_depth_threshold",&
         "snow depth threshold",rc)
    if ( rc == 1 ) then
       call register_dataEntry(LIS_MOC_LSM_COUNT, LIS_MOC_SNOWTHRESH,&
            LIS_histData(n)%head_lsm_list,&
            n,1,ntiles,(/"m"/),1,(/"-"/),1,1,1,&
            model_patch=.true.)
    endif
    
    !<- end of RUC addition -> 
    
    call ESMF_ConfigDestroy(modelSpecConfig,rc=rc)

    allocate(LIS_histData(n)%ptr_into_lsm_list(LIS_MOC_LSM_COUNT))
    allocate(LIS_histData(n)%ptr_into_routing_list(LIS_MOC_ROUTING_COUNT))
    allocate(LIS_histData(n)%ptr_into_rtm_list(LIS_MOC_RTM_COUNT))
    allocate(LIS_histData(n)%ptr_into_irrig_list(LIS_MOC_IRRIG_COUNT))

    call set_ptr_into_list(LIS_MOC_LSM_COUNT, &
         LIS_histData(n)%head_lsm_list, &
         LIS_histData(n)%ptr_into_lsm_list)

    call set_ptr_into_list(LIS_MOC_ROUTING_COUNT, &
         LIS_histData(n)%head_routing_list, &
         LIS_histData(n)%ptr_into_routing_list)

    call set_ptr_into_list(LIS_MOC_RTM_COUNT, &
         LIS_histData(n)%head_rtm_list, &
         LIS_histData(n)%ptr_into_rtm_list)

    call set_ptr_into_list(LIS_MOC_IRRIG_COUNT, &
         LIS_histData(n)%head_irrig_list, &
         LIS_histData(n)%ptr_into_irrig_list)

    call LIS_resetOutputVars(n,1) !for LSM
    call LIS_resetOutputVars(n,2) !for ROUTING
    call LIS_resetOutputVars(n,3) !for RTM
    call LIS_resetOutputVars(n,4) !for Irrigation


end subroutine LIS_histDataInit

!BOP
!
! !ROUTINE: get_moc_attributes
!  \label{get_moc_attributes}
!
! !INTERFACE:
subroutine get_moc_attributes(modelSpecConfig, head_dataEntry, short_name, &
                              standard_name, long_name, status)
!
! !DESCRIPTION:
! This routine reads the model output configuration attributes for
! the given history output variable.
!
! For each user-selected history output variable, this routine allocates
! the history output object for the given history output linked list, and it
! then assigns the user-definable attributes.
!
!EOP

! !USES:
   use ESMF 

   implicit none

! !ARGUMENTS:
! \begin{description}
!    \item[modelSpecConfig]
!       ESMF handle to the Model Output Configuration file
!    \item[head\_dataEntry]
!       head of the history output linked list to configure
!    \item[short\_name]
!       short name for the output variable (will be used in the
!       stats file)
!    \item[standard\_name]
!       CF-style name for the output variable
!    \item[long\_name]
!       descriptive long name for the output variable
!    \item[status]
!       input: flag indicating if the entry for the output variable \newline
!              exists in the configuration (0-yes, 1-no)
!       output: flag indicating whether the output variable was selected \newline
!               in the configuration (0-no, 1-yes)
! \end{description}
   type(ESMF_Config), intent(inout)     :: modelSpecConfig
   type(LIS_metadataEntry), pointer, intent(out) :: head_dataEntry
   character(len=*), intent(in)         :: short_name
   character(len=*), intent(in)         :: standard_name
   character(len=*), intent(in)         :: long_name
   integer, intent(inout)               :: status
   
   integer                              :: selectOpt
   character(len=20)                    :: cfunit
   integer                              :: rc

   type(LIS_metadataEntry), pointer     :: current, dataEntry
   
   ! status is an in/out argument.
   !
   ! Its input value is set by a call to ESMF_ConfigFindLabel, and it
   ! indicates whether a label was found in the modelSpecConfig.
   !
   ! Its output value indicates whether the output variable associated
   ! with the label was selected for output. 0 = no; 1 = yes.

   if ( status == ESMF_SUCCESS ) then ! found label in modelSpecConfig

      call ESMF_ConfigGetAttribute(modelSpecConfig,selectOpt,&
           default=0,rc=rc)

      if ( selectOpt == 1 ) then

         status = 1 ! output variable was selected

         allocate(dataEntry)
         if ( .not. associated(head_dataEntry) ) then
            head_dataEntry => dataEntry
         else
            current => head_dataEntry
            do while ( associated(current%next) )
               current => current%next
            enddo
            current%next => dataEntry
         endif
         dataEntry%next => null()

         dataEntry%selectOpt = selectOpt

         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%units,rc=rc)
         call convertToCFunits(dataEntry%units, cfunit)

         dataEntry%units = trim(cfunit)

         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%dir,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%timeAvgOpt,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%minMaxOpt,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%stdOpt,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%vlevels,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%varId_def,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%gribSF,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%gribDis,rc=rc)
         call ESMF_ConfigGetAttribute(modelSpecConfig,dataEntry%gribCat,rc=rc)
      
         dataEntry%short_name=trim(short_name)
         dataEntry%standard_name=trim(standard_name)
         dataEntry%long_name=trim(long_name)
      else
         status = 0 ! output variable was not selected
      endif
   else
      status = 0 ! label was not found in modelSpecConfig, therefore
                 ! the output variable was not selected
   endif
   
end subroutine get_moc_attributes


!BOP
! 
! !ROUTINE: allocate_dataEntry
! \label{allocate_dataEntry}
! 
! !INTERFACE: 
  subroutine allocate_dataEntry(dataEntry,nunits,ntiles,unittypes,ndirs,&
       dirtypes, model_patch)

    implicit none
! !ARGUMENTS: 
    type(LIS_metadataEntry), pointer :: dataEntry
    integer                 :: nunits
    integer                 :: ntiles
    character(len=*)        :: unittypes(nunits)
    integer                 :: ndirs
    character(len=*)        :: dirtypes(ndirs)
    logical                 :: model_patch
! 
! !DESCRIPTION: 
!  This routine initializes the datastructures required for the 
!  specified output variable. 
!EOP
    integer                 :: i 
    character(len=20)       :: cfunit

    if(dataEntry%selectOpt.ne.0) then 
       if(dataEntry%timeAvgOpt.eq.2) then 
          allocate(dataEntry%modelOutput(2,ntiles,dataEntry%vlevels))
       else
          allocate(dataEntry%modelOutput(1,ntiles,dataEntry%vlevels))
       endif
       allocate(dataEntry%count(ntiles,dataEntry%vlevels))
       dataEntry%modelOutput = 0 
       dataEntry%count = 0
       dataEntry%diagFlag = 0 

       dataEntry%nunits = nunits
       allocate(dataEntry%unittypes(nunits))       

       do i=1,nunits
          call convertToCFunits(unittypes(i),cfunit)
          dataEntry%unittypes(i) = cfunit
       enddo

       dataEntry%ndirs = ndirs
       allocate(dataEntry%dirtypes(ndirs))       
       do i=1,ndirs
          dataEntry%dirtypes(i) = dirtypes(i)
       enddo

       if(dataEntry%minMaxOpt.ne.0) then 
          allocate(dataEntry%minimum(ntiles,dataEntry%vlevels))
          allocate(dataEntry%maximum(ntiles,dataEntry%vlevels))
          ! Initialize the minimux and maximum fields to implausible values.
          dataEntry%minimum = LIS_MOC_MAX_NUM
          dataEntry%maximum = LIS_MOC_MIN_NUM
       endif
    endif
  end subroutine allocate_dataEntry

!BOP
! 
! !ROUTINE: convertToCFunits
! \label{convertToCFunits}
!
! !INTERFACE: 
  subroutine convertToCFunits(unit,cfunit)
!
! !DESCRIPTION: 
!   This routine converts the LIS specified units to CF compliant
!   unit specifications. 
! 
!EOP

    implicit none

    character(len=*)  :: unit
    character(len=*), intent(out)  :: cfunit
    cfunit = ""
    if(unit.eq."W/m2") then 
       cfunit = "W m-2"
    elseif(unit.eq."J/m2") then 
       cfunit = "J m-2"
    elseif(unit.eq."kg/m2") then 
       cfunit = "kg m-2"
    elseif(unit.eq."kg/m2s") then 
       cfunit = "kg m-2 s-1"
    elseif(unit.eq."kg/m2s2") then 
       cfunit = "kg m-2 s-2"
    elseif(unit.eq."kg/m3") then 
       cfunit = "kg m-3"
    elseif(unit.eq."m3/m3") then 
       cfunit = "m^3 m-3"
    elseif(unit.eq."m/s") then 
       cfunit = "m s-1"
    elseif(unit.eq."s/m") then 
       cfunit = "s m-1"
    elseif(unit.eq."m3/s") then 
       cfunit = "m3 s-1"
    elseif(unit.eq."N/m2") then 
       cfunit = "N m-2"
    elseif(unit.eq."g/m2") then 
       cfunit = "g m-2"
    elseif(unit.eq."g/m2s") then 
       cfunit = "g m-2 s-1"
    elseif(unit.eq."kg/kg") then 
       cfunit = "kg kg-1"
    elseif(unit.eq."g/g") then 
       cfunit = "g g-1"
    elseif(unit.eq."umol/m2s") then 
       cfunit = "umol m-2 s-1"
    elseif(unit.eq."mm/hr") then 
       cfunit = "mm hr-1"
    elseif(unit.eq."km/day") then 
       cfunit = "km day-1"
    elseif(unit.eq."J/kg") then 
       cfunit = "J kg-1"
    else
       cfunit = unit
    endif
  end subroutine convertToCFunits

!BOP
! !ROUTINE: LIS_diagnoseSurfaceOutputVar
! \label{LIS_diagnoseSurfaceOutputVar}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseSurfaceOutputVar(n, t, index, vlevel, value, unit,&
                                          direction, valid_min, valid_max, &
                                          surface_type)
    use  LIS_coreMod, only : LIS_domain
    implicit none
! !ARGUMENTS:
    integer, intent(in)           :: n
    integer, intent(in)           :: t
    integer, intent(in)           :: index
    integer, intent(in)           :: vlevel
    real   , intent(in)           :: value
    character(len=*), intent(in)  :: unit
    character(len=*)              :: direction
    real,    intent(in), optional :: valid_min
    real,    intent(in), optional :: valid_max
    integer, intent(in), optional :: surface_type
! 
! !DESCRIPTION: 
!  This routine is the user callable interface to the LIS\_diagnoseOutputVar
!  routine.
!
!  See LIS\_diagnoseOutputVar for more details.
!EOP    
    integer                     :: tid
    logical                     :: model_patch

    if(present(surface_type)) then 
       tid = LIS_surface(n,surface_type)%tile(t)%tile_id
       model_patch = .true. 
    else
       tid = t
       model_patch = .false. 
    endif
    call LIS_diagnoseOutputVar(LIS_histData(n)%head_lsm_list,   &
         LIS_MOC_LSM_COUNT, LIS_histData(n)%ptr_into_lsm_list,&
         n, tid, index, vlevel, value, unit,&
         direction,valid_min,valid_max,model_patch)
    
  end subroutine LIS_diagnoseSurfaceOutputVar

!BOP
! !ROUTINE: LIS_diagnoseRoutingOutputVar
! \label{LIS_diagnoseRoutingOutputVar}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseRoutingOutputVar(n, t, index, vlevel, value, unit,&
                                          direction,valid_min,valid_max)
    implicit none
! !ARGUMENTS:
    integer, intent(in)           :: n    
    integer, intent(in)           :: t
    integer, intent(in)           :: index
    integer, intent(in)           :: vlevel
    real   , intent(in)           :: value
    character(len=*), intent(in)  :: unit
    character(len=*)              :: direction
    real,    intent(in), optional :: valid_min
    real,    intent(in), optional :: valid_max
! 
! !DESCRIPTION: 
!  This routine is the user callable interface to the LIS\_diagnoseOutputVar
!  routine.
!
!  See LIS\_diagnoseOutputVar for more details.
!EOP    
!index 1 is hardcoded since we know that routing models are only over land. 
    integer         :: tid
    logical         :: model_patch

!    tid = LIS_surface(n,1)%tile(t)%tile_id
    model_patch = .true.
    call LIS_diagnoseOutputVar(LIS_histData(n)%head_routing_list, &
         LIS_MOC_ROUTING_COUNT, LIS_histData(n)%ptr_into_routing_list,&
         n, t, index, vlevel, value, unit,  &
         direction,valid_min,valid_max,model_patch)
    
  end subroutine LIS_diagnoseRoutingOutputVar

!BOP
! !ROUTINE: LIS_diagnoseRTMOutputVar
! \label{LIS_diagnoseRTMOutputVar}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseRTMOutputVar(n, t, index, vlevel, value, unit,&
                                      direction,valid_min,valid_max)
    use  LIS_coreMod, only : LIS_domain
    implicit none
! !ARGUMENTS:
    integer, intent(in)           :: n    
    integer, intent(in)           :: t
    integer, intent(in)           :: index
    integer, intent(in)           :: vlevel
    real   , intent(in)           :: value
    character(len=*), intent(in)  :: unit
    character(len=*)              :: direction
    real,    intent(in), optional :: valid_min
    real,    intent(in), optional :: valid_max
! 
! !DESCRIPTION: 
!  This routine is the user callable interface to the LIS\_diagnoseOutputVar
!  routine.
!
!  See LIS\_diagnoseOutputVar for more details.
!EOP    
  call LIS_diagnoseOutputVar(LIS_histData(n)%head_rtm_list, &
       LIS_MOC_RTM_COUNT, LIS_histData(n)%ptr_into_rtm_list,&
       n, t, index, vlevel, value, unit,  &
       direction,valid_min,valid_max)
    
  end subroutine LIS_diagnoseRTMOutputVar

!BOP
! !ROUTINE: LIS_diagnoseIrrigationOutputVar
! \label{LIS_diagnoseIrrigationOutputVar}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseIrrigationOutputVar(n, t, index, vlevel, value, unit,&
                                      direction,valid_min,valid_max)
    use  LIS_coreMod, only : LIS_domain
    implicit none
! !ARGUMENTS:
    integer, intent(in)           :: n    
    integer, intent(in)           :: t
    integer, intent(in)           :: index
    integer, intent(in)           :: vlevel
    real   , intent(in)           :: value
    character(len=*), intent(in)  :: unit
    character(len=*)              :: direction
    real,    intent(in), optional :: valid_min
    real,    intent(in), optional :: valid_max
! 
! !DESCRIPTION: 
!  This routine is the user callable interface to the LIS\_diagnoseOutputVar
!  routine.
!
!  See LIS\_diagnoseOutputVar for more details.
!EOP    
    call LIS_diagnoseOutputVar(LIS_histData(n)%head_irrig_list, &
         LIS_MOC_IRRIG_COUNT, LIS_histData(n)%ptr_into_irrig_list,&
         n, t, index, vlevel, value, unit,  &
         direction,valid_min,valid_max)
    
end subroutine LIS_diagnoseIrrigationOutputVar

!BOP
! !ROUTINE: LIS_diagnoseOutputVar
! \label{LIS_diagnoseOutputVar}
! 
! !INTERFACE: 
  subroutine LIS_diagnoseOutputVar(head_dataEntry,&
                                   count, ptr_into_list,&
                                   n, t, index, vlevel, value, unit,&
                                   direction,valid_min,valid_max,&
                                   model_patch)
    use  LIS_coreMod, only : LIS_domain
    use  LIS_logMod, only : LIS_logunit, LIS_endrun
    implicit none
! !ARGUMENTS:
    type(LIS_metadataEntry), pointer, intent(in) :: head_dataEntry
    integer, intent(in)           :: count
    type(dep), dimension(count), intent(in) :: ptr_into_list
    integer, intent(in)           :: n    
    integer, intent(in)           :: t
    integer, intent(in)           :: index
    integer, intent(in)           :: vlevel
    real   , intent(in)           :: value
    character(len=*), intent(in)  :: unit
    character(len=*)              :: direction
    real,    intent(in), optional :: valid_min
    real,    intent(in), optional :: valid_max
    logical, intent(in), optional :: model_patch
! 
! !DESCRIPTION: 
!  This routine maps the LSM, routing, or RTM specific variable to the 
!  appropriate variable in the given history output linked list, and it also
!  enables time averaging of the variable, if specified. 
! 
!   The arguments are: 
!   \begin{description}
!   \item[head\_dataEntry]  head of the given history linked list \newline
!   \item[n]  index of the nest \newline
!   \item[t]  index of the specified variable in the tile space \newline
!   \item[index]  the index of the specified variable in the generic (ALMA)
!                 variable list \newline
!   \item[vlevel] vertical dimension of the variable, if any \newline
!   \item[value]  value of the specified variable \newline
!   \end{description}
!EOP    
    real    :: vmin, vmax
    integer :: gindex
    type(LIS_metadataEntry), pointer :: dataEntry
    logical :: mpatch

    if ( index /= -9999 ) then

       if(PRESENT(valid_min)) then 
          vmin = valid_min
       else
          vmin = -1.0E+15
       endif
       
       if(PRESENT(valid_max)) then 
          vmax = valid_max
       else
          vmax = 1.0E+15
       endif
       
       if(PRESENT(model_patch)) then 
          mpatch = model_patch
       else
          mpatch = .false.
       endif

!#if 0
!    dataEntry => head_dataEntry
!
!    do while ( associated(dataEntry) )
!       if ( dataEntry%index == index ) then
!          exit
!       endif
!       dataEntry => dataEntry%next
!    enddo
!
!    if ( associated(dataEntry) ) then
!       if ( dataEntry%selectOpt /= 0 ) then
!          gindex = LIS_domain(n)%tile(t)%index
!          call diagnoseDataEntry(n,dataEntry,                               &
!                                 t,vlevel,value,unit,direction,vmin,vmax, &
!                                 LIS_domain(n)%grid(gindex)%ntiles,       &
!                                 LIS_domain(n)%grid(gindex)%subgrid_tiles,&
!                                 mpatch)
!       endif
!    endif
!#else
       dataEntry => ptr_into_list(index)%dataEntryPtr
       gindex = LIS_domain(n)%tile(t)%index
       call diagnoseDataEntry(n,dataEntry,                             &
            t,vlevel,value,unit,direction,vmin,vmax, &
            LIS_domain(n)%grid(gindex)%ntiles,       &
            LIS_domain(n)%grid(gindex)%subgrid_tiles,&
            mpatch)
!#endif

    endif
    
  end subroutine LIS_diagnoseOutputVar


!BOP
! 
! !ROUTINE: diagnoseDataEntry
! \label{diagnoseDataEntry}
! 
! !INTERFACE:
  subroutine diagnoseDataEntry(n,dataentry, t, vlevel,in_value, unit, &
                               direction, vmin, vmax, nsiblings, &
                               siblings,model_patch)
! !USES: 

    implicit none
! !ARGUMENTS:     
    integer                 :: n 
    type(LIS_metadataEntry) :: dataEntry
    integer                 :: t
    integer                 :: vlevel
    real, intent(in)        :: in_value
    character(len=*)        :: unit
    character(len=*)        :: direction
    real                    :: vmin
    real                    :: vmax
    integer, intent(in)     :: nsiblings
    integer, intent(in), dimension(nsiblings) :: siblings
    logical                 :: model_patch

! 
! !DESCRIPTION: 
!  This routine maps a single output variable to the appropriate variable 
!  in the generic list of the LIS history writer. 
!EOP
    integer                 :: i
    logical                 :: unit_status
    logical                 :: dir_status
    real                    :: mfactor
    real                    :: value
    integer                 :: sftype
!    character(len=20)            :: cfunit
    
    sftype = LIS_domain(n)%tile(t)%sftype
       
    unit_status = .false.
    do i=1,dataEntry%nunits
       if(unit.eq.dataEntry%unittypes(i)) then 
          unit_status = .true. 
          exit
       endif
    enddo
    
    dir_status = .false. 
    do i=1,dataEntry%ndirs
       if(direction.eq.dataEntry%dirtypes(i)) then 
          dir_status = .true. 
          exit
       endif
    enddo
    
    if(unit_status.and.dir_status) then 
       if(unit.eq.dataEntry%units) then
          ! it is assumed that there will be only two 
          ! directions. 
          if(direction.eq.dataEntry%dir) then 
             mfactor = 1.0
          else
             mfactor = -1.0
          endif
          
          if(in_value.ne.LIS_rc%udef) then 
             ! Correct the direction of value
             value = in_value * mfactor
          else
             value = in_value
          endif
          
          if(mfactor.eq.1) then 
             dataEntry%valid_min = vmin
             dataEntry%valid_max = vmax
          else
             dataEntry%valid_min = vmax
             dataEntry%valid_max = vmin
          endif
          if(value.ne.LIS_rc%udef) then 
             ! accumulate values and record instantaneous values
             if(dataEntry%timeAvgOpt.eq.2) then 
                dataEntry%modelOutput(1,t,vlevel) = &
                     dataEntry%modelOutput(1,t,vlevel) + value
                dataEntry%modelOutput(2,t,vlevel) = value
                !$OMP CRITICAL 
                dataEntry%count(t,vlevel) = &
                     dataEntry%count(t,vlevel)+1
                !$OMP END CRITICAL 
                ! accumulate values
             elseif(dataEntry%timeAvgOpt.eq.1 .or. &
                  dataEntry%timeAvgOpt.eq.3) then 
                dataEntry%modelOutput(1,t,vlevel) = &
                     dataEntry%modelOutput(1,t,vlevel) + value
                !$OMP CRITICAL 
                dataEntry%count(t,vlevel) = &
                     dataEntry%count(t,vlevel)+1
                !$OMP END CRITICAL 
                ! record instantaneous values
             else 
                dataEntry%modelOutput(1,t,vlevel) = value
                dataEntry%count(t,vlevel) = 1
             endif
             
             if ( dataEntry%minMaxOpt /= 0 ) then
                if ( value < dataEntry%minimum(t,vlevel) ) then
                   do i = 1, nsiblings
                      dataEntry%minimum(siblings(i),vlevel) = value
                   enddo
                endif
                if ( value > dataEntry%maximum(t,vlevel) ) then
                   do i = 1, nsiblings
                      dataEntry%maximum(siblings(i),vlevel) = value
                   enddo
                endif
             endif
          endif
          dataEntry%diagflag = 1 
          
       endif
    endif
    if(.not.unit_status) then 
       write(LIS_logunit,*) '[ERR] ',trim(dataEntry%units),&
            ' for field ',trim(dataEntry%standard_name),' is not defined '
       write(LIS_logunit,*) '[ERR] for diagnostic output...'
       write(LIS_logunit,*) '[ERR] supported unit types: ',dataEntry%unittypes
       write(LIS_logunit,*) '[ERR] Program stopping ..'
       call LIS_endrun()       
    endif
    if(.not.dir_status) then 
       write(LIS_logunit,*) '[ERR] ',trim(dataEntry%dir),&
            ' for field ',trim(dataEntry%standard_name),' is not defined '
       write(LIS_logunit,*) '[ERR] for diagnostic output...'
       write(LIS_logunit,*) '[ERR] supported direction types: ',&
            dataEntry%dir
       write(LIS_logunit,*) '[ERR] Program stopping ..'
       call LIS_endrun()       
    endif
  end subroutine diagnoseDataEntry
!BOP
! !ROUTINE: LIS_resetOutputVars
! \label{LIS_resetOutputVars}
! 
! !INTERFACE: 
  subroutine LIS_resetOutputVars(n,group)
! !ARGUMENTS:
    integer, intent(in) :: n 
    integer, intent(in) :: group
! 
! !DESCRIPTION: 
!   This routine resets the specified variables to enable time averaging 
!    for the next history output step. 
!
!   It also resets the minimum and maximum fields.
!   
!   The arguments are: 
!   \begin{description}
!   \item[n]  index of the nest \newline
!   \item[group]  output group (1- LSM, 2-ROUTING, 3-RTM) \newline
!   \end{description}
!EOP
    integer :: index
    type(LIS_metadataEntry), pointer :: dataEntry 

    if(group.eq.1) then !LSM output
       dataEntry => LIS_histData(n)%head_lsm_list
    elseif(group.eq.2) then !ROUTING
       dataEntry => LIS_histData(n)%head_routing_list
    elseif(group.eq.3) then !RTM
       dataEntry => LIS_histData(n)%head_rtm_list
    elseif(group.eq.4) then !Irrig
       dataEntry => LIS_histData(n)%head_irrig_list
    endif

    do while ( associated(dataEntry) )
       if(dataEntry%selectOpt.ne.0) then 
          call resetOutputVar(dataEntry)
       endif
       dataEntry => dataEntry%next
    enddo

  end subroutine LIS_resetOutputVars


!BOP
! !ROUTINE: resetOutputVar
! \label{resetOutputVar}
! 
! !INTERFACE: 
  subroutine resetOutputVar(dataEntry)

  implicit none

! !ARGUMENTS:
   type(LIS_metadataEntry), pointer :: dataEntry

! !DESCRIPTION: 
!   This routine resets the specified variable to enable time averaging 
!   for the next history output step. 
!
!   It also resets the minimum and maximum fields.
!
! The arguments are:
!   \begin{description}
!   \item[dataEntry] data structure to reset
!   \end{description}
!EOP
      dataEntry%modelOutput = 0.0
      dataEntry%count = 0 
      dataEntry%diagFlag = 0 

      if ( dataEntry%minMaxOpt .ne. 0 ) then
         dataEntry%minimum = LIS_MOC_MAX_NUM
         dataEntry%maximum = LIS_MOC_MIN_NUM
      endif

  end subroutine resetOutputVar

          
!BOP
! !ROUTINE: LIS_rescaleCount
! \label{LIS_rescaleCount}
! 
! !INTERFACE: 
  subroutine LIS_rescaleCount(n,group)
! 
    use LIS_mpiMod

     implicit none
! !ARGUMENTS:
     integer, intent(in) :: n, group

! !DESCRIPTION: 
!   This routine rescales the count.  dataEntry%count records
!   how many time the corresponding variable has been diagnosed.
!   This count is incremented each time the diagnose routine is called
!   in terms of both time and space.  Say variable var has 10 tiles, and
!   it is diagnosed for 3 model time-steps, then its count will be 30.
!   Thus count must be rescaled by the number of tiles to get the desired
!   number of model time-step calls.
!
!   count is used to time-average output.  count is not used when writing
!   accumulated output.  count is 1 when diagnosing instantaneous output.
!
!   The arguments are: 
!   \begin{description}
!   \item[n]  index of the nest \newline
!   \item[group]  output group (1 = LSM; 2 = ROUTING; 3 = RTM) \newline
!   \end{description}
!EOP

     type(LIS_metadataEntry), pointer :: dataEntry 
     integer :: count
     integer :: k, m
     integer :: ierr
     
     if(group.eq.1) then !LSM output
        dataEntry => LIS_histData(n)%head_lsm_list
     elseif(group.eq.2) then !ROUTING
        dataEntry => LIS_histData(n)%head_routing_list
     elseif(group.eq.3) then !RTM
        dataEntry => LIS_histData(n)%head_rtm_list
     elseif(group.eq.4) then !Irrigation
        dataEntry => LIS_histData(n)%head_irrig_list
     endif

     do while ( associated(dataEntry) )
        if(dataEntry%selectOpt.ne.0) then 
           do k=1,dataEntry%vlevels
!#if (defined SPMD)
!              call mpi_reduce(sum(dataEntry%count(:,k)),&
!                   count, 1, MPI_INTEGER, MPI_SUM, 0, &
!                   LIS_mpi_comm,ierr)
!#else
!              count = sum(dataEntry%count(:,k))
!#endif
#if (defined SPMD)
              call mpi_reduce(dataEntry%diagFlag,&
                   count, 1, MPI_INTEGER, MPI_SUM, 0, &
                   LIS_mpi_comm,ierr)
#else
              count = dataEntry%diagflag
#endif
              if(LIS_masterproc) then 
                 if(count.eq.0) then 
                    write(LIS_logunit,*) '[ERR] ',dataEntry%short_name,&
                         ' field is not defined'
                    write(LIS_logunit,*) '[ERR] for diagnostic output...'
                    write(LIS_logunit,*) '[ERR] Please exclude it from the model output attributes table'
                    write(LIS_logunit,*) '[ERR] Program stopping ..'
                    call LIS_endrun()
                 endif
              endif
           enddo

#if 0        
           ! Rescale the count when the timeAvgOpt indicates to use
           ! time averaging.
           ! If timeAvgOpt indicates instantaneous-only or accumulate
           ! then there is no need to rescale the count

           if ( (dataEntry%timeAvgOpt == 1) .or. &
                (dataEntry%timeAvgOpt == 2) ) then 
              do k=1,dataEntry%vlevels
                 do m=1,LIS_rc%max_model_types
                    if(dataEntry%count(m,k).gt.0) then 
                       dataEntry%count(m,k) = dataEntry%count(m,k) / &
                                              LIS_rc%npatch(n,m)
                   endif
                 enddo
              enddo
           endif
#endif
           
        endif
        dataEntry => dataEntry%next
     enddo

  end subroutine LIS_rescaleCount


!BOP
! !ROUTINE: register_dataEntry
! \label{register_dataEntry}
! 
! !INTERFACE: 
  subroutine register_dataEntry(count,lis_moc_index,head_dataEntry,&
                                n,nunits,ntiles,&
                                unittypes,ndirs,dirtypes,form,gribsfc,griblvl,&
                                model_patch)
  implicit none
! !ARGUMENTS:
    integer                 :: count
    integer                 :: lis_moc_index
    type(LIS_metadataEntry), pointer :: head_dataEntry
    integer                 :: n
    integer                 :: nunits
    integer                 :: ntiles
    character(len=*)        :: unittypes(nunits)
    integer                 :: ndirs
    character(len=*)        :: dirtypes(ndirs)
    integer                 :: form
    integer                 :: gribsfc
    integer                 :: griblvl
    logical, optional       :: model_patch
! 
! !DESCRIPTION: 
!  This routine completes the assignment of elements of the history
!  output linked list object corresponding to the variable given by
!  lis\_moc\_index.
!
!   The arguments are: 
!   \begin{description}
!   \item[lis\_moc\_index] variable to process
!   \item[head\_dataEntry] head of the given history output linked list
!   \item[n] nest
!   \item[nunits] number of units associated with this variable
!   \item[ntiles] number of tiles
!   \item[unittypes] types of units associated with this variable
!   \item[ndirs] number of direction types associated with this variable
!   \item[dirtypes] types of direction types associated with this variable
!   \item[form] form for stats file
!   \item[gribsfc] type of GRIB surface
!   \item[griblvl] type of GRIB level
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[allocate\_dataEntry] (\ref{allocate_dataEntry}) \newline
!      allocates the data structures related to the variable being output
!   \end{description}
!EOP
!   
      type(LIS_metadataEntry), pointer :: dataEntry
      logical                          :: model_patch_tmp
      integer                          :: i, ierr

      if(present(model_patch)) then 
         model_patch_tmp = model_patch
      else
         model_patch_tmp = .false.
      endif

      count = count + 1
      lis_moc_index = count

      dataEntry => head_dataEntry
      do while ( associated(dataEntry%next) )
         dataEntry => dataEntry%next
      enddo

      dataEntry%index   = lis_moc_index
      dataEntry%form    = form
      dataEntry%gribSfc = gribsfc
      dataEntry%gribLvl = griblvl

      call allocate_dataEntry(dataEntry,nunits,ntiles,unittypes, &
                              ndirs,dirtypes, model_patch_tmp)

      ierr = -1
      do i = 1,dataEntry%ndirs
         if ( dataEntry%dir == dataEntry%dirtypes(i) ) then 
            ierr = i
            exit
         endif
      enddo

      if ( ierr == -1 ) then 
         write(LIS_logunit,*) '[ERR] register_dataEntry: ', &
                              dataEntry%standard_name,       &
                              ' with direction type of ',    &
                              dataEntry%dir,                 &
                              ' is not defined'
         write(LIS_logunit,*) '[ERR] for diagnostic output.'
         write(LIS_logunit,*) '[ERR] Acceptable options are:'
         do i = 1,dataEntry%ndirs
            write(LIS_logunit,*) '[ERR] ',trim(dataEntry%dirtypes(i))
         enddo
         write(LIS_logunit,*) '[ERR] Program stopping.'
         call LIS_endrun()
      endif          

      ierr = -1
      do i = 1,dataEntry%nunits
         if ( dataEntry%units == dataEntry%unittypes(i) ) then 
            ierr = i
            exit
         endif
      enddo    

      if ( ierr == -1 ) then 
         write(LIS_logunit,*) '[ERR] register_dataEntry: ', &
                              dataEntry%standard_name,       &
                              ' in units of ',               &
                              dataEntry%units,               &
                              ' is not defined'
         write(LIS_logunit,*) '[ERR] for diagnostic output.'
         write(LIS_logunit,*) '[ERR] Program stopping.'
         do i = 1,dataEntry%nunits
            write(LIS_logunit,*) '[ERR] Acceptable options are:'
            write(LIS_logunit,*) '[ERR] ',trim(dataEntry%unittypes(i))
         enddo    
         call LIS_endrun()
      endif

  end subroutine register_dataEntry

!BOP
! !ROUTINE: set_ptr_into_list
! \label{set_ptr_into_list}
! 
! !INTERFACE: 
  subroutine set_ptr_into_list(count, head, array)
     implicit none
! !ARGUMENTS:
     integer                          :: count
     type(LIS_metadataEntry), pointer :: head
     type(dep), dimension(count)      :: array
! 
! !DESCRIPTION: 
!  This routine takes an array of pointers and sets each element of the
!  array to point directly to its corresponding element in the given
!  history output linked list.  This allows for direct access to the
!  elements in the history output linked list.
!
!   The arguments are: 
!   \begin{description}
!   \item[count] number of elements in the given history output linked list
!   \item[head] head of the given history output linked list
!   \item[array] array of pointers to point directly to the elements of \newline
!                the given history output linked list
!   \end{description}
!EOP
!   

     type(LIS_metadataEntry), pointer :: dataEntry
     integer :: i

     dataEntry => head
     do i = 1, count
        array(i)%dataEntryPtr => dataEntry
        dataEntry => dataEntry%next
     enddo

  end subroutine set_ptr_into_list

end module LIS_histDataMod
