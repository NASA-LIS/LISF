!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_histDataMod
! \label(LVT_histDataMod)
!
! !INTERFACE:
module LVT_histDataMod
! 
! !USES:
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_statsDataMod
  use LVT_pluginIndices

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module is used by the user to define the metadata associated with 
!   the variables that are handled by LVT. 
!   The module lists an exhaustive list of the the land relevant
!   variables. The user can choose a subset from this list through the 
!   lis configuration file for model output. The list is based on the 
!   the variable definitions from ALMA specification. 
!   \textsl{http://www.lmd.jussieu.fr/ALMA/}
! 
!   LVT defines two separate lists (datastream 1 and datastream 2) 
!   of variables chosen in an analysis instance. The variables defined
!   in the datastream 1 are compared against the variables in datastream 2
!   using the specified set of analysis metrics. 
!   
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!  17 Oct 2018  Mahdi Navari  Enhanced the LVT reader to read the 
!               Veg. Water Content (VWC) from SMAP SM dataset 
!  19 Nov 2018  Mahdi Navari added suport to read SMAP_L3 brightness temperature
!
! 
!EOP
!BOP
  public :: LVT_histDataInit
  public :: LVT_logSingleDataStreamVar
  public :: LVT_checkDatastreamSetup
  public :: LVT_getDataStream1Ptr
  public :: LVT_getDataStream2Ptr
  public :: LVT_getDataStream3Ptr
  public :: LVT_getstatsEntryPtr
  public :: LVT_getdataEntryUnits
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
! MOC - Model Output Convention
!-----------------------------------------------------------------------------
  public :: LVT_histData

  public :: LVT_MOC_SWNET   
  public :: LVT_MOC_LWNET   
  public :: LVT_MOC_QLE     
  public :: LVT_MOC_QH      
  public :: LVT_MOC_QG      
  public :: LVT_MOC_QF      
  public :: LVT_MOC_QV      
  public :: LVT_MOC_QTAU    
  public :: LVT_MOC_QA         
  public :: LVT_MOC_DELSURFHEAT
  public :: LVT_MOC_DELCOLDCONT
  public :: LVT_MOC_BR
  public :: LVT_MOC_EF
  public :: LVT_MOC_SNOWF     
  public :: LVT_MOC_RAINF     
  public :: LVT_MOC_EVAP      
  public :: LVT_MOC_QS        
  public :: LVT_MOC_QREC      
  public :: LVT_MOC_QSB       
  public :: LVT_MOC_QSM       
  public :: LVT_MOC_QFZ       
  public :: LVT_MOC_QST       
  public :: LVT_MOC_DELSOILMOIST
  public :: LVT_MOC_DELSWE    
  public :: LVT_MOC_DELSURFSTOR
  public :: LVT_MOC_DELINTERCEPT
  public :: LVT_MOC_SNOWT     
  public :: LVT_MOC_VEGT      
  public :: LVT_MOC_BARESOILT 
  public :: LVT_MOC_AVGSURFT  
  public :: LVT_MOC_GROUNDAVGT
  public :: LVT_MOC_GROUNDVEGT
  public :: LVT_MOC_RADT      
  public :: LVT_MOC_ALBEDO 
  public :: LVT_MOC_VISDIRALBEDO 
  public :: LVT_MOC_VISDIFALBEDO 
  public :: LVT_MOC_NIRDIRALBEDO 
  public :: LVT_MOC_NIRDIFALBEDO 
  public :: LVT_MOC_SWDIRALBEDO 
  public :: LVT_MOC_SWDIFALBEDO 
  public :: LVT_MOC_SWE
  public :: LVT_MOC_SNOWICE
  public :: LVT_MOC_SWEVEG    
  public :: LVT_MOC_SNOWAGE
  public :: LVT_MOC_SURFSTOR  
  public :: LVT_MOC_SOILMOIST 
  public :: LVT_MOC_VEGWATERCONTENT 
  public :: LVT_MOC_VOD
  public :: LVT_MOC_L3TBv_D 
  public :: LVT_MOC_L3TBv_A 
  public :: LVT_MOC_L3TBh_D 
  public :: LVT_MOC_L3TBh_A 
  public :: LVT_MOC_SOILTEMP  
  public :: LVT_MOC_SMLIQFRAC
  public :: LVT_MOC_SMFROZFRAC
  public :: LVT_MOC_SOILWET   
  public :: LVT_MOC_MATRICPOTENTIAL
  public :: LVT_MOC_POTEVAP   
  public :: LVT_MOC_ECANOP    
  public :: LVT_MOC_TVEG      
  public :: LVT_MOC_ESOIL     
  public :: LVT_MOC_EWATER    
  public :: LVT_MOC_ROOTMOIST 
  public :: LVT_MOC_ROOTTEMP
  public :: LVT_MOC_CANOPINT  
  public :: LVT_MOC_EVAPSNOW     
  public :: LVT_MOC_SUBSNOW   
  public :: LVT_MOC_SUBSURF   
  public :: LVT_MOC_ACOND   
  public :: LVT_MOC_WATERTABLED
  public :: LVT_MOC_TWS
  public :: LVT_MOC_GWS
  public :: LVT_MOC_WT
  public :: LVT_MOC_SNOWCOVER
  public :: LVT_MOC_SALBEDO
  public :: LVT_MOC_SNOWTPROF
  public :: LVT_MOC_SNOWDEPTH
  public :: LVT_MOC_SLIQFRAC
  public :: LVT_MOC_LWUP
  public :: LVT_MOC_GPP
  public :: LVT_MOC_NPP
  public :: LVT_MOC_NEE
  public :: LVT_MOC_AUTORESP
  public :: LVT_MOC_HETERORESP
  public :: LVT_MOC_LEAFRESP
  public :: LVT_MOC_TOTSOILCARB
  public :: LVT_MOC_TOTLIVBIOM
  
  public :: LVT_MOC_WINDFORC  
  public :: LVT_MOC_RAINFFORC 
  public :: LVT_MOC_SNOWFFORC 
  public :: LVT_MOC_CRAINFFORC 
  public :: LVT_MOC_TAIRFORC  
  public :: LVT_MOC_QAIRFORC  
  public :: LVT_MOC_PSURFFORC 
  public :: LVT_MOC_SWDOWNFORC
  public :: LVT_MOC_LWDOWNFORC
  public :: LVT_MOC_DIRECTSWFORC
  public :: LVT_MOC_DIFFUSESWFORC
  public :: LVT_MOC_NWINDFORC  
  public :: LVT_MOC_EWINDFORC  
  public :: LVT_MOC_FHEIGHTFORC  
  public :: LVT_MOC_CHFORC  
  public :: LVT_MOC_CMFORC  
  public :: LVT_MOC_EMISSFORC  
  public :: LVT_MOC_MIXRATIOFORC  
  public :: LVT_MOC_COSZENFORC  
  public :: LVT_MOC_ALBEDOFORC  
  public :: LVT_MOC_PARDRFORC  
  public :: LVT_MOC_PARDFFORC  
!<for vic>
  public :: LVT_MOC_SNOWFLAGFORC
  public :: LVT_MOC_DENSITYFORC
  public :: LVT_MOC_VAPORPRESSFORC
  public :: LVT_MOC_VAPORPRESSDEFICITFORC
  public :: LVT_MOC_ARESIST
!</for vic>

  public :: LVT_MOC_LANDMASK  
  public :: LVT_MOC_LANDCOVER 
  public :: LVT_MOC_SOILTYPE  
  public :: LVT_MOC_SANDFRAC
  public :: LVT_MOC_CLAYFRAC
  public :: LVT_MOC_SILTFRAC
  public :: LVT_MOC_POROSITY
  public :: LVT_MOC_SOILCOLOR 
  public :: LVT_MOC_ELEVATION
  public :: LVT_MOC_SLOPE
  public :: LVT_MOC_ASPECT
  public :: LVT_MOC_LAI       
  public :: LVT_MOC_SAI       
  public :: LVT_MOC_SNFRALBEDO
  public :: LVT_MOC_MXSNALBEDO
  public :: LVT_MOC_GREENNESS 
  public :: LVT_MOC_NDVI
  public :: LVT_MOC_SIF
  public :: LVT_MOC_TEMPBOT  

  public :: LVT_MOC_CCOND
  public :: LVT_MOC_VPD
  public :: LVT_MOC_RELSMC
  public :: LVT_MOC_RHMIN
  public :: LVT_MOC_TOTALPRECIP     
  public :: LVT_MOC_RAINFCONV
  public :: LVT_MOC_SOILET
  public :: LVT_MOC_Z0BRD
  public :: LVT_MOC_ROUGHNESS
  public :: LVT_MOC_T2DIAG
  public :: LVT_MOC_Q2DIAG
  public :: LVT_MOC_RNET   
  public :: LVT_MOC_CH
  public :: LVT_MOC_CM
  public :: LVT_MOC_MIXRATIO

  public :: LVT_MOC_DR_CATEGORY
  public :: LVT_MOC_PERCENTILE

  !SAC/SNOW17
  public :: LVT_MOC_SACQS
  public :: LVT_MOC_SACQSB
  public :: LVT_MOC_SACSWE
  public :: LVT_MOC_SACPOTEVAP
  public :: LVT_MOC_SACEVAP
  public :: LVT_MOC_SACUZTWC
  public :: LVT_MOC_SACUZFWC
  public :: LVT_MOC_SACLZTWC
  public :: LVT_MOC_SACLZFSC
  public :: LVT_MOC_SACLZFPC
  public :: LVT_MOC_SACADIMPC
  !public :: LVT_MOC_SACSOILMOIST1
  public :: LVT_MOC_SACSOILMOIST2
  !public :: LVT_MOC_SACSOILMOIST3
  public :: LVT_MOC_SACSOILMOIST4
  public :: LVT_MOC_SACSOILMOIST5
  !public :: LVT_MOC_SACSOILMOIST6
  public :: LVT_MOC_SACSOILTEMP1
  public :: LVT_MOC_SACSOILTEMP2
  public :: LVT_MOC_SACSFRAC
  public :: LVT_MOC_SACSDEPTH
  public :: LVT_MOC_SACUZTWM
  public :: LVT_MOC_SACUZFWM
  public :: LVT_MOC_SACLZTWM
  public :: LVT_MOC_SACLZFSM
  public :: LVT_MOC_SACLZFPM
  public :: LVT_MOC_SACADIMP
  public :: LVT_MOC_SACRAINF
  public :: LVT_MOC_SACTAIR

  public :: LVT_MOC_SNOW17SWE
  public :: LVT_MOC_SNOW17TWE
  public :: LVT_MOC_SNOW17AEADJ
  public :: LVT_MOC_SNOW17SWEINCR
  public :: LVT_MOC_SNOW17WEINCR
  public :: LVT_MOC_SNOW17TWEINCR
  public :: LVT_MOC_SNOW17SFRAC1
  public :: LVT_MOC_SNOW17NEGHS
  public :: LVT_MOC_SNOW17LIQW
  public :: LVT_MOC_SNOW17DS
  public :: LVT_MOC_SNOW17AESCINCR
  public :: LVT_MOC_SNOW17CNHS
  public :: LVT_MOC_SNOW17SFRAC4
  public :: LVT_MOC_SNOW17ACCMAX
  public :: LVT_MOC_SNOW17SFRAC6
  public :: LVT_MOC_SNOW17SFRAC7
  public :: LVT_MOC_SNOW17SFRAC8
  public :: LVT_MOC_SNOW17SFRAC9
  public :: LVT_MOC_SNOW17SFRAC10
  public :: LVT_MOC_SNOW17SFRAC11
  public :: LVT_MOC_SNOW17SFRAC12
  public :: LVT_MOC_SNOW17SFRAC13
  public :: LVT_MOC_SNOW17SFRAC14
  public :: LVT_MOC_SNOW17SFRAC15
  public :: LVT_MOC_SNOW17RAINF
  public :: LVT_MOC_SNOW17SNOWF
  public :: LVT_MOC_SNOW17TAIR
  public :: LVT_MOC_SNOW17RMLT
  public :: LVT_MOC_SNOW17RMINCR

  public ::   LVT_MOC_SACUZTWH
  public ::   LVT_MOC_SACUZFWH
  public ::   LVT_MOC_SACLZTWH
  public ::   LVT_MOC_SACLZFSH
  public ::   LVT_MOC_SACLZFPH
  
  public ::   LVT_MOC_SACSWINT
  public ::   LVT_MOC_SACTSINT
  public ::   LVT_MOC_SACSWHINT
  public ::   LVT_MOC_SACFROST

!</for vic>
  public :: LVT_MOC_VIC_PET_SATSOIL
  public :: LVT_MOC_VIC_PET_H2OSURF
  public :: LVT_MOC_VIC_PET_SHORT
  public :: LVT_MOC_VIC_PET_TALL
  public :: LVT_MOC_VIC_PET_NATVEG
  public :: LVT_MOC_VIC_PET_VEGNOCR
!</for vic>

!FLDAS
  public :: LVT_MOC_PETFORC
  public :: LVT_MOC_REFETFORC
  ! FLDAS-WRSI OUTPUTS LIST
  public :: LVT_MOC_REFET
  public :: LVT_MOC_ETa
  public :: LVT_MOC_SOS
  public :: LVT_MOC_WRSI
  public :: LVT_MOC_KF2
  public :: LVT_MOC_SumWR
  public :: LVT_MOC_SumET
  public :: LVT_MOC_SWI
  public :: LVT_MOC_SOSa
  public :: LVT_MOC_TotalSurplusWater
  public :: LVT_MOC_MaxSurplusWater
  public :: LVT_MOC_TotalWaterDeficit
  public :: LVT_MOC_MaxWaterDeficit
  public :: LVT_MOC_TotalAETInitial
  public :: LVT_MOC_TotalWRInitial
  public :: LVT_MOC_TotalSurplusWaterInitial
  public :: LVT_MOC_TotalWaterDeficitInitial
  public :: LVT_MOC_TotalAETVeg
  public :: LVT_MOC_TotalWRVeg
  public :: LVT_MOC_TotalSurplusWaterVeg
  public :: LVT_MOC_TotalWaterDeficitVeg
  public :: LVT_MOC_TotalAETFlower
  public :: LVT_MOC_TotalWRFlower
  public :: LVT_MOC_TotalSurplusWaterFlower
  public :: LVT_MOC_TotalWaterDeficitFlower
  public :: LVT_MOC_TotalAETRipe
  public :: LVT_MOC_TotalWRRipe
  public :: LVT_MOC_TotalSurplusWaterRipe
  public :: LVT_MOC_TotalWaterDeficitRipe
  public :: LVT_MOC_PermWiltDate
  public :: LVT_MOC_Wilting1
  public :: LVT_MOC_Wilting2
  public :: LVT_MOC_WRSIa
  public :: LVT_MOC_growing_season
  public :: LVT_MOC_WHC
  public :: LVT_MOC_LGP
  public :: LVT_MOC_WR_TimeStep ! SY
  public :: LVT_MOC_AET_TimeStep ! SY
  public :: LVT_MOC_WRSI_TimeStep ! SY
  public :: LVT_MOC_SurplusWater_TimeStep ! SY

! NLDAS
  public :: LVT_MOC_CAPEFORC

  public ::   LVT_MOC_LAKE_T_SNOW
  public ::   LVT_MOC_LAKE_T_ICE
  public ::   LVT_MOC_LAKE_T_MNW
  public ::   LVT_MOC_LAKE_T_WML
  public ::   LVT_MOC_LAKE_T_BOT
  public ::   LVT_MOC_LAKE_T_B1
  public ::   LVT_MOC_LAKE_C_T
  public ::   LVT_MOC_LAKE_H_SNOW
  public ::   LVT_MOC_LAKE_H_ICE
  public ::   LVT_MOC_LAKE_H_ML
  public ::   LVT_MOC_LAKE_H_B1
  public ::   LVT_MOC_LAKE_T_SFC
  public ::   LVT_MOC_LAKE_ALBEDO_WATER
  public ::   LVT_MOC_LAKE_ALBEDO_ICE
  public ::   LVT_MOC_LAKE_ALBEDO_SNOW
  public ::   LVT_MOC_LAKE_UFR_A
  public ::   LVT_MOC_LAKE_UFR_W
  public ::   LVT_MOC_LAKE_WCONV
  public ::   LVT_MOC_LAKE_Q_SE
  public ::   LVT_MOC_LAKE_Q_LA
  public ::   LVT_MOC_LAKE_I_W
  public ::   LVT_MOC_LAKE_Q_LWA
  public ::   LVT_MOC_LAKE_Q_LWW
  public ::   LVT_MOC_LAKE_Q_BOT 
!derived
  public :: LVT_MOC_EBAL
  public :: LVT_MOC_WBAL
  public :: LVT_MOC_EVAPBAL
  public :: LVT_MOC_SWEOVERP
  public :: LVT_MOC_ETOVERP
  public :: LVT_MOC_QSOVERP
  public :: LVT_MOC_QSBOVERP
  public :: LVT_MOC_RUNOFF
  public :: LVT_MOC_dS ! P-ET-Q

  public :: LVT_MOC_ECANOPOVERQLE
  public :: LVT_MOC_TVEGOVERQLE
  public :: LVT_MOC_ESOILOVERQLE

  public :: LVT_MOC_RIVSTO
  public :: LVT_MOC_RIVDPH
  public :: LVT_MOC_RIVVEL
  public :: LVT_MOC_STREAMFLOW
  public :: LVT_MOC_FLDOUT
  public :: LVT_MOC_FLDEVAP
  public :: LVT_MOC_FLDSTO
  public :: LVT_MOC_FLDDPH
  public :: LVT_MOC_FLDVEL
  public :: LVT_MOC_FLDFRC
  public :: LVT_MOC_FLDARE  
  public :: LVT_MOC_SFCELV
  public :: LVT_MOC_RNFSTO
  public :: LVT_MOC_BSFSTO

  public :: LVT_MOC_RTM_EMISSIVITY
  public :: LVT_MOC_RTM_TB
  public :: LVT_MOC_RTM_SM

  public :: LVT_MOC_IRRIGATEDWATER

  public :: LVT_MOC_COUNT

  public :: LVT_MOC_DA_ENSSPREAD
  public :: LVT_MOC_DA_INCR
  public :: LVT_MOC_DA_NINNOV
  public :: LVT_MOC_DA_OBSCOUNT

  public :: LVT_MOC_TAIRFORC_MIN
  public :: LVT_MOC_TAIRFORC_MAX

  public :: LVT_MOC_ESI

  public :: LVT_temp_maxvEntry
  public :: LVT_temp_minvEntry

   ! ALMA ENERGY BALANCE COMPONENTS
  integer :: LVT_MOC_SWNET(3)      = -9999
  integer :: LVT_MOC_LWNET(3)      = -9999
  integer :: LVT_MOC_QLE(3)        = -9999
  integer :: LVT_MOC_QH(3)         = -9999
  integer :: LVT_MOC_QG(3)         = -9999
  integer :: LVT_MOC_QF(3)         = -9999
  integer :: LVT_MOC_QV(3)         = -9999
  integer :: LVT_MOC_QTAU(3)       = -9999
  integer :: LVT_MOC_QA(3)         = -9999
  integer :: LVT_MOC_DELSURFHEAT(3) = -9999
  integer :: LVT_MOC_DELCOLDCONT(3) = -9999

   ! ALMA WATER BALANCE COMPONENTS
  integer :: LVT_MOC_SNOWF(3)      = -9999
  integer :: LVT_MOC_RAINF(3)      = -9999
  integer :: LVT_MOC_EVAP(3)       = -9999
  integer :: LVT_MOC_QS(3)         = -9999
  integer :: LVT_MOC_QREC(3)       = -9999
  integer :: LVT_MOC_QSB(3)        = -9999
  integer :: LVT_MOC_QSM(3)        = -9999
  integer :: LVT_MOC_QFZ(3)        = -9999
  integer :: LVT_MOC_QST(3)        = -9999
  integer :: LVT_MOC_DELSOILMOIST(3) = -9999
  integer :: LVT_MOC_DELSWE(3)     = -9999
  integer :: LVT_MOC_DELSURFSTOR(3)  = -9999
  integer :: LVT_MOC_DELINTERCEPT(3) = -9999

   ! ALMA SURFACE STATE VARIABLES
  integer :: LVT_MOC_SNOWT(3)      = -9999
  integer :: LVT_MOC_VEGT(3)       = -9999
  integer :: LVT_MOC_BARESOILT(3)  = -9999
  integer :: LVT_MOC_AVGSURFT(3)   = -9999
  integer :: LVT_MOC_GROUNDAVGT(3) = -9999
  integer :: LVT_MOC_GROUNDVEGT(3) = -9999
  integer :: LVT_MOC_RADT(3)       = -9999
  integer :: LVT_MOC_ALBEDO(3)     = -9999
  integer :: LVT_MOC_VISDIRALBEDO(3) = -9999
  integer :: LVT_MOC_VISDIFALBEDO(3) = -9999
  integer :: LVT_MOC_NIRDIRALBEDO(3) = -9999
  integer :: LVT_MOC_NIRDIFALBEDO(3) = -9999
  integer :: LVT_MOC_SWDIRALBEDO(3) = -9999
  integer :: LVT_MOC_SWDIFALBEDO(3) = -9999
  integer :: LVT_MOC_SWE(3)        = -9999
  integer :: LVT_MOC_SNOWICE(3)    = -9999
  integer :: LVT_MOC_SWEVEG(3)     = -9999
  integer :: LVT_MOC_SNOWAGE(3)    = -9999 
  integer :: LVT_MOC_SURFSTOR(3)   = -9999
  integer :: LVT_MOC_VEGWATERCONTENT(3)   = -9999 
  integer :: LVT_MOC_VOD(3)   = -9999 
! integer :: LVT_MOC_L3TB(3)   = -9999 !MN
 integer :: LVT_MOC_L3TBv_D(3)   = -9999 !MN
 integer :: LVT_MOC_L3TBv_A(3)   = -9999 !MN
 integer :: LVT_MOC_L3TBh_D(3)   = -9999 !MN
 integer :: LVT_MOC_L3TBh_A(3)   = -9999 !MN

   ! ALMA SUBSURFACE STATE VARIABLES
   integer :: LVT_MOC_SOILMOIST(3)  = -9999
   integer :: LVT_MOC_SOILTEMP(3)   = -9999
   integer :: LVT_MOC_SMLIQFRAC(3)  = -9999
   integer :: LVT_MOC_SMFROZFRAC(3) = -9999
   integer :: LVT_MOC_SOILWET(3)    = -9999
   integer :: LVT_MOC_MATRICPOTENTIAL(3)    = -9999

   ! ALMA EVAPORATION COMPONENTS
   integer :: LVT_MOC_POTEVAP(3)    = -9999
   integer :: LVT_MOC_ECANOP(3)     = -9999
   integer :: LVT_MOC_TVEG(3)       = -9999
   integer :: LVT_MOC_ESOIL(3)      = -9999
   integer :: LVT_MOC_EWATER(3)     = -9999
   integer :: LVT_MOC_ROOTMOIST(3)  = -9999
   integer :: LVT_MOC_CANOPINT(3)   = -9999
   integer :: LVT_MOC_EVAPSNOW(3)   = -9999
   integer :: LVT_MOC_SUBSNOW(3)    = -9999
   integer :: LVT_MOC_SUBSURF(3)    = -9999
   integer :: LVT_MOC_ACOND(3)      = -9999

   ! ALMA OTHER HYDROLOGIC VARIABLES
  integer :: LVT_MOC_WATERTABLED(3)= -9999
  integer :: LVT_MOC_TWS(3)        = -9999
  integer :: LVT_MOC_GWS(3)        = -9999
  integer :: LVT_MOC_WT(3)         = -9999

   ! ALMA COLD SEASON PROCESSES
  integer :: LVT_MOC_SNOWCOVER(3)  = -9999
  integer :: LVT_MOC_SALBEDO(3)    = -9999
  integer :: LVT_MOC_SNOWTPROF(3)  = -9999
  integer :: LVT_MOC_SNOWDEPTH(3)  = -9999
  integer :: LVT_MOC_SLIQFRAC(3)   = -9999

   ! ALMA VARIABLES TO BE COMPARED WITH REMOTE SENSED DATA
   integer :: LVT_MOC_LWUP(3)       = -9999

   ! ALMA CARBON VARIABLES
   integer :: LVT_MOC_GPP(3)        = -9999
   integer :: LVT_MOC_NPP(3)        = -9999
   integer :: LVT_MOC_NEE(3)        = -9999
   integer :: LVT_MOC_AUTORESP(3)   = -9999
   integer :: LVT_MOC_HETERORESP(3) = -9999
   integer :: LVT_MOC_LEAFRESP(3)   = -9999
   integer :: LVT_MOC_TOTSOILCARB(3)= -9999
   integer :: LVT_MOC_TOTLIVBIOM(3) = -9999

   ! ALMA FORCING VARIABLES
   integer :: LVT_MOC_WINDFORC(3)   = -9999
   integer :: LVT_MOC_RAINFFORC(3)  = -9999
   integer :: LVT_MOC_SNOWFFORC(3)  = -9999
   integer :: LVT_MOC_CRAINFFORC(3) = -9999
   integer :: LVT_MOC_TAIRFORC(3)   = -9999
   integer :: LVT_MOC_QAIRFORC(3)   = -9999
   integer :: LVT_MOC_PSURFFORC(3)  = -9999
   integer :: LVT_MOC_SWDOWNFORC(3) = -9999
   integer :: LVT_MOC_LWDOWNFORC(3) = -9999

   ! CLSM FORCING VARIABLES
   integer :: LVT_MOC_PARDRFORC(3)  = -9999
   integer :: LVT_MOC_PARDFFORC(3)  = -9999

   ! PARAMETER OUTPUT - EXPERIMENTAL (USE W/WRF-WPS)
   integer :: LVT_MOC_LANDMASK(3)   = -9999
   integer :: LVT_MOC_LANDCOVER(3)  = -9999
   integer :: LVT_MOC_SOILTYPE(3)   = -9999
   integer :: LVT_MOC_SANDFRAC(3)   = -9999
   integer :: LVT_MOC_CLAYFRAC(3)   = -9999
   integer :: LVT_MOC_SILTFRAC(3)   = -9999
   integer :: LVT_MOC_POROSITY(3)   = -9999
   integer :: LVT_MOC_SOILCOLOR(3)  = -9999
   integer :: LVT_MOC_ELEVATION(3)  = -9999
   integer :: LVT_MOC_SLOPE(3)      = -9999
   integer :: LVT_MOC_ASPECT(3)     = -9999
   integer :: LVT_MOC_LAI(3)        = -9999
   integer :: LVT_MOC_SAI(3)        = -9999
   integer :: LVT_MOC_SNFRALBEDO(3) = -9999
   integer :: LVT_MOC_MXSNALBEDO(3) = -9999
   integer :: LVT_MOC_GREENNESS(3)  = -9999
   integer :: LVT_MOC_NDVI(3)       = -9999
   integer :: LVT_MOC_SIF(3)       = -9999
   integer :: LVT_MOC_TEMPBOT(3)   = -9999

   ! NLDAS OUTPUT
   integer :: LVT_MOC_CCOND(3)    = -9999
   integer :: LVT_MOC_VPD(3)    = -9999

   ! ADDITIONAL AFWA VARIABLES
   integer :: LVT_MOC_RELSMC(3)       = -9999
   integer :: LVT_MOC_RHMIN(3)        = -9999
   integer :: LVT_MOC_ROOTTEMP(3)  = -9999
   integer :: LVT_MOC_TOTALPRECIP(3) = -9999
   integer :: LVT_MOC_RAINFCONV(3) = -9999

   ! multivariate diagnostics (Bowen Ratio, Evaporative fraction)
   integer :: LVT_MOC_BR(3) = -9999
   integer :: LVT_MOC_EF(3) = -9999

   ! ADDITIONAL COUPLING FORCING VARIABLES
   integer :: LVT_MOC_DIRECTSWFORC(3)  = -9999
   integer :: LVT_MOC_DIFFUSESWFORC(3) = -9999
   integer :: LVT_MOC_NWINDFORC(3)     = -9999
   integer :: LVT_MOC_EWINDFORC(3)     = -9999
   integer :: LVT_MOC_FHEIGHTFORC(3)   = -9999
   integer :: LVT_MOC_CHFORC(3)        = -9999
   integer :: LVT_MOC_CMFORC(3)        = -9999
   integer :: LVT_MOC_EMISSFORC(3)     = -9999
   integer :: LVT_MOC_MIXRATIOFORC(3)  = -9999
   integer :: LVT_MOC_COSZENFORC(3)    = -9999
   integer :: LVT_MOC_ALBEDOFORC(3)    = -9999

   integer :: LVT_MOC_DR_CATEGORY(3)   = -9999
   integer :: LVT_MOC_PERCENTILE(3)   = -9999

   ! ADDITIONAL Noah3.x variables
   integer :: LVT_MOC_SOILET(3)  = -9999
   integer :: LVT_MOC_Z0BRD(3)   = -9999
   integer :: LVT_MOC_ROUGHNESS(3)   = -9999

   !t2,q2 diagnostics
   integer :: LVT_MOC_T2DIAG(3) = -9999
   integer :: LVT_MOC_Q2DIAG(3) = -9999
   integer :: LVT_MOC_RNET(3) = -9999
   integer :: LVT_MOC_CH(3)     = -9999
   integer :: LVT_MOC_CM(3)     = -9999
   integer :: LVT_MOC_MIXRATIO(3) = -9999

!<for vic>
   !Additional VIC forcing variables
   integer :: LVT_MOC_SNOWFLAGFORC(3)          = -9999
   integer :: LVT_MOC_DENSITYFORC(3)           = -9999
   integer :: LVT_MOC_VAPORPRESSFORC(3)        = -9999
   integer :: LVT_MOC_VAPORPRESSDEFICITFORC(3) = -9999
   integer :: LVT_MOC_ARESIST(3)               = -9999
!</for vic>

   ! Sacrament/Snow17 variables
   integer :: LVT_MOC_SACQS(3)             = -9999
   integer :: LVT_MOC_SACQSB(3)            = -9999
   integer :: LVT_MOC_SACSWE(3)            = -9999
   integer :: LVT_MOC_SACPOTEVAP(3)        = -9999
   integer :: LVT_MOC_SACEVAP(3)           = -9999
   integer :: LVT_MOC_SACUZTWC(3)          = -9999
   integer :: LVT_MOC_SACUZFWC(3)          = -9999
   integer :: LVT_MOC_SACLZTWC(3)          = -9999
   integer :: LVT_MOC_SACLZFSC(3)          = -9999
   integer :: LVT_MOC_SACLZFPC(3)          = -9999
   integer :: LVT_MOC_SACADIMPC(3)         = -9999
   !integer :: LVT_MOC_SACSOILMOIST1     = -9999
   integer :: LVT_MOC_SACSOILMOIST2(3)     = -9999
   !integer :: LVT_MOC_SACSOILMOIST3     = -9999
   integer :: LVT_MOC_SACSOILMOIST4(3)     = -9999
   integer :: LVT_MOC_SACSOILMOIST5(3)     = -9999
   !integer :: LVT_MOC_SACSOILMOIST6     = -9999
   integer :: LVT_MOC_SACSOILTEMP1(3)      = -9999
   integer :: LVT_MOC_SACSOILTEMP2(3)      = -9999
   integer :: LVT_MOC_SACSFRAC(3)          = -9999
   integer :: LVT_MOC_SACSDEPTH(3)         = -9999
   integer :: LVT_MOC_SACUZTWM(3)          = -9999
   integer :: LVT_MOC_SACUZFWM(3)          = -9999
   integer :: LVT_MOC_SACLZTWM(3)          = -9999
   integer :: LVT_MOC_SACLZFSM(3)          = -9999
   integer :: LVT_MOC_SACLZFPM(3)          = -9999
   integer :: LVT_MOC_SACADIMP(3)          = -9999
   integer :: LVT_MOC_SACRAINF(3)          = -9999
   integer :: LVT_MOC_SACTAIR(3)           = -9999
   integer :: LVT_MOC_SNOW17SWE(3)         = -9999
   integer :: LVT_MOC_SNOW17TWE(3)         = -9999
   integer :: LVT_MOC_SNOW17AEADJ(3)       = -9999
   integer :: LVT_MOC_SNOW17SWEINCR(3)     = -9999
   integer :: LVT_MOC_SNOW17WEINCR(3)      = -9999
   integer :: LVT_MOC_SNOW17TWEINCR(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC1(3)      = -9999
   integer :: LVT_MOC_SNOW17NEGHS(3)       = -9999
   integer :: LVT_MOC_SNOW17LIQW(3)        = -9999
   integer :: LVT_MOC_SNOW17DS(3)          = -9999
   integer :: LVT_MOC_SNOW17AESCINCR(3)    = -9999
   integer :: LVT_MOC_SNOW17CNHS(3)        = -9999
   integer :: LVT_MOC_SNOW17SFRAC4(3)      = -9999
   integer :: LVT_MOC_SNOW17ACCMAX(3)      = -9999
   integer :: LVT_MOC_SNOW17SFRAC6(3)      = -9999
   integer :: LVT_MOC_SNOW17SFRAC7(3)      = -9999
   integer :: LVT_MOC_SNOW17SFRAC8(3)      = -9999
   integer :: LVT_MOC_SNOW17SFRAC9(3)      = -9999
   integer :: LVT_MOC_SNOW17SFRAC10(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC11(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC12(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC13(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC14(3)     = -9999
   integer :: LVT_MOC_SNOW17SFRAC15(3)     = -9999
   integer :: LVT_MOC_SNOW17RAINF(3)       = -9999
   integer :: LVT_MOC_SNOW17SNOWF(3)       = -9999
   integer :: LVT_MOC_SNOW17TAIR(3)        = -9999
   integer :: LVT_MOC_SNOW17RMLT(3)        = -9999
   integer :: LVT_MOC_SNOW17RMINCR(3)      = -9999
   
   integer ::   LVT_MOC_SACUZTWH(3)  = -9999
   integer ::   LVT_MOC_SACUZFWH(3)  = -9999
   integer ::   LVT_MOC_SACLZTWH(3)  = -9999
   integer ::   LVT_MOC_SACLZFSH(3)  = -9999
   integer ::   LVT_MOC_SACLZFPH(3)  = -9999
   
   integer ::   LVT_MOC_SACSWINT(3)  = -9999
   integer ::   LVT_MOC_SACTSINT(3)  = -9999
   integer ::   LVT_MOC_SACSWHINT(3)  = -9999
   integer ::   LVT_MOC_SACFROST(3)  = -9999

!Noah-MP variables
   integer ::   LVT_MOC_LEAFMASS(3)  = -9999
   integer ::   LVT_MOC_ROOTMASS(3)  = -9999
   integer ::   LVT_MOC_STEMMASS(3)  = -9999
   integer ::   LVT_MOC_WOODMASS(3)  = -9999
   integer ::   LVT_MOC_CARBON_DEEPSOIL(3)  = -9999
   integer ::   LVT_MOC_CARBON_SHALLOWSOIL(3)  = -9999

!<for vic>
   integer :: LVT_MOC_VIC_PET_SATSOIL(3)   = -9999
   integer :: LVT_MOC_VIC_PET_H2OSURF(3)   = -9999
   integer :: LVT_MOC_VIC_PET_SHORT(3)     = -9999
   integer :: LVT_MOC_VIC_PET_TALL(3)      = -9999
   integer :: LVT_MOC_VIC_PET_NATVEG(3)    = -9999
   integer :: LVT_MOC_VIC_PET_VEGNOCR(3)   = -9999
!</for vic>

   !FLDAS
   integer :: LVT_MOC_PETFORC(3)         = -9999
   integer :: LVT_MOC_REFETFORC(3)       = -9999

   ! FLDAS-WRSI OUTPUTS LIST
   integer :: LVT_MOC_REFET(3) = -9999
   integer :: LVT_MOC_ETa(3) = -9999
   integer :: LVT_MOC_SOS(3) = -9999
   integer :: LVT_MOC_WRSI(3) = -9999
   integer :: LVT_MOC_KF2(3) = -9999
   integer :: LVT_MOC_SumWR(3) = -9999
   integer :: LVT_MOC_SumET(3) = -9999
   integer :: LVT_MOC_SWI(3) = -9999
   integer :: LVT_MOC_SOSa(3) = -9999
   integer :: LVT_MOC_TotalSurplusWater(3) = -9999
   integer :: LVT_MOC_MaxSurplusWater(3) = -9999
   integer :: LVT_MOC_TotalWaterDeficit(3) = -9999
   integer :: LVT_MOC_MaxWaterDeficit(3) = -9999
   integer :: LVT_MOC_TotalAETInitial(3) = -9999
   integer :: LVT_MOC_TotalWRInitial(3) = -9999
   integer :: LVT_MOC_TotalSurplusWaterInitial(3) = -9999
   integer :: LVT_MOC_TotalWaterDeficitInitial(3) = -9999
   integer :: LVT_MOC_TotalAETVeg(3) = -9999
   integer :: LVT_MOC_TotalWRVeg(3) = -9999
   integer :: LVT_MOC_TotalSurplusWaterVeg(3) = -9999
   integer :: LVT_MOC_TotalWaterDeficitVeg(3) = -9999
   integer :: LVT_MOC_TotalAETFlower(3) = -9999
   integer :: LVT_MOC_TotalWRFlower(3) = -9999
   integer :: LVT_MOC_TotalSurplusWaterFlower(3) = -9999
   integer :: LVT_MOC_TotalWaterDeficitFlower(3) = -9999
   integer :: LVT_MOC_TotalAETRipe(3) = -9999
   integer :: LVT_MOC_TotalWRRipe(3) = -9999
   integer :: LVT_MOC_TotalSurplusWaterRipe(3) = -9999
   integer :: LVT_MOC_TotalWaterDeficitRipe(3) = -9999
   integer :: LVT_MOC_PermWiltDate(3) = -9999
   integer :: LVT_MOC_Wilting1(3) = -9999
   integer :: LVT_MOC_Wilting2(3) = -9999
   integer :: LVT_MOC_WRSIa(3) = -9999
   integer :: LVT_MOC_growing_season(3) = -9999
   integer :: LVT_MOC_WHC(3)                      = -9999
   integer :: LVT_MOC_LGP(3)                      = -9999
   integer :: LVT_MOC_WR_TimeStep(3)              = -9999 ! SY
   integer :: LVT_MOC_AET_TimeStep(3)             = -9999 ! SY
   integer :: LVT_MOC_WRSI_TimeStep(3)            = -9999 ! SY
   integer :: LVT_MOC_SurplusWater_TimeStep(3)    = -9999 ! SY

   !NLDAS
   integer :: LVT_MOC_CAPEFORC(3)         = -9999

   integer :: LVT_MOC_LAKE_T_SNOW(3)    =   -9999
   integer :: LVT_MOC_LAKE_T_ICE(3) =   -9999
   integer :: LVT_MOC_LAKE_T_MNW(3) =   -9999
   integer :: LVT_MOC_LAKE_T_WML(3) =   -9999
   integer :: LVT_MOC_LAKE_T_BOT(3) =   -9999
   integer :: LVT_MOC_LAKE_T_B1 (3) =   -9999
   integer :: LVT_MOC_LAKE_C_T(3)   =   -9999
   integer :: LVT_MOC_LAKE_H_SNOW(3)    =   -9999
   integer :: LVT_MOC_LAKE_H_ICE(3) =   -9999
   integer :: LVT_MOC_LAKE_H_ML(3)  =   -9999
   integer :: LVT_MOC_LAKE_H_B1(3)  =   -9999
   integer :: LVT_MOC_LAKE_T_SFC(3) =   -9999
   integer :: LVT_MOC_LAKE_ALBEDO_WATER(3)  =   -9999
   integer :: LVT_MOC_LAKE_ALBEDO_ICE(3)    =   -9999
   integer :: LVT_MOC_LAKE_ALBEDO_SNOW(3)   =   -9999
   integer :: LVT_MOC_LAKE_UFR_A(3) =   -9999
   integer :: LVT_MOC_LAKE_UFR_W(3) =   -9999
   integer :: LVT_MOC_LAKE_WCONV(3) =   -9999
   integer :: LVT_MOC_LAKE_Q_SE(3)  =   -9999
   integer :: LVT_MOC_LAKE_Q_LA(3)  =   -9999
   integer :: LVT_MOC_LAKE_I_W(3)   =   -9999
   integer :: LVT_MOC_LAKE_Q_LWA(3) =   -9999
   integer :: LVT_MOC_LAKE_Q_LWW(3) =   -9999
   integer :: LVT_MOC_LAKE_Q_BOT(3) =   -9999

   integer :: LVT_MOC_EBAL(3)                     = -9999 
   integer :: LVT_MOC_WBAL(3)                     = -9999 
   integer :: LVT_MOC_EVAPBAL(3)                  = -9999 
   integer :: LVT_MOC_SWEOVERP(3)                 = -9999
   integer :: LVT_MOC_ETOVERP(3)                  = -9999
   integer :: LVT_MOC_QSOVERP(3)                  = -9999
   integer :: LVT_MOC_QSBOVERP(3)                 = -9999
   integer :: LVT_MOC_RUNOFF(3)                   = -9999
   integer :: LVT_MOC_dS(3)                   = -9999

   integer :: LVT_MOC_ECANOPOVERQLE(3)     = -9999
   integer :: LVT_MOC_TVEGOVERQLE(3)     = -9999
   integer :: LVT_MOC_ESOILOVERQLE(3)     = -9999

   integer :: LVT_MOC_STREAMFLOW(3) = -9999
   integer :: LVT_MOC_RIVSTO(3) = -9999
   integer :: LVT_MOC_RIVDPH(3) = -9999
   integer :: LVT_MOC_RIVVEL(3) = -9999
   integer :: LVT_MOC_FLDOUT(3) = -9999
   integer :: LVT_MOC_FLDEVAP(3) = -9999
   integer :: LVT_MOC_FLDSTO(3) = -9999
   integer :: LVT_MOC_FLDDPH(3) = -9999
   integer :: LVT_MOC_FLDVEL(3) = -9999
   integer :: LVT_MOC_FLDFRC(3) = -9999
   integer :: LVT_MOC_FLDARE(3) = -9999
   integer :: LVT_MOC_SFCELV(3) = -9999
   integer :: LVT_MOC_RNFSTO(3) = -9999
   integer :: LVT_MOC_BSFSTO(3) = -9999

   !variables related to RTMs
   integer :: LVT_MOC_RTM_EMISSIVITY(3)           = -9999 
   integer :: LVT_MOC_RTM_TB(3)                   = -9999 
   integer :: LVT_MOC_RTM_SM(3)                   = -9999 

   !DA diagnostics
   integer :: LVT_MOC_DA_ENSSPREAD(3)  = -9999
   integer :: LVT_MOC_DA_INCR(3)  = -9999
   integer :: LVT_MOC_DA_NINNOV(3)  = -9999
   integer :: LVT_MOC_DA_OBSCOUNT(3)  = -9999

   integer :: LVT_MOC_TAIRFORC_MIN(3)
   integer :: LVT_MOC_TAIRFORC_MAX(3)

   integer :: LVT_MOC_IRRIGATEDWATER(3)           = -9999 
   
   integer :: LVT_MOC_ESI(3) = -9999

   integer :: LVT_MOC_COUNT(3)

#if 0 
   ! SPECIAL CASE INDICES
   ! These are required because Min/Max support cannot be generically
   ! handled for GRIB output.  The routine writeSingleGrib1Var maps
   ! these two entries to LVT_MOC_TAIRFORC.
   ! They should not be counted in the LVT_MOC_COUNT total count.
   integer :: LVT_MOC_TAIRFORC_MIN = 142
   integer :: LVT_MOC_TAIRFORC_MAX = 143

   ! READ ABOVE NOTE ABOUT SPECIAL CASE INDICES
   integer :: LVT_MOC_COUNT      = 148
   ! Add the special cases.  LVT_MOC_GRIB_COUNT should be used only in
   ! LVT_gribMod.F90.
   integer :: LVT_MOC_GRIB_COUNT = 148 
#endif
!EOP

   type, public :: LVT_metadataEntry
      !index among all selected variables
      integer                   :: index     
      !long name (CF-based)
      character*100             :: long_name 
      !standard name (CF-based)
      character*100             :: standard_name
      !short name (ALMA based)
      character*100             :: short_name  
      !specified unit of the variable
      character*20              :: units  
      !number of unit types supported
      integer                   :: nunits 
      !supported unit types
      character*20, allocatable :: unittypes(:) 
      !number of supported direction types
      integer                   :: ndirs  
      !specified direction of the variable
      character*20              :: dir    
      !supported direction types
      character*20, allocatable :: dirtypes(:) 
      ! (scientific - E, else - F)
      character*1               :: format  
      !valid min values (for each unit type)
      real,         allocatable :: valid_min(:)
      !valid max values (for each unit type) 
      real,         allocatable :: valid_max(:)
      !Number of vertical levels
      integer                   :: vlevels    
      !time averging used (0-instantaneous, 1- time averaged, 2-accumulated)
      integer                   :: timeAvgOpt 
      !Number of vertical levels chosen in the analysis
      integer                   :: selectNlevs

      integer                   :: startNlevs
      integer                   :: endNlevs
      !default variable id  
      integer                   :: varid_def  
      !variable used to store the number of counts of the variable between
      !each time averaging window
      integer, allocatable      :: count(:,:,:) 
      !variable used to store the values of the variable
      integer, allocatable      :: count_status(:,:,:) 

      real, allocatable         :: value(:,:,:) 
      !whether stdev of the measurement is specified
      logical                   :: stdev_flag 
      !standard deviation of the measurement
      real,        allocatable  :: stdev(:,:)  
      !variable used to store the number of counts for the standard
      !deviation of the variable. 
      integer,     allocatable  :: count_stdev(:,:) 
      !variable used to determine if it needs to be derived (as opposed to
      !reading directly from the data files
      integer                   :: computeVar

      type(LVT_metadataEntry), pointer :: next
  end type LVT_metadataEntry

  ! To create an array of allocatables, you must create a derived
  ! datatype containing the allocatable.  Then you create
  ! and array of this datatype.
  type, public :: dep
     type(LVT_metadataEntry), pointer :: dataEntryPtr
  end type dep

  type, public :: statsdep
     type(LVT_statsEntry), pointer :: dataEntryPtr
  end type statsdep


  type, public :: output_meta
     integer                 :: xtimeID
     type(LVT_metadataEntry) :: xlat
     type(LVT_metadataEntry) :: xlon
!hack to get the HYCOM SST fields included
     type(LVT_metadataEntry) :: watertemp
!EMK...Add sea ice fraction (AICE)
     type(LVT_metadataEntry) :: aice
!EMK...Add sea ice thickness (HI)
     type(LVT_metadataEntry) :: hi

     !list of LSM variables for datastream 1 
     type(LVT_metadataEntry), pointer :: head_ds1_list 
     type(dep), allocatable, dimension(:) :: ptr_into_ds1_list

     !list of LSM variables for datastream 2 
     type(LVT_metadataEntry), pointer :: head_ds2_list 
     type(dep), allocatable, dimension(:) :: ptr_into_ds2_list

     !list of LSM variables for datastream 3
     type(LVT_metadataEntry), pointer :: head_ds3_list 
     type(dep), allocatable, dimension(:) :: ptr_into_ds3_list

     !list of lsm related stats variables
     type(LVT_statsEntry), pointer :: head_stats_list

     type(statsdep), allocatable, dimension(:) :: ptr_into_stats_list

     type(LVT_metadataEntry)                   :: strat_varEntry
  end type output_meta

  type(output_meta)     :: LVT_histData
  type(LVT_metadataEntry)   :: LVT_temp_maxvEntry
  type(LVT_metadataEntry)   :: LVT_temp_minvEntry


contains
 
!BOP
! 
! !ROUTINE: LVT_histDataInit
! \label{LVT_histDataInit}
!
! !INTERFACE: 
  subroutine LVT_histDataInit()
! 
! !USES:   

    implicit none
!
! !DESCRIPTION: 
!  This routine initializes the required data structures to hold the 
!  selected list of variables
!
!   The routines invoked are: 
!   \begin{description}
!    \item[register\_dataEntry] (\ref{register_dataEntry}) \newline
!      allocates the data structures for the selected variables
!    \item[set\_ptr\_into\_list] (\ref{set_ptr_into_list}) \newline
!      assigns the pointers for the selected variables. 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    type(ESMF_Config) :: modelSpecConfig
    integer           :: k
    integer           :: ftn
    integer           :: c,r,i
    integer           :: nsize
    integer           :: rc
    integer           :: iloc
    integer           :: status_m, status_o
    integer           :: max_count
    character*200     :: currentLine
    integer           :: file_lines
    logical           :: file_exists
    logical           :: table_found
    character*20, allocatable :: ds_name(:,:), ds_unit(:,:), ds_dir(:,:)
    integer,      allocatable :: ds_startNlevs(:,:),ds_endNlevs(:,:), &
         ds_vlevels(:,:),ds_timeAvgOpt(:,:)

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%sp_avg_mode, &
         label="Spatial averaging mode:",&
         default="pixel-by-pixel",&
         rc=rc)
    if(rc.ne.0) then 
       write(LVT_logunit,*) "[ERR] Spatial averaging mode: not defined"
       write(LVT_logunit,*) "[ERR] Supported options are: 'pixel-by-pixel' (default)"
       write(LVT_logunit,*) "[ERR] or 'region-based'"
       call LVT_endrun()
    endif
    
    if(LVT_rc%sp_avg_mode.eq."region-based") then 
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%reg_maskfile, &
            label="Regional mask file for spatial averaging:",&
            rc=rc)
       call LVT_verify(rc,&
            'Regional mask file for spatial averaging: not defined')
       
!--------------------------------------------------------------------------
!  The regional mask file must be binary, direct access and in the
!  same domain, resolution, and map projection of the LVT running domain
!--------------------------------------------------------------------------

       allocate(LVT_rc%regmask(LVT_rc%lnc,LVT_rc%lnr))
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=LVT_rc%reg_maskfile,form='unformatted',access='direct',&
            recl=LVT_rc%lnc*LVT_rc%lnr*4)
       read(ftn,rec=1) LVT_rc%regmask
       call LVT_releaseUnitNumber(ftn)
       
       LVT_rc%regmask_max = -1000.0
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             if(LVT_rc%regmask(c,r).gt.LVT_rc%regmask_max) then 
                LVT_rc%regmask_max = LVT_rc%regmask(c,r)
             endif
          enddo
       enddo

    endif

    if(LVT_rc%computeEnsMetrics.eq.1) then 
       do k=1,LVT_rc%nDataStreams
          if(LVT_rc%obssource(k).eq."LIS output") then 
             LVT_rc%npts = LVT_LIS_rc(k)%ntiles
          endif
       enddo
    else
       LVT_rc%npts = LVT_rc%ngrid
    endif
    
    nsize = LVT_rc%npts
    ! if the LIS output is structured in tile space, then the arrays are 
    ! sized to be ntiles. Else the arrays are sized to be ngrid. If the ouput i 
    ! written as 2d grid, then then it will be stored in the 1-d land space. 

    LVT_rc%prev_mo_tavg = LVT_rc%mo

    LVT_histData%head_ds1_list     => null()

    LVT_histData%head_ds2_list     => null()

    LVT_histData%head_stats_list     => null()

    LVT_MOC_COUNT     = 0

    
    ftn = LVT_getNextUnitNumber()   
    open(ftn,file=LVT_rc%configFile, form='formatted')
    file_lines = 0 
    rc = 0 
    table_found = .false. 

    do while(rc.eq.0) 
       read(ftn,'(a)',iostat=rc) currentLine
       if(currentLine.eq."LVT datastream attributes table::") then 
          table_found = .true. 
       endif
       if(table_found) then 
          file_lines = file_lines+1
       endif
       if(currentLine.eq."::") then 
          table_found = .false. 
          exit
       endif
    enddo
    call LVT_releaseUnitNumber(ftn)
    file_lines = file_lines - 2 !minus the top and bottom row

    allocate(ds_name(LVT_rc%nDataStreams,file_lines))
    allocate(ds_startNlevs(LVT_rc%nDataStreams,file_lines))
    allocate(ds_endNlevs(LVT_rc%nDataStreams,file_lines))
    allocate(ds_unit(LVT_rc%nDataStreams,file_lines))
    allocate(ds_dir(LVT_rc%nDataStreams,file_lines))
    allocate(ds_timeAvgOpt(LVT_rc%nDataStreams,file_lines))
    allocate(ds_vlevels(LVT_rc%nDataStreams,file_lines))

    call ESMF_ConfigFindLabel(LVT_config,&
         "LVT datastream attributes table::",rc=rc)
    do i=1,file_lines
       call ESMF_ConfigNextLine(LVT_config, rc=rc)
       do k=1,LVT_rc%nDataStreams
          call ESMF_ConfigGetAttribute(LVT_config, ds_name(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_startNlevs(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_endNlevs(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_unit(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_dir(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_timeAvgOpt(k,i),rc=rc)
          call ESMF_ConfigGetAttribute(LVT_config, ds_vlevels(k,i),rc=rc)
       enddo

       call register_dataEntry(LVT_MOC_COUNT, nsize,&
            LVT_histData%head_ds1_list,&
            LVT_histData%head_ds2_list,&
            LVT_histData%head_ds3_list,&
            LVT_histData%head_stats_list,&
            ds_name(:,i), &
            ds_startNlevs(:,i), &
            ds_endNlevs(:,i),&
            ds_unit(:,i),&
            ds_dir(:,i), &
            ds_timeAvgOpt(:,i), &
            ds_vlevels(:,i))
    enddo

    deallocate(ds_name)
    deallocate(ds_startNlevs)
    deallocate(ds_endNlevs)
    deallocate(ds_unit)
    deallocate(ds_dir)
    deallocate(ds_timeAvgOpt)
    deallocate(ds_vlevels)


    allocate(LVT_histData%ptr_into_ds1_list(LVT_MOC_COUNT(1)))
    allocate(LVT_histData%ptr_into_ds2_list(LVT_MOC_COUNT(2)))
    allocate(LVT_histData%ptr_into_ds3_list(LVT_MOC_COUNT(3)))

    max_count = max(LVT_MOC_COUNT(1),LVT_MOC_COUNT(2),LVT_MOC_COUNT(3))
    allocate(LVT_histData%ptr_into_stats_list(max_count))

    call set_ptr_into_list(LVT_MOC_COUNT, &
                           LVT_histData%head_ds1_list, &
                           LVT_histData%ptr_into_ds1_list,&
                           LVT_histData%head_ds2_list,&
                           LVT_histData%ptr_into_ds2_list,&
                           LVT_histData%head_ds3_list, &
                           LVT_histData%ptr_into_ds3_list)


    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lis_sf_d,&
         label="LVT surface soil layer thickness:",rc=rc)
    call LVT_verify(rc, 'LVT surface soil layer thickness: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lis_rz_d,&
         label="LVT root zone soil layer thickness:",rc=rc)
    call LVT_verify(rc, 'LVT root zone soil layer thickness: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%vinterp_option,&
         label="LVT vertical interpolation option:",default=1,rc=rc)
    call LVT_warning(rc, 'LVT vertical interpolation option: not defined')

    LVT_rc%nsmlayers = 0 
    LVT_rc%nstlayers = 0 
    if(LVT_MOC_SOILMOIST(1).gt.0) then 
       LVT_rc%nsmlayers =  &
            LVT_histData%ptr_into_ds1_list(&
            LVT_MOC_SOILMOIST(1))%dataEntryPtr%vlevels
    endif

    if(LVT_MOC_SMLIQFRAC(1).gt.0) then 
       LVT_rc%nsmlayers =  &
            LVT_histData%ptr_into_ds1_list(&
            LVT_MOC_SMLIQFRAC(1))%dataEntryPtr%vlevels
    endif
    if(LVT_MOC_SOILTEMP(1).gt.0) then 
       LVT_rc%nstlayers =  &
            LVT_histData%ptr_into_ds1_list(&
            LVT_MOC_SOILTEMP(1))%dataEntryPtr%vlevels
    endif

  end subroutine LVT_histDataInit
     
!BOP
! 
! !ROUTINE: register_dataEntry
! \label{register_dataEntry}
!
! !INTERFACE: 
 subroutine register_dataEntry(var_count, nsize, &
       head_ds1Entry, head_ds2Entry, head_ds3Entry, head_statsEntry,&
       ds_name, ds_startNlevs, ds_endNlevs, ds_unit, &
       ds_dir, ds_timeAvgOpt, ds_vlevels)   
   implicit none
!
! !ARGUMENTS: 
!    
   integer                          :: var_count(3)
   integer                          :: nsize
   type(LVT_metadataEntry), pointer :: head_ds1Entry
   type(LVT_metadataEntry), pointer :: head_ds2Entry
   type(LVT_metadataEntry), pointer :: head_ds3Entry
   type(LVT_statsEntry),    pointer :: head_statsEntry
   character(len=*)                 :: ds_name(LVT_rc%nDataStreams)
   integer                          :: ds_startNlevs(LVT_rc%nDataStreams)
   integer                          :: ds_endNlevs(LVT_rc%nDataStreams)
   character(len=*)                 :: ds_unit(LVT_rc%nDataStreams)
   character(len=*)                 :: ds_dir(LVT_rc%nDataStreams)
   integer                          :: ds_timeAvgOpt(LVT_rc%nDataStreams)
   integer                          :: ds_vlevels(LVT_rc%nDataStreams)

!
!
! !DESCRIPTION: 
!  This routine initializes the datastructures required for a 
!  specified output variable. 
! 
!  The arguments are: 
!  \begin{description}
!    \item[var\_count]
!       current variable count (assigned to global index for the
!       data entry)
!    \item[nsize]
!       size of the LVT domain (tile space)
!    \item[head\_ds1Entry]
!       datastream 1 object representing the variable being processed
!    \item[head\_ds2Entry]
!       datastream 2 object representing the variable being processed
!    \item[head\_ds3Entry]
!       datastream 3 object representing the variable being processed
!    \item[head\_statsEntry]
!       stats object representing the variable being processed
!    \item[ds\_name]
!       name of the datastream variables
!    \item[ds\_selectNlevs]
!       number of selected levels of the datastream variables
!    \item[ds\_unit]
!       unit of the datastream variables
!    \item[ds\_dir]
!       direction of the datastream variables
!    \item[ds\_timeAvgOpt]
!       time averaging specification of the datastream variables
!    \item[ds\_vlevels]
!       total number of levels of the datastream variables
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[updateDataEntryMetaData] (\ref{updateDataEntryMetaData}) \newline
!      updates the metadata for each data entry
!   \end{description}
!EOP

   type(LVT_metadataEntry),  pointer :: current
   type(LVT_metadataEntry),  pointer :: ds1dataEntry
   type(LVT_metadataEntry),  pointer :: ds2dataEntry
   type(LVT_metadataEntry),  pointer :: ds3dataEntry
   type(LVT_statsEntry),     pointer :: currentStats
   type(LVT_statsEntry),     pointer :: statsEntry
   logical                           :: enabled
   integer                           :: k

   if(ds_startNlevs(1).ne.0) then 
      var_count(1) = var_count(1) + 1
      allocate(ds1dataEntry)
      if ( .not. associated(head_ds1Entry) ) then
         head_ds1Entry => ds1dataEntry
      else
         current => head_ds1Entry
         do while ( associated(current%next) )
            current => current%next
         enddo
         current%next => ds1dataEntry
      endif
      ds1dataEntry%short_name    = ds_name(1)
      ds1dataEntry%startNlevs    = ds_startNlevs(1)
      ds1dataEntry%endNlevs      = ds_endNlevs(1)
      ds1dataEntry%selectNlevs   = ds_endNlevs(1) - ds_startNlevs(1) + 1
      ds1dataEntry%units         = ds_unit(1)
      ds1dataEntry%dir           = ds_dir(1)
      ds1dataEntry%timeAvgOpt    = ds_timeAvgOpt(1)
      ds1dataEntry%vlevels       = ds_vlevels(1)
      ds1dataEntry%computeVar    = 0 

      call updateDataEntryMetaData(1, ds_name(1), var_count(1), &
           ds1dataEntry)
      
      ds1dataEntry%next => null()

      ds1dataEntry => head_ds1Entry
      do while(associated(ds1dataEntry%next)) 
         ds1dataEntry =>ds1dataEntry%next
      enddo
      ds1dataEntry%index = var_count(1)
      
      allocate(ds1dataEntry%value(nsize,LVT_rc%nensem,ds1dataEntry%vlevels))
      allocate(ds1dataEntry%count(nsize,LVT_rc%nensem,ds1dataEntry%vlevels))
      allocate(ds1dataEntry%count_status(nsize,LVT_rc%nensem,ds1dataEntry%vlevels))
      
      ds1dataEntry%value = 0 
      ds1dataEntry%count = 0 
      ds1dataEntry%count_status = 0 
      
      ds1dataEntry%stdev_flag = .false.
   endif

   if(ds_startNlevs(2).ge.1) then 
      var_count(2) = var_count(2) + 1
      allocate(ds2dataEntry)
      if ( .not. associated(head_ds2Entry) ) then
         head_ds2Entry => ds2dataEntry
      else
         current => head_ds2Entry
         do while ( associated(current%next) )
            current => current%next
         enddo
         current%next => ds2dataEntry
      endif
      ds2dataEntry%short_name    = ds_name(2)
      ds2dataEntry%startNlevs    = ds_startNlevs(2)
      ds2dataEntry%endNlevs      = ds_endNlevs(2)
      ds2dataEntry%selectNlevs   = ds_endNlevs(2) - ds_startNlevs(2) + 1
      ds2dataEntry%units         = ds_unit(2)
      ds2dataEntry%dir           = ds_dir(2)
      ds2dataEntry%timeAvgOpt    = ds_timeAvgOpt(2)
      ds2dataEntry%vlevels       = ds_vlevels(2)
      ds2dataEntry%computeVar    = 0 

      call updateDataEntryMetaData(2, ds_name(2), var_count(2), &
           ds2dataEntry)

      ds2dataEntry%next => null()

      ds2dataEntry => head_ds2Entry
      do while(associated(ds2dataEntry%next)) 
         ds2dataEntry =>ds2dataEntry%next
      enddo
      ds2dataEntry%index = var_count(2)
      
      allocate(ds2dataEntry%value(nsize,LVT_rc%nensem,ds2dataEntry%vlevels))
      allocate(ds2dataEntry%count(nsize,LVT_rc%nensem,ds2dataEntry%vlevels))
      allocate(ds2dataEntry%count_status(nsize,LVT_rc%nensem,ds2dataEntry%vlevels))
      
      ds2dataEntry%value = 0 
      ds2dataEntry%count = 0 
      ds2dataEntry%count_status = 0 
   
      ds2dataEntry%stdev_flag = .false.
   endif

   if(LVT_rc%nDataStreams.gt.2) then 
      if(ds_startNlevs(3).ge.1) then 
         var_count(3) = var_count(3) + 1
         allocate(ds3dataEntry)
         if ( .not. associated(head_ds3Entry) ) then
            head_ds3Entry => ds3dataEntry
         else
            current => head_ds3Entry
            do while ( associated(current%next) )
               current => current%next
            enddo
            current%next => ds3dataEntry
         endif
         ds3dataEntry%short_name    = ds_name(3)
         ds3dataEntry%startNlevs    = ds_startNlevs(3)
         ds3dataEntry%endNlevs      = ds_endNlevs(3)
         ds3dataEntry%selectNlevs   = ds_endNlevs(3) - ds_startNlevs(3) + 1
         ds3dataEntry%units         = ds_unit(3)
         ds3dataEntry%dir           = ds_dir(3)
         ds3dataEntry%timeAvgOpt    = ds_timeAvgOpt(3)
         ds3dataEntry%vlevels       = ds_vlevels(3)
         ds3dataEntry%computeVar    = 0 
         
         call updateDataEntryMetaData(3, ds_name(3), var_count(3), &
              ds3dataEntry)
         
         ds3dataEntry%next => null()
         
         ds3dataEntry => head_ds3Entry
         do while(associated(ds3dataEntry%next)) 
            ds3dataEntry =>ds3dataEntry%next
         enddo
         ds3dataEntry%index = var_count(3)
      
         allocate(ds3dataEntry%value(nsize,LVT_rc%nensem,ds3dataEntry%vlevels))
         allocate(ds3dataEntry%count(nsize,LVT_rc%nensem,ds3dataEntry%vlevels))
         allocate(ds3dataEntry%count_status(nsize,LVT_rc%nensem,ds3dataEntry%vlevels))
         
         ds3dataEntry%value = 0 
         ds3dataEntry%count = 0 
         ds3dataEntry%count_status = 0 
         
         ds3dataEntry%stdev_flag = .false.
      endif
   endif
!stats data structures
   enabled = .false. 
   do k=1,LVT_rc%nDataStreams
      if(ds_startNlevs(k).ge.1) then 
         enabled = .true.
      endif
   enddo

   if(enabled) then 
      allocate(statsEntry)
      if(.not.associated(head_statsEntry)) then 
         head_statsEntry => statsEntry
      else
         currentStats => head_statsEntry
         do while(associated(currentStats%next))
            currentStats => currentStats%next
         enddo
         currentStats%next => statsEntry
      endif
      statsEntry%next => null()
      
      if (LVT_rc%obssource(2).eq."none") then
         statsEntry%standard_name = trim(ds1dataEntry%standard_name)
         statsEntry%long_name     = trim(ds1dataEntry%long_name)
         statsEntry%short_name    = trim(ds_name(1))
         statsEntry%units         = trim(ds_unit(1))
         statsEntry => head_statsEntry
      elseif(LVT_rc%nDataStreams.eq.2) then 
         statsEntry%standard_name = trim(ds1dataEntry%standard_name)//&
              "_v_"//trim(ds2dataEntry%standard_name)
         statsEntry%long_name     = trim(ds1dataEntry%long_name)//&
              "_v_"//trim(ds2dataEntry%long_name)
         statsEntry%short_name    = trim(ds_name(1))//"_v_"//trim(ds_name(2))
         statsEntry%units         = trim(ds_unit(1))//"_v_"//trim(ds_unit(2))
         statsEntry => head_statsEntry
      elseif(LVT_rc%nDataStreams.eq.3) then 
         statsEntry%standard_name = trim(ds1dataEntry%standard_name)//&
              "_v_"//trim(ds2dataEntry%standard_name)//&
              "_v_"//trim(ds3dataEntry%standard_name)
         statsEntry%long_name     = trim(ds1dataEntry%long_name)//&
              "_v_"//trim(ds2dataEntry%long_name)//&
              "_v_"//trim(ds3dataEntry%long_name)
         statsEntry%short_name    = trim(ds_name(1))//"_v_"//trim(ds_name(2))//&
              "_v_"//trim(ds_name(3))
         statsEntry%units         = trim(ds_unit(1))//"_v_"//trim(ds_unit(1))//&
              "_v_"//trim(ds_unit(3))
         statsEntry => head_statsEntry
      endif

      do while(associated(statsEntry%next))
         statsEntry=>statsEntry%next
      enddo
      
      statsEntry%selectOpt = 1
   endif

 end subroutine register_dataEntry

!BOP
! 
! !ROUTINE: updateDataEntryMetaData
! \label{updateDataEntryMetaData}
! 
! !INTERFACE: 
 subroutine updateDataEntryMetaData(source, name, var_count,dataEntry)
! !ARGUMENTS:
   integer                 :: source
   character(len=*)        :: name 
   integer                 :: var_count
   type(LVT_metadataEntry) :: dataEntry
! 
! !DESCRIPTION:
!  This routine populates the dataEntry with the relevant metadata such as
!  the standard name, long name, supported units and supported directions. 
!  The routine also updates the global indices used to track the 
!  variables that are enabled in the LVT instance. 
! 
!  The arguments are: 
!  \begin{description}
!    \item[source]
!       index of the data source (1 or 2)
!    \item[name]
!       short variable name of the data entry
!    \item[var\_count]
!       current variable count (assigned to global index for the
!       data entry)
!    \item[dataEntry]
!       object representing the variable being processed
!  \end{description}
!
!EOP
   if(name.eq."Swnet") then 
      if(LVT_MOC_SWNET(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWNET(source) = var_count
         dataEntry%standard_name = "surface_net_downward_shortwave_flux"
         dataEntry%long_name = "net downward shortwave radiation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Lwnet") then 
      if(LVT_MOC_LWNET(source).eq.LVT_rc%udef) then 
         LVT_MOC_LWNET(source) = var_count
         dataEntry%standard_name = "surface_net_downward_longwave_flux"
         dataEntry%long_name = "net downward longwave radiation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-500.0/)
         dataEntry%valid_max = (/500.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Rnet") then 
      if(LVT_MOC_RNET(source).eq.LVT_rc%udef) then 
         LVT_MOC_RNET(source) = var_count
         dataEntry%standard_name = "net_radiation_flux"
         dataEntry%long_name = "total net radiation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Qle") then 
      if(LVT_MOC_QLE(source).eq.LVT_rc%udef) then 
         LVT_MOC_QLE(source) = var_count
         dataEntry%standard_name = "surface_upward_latent_heat_flux"
         dataEntry%long_name = "latent heat flux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-700.0/)
         dataEntry%valid_max = (/700.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Qh") then 
      if(LVT_MOC_QH(source).eq.LVT_rc%udef) then 
         LVT_MOC_QH(source) = var_count
         dataEntry%standard_name = "surface_upward_sensible_heat_flux"
         dataEntry%long_name = "sensible heat flux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-600.0/)
         dataEntry%valid_max = (/600.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Qg") then 
      if(LVT_MOC_QG(source).eq.LVT_rc%udef) then 
         LVT_MOC_QG(source) = var_count
         dataEntry%standard_name = "downward_heat_flux_in_soil"
         dataEntry%long_name = "soil heat flux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-500.0/)
         dataEntry%valid_max = (/500.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Qf") then 
      if(LVT_MOC_QF(source).eq.LVT_rc%udef) then 
         LVT_MOC_QF(source) = var_count
         dataEntry%standard_name ="energy_of_fusion" 
         dataEntry%long_name = "energy of fusion"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-1200.0/)
         dataEntry%valid_max = (/1200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"S2L","L2S"/)
      endif
   elseif(name.eq."Qv") then 
      if(LVT_MOC_QV(source).eq.LVT_rc%udef) then 
         LVT_MOC_QV(source) = var_count
         dataEntry%standard_name = "surface_snow_sublimation_heat_flux"
         dataEntry%long_name = "energy of sublimation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-600.0/)
         dataEntry%valid_max = (/600.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"S2V","V2S"/)
      endif
   elseif(name.eq."Qtau") then 
      if(LVT_MOC_QTAU(source).eq.LVT_rc%udef) then 
         LVT_MOC_QTAU(source) = var_count
         dataEntry%standard_name = "momentum_flux"
         dataEntry%long_name = "momentum flux"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2   ","m/s    ","kg/m/s2"/)
         dataEntry%valid_min = (/-100.0,-100.0,-100.0/)
         dataEntry%valid_max = (/100.0,100.0,100.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Qa") then 
      if(LVT_MOC_QA(source).eq.LVT_rc%udef) then 
         LVT_MOC_QA(source) = var_count
         dataEntry%standard_name ="advective_energy" 
         dataEntry%long_name = "advective energy"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-50.0/)
         dataEntry%valid_max = (/50.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."DelSurfHeat") then 
      if(LVT_MOC_DELSURFHEAT(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELSURFHEAT(source) = var_count
         dataEntry%standard_name = "change_in_heat_storage"
         dataEntry%long_name = "change in heat storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"J/m2"/)
         dataEntry%valid_min = (/-500.0/)
         dataEntry%valid_max = (/500.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC","DEC"/)
      endif
   elseif(name.eq."DelColdCont") then 
      if(LVT_MOC_DELCOLDCONT(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELCOLDCONT(source) = var_count
         dataEntry%standard_name = "change_in_cold_content"
         dataEntry%long_name = "change in cold content"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"J/m2"/)
         dataEntry%valid_min = (/-600.0/)
         dataEntry%valid_max = (/1000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC","DEC"/)
      endif
   elseif(name.eq."BR") then 
      if(LVT_MOC_BR(source).eq.LVT_rc%udef) then 
         LVT_MOC_BR(source) = var_count
         dataEntry%standard_name ="bowen_ratio" 
         dataEntry%long_name = "bowen ratio"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."EF") then 
      if(LVT_MOC_EF(source).eq.LVT_rc%udef) then 
         LVT_MOC_EF(source) = var_count
         dataEntry%standard_name = "evaporative_fraction"
         dataEntry%long_name = "evaporative fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Snowf") then 
      if(LVT_MOC_SNOWF(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWF(source) = var_count
         dataEntry%standard_name ="snowfall_rate" 
         dataEntry%long_name = "snowfall rate"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0,0.0/)
         dataEntry%valid_max = (/0.0085, 750.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."Rainf") then 
      if(LVT_MOC_RAINF(source).eq.LVT_rc%udef) then 
         LVT_MOC_RAINF(source) = var_count
         dataEntry%standard_name = "rainfall_rate"
         dataEntry%long_name = "rainfall rate"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0,0.0/)
         dataEntry%valid_max = (/0.02, 2000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."RainfConv") then 
      if(LVT_MOC_RAINFCONV(source).eq.LVT_rc%udef) then 
         LVT_MOC_RAINFCONV(source) = var_count
         dataEntry%standard_name = "convective_rainfall_rate"
         dataEntry%long_name = "convective rainfall"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0,0.0/)
         dataEntry%valid_max = (/0.02, 2000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."Evap") then 
      if(LVT_MOC_EVAP(source).eq.LVT_rc%udef) then 
         LVT_MOC_EVAP(source) = var_count
         dataEntry%standard_name ="total_evapotranspiration" 
         dataEntry%long_name = "total evapotranspiration"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 ", "mm/hr "/)
         dataEntry%valid_min = (/-0.0003, -9999.0, -9999.0/)
         dataEntry%valid_max = (/0.0003, -9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."Qs") then 
      if(LVT_MOC_QS(source).eq.LVT_rc%udef) then 
         LVT_MOC_QS(source) = var_count
         dataEntry%standard_name = "surface_runoff_amount"
         dataEntry%long_name = "surface runoff"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s  ","kg/m2   ","mm/month"/)
         dataEntry%valid_min = (/0.0, 0.0, 0.0/)
         dataEntry%valid_max = (/5.0, 43200.0, 43200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."Qsb") then 
      if(LVT_MOC_QSB(source).eq.LVT_rc%udef) then 
         LVT_MOC_QSB(source) = var_count
         dataEntry%standard_name = "subsurface_runoff_amount"
         dataEntry%long_name = "subsurface runoff amount"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s  ","kg/m2   ","mm/month"/)
         dataEntry%valid_min = (/0.0, 0.0, 0.0/)
         dataEntry%valid_max = (/5.0, 43200.0, 43200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."Qrec") then 
      if(LVT_MOC_QREC(source).eq.LVT_rc%udef) then 
         LVT_MOC_QREC(source) = var_count
         dataEntry%standard_name ="recharge" 
         dataEntry%long_name = "recharge"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/5.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."Qsm") then 
      if(LVT_MOC_QSM(source).eq.LVT_rc%udef) then 
         LVT_MOC_QSM(source) = var_count
         dataEntry%standard_name = "snowmelt"
         dataEntry%long_name = "snowmelt"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/0.005, 450.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"S2L", "L2S"/)
      endif
   elseif(name.eq."Qfz") then 
      if(LVT_MOC_QFZ(source).eq.LVT_rc%udef) then 
         LVT_MOC_QFZ(source) = var_count
         dataEntry%standard_name = "refreezing_of_water" 
         dataEntry%long_name = "refreezing water"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/0.005/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"L2S", "S2L"/)
      endif
   elseif(name.eq."Qst") then 
      if(LVT_MOC_QST(source).eq.LVT_rc%udef) then 
         LVT_MOC_QST(source) = var_count
         dataEntry%standard_name = "snow_throughfall" 
         dataEntry%long_name = "snow throughfall"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/0.005/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DelSoilMoist") then 
      if(LVT_MOC_DELSOILMOIST(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELSOILMOIST(source) = var_count
         dataEntry%standard_name = "change_in_soil_moisture"
         dataEntry%long_name = "change in soil moisture"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 "/)
         dataEntry%valid_min = (/-2000.0/)
         dataEntry%valid_max = (/2000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC", "DEC"/)
      endif
   elseif(name.eq."DelSWE") then 
      if(LVT_MOC_DELSWE(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELSWE(source) = var_count
         dataEntry%standard_name = "change_in_swe" 
         dataEntry%long_name = "change in snow water equivalent"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 "/)
         dataEntry%valid_min = (/-2000.0/)
         dataEntry%valid_max = (/2000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC", "DEC"/)
      endif
   elseif(name.eq."DelSurfStor") then 
      if(LVT_MOC_DELSURFSTOR(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELSURFSTOR(source) = var_count
         dataEntry%standard_name = "change_in_surface_water_storage"
         dataEntry%long_name = "change in surface water storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 "/)
         dataEntry%valid_min = (/-2000.0/)
         dataEntry%valid_max = (/2000.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC", "DEC"/)
      endif
   elseif(name.eq."DelIntercept") then 
      if(LVT_MOC_DELINTERCEPT(source).eq.LVT_rc%udef) then 
         LVT_MOC_DELINTERCEPT(source) = var_count
         dataEntry%standard_name = "change_in_interception_storage"
         dataEntry%long_name = "change in interception storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 "/)
         dataEntry%valid_min = (/-100.0/)
         dataEntry%valid_max = (/100.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"INC", "DEC"/)
      endif
   elseif(name.eq."SnowT") then 
      if(LVT_MOC_QSM(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWT(source) = var_count
         dataEntry%standard_name = "temperature_in_surface_snow"
         dataEntry%long_name = "temperature in surface snow"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/280.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VegT") then 
      if(LVT_MOC_VEGT(source).eq.LVT_rc%udef) then 
         LVT_MOC_VEGT(source) = var_count
         dataEntry%standard_name ="canopy_temperature" 
         dataEntry%long_name = "canopy temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."BareSoilT") then 
      if(LVT_MOC_BARESOILT(source).eq.LVT_rc%udef) then 
         LVT_MOC_BARESOILT(source) = var_count
         dataEntry%standard_name = "bare_soil_temperature"
         dataEntry%long_name = "bare soil temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/343.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AvgSurfT") then 
      if(LVT_MOC_AVGSURFT(source).eq.LVT_rc%udef) then 
         LVT_MOC_AVGSURFT(source) = var_count
         dataEntry%standard_name = "surface_temperature"
         dataEntry%long_name = "surface temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AvgGrndT") then 
      if(LVT_MOC_GROUNDAVGT(source).eq.LVT_rc%udef) then 
         LVT_MOC_GROUNDAVGT(source) = var_count
         dataEntry%standard_name = "ground_surface_temperature"
         dataEntry%long_name = "ground surface temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VegGrndT") then 
      if(LVT_MOC_GROUNDVEGT(source).eq.LVT_rc%udef) then 
         LVT_MOC_GROUNDVEGT(source) = var_count
         dataEntry%standard_name = "vegetated_ground_surface_temperature"
         dataEntry%long_name = "vegetated ground surface temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RadT") then 
      if(LVT_MOC_RADT(source).eq.LVT_rc%udef) then 
         LVT_MOC_RADT(source) = var_count
         dataEntry%standard_name = "surface_radiative_temperature"
         dataEntry%long_name = "surface radiative temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/353.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Albedo") then 
      if(LVT_MOC_ALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_ALBEDO(source) = var_count
         dataEntry%standard_name = "surface_albedo" 
         dataEntry%long_name = "surface albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoVisDir") then 
      if(LVT_MOC_VISDIRALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_VISDIRALBEDO(source) = var_count
         dataEntry%standard_name = "direct_visible_albedo" 
         dataEntry%long_name = "direct visible albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoVisDif") then 
      if(LVT_MOC_VISDIFALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_VISDIFALBEDO(source) = var_count
         dataEntry%standard_name = "diffused_visible_albedo" 
         dataEntry%long_name = "diffused visible albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoNIRDir") then 
      if(LVT_MOC_NIRDIRALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_NIRDIRALBEDO(source) = var_count
         dataEntry%standard_name = "direct_NIR_albedo" 
         dataEntry%long_name = "direct NIR albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoNIRDif") then 
      if(LVT_MOC_NIRDIFALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_NIRDIFALBEDO(source) = var_count
         dataEntry%standard_name = "diffused_NIR_albedo" 
         dataEntry%long_name = "diffused NIR albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoSWDir") then 
      if(LVT_MOC_SWDIRALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWDIRALBEDO(source) = var_count
         dataEntry%standard_name = "direct_SW_albedo" 
         dataEntry%long_name = "direct SW albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AlbedoSWDif") then 
      if(LVT_MOC_SWDIFALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWDIFALBEDO(source) = var_count
         dataEntry%standard_name = "diffused_SW_albedo" 
         dataEntry%long_name = "diffused SW albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SWE") then 
      if(LVT_MOC_SWE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWE(source) = var_count
         dataEntry%standard_name = "liquid_water_content_of_surface_snow" 
         dataEntry%long_name ="snow water equivalent" 
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2","m    "/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/2000.0, 2.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SnowIce") then 
      if(LVT_MOC_SNOWICE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWICE(source) = var_count
         dataEntry%standard_name = "snow_ice"
         dataEntry%long_name ="snow ice"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2","mm   "/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/2000.0, 2.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SnowDepth") then 
      if(LVT_MOC_SNOWDEPTH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWDEPTH(source) = var_count
         dataEntry%standard_name = "snow_depth"
         dataEntry%long_name ="snow depth" 
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m ","cm","mm"/)
         dataEntry%valid_min = (/0.0, 0.0, 0.0 /)
         dataEntry%valid_max = (/10.0, 1000.0, 10000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SWEVeg") then 
      if(LVT_MOC_SWEVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWEVEG(source) = var_count
         dataEntry%standard_name = "swe_intercepted_by_vegetation" 
         dataEntry%long_name = "swe intercepted by vegetation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SnowAge") then 
      if(LVT_MOC_SNOWAGE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWAGE(source) = var_count
         dataEntry%standard_name = "snow_age"
         dataEntry%long_name = "snow age"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SurfStor") then 
      if(LVT_MOC_SURFSTOR(source).eq.LVT_rc%udef) then 
         LVT_MOC_SURFSTOR(source) = var_count
         dataEntry%standard_name = "surface_water_storage"
         dataEntry%long_name = "surface water storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/2000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SoilMoist") then 
      if(LVT_MOC_SOILMOIST(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILMOIST(source) = var_count
         dataEntry%standard_name = "soil_moisture_content"
         dataEntry%long_name = "soil moisture content"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2","m3/m3"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/2000.0, 0.5/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SoilTemp") then 
      if(LVT_MOC_SOILTEMP(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILTEMP(source) = var_count
         dataEntry%standard_name = "soil_temperature"
         dataEntry%long_name = "soil temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SmLiqFrac") then 
      if(LVT_MOC_SMLIQFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SMLIQFRAC(source) = var_count
         dataEntry%standard_name = "liquid_fraction_of_soil_moisture"
         dataEntry%long_name = "average layer fraction of liquid moisture"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-    ","m3/m3"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 0.5/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SmFrozFrac") then 
      if(LVT_MOC_SMFROZFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SMFROZFRAC(source) = var_count
         dataEntry%standard_name = "frozen_fraction_of_soil_moisture"
         dataEntry%long_name = "average layer fraction of frozen moisture"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-    ","m3/m3"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 0.5/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SoilWet") then 
      if(LVT_MOC_SOILWET(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILWET(source) = var_count
         dataEntry%standard_name = "total_soil_wetness" 
         dataEntry%long_name = "total soil wetness"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"- ","% ","mm"/)
         dataEntry%valid_min = (/-0.2, -9999.0, -9999.0/)
         dataEntry%valid_max = (/1.2, -9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."MatricPotential") then 
      if(LVT_MOC_MATRICPOTENTIAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_MATRICPOTENTIAL(source) = var_count
         dataEntry%standard_name = "matric_potential" 
         dataEntry%long_name = "matric potential"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m ","mm"/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."PotEvap") then 
      if(LVT_MOC_POTEVAP(source).eq.LVT_rc%udef) then 
         LVT_MOC_POTEVAP(source) = var_count
         dataEntry%standard_name ="potential_evapotranspiration" 
         dataEntry%long_name = "potential evapotranspiration"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s", "mm/hr ","W/m2  "/)
         dataEntry%valid_min = (/-0.0006, -9999.0, -9999.0/)
         dataEntry%valid_max = (/0.0006, -9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."ECanop") then 
      if(LVT_MOC_ECANOP(source).eq.LVT_rc%udef) then 
         LVT_MOC_ECANOP(source) = var_count
         dataEntry%standard_name = "interception_evaporation"
         dataEntry%long_name = "interception evaporation"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s", "W/m2  "/)
         dataEntry%valid_min = (/-0.0003, -9999.0/)
         dataEntry%valid_max = (/0.0003, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."TVeg") then 
      if(LVT_MOC_TVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_TVEG(source) = var_count
         dataEntry%standard_name = "vegetation_transpiration"
         dataEntry%long_name = "vegetation transpiration"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s", "W/m2  "/)
         dataEntry%valid_min = (/-0.0003,-9999.0/)
         dataEntry%valid_max = (/0.0003, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."ESoil") then 
      if(LVT_MOC_ESOIL(source).eq.LVT_rc%udef) then 
         LVT_MOC_ESOIL(source) = var_count
         dataEntry%standard_name = "bare_soil_evaporation" 
         dataEntry%long_name = "bare soil evaporation"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s", "W/m2  "/)
         dataEntry%valid_min = (/-0.0003,-9999.0/)
         dataEntry%valid_max = (/0.0003, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."EWater") then 
      if(LVT_MOC_EWATER(source).eq.LVT_rc%udef) then 
         LVT_MOC_EWATER(source) = var_count
         dataEntry%standard_name ="open_water_evaporation" 
         dataEntry%long_name = "open water evaporation"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s", "W/m2  "/)
         dataEntry%valid_min = (/-0.0003,-9999.0/)
         dataEntry%valid_max = (/0.0003, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP", "DN"/)
      endif
   elseif(name.eq."RootMoist") then 
      if(LVT_MOC_ROOTMOIST(source).eq.LVT_rc%udef) then 
         LVT_MOC_ROOTMOIST(source) = var_count
         dataEntry%standard_name = "root_zone_soil_moisture"
         dataEntry%long_name = "root zone soil moisture"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/m3","kg/m2"/)
         dataEntry%valid_min = (/0.0,0.0/)
         dataEntry%valid_max = (/0.5, 2000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."CanopInt") then 
      if(LVT_MOC_CANOPINT(source).eq.LVT_rc%udef) then 
         LVT_MOC_CANOPINT(source) = var_count
         dataEntry%standard_name = "total_canopy_water_storage"
         dataEntry%long_name =  "total canopy water storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."EvapSnow") then 
      if(LVT_MOC_EVAPSNOW(source).eq.LVT_rc%udef) then 
         LVT_MOC_EVAPSNOW(source) = var_count
         dataEntry%standard_name = "snow_evaporation"
         dataEntry%long_name = "snow evaporation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-0.003/)
         dataEntry%valid_max = (/0.003/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SubSnow") then 
      if(LVT_MOC_SUBSNOW(source).eq.LVT_rc%udef) then 
         LVT_MOC_SUBSNOW(source) = var_count
         dataEntry%standard_name ="snow_sublimation" 
         dataEntry%long_name = "snow sublimation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-0.003/)
         dataEntry%valid_max = (/0.003/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SubSurf") then 
      if(LVT_MOC_SUBSURF(source).eq.LVT_rc%udef) then 
         LVT_MOC_SUBSURF(source) = var_count
         dataEntry%standard_name ="sublimation_of_the_snow_free_area" 
         dataEntry%long_name =  "sublimation of the snow free area"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-0.003/)
         dataEntry%valid_max = (/0.003/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ACond") then 
      if(LVT_MOC_ACOND(source).eq.LVT_rc%udef) then 
         LVT_MOC_ACOND(source) = var_count
         dataEntry%standard_name = "aerodynamic_conductance"
         dataEntry%long_name = "aerodynamic conductance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."CCond") then 
      if(LVT_MOC_CCOND(source).eq.LVT_rc%udef) then 
         LVT_MOC_CCOND(source) = var_count
         dataEntry%standard_name = "canopy_conductance"
         dataEntry%long_name = "canopy conductance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"s/m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VPD") then 
      if(LVT_MOC_VPD(source).eq.LVT_rc%udef) then 
         LVT_MOC_VPD(source) = var_count
         dataEntry%standard_name = "vapor_pressure_deficit"
         dataEntry%long_name = "vapor pressure deficit"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"Pa"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WaterTableD") then 
      if(LVT_MOC_WATERTABLED(source).eq.LVT_rc%udef) then 
         LVT_MOC_WATERTABLED(source) = var_count
         dataEntry%standard_name = "water_table_depth"
         dataEntry%long_name = "water table depth"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm", "m "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TWS") then 
      if(LVT_MOC_TWS(source).eq.LVT_rc%udef) then 
         LVT_MOC_TWS(source) = var_count
         dataEntry%standard_name = "terrestrial_water_storage" 
         dataEntry%long_name = "terrestrial water storage"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm", "m "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."GWS") then 
      if(LVT_MOC_GWS(source).eq.LVT_rc%udef) then 
         LVT_MOC_GWS(source) = var_count
         dataEntry%standard_name = "ground_water_storage" 
         dataEntry%long_name = "ground water storage"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm", "m "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Snowcover") then 
      if(LVT_MOC_SNOWCOVER(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWCOVER(source) = var_count
         dataEntry%standard_name = "surface_snow_area_fraction"
         dataEntry%long_name = "surface snow area fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SAlbedo") then 
      if(LVT_MOC_SALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_SALBEDO(source) = var_count
         dataEntry%standard_name = "snow_albedo"
         dataEntry%long_name = "snow albedo"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SnowTProf") then 
      if(LVT_MOC_SNOWTPROF(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWTPROF(source) = var_count
         dataEntry%standard_name ="snow_temperature_profile" 
         dataEntry%long_name = "snow temperature profile"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SLiqFrac") then 
      if(LVT_MOC_SLIQFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SLIQFRAC(source) = var_count
         dataEntry%standard_name ="snow_liquid_fraction_on_ground" 
         dataEntry%long_name = "snow liquid fraction on ground"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."LWup") then 
      if(LVT_MOC_LWUP(source).eq.LVT_rc%udef) then 
         LVT_MOC_LWUP(source) = var_count
         dataEntry%standard_name ="longwave_radiation_up" 
         dataEntry%long_name = "longwave radiation up"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."GPP") then 
      if(LVT_MOC_GPP(source).eq.LVT_rc%udef) then 
         LVT_MOC_GPP(source) = var_count
         dataEntry%standard_name ="gross_primary_production" 
         dataEntry%long_name = "gross primary production"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ", "umol/m2s", "g/m2s   "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."NPP") then 
      if(LVT_MOC_NPP(source).eq.LVT_rc%udef) then 
         LVT_MOC_NPP(source) = var_count
         dataEntry%standard_name ="net_primary_production" 
         dataEntry%long_name = "net primary production"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ", "umol/m2s", "g/m2s   "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."NEE") then 
      if(LVT_MOC_NEE(source).eq.LVT_rc%udef) then 
         LVT_MOC_NEE(source) = var_count
         dataEntry%standard_name ="net_ecosystem_exchange" 
         dataEntry%long_name = "net ecosystem exchange"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ", "umol/m2s", "g/m2s   "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."AutoResp") then 
      if(LVT_MOC_AUTORESP(source).eq.LVT_rc%udef) then 
         LVT_MOC_AUTORESP(source) = var_count
         dataEntry%standard_name ="autoresp" 
         dataEntry%long_name = "autoresp"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ","umol/m2s"/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
 elseif(name.eq."HeteroResp") then 
      if(LVT_MOC_HETERORESP(source).eq.LVT_rc%udef) then 
         LVT_MOC_HETERORESP(source) = var_count
         dataEntry%standard_name ="heteroresp" 
         dataEntry%long_name = "heteroresp"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ","umol/m2s"/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
 elseif(name.eq."LeafResp") then 
      if(LVT_MOC_LEAFRESP(source).eq.LVT_rc%udef) then 
         LVT_MOC_LEAFRESP(source) = var_count
         dataEntry%standard_name ="leafresp" 
         dataEntry%long_name = "leafresp"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s2 ","umol/m2s"/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."TotSoilCarb") then 
      if(LVT_MOC_TOTSOILCARB(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTSOILCARB(source) = var_count
         dataEntry%standard_name ="total_soil_and_litter_carbon_content" 
         dataEntry%long_name = "total soil and litter carbon content" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotLivBiom") then 
      if(LVT_MOC_TOTLIVBIOM(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTLIVBIOM(source) = var_count
         dataEntry%standard_name ="total_living_biomass_carbon_content"
         dataEntry%long_name = "total living biomass carbon content"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SoilET") then 
      if(LVT_MOC_SOILET(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILET(source) = var_count
         dataEntry%standard_name ="soil_evaporation" 
         dataEntry%long_name = "soil evaporation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Z0brd") then 
      if(LVT_MOC_Z0BRD(source).eq.LVT_rc%udef) then 
         LVT_MOC_Z0BRD(source) = var_count
         dataEntry%standard_name ="z0brd" 
         dataEntry%long_name = "z0brd"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Ch") then 
      if(LVT_MOC_CH(source).eq.LVT_rc%udef) then 
         LVT_MOC_CH(source) = var_count
         dataEntry%standard_name ="heat_exchange_coefficient" 
         dataEntry%long_name = "heat exchange coefficient"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Cm") then 
      if(LVT_MOC_CM(source).eq.LVT_rc%udef) then 
         LVT_MOC_CM(source) = var_count
         dataEntry%standard_name ="momentum_exchange_coefficient" 
         dataEntry%long_name = "momentum exchange coefficient"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."T2diag") then 
      if(LVT_MOC_T2DIAG(source).eq.LVT_rc%udef) then 
         LVT_MOC_T2DIAG(source) = var_count
         dataEntry%standard_name ="diagnostic t2" 
         dataEntry%long_name = "diagnostic t2"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Q2diag") then 
      if(LVT_MOC_Q2DIAG(source).eq.LVT_rc%udef) then 
         LVT_MOC_Q2DIAG(source) = var_count
         dataEntry%standard_name ="diagnostic q2" 
         dataEntry%long_name = "diagnostic q2"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/kg"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RootTemp") then 
      if(LVT_MOC_ROOTTEMP(source).eq.LVT_rc%udef) then 
         LVT_MOC_ROOTTEMP(source) = var_count
         dataEntry%standard_name ="root_zone_temperature" 
         dataEntry%long_name = "root zone temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Wind_f") then 
      if(LVT_MOC_WINDFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_WINDFORC(source) = var_count
         dataEntry%standard_name ="wind_speed" 
         dataEntry%long_name = "wind speed"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s   ","km/day"/)
         dataEntry%valid_min = (/-75.0, -6500.0/)
         dataEntry%valid_max = (/75.0, 6500.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Rainf_f") then 
      if(LVT_MOC_RAINFFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_RAINFFORC(source) = var_count
         dataEntry%standard_name ="rainfall_flux" 
         dataEntry%long_name = "rainfall flux"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0, -9999.0/)
         dataEntry%valid_max = (/0.02,-9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Snowf_f") then 
      if(LVT_MOC_SNOWFFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWFFORC(source) = var_count
         dataEntry%standard_name ="snowfall_flux" 
         dataEntry%long_name = "snowfall flux"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0, -9999.0/)
         dataEntry%valid_max = (/0.02,-9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."CRainf_f") then 
      if(LVT_MOC_CRAINFFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_CRAINFFORC(source) = var_count
         dataEntry%standard_name ="convective_rainfall_flux" 
         dataEntry%long_name = "convective rainfall flux"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/0.0, -9999.0/)
         dataEntry%valid_max = (/0.02,-9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Tair_f") then 
      if(LVT_MOC_TAIRFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_TAIRFORC(source) = var_count
         dataEntry%standard_name ="air_temperature" 
         dataEntry%long_name = "air temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Qair_f") then 
      if(LVT_MOC_QAIRFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_QAIRFORC(source) = var_count
         dataEntry%standard_name ="specific_humidity" 
         dataEntry%long_name = "specific humidity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/kg"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/0.03/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Psurf_f") then 
      if(LVT_MOC_PSURFFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_PSURFFORC(source) = var_count
         dataEntry%standard_name ="surface_air_pressure" 
         dataEntry%long_name = "surface pressure"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"Pa"/)
         dataEntry%valid_min = (/5000.0/)
         dataEntry%valid_max = (/11000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SWdown_f") then 
      if(LVT_MOC_SWDOWNFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWDOWNFORC(source) = var_count
         dataEntry%standard_name ="surface_downwelling_shortwave_flux_in_air" 
         dataEntry%long_name = "downward shortwave radiation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1360.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."LWdown_f") then 
      if(LVT_MOC_LWDOWNFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_LWDOWNFORC(source) = var_count
         dataEntry%standard_name ="surface_downwelling_longwave_flux_in_air" 
         dataEntry%long_name = "downward longwave radiation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/750.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."DirectSW_f") then 
      if(LVT_MOC_DIRECTSWFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_DIRECTSWFORC(source) = var_count
         dataEntry%standard_name ="surface_direct_downwelling_shortwave_flux_in_air" 
         dataEntry%long_name = "direct shortwave flux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1360.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DiffuseSW_f") then 
      if(LVT_MOC_DIFFUSESWFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_DIFFUSESWFORC(source) = var_count
         dataEntry%standard_name ="surface_diffuse_downwelling_shortwave_flux_in_air" 
         dataEntry%long_name = "diffuse shortwave flux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1360.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."NWind_f") then 
      if(LVT_MOC_NWINDFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_NWINDFORC(source) = var_count
         dataEntry%standard_name ="northward_wind" 
         dataEntry%long_name = "northward wind"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/-75.0/)
         dataEntry%valid_max = (/75.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"N","E"/)
      endif
   elseif(name.eq."EWind_f") then 
      if(LVT_MOC_EWINDFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_EWINDFORC(source) = var_count
         dataEntry%standard_name ="eastward_wind" 
         dataEntry%long_name = "eastward wind"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/-75.0/)
         dataEntry%valid_max = (/75.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"N","E"/)
      endif
   elseif(name.eq."FHeight_f") then 
      if(LVT_MOC_FHEIGHTFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_FHEIGHTFORC(source) = var_count
         dataEntry%standard_name ="forcing_height" 
         dataEntry%long_name = "forcing height"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Ch_f") then 
      if(LVT_MOC_CHFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_CHFORC(source) = var_count
         dataEntry%standard_name ="heat_exchange_coefficient" 
         dataEntry%long_name = "heat exchange coefficient"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Cm_f") then 
      if(LVT_MOC_CMFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_CMFORC(source) = var_count
         dataEntry%standard_name ="momentum_exchange_coefficient" 
         dataEntry%long_name = "momentum exchange coefficient"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."MixRatio_f") then 
      if(LVT_MOC_MIXRATIOFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_MIXRATIOFORC(source) = var_count
         dataEntry%standard_name ="mixing_ratio" 
         dataEntry%long_name = "mixing ratio"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."CosZenith_f") then 
      if(LVT_MOC_COSZENFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_COSZENFORC(source) = var_count
         dataEntry%standard_name ="cosine_solar_zenith_angle" 
         dataEntry%long_name = "cosine of solar zenith angle"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Albedo_f") then 
      if(LVT_MOC_ALBEDOFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_ALBEDOFORC(source) = var_count
         dataEntry%standard_name = "albedo" 
         dataEntry%long_name = "albedo"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."PARDR_f") then 
      if(LVT_MOC_PARDRFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_PARDRFORC(source) = var_count
         dataEntry%standard_name ="surface_downward_PAR_direct" 
         dataEntry%long_name = "surface downward PAR direct"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."PARDF_f") then 
      if(LVT_MOC_PARDFFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_PARDFFORC(source) = var_count
         dataEntry%standard_name ="surface_downward_PAR_diffuse" 
         dataEntry%long_name = "surface downward PAR diffuse"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Drought_Category") then 
      if(LVT_MOC_DR_CATEGORY(source).eq.LVT_rc%udef) then 
         LVT_MOC_DR_CATEGORY(source) = var_count
         dataEntry%standard_name ="drought_category" 
         dataEntry%long_name = "drought category"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Percentile") then 
      if(LVT_MOC_PERCENTILE(source).eq.LVT_rc%udef) then 
         LVT_MOC_PERCENTILE(source) = var_count
         dataEntry%standard_name ="percentile"
         dataEntry%long_name = "percentile"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Landmask") then 
      if(LVT_MOC_LANDMASK(source).eq.LVT_rc%udef) then 
         LVT_MOC_LANDMASK(source) = var_count
         dataEntry%standard_name ="landmask" 
         dataEntry%long_name = "landmask"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Landcover") then 
      if(LVT_MOC_LANDCOVER(source).eq.LVT_rc%udef) then 
         LVT_MOC_LANDCOVER(source) = var_count
         dataEntry%standard_name ="landcover" 
         dataEntry%long_name = "landcover"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Soiltype") then 
      if(LVT_MOC_SOILTYPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILTYPE(source) = var_count
         dataEntry%standard_name ="soiltype" 
         dataEntry%long_name = "soiltype"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SandFrac") then 
      if(LVT_MOC_SANDFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SANDFRAC(source) = var_count
         dataEntry%standard_name ="sand_fraction" 
         dataEntry%long_name = "sand fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ClayFrac") then 
      if(LVT_MOC_CLAYFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_CLAYFRAC(source) = var_count
         dataEntry%standard_name ="clay_fraction" 
         dataEntry%long_name = "clay fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SiltFrac") then 
      if(LVT_MOC_SILTFRAC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SILTFRAC(source) = var_count
         dataEntry%standard_name ="silt_fraction" 
         dataEntry%long_name = "silt fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Porosity") then 
      if(LVT_MOC_POROSITY(source).eq.LVT_rc%udef) then 
         LVT_MOC_POROSITY(source) = var_count
         dataEntry%standard_name ="porosity" 
         dataEntry%long_name = "porosity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Soilcolor") then 
      if(LVT_MOC_SOILCOLOR(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOILCOLOR(source) = var_count
         dataEntry%standard_name ="soil_color" 
         dataEntry%long_name = "soil color"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Elevation") then 
      if(LVT_MOC_ELEVATION(source).eq.LVT_rc%udef) then 
         LVT_MOC_ELEVATION(source) = var_count
         dataEntry%standard_name ="elevation" 
         dataEntry%long_name = "elevation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Slope") then 
      if(LVT_MOC_SLOPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SLOPE(source) = var_count
         dataEntry%standard_name ="slope" 
         dataEntry%long_name = "slope"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Aspect") then 
      if(LVT_MOC_ASPECT(source).eq.LVT_rc%udef) then 
         LVT_MOC_ASPECT(source) = var_count
         dataEntry%standard_name ="aspect" 
         dataEntry%long_name = "aspect"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."LAI") then 
      if(LVT_MOC_LAI(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAI(source) = var_count
         dataEntry%standard_name ="leaf_area_index" 
         dataEntry%long_name = "leaf area index"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SAI") then 
      if(LVT_MOC_SAI(source).eq.LVT_rc%udef) then 
         LVT_MOC_SAI(source) = var_count
         dataEntry%standard_name ="stem_area_index" 
         dataEntry%long_name = "stem area index"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Snfralbedo") then 
      if(LVT_MOC_SNFRALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNFRALBEDO(source) = var_count
         dataEntry%standard_name ="snow_free_albedo"
         dataEntry%long_name = "snow free albedo"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Mxsnalbedo") then 
      if(LVT_MOC_MXSNALBEDO(source).eq.LVT_rc%udef) then 
         LVT_MOC_MXSNALBEDO(source) = var_count
         dataEntry%standard_name ="maximum_snow_free_albedo"
         dataEntry%long_name = "maximum snow free albedo"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Greenness") then
      if(LVT_MOC_GREENNESS(source).eq.LVT_rc%udef) then 
         LVT_MOC_GREENNESS(source) = var_count
         dataEntry%standard_name ="green_vegetation_fraction"
         dataEntry%long_name = "green vegetation fraction"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/0.0,0.0/)
         dataEntry%valid_max = (/1.0,100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."NDVI") then 
      if(LVT_MOC_NDVI(source).eq.LVT_rc%udef) then 
         LVT_MOC_NDVI(source) = var_count
         dataEntry%standard_name ="normalized_difference_vegetation_index"
         dataEntry%long_name = "normalized difference vegetation index"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VWC") then ! MN 
      if(LVT_MOC_VEGWATERCONTENT(source).eq.LVT_rc%udef) then 
         LVT_MOC_VEGWATERCONTENT(source) = var_count
         dataEntry%standard_name = "vegetation_water_content" 
         dataEntry%long_name ="vegetation water content" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VOD") then 
      if(LVT_MOC_VOD(source).eq.LVT_rc%udef) then 
         LVT_MOC_VOD(source) = var_count
         dataEntry%standard_name = "vegetation_optical_depth" 
         dataEntry%long_name ="vegetation optical depth" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/100.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
#if 0
elseif(name.eq."SMAPL3TB") then ! MN 
      if(LVT_MOC_L3TB(source).eq.LVT_rc%udef) then 
         LVT_MOC_L3TB(source) = var_count
         dataEntry%standard_name = "brightness_temperature" 
         dataEntry%long_name ="brightness temperature" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/) 
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
#endif
elseif(name.eq."SMAPL3TBv_D") then ! MN 
      if(LVT_MOC_L3TBv_D(source).eq.LVT_rc%udef) then 
         LVT_MOC_L3TBv_D(source) = var_count
         dataEntry%standard_name = "brightness_temperature" 
         dataEntry%long_name ="brightness temperature" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/) 
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
elseif(name.eq."SMAPL3TBv_A") then ! MN 
      if(LVT_MOC_L3TBv_A(source).eq.LVT_rc%udef) then 
         LVT_MOC_L3TBv_A(source) = var_count
         dataEntry%standard_name = "brightness_temperature" 
         dataEntry%long_name ="brightness temperature" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/) 
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
elseif(name.eq."SMAPL3TBh_D") then ! MN 
      if(LVT_MOC_L3TBh_D(source).eq.LVT_rc%udef) then 
         LVT_MOC_L3TBh_D(source) = var_count
         dataEntry%standard_name = "brightness_temperature" 
         dataEntry%long_name ="brightness temperature" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/) 
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif

elseif(name.eq."SMAPL3TBh_A") then ! MN 
      if(LVT_MOC_L3TBh_A(source).eq.LVT_rc%udef) then 
         LVT_MOC_L3TBh_A(source) = var_count
         dataEntry%standard_name = "brightness_temperature" 
         dataEntry%long_name ="brightness temperature" 
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/0.0/) 
         dataEntry%valid_max = (/350.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif

   elseif(name.eq."SIF") then 
      if(LVT_MOC_SIF(source).eq.LVT_rc%udef) then 
         LVT_MOC_SIF(source) = var_count
         dataEntry%standard_name ="solar_induced_fluorescence"
         dataEntry%long_name = "solar induced fluorescence"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mW/m^2/nm/sr"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Tempbot") then 
      if(LVT_MOC_TEMPBOT(source).eq.LVT_rc%udef) then 
         LVT_MOC_TEMPBOT(source) = var_count
         dataEntry%standard_name ="bottom_temperature"
         dataEntry%long_name = "bottom temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/333.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."PET_f") then 
      if(LVT_MOC_PETFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_PETFORC(source) = var_count
         dataEntry%standard_name ="potential evaporation"
         dataEntry%long_name = "potential evaporation"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RefET_f") then 
      if(LVT_MOC_REFETFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_REFETFORC(source) = var_count
         dataEntry%standard_name ="reference_ET_forcing"
         dataEntry%long_name = "reference ET forcing"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RefET") then 
      if(LVT_MOC_REFET(source).eq.LVT_rc%udef) then 
         LVT_MOC_REFET(source) = var_count
         dataEntry%standard_name ="reference_ET"
         dataEntry%long_name = "reference ET"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."CAPE_f") then 
      if(LVT_MOC_CAPEFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_CAPEFORC(source) = var_count
         dataEntry%standard_name ="convective_available_potential_energy"
         dataEntry%long_name = "convective available potential energy"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"J/kg"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ETa") then 
      if(LVT_MOC_ETa(source).eq.LVT_rc%udef) then 
         LVT_MOC_ETa(source) = var_count
         dataEntry%standard_name ="ET_anomaly"
         dataEntry%long_name = "ET anomaly"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SOS") then 
      if(LVT_MOC_SOS(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOS(source) = var_count
         dataEntry%standard_name ="start_of_season"
         dataEntry%long_name = "start of season"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WRSI") then 
      if(LVT_MOC_WRSI(source).eq.LVT_rc%udef) then 
         LVT_MOC_WRSI(source) = var_count
         dataEntry%standard_name ="water_requirements_satisfaction_index"
         dataEntry%long_name = "water requirements satisfaction index"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."KF2") then 
      if(LVT_MOC_KF2(source).eq.LVT_rc%udef) then 
         LVT_MOC_KF2(source) = var_count
         dataEntry%standard_name ="kf2"
         dataEntry%long_name = "kf2"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SumWR") then 
      if(LVT_MOC_SUMWR(source).eq.LVT_rc%udef) then 
         LVT_MOC_SUMWR(source) = var_count
         dataEntry%standard_name ="SumWR"
         dataEntry%long_name = "SumWR"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SumET") then 
      if(LVT_MOC_SUMET(source).eq.LVT_rc%udef) then 
         LVT_MOC_SUMET(source) = var_count
         dataEntry%standard_name ="SumET"
         dataEntry%long_name = "SumET"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SWI") then 
      if(LVT_MOC_SWI(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWI(source) = var_count
         dataEntry%standard_name ="soil_water_index"
         dataEntry%long_name = "soil_water_index"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SOSa") then 
      if(LVT_MOC_SOSA(source).eq.LVT_rc%udef) then 
         LVT_MOC_SOSA(source) = var_count
         dataEntry%standard_name ="SOSa"
         dataEntry%long_name = "SOSa"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWater") then 
      if(LVT_MOC_TOTALSURPLUSWATER(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALSURPLUSWATER(source) = var_count
         dataEntry%standard_name ="total_surplus_water"
         dataEntry%long_name = "total surplus water"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."MaxSurplusWater") then 
      if(LVT_MOC_MAXSURPLUSWATER(source).eq.LVT_rc%udef) then 
         LVT_MOC_MAXSURPLUSWATER(source) = var_count
         dataEntry%standard_name ="max_surplus_water"
         dataEntry%long_name = "max surplus water"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWaterDeficit") then 
      if(LVT_MOC_TOTALWATERDEFICIT(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWATERDEFICIT(source) = var_count
         dataEntry%standard_name ="total_water_deficit"
         dataEntry%long_name = "total water deficit"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."MaxWaterDeficit") then 
      if(LVT_MOC_MAXWATERDEFICIT(source).eq.LVT_rc%udef) then 
         LVT_MOC_MAXWATERDEFICIT(source) = var_count
         dataEntry%standard_name ="max_water_deficit"
         dataEntry%long_name = "max water deficit"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalAETInitial") then 
      if(LVT_MOC_TOTALAETINITIAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALAETINITIAL(source) = var_count
         dataEntry%standard_name ="total_AET_initial"
         dataEntry%long_name = "total AET initial"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWRInitial") then 
      if(LVT_MOC_TOTALWRINITIAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWRINITIAL(source) = var_count
         dataEntry%standard_name ="total_WR_initial"
         dataEntry%long_name = "total WR initial"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWaterInitial") then 
      if(LVT_MOC_TOTALSURPLUSWATERINITIAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALSURPLUSWATERINITIAL(source) = var_count
         dataEntry%standard_name ="total_surplus_water_initial"
         dataEntry%long_name = "total surplus water initial"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWaterDeficitInitial") then 
      if(LVT_MOC_TOTALWATERDEFICITINITIAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWATERDEFICITINITIAL(source) = var_count
         dataEntry%standard_name ="total_water_deficit_initial"
         dataEntry%long_name = "total water deficit initial"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalAETVeg") then 
      if(LVT_MOC_TOTALAETVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALAETVEG(source) = var_count
         dataEntry%standard_name ="total_AET_veg"
         dataEntry%long_name = "total AET veg"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWRVeg") then 
      if(LVT_MOC_TOTALWRVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWRVEG(source) = var_count
         dataEntry%standard_name ="total_WR_veg"
         dataEntry%long_name = "total WR veg"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWaterVeg") then 
      if(LVT_MOC_TOTALSURPLUSWATERVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALSURPLUSWATERVEG(source) = var_count
         dataEntry%standard_name ="total_surplus_water_veg"
         dataEntry%long_name = "total surplus water veg"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWaterDeficitVeg") then 
      if(LVT_MOC_TOTALWATERDEFICITVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWATERDEFICITVEG(source) = var_count
         dataEntry%standard_name ="total_water_deficit_veg"
         dataEntry%long_name = "total water deficit veg"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalAETFlower") then 
      if(LVT_MOC_TOTALAETFLOWER(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALAETFLOWER(source) = var_count
         dataEntry%standard_name ="total_AET_flower"
         dataEntry%long_name = "total AET flower"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWRFlower") then 
      if(LVT_MOC_TOTALWRFLOWER(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWRFLOWER(source) = var_count
         dataEntry%standard_name ="total_WR_flower"
         dataEntry%long_name = "total WR flower"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWaterFlower") then 
      if(LVT_MOC_TOTALSURPLUSWATERFLOWER(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALSURPLUSWATERFLOWER(source) = var_count
         dataEntry%standard_name ="total_surplus_water_flower"
         dataEntry%long_name = "total surplus water flower"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWaterDeficitFlower") then 
      if(LVT_MOC_TOTALWATERDEFICITFLOWER(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWATERDEFICITFLOWER(source) = var_count
         dataEntry%standard_name ="total_surplus_water_deficit_flower"
         dataEntry%long_name = "total surplus water deficit flower"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalAETRipe") then 
      if(LVT_MOC_TOTALAETRIPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALAETRIPE(source) = var_count
         dataEntry%standard_name ="total_AET_ripe"
         dataEntry%long_name = "total AET ripe"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWRRipe") then 
      if(LVT_MOC_TOTALWRRIPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWRRIPE(source) = var_count
         dataEntry%standard_name ="total_WR_ripe"
         dataEntry%long_name = "total WR ripe"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalSurplusWaterRipe") then 
      if(LVT_MOC_TOTALSURPLUSWATERRIPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALSURPLUSWATERRIPE(source) = var_count
         dataEntry%standard_name ="total_surplus_water_ripe"
         dataEntry%long_name = "total surplus water ripe"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalWaterDeficitRipe") then 
      if(LVT_MOC_TOTALWATERDEFICITRIPE(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALWATERDEFICITRIPE(source) = var_count
         dataEntry%standard_name ="total_water_deficit_ripe"
         dataEntry%long_name = "total water deficit ripe"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."PermWiltDate") then 
      if(LVT_MOC_PERMWILTDATE(source).eq.LVT_rc%udef) then 
         LVT_MOC_PERMWILTDATE(source) = var_count
         dataEntry%standard_name ="perm_wilt_date"
         dataEntry%long_name = "total wilt date"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Wilting1") then 
      if(LVT_MOC_WILTING1(source).eq.LVT_rc%udef) then 
         LVT_MOC_WILTING1(source) = var_count
         dataEntry%standard_name ="wilting1"
         dataEntry%long_name = "wilting1"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."wilting2") then 
      if(LVT_MOC_WILTING2(source).eq.LVT_rc%udef) then 
         LVT_MOC_wilting2(source) = var_count
         dataEntry%standard_name ="wilting2"
         dataEntry%long_name = "wilting2"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WRSIa") then 
      if(LVT_MOC_WRSIa(source).eq.LVT_rc%udef) then 
         LVT_MOC_WRSIa(source) = var_count
         dataEntry%standard_name ="water_requirements_satisfaction_index_anomaly"
         dataEntry%long_name = "water requirements satisfaction index anomaly"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)   
      endif
   elseif(name.eq."growing_season") then 
      if(LVT_MOC_GROWING_SEASON(source).eq.LVT_rc%udef) then 
         LVT_MOC_GROWING_SEASON(source) = var_count
         dataEntry%standard_name ="growing_season"
         dataEntry%long_name = "growing season"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WHC") then 
      if(LVT_MOC_WHC(source).eq.LVT_rc%udef) then 
         LVT_MOC_WHC(source) = var_count
         dataEntry%standard_name ="WHC"
         dataEntry%long_name = "WHC"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."LGP") then 
      if(LVT_MOC_LGP(source).eq.LVT_rc%udef) then 
         LVT_MOC_LGP(source) = var_count
         dataEntry%standard_name ="LGP"
         dataEntry%long_name = "LGP"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WR_TimeStep") then 
      if(LVT_MOC_WR_TIMESTEP(source).eq.LVT_rc%udef) then 
         LVT_MOC_WR_TIMESTEP(source) = var_count
         dataEntry%standard_name ="WR_TimeStep"
         dataEntry%long_name = "WR_TimeStep"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AET_TimeStep") then 
      if(LVT_MOC_AET_TIMESTEP(source).eq.LVT_rc%udef) then 
         LVT_MOC_AET_TIMESTEP(source) = var_count
         dataEntry%standard_name ="AET_TimeStep"
         dataEntry%long_name = "AET_TimeStep"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WRSI_TimeStep") then 
      if(LVT_MOC_WRSI_TIMESTEP(source).eq.LVT_rc%udef) then 
         LVT_MOC_WRSI_TIMESTEP(source) = var_count
         dataEntry%standard_name ="WRSI_TimeStep"
         dataEntry%long_name = "WRSI_TimeStep"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Snowflag_f") then 
      if(LVT_MOC_SNOWFLAGFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOWFLAGFORC(source) = var_count
         dataEntry%standard_name ="Snowflag"
         dataEntry%long_name = "Snowflag"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Density_f") then 
      if(LVT_MOC_DENSITYFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_DENSITYFORC(source) = var_count
         dataEntry%standard_name ="Atmospheric_Density"
         dataEntry%long_name = "Atmospheric_Density"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m3"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VaporPress_f") then 
      if(LVT_MOC_VAPORPRESSFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_VAPORPRESSFORC(source) = var_count
         dataEntry%standard_name ="Vapor_Pressure"
         dataEntry%long_name = "Vapor_Pressure"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."VaporPressDeficit_f") then 
      if(LVT_MOC_VAPORPRESSDEFICITFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_VAPORPRESSDEFICITFORC(source) = var_count
         dataEntry%standard_name ="Vapor_Pressure_Deficit"
         dataEntry%long_name = "Vapor_Pressure_Deficit"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."AResist") then 
      if(LVT_MOC_ARESIST(source).eq.LVT_rc%udef) then 
         LVT_MOC_ARESIST(source) = var_count
         dataEntry%standard_name ="Aerodynamic_Resistance"
         dataEntry%long_name = "Aerodynamic_Resistance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"s/m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_tsint") then 
      if(LVT_MOC_SACTSINT(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACTSINT(source) = var_count
         dataEntry%standard_name ="sac_soil_temperature_of_intended_layer"
         dataEntry%long_name = "sac soil temperature of intended layer"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_swint") then 
      if(LVT_MOC_SACSWINT(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACSWINT(source) = var_count
         dataEntry%standard_name ="sac_total_volumetric_soil_moisture_content_of_intended_layer"
         dataEntry%long_name = "sac total volumetric soil moisture content of intended layer"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/m3"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_swhint") then 
      if(LVT_MOC_SACSWHINT(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACSWHINT(source) = var_count
         dataEntry%standard_name ="sac_liquid_volumetric_soil_moisture_content_of_intended_layer"
         dataEntry%long_name = "sac liquid volumetric soil moisture content of intended layer"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/m3"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_frost") then 
      if(LVT_MOC_SACFROST(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACFROST(source) = var_count
         dataEntry%standard_name ="sac_frost"
         dataEntry%long_name = "sac frost"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_uztwc") then 
      if(LVT_MOC_SACUZTWC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACUZTWC(source) = var_count
         dataEntry%standard_name ="sac_uztwc"
         dataEntry%long_name = "sac uztwc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_uzfwc") then 
      if(LVT_MOC_SACUZFWC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACUZFWC(source) = var_count
         dataEntry%standard_name ="sac_uzfwc"
         dataEntry%long_name = "sac uzfwc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lztwc") then 
      if(LVT_MOC_SACLZTWC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZTWC(source) = var_count
         dataEntry%standard_name ="sac_lztwc"
         dataEntry%long_name = "sac lztwc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lzfsc") then 
      if(LVT_MOC_SACLZFSC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZFSC(source) = var_count
         dataEntry%standard_name ="sac_lzfsc"
         dataEntry%long_name = "sac lzfsc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lzfpc") then 
      if(LVT_MOC_SACLZFPC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZFPC(source) = var_count
         dataEntry%standard_name ="sac_lzfpc"
         dataEntry%long_name = "sac lzfpc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_adimpc") then 
      if(LVT_MOC_SACADIMPC(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACADIMPC(source) = var_count
         dataEntry%standard_name ="sac_adimpc"
         dataEntry%long_name = "sac adimpc"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm","- "/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_uztwh") then 
      if(LVT_MOC_SACUZTWH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACUZTWH(source) = var_count
         dataEntry%standard_name ="sac_uztwh"
         dataEntry%long_name = "sac uztwh"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_uzfwh") then 
      if(LVT_MOC_SACUZFWH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACUZFWH(source) = var_count
         dataEntry%standard_name ="sac_uzfwh"
         dataEntry%long_name = "sac uzfwh"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lztwh") then 
      if(LVT_MOC_SACLZTWH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZTWH(source) = var_count
         dataEntry%standard_name ="sac_lztwh"
         dataEntry%long_name = "sac lztwh"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lzfsh") then 
      if(LVT_MOC_SACLZFSH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZFSH(source) = var_count
         dataEntry%standard_name ="sac_lzfsh"
         dataEntry%long_name = "sac lzfsh"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."sac_lzfph") then 
      if(LVT_MOC_SACLZFPH(source).eq.LVT_rc%udef) then 
         LVT_MOC_SACLZFPH(source) = var_count
         dataEntry%standard_name ="sac_lzfph"
         dataEntry%long_name = "sac lzfph"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_swe") then 
      if(LVT_MOC_SNOW17SWE(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17SWE(source) = var_count
         dataEntry%standard_name ="snow17_swe"
         dataEntry%long_name = "snow17 swe"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_aeadj") then 
      if(LVT_MOC_SNOW17AEADJ(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17AEADJ(source) = var_count
         dataEntry%standard_name ="snow17_aeadj"
         dataEntry%long_name = "snow17 aeadj"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_neghs") then 
      if(LVT_MOC_SNOW17NEGHS(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17NEGHS(source) = var_count
         dataEntry%standard_name ="snow17_neghs"
         dataEntry%long_name = "snow17 neghs"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_liqw") then 
      if(LVT_MOC_SNOW17LIQW(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17LIQW(source) = var_count
         dataEntry%standard_name ="snow17_liqw"
         dataEntry%long_name = "snow17 liqw"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_accmax") then 
      if(LVT_MOC_SNOW17ACCMAX(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17ACCMAX(source) = var_count
         dataEntry%standard_name ="snow17_accmax"
         dataEntry%long_name = "snow17 accmax"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."snow17_rmlt") then 
      if(LVT_MOC_SNOW17RMLT(source).eq.LVT_rc%udef) then 
         LVT_MOC_SNOW17RMLT(source) = var_count
         dataEntry%standard_name ="snow17_rmlt"
         dataEntry%long_name = "snow17 rmlt"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_satsoil") then 
      if(LVT_MOC_VIC_PET_SATSOIL(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_SATSOIL(source) = var_count
         dataEntry%standard_name ="vic_pet_satsoil"
         dataEntry%long_name = "vic pet satsoil"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_h2osurf") then 
      if(LVT_MOC_VIC_PET_H2OSURF(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_H2OSURF(source) = var_count
         dataEntry%standard_name ="vic_pet_h2osurf"
         dataEntry%long_name = "vic pet h2osurf"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_short") then 
      if(LVT_MOC_VIC_PET_SHORT(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_SHORT(source) = var_count
         dataEntry%standard_name ="vic_pet_short"
         dataEntry%long_name = "vic pet short"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_tall") then 
      if(LVT_MOC_VIC_PET_TALL(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_TALL(source) = var_count
         dataEntry%standard_name ="vic_pet_tall"
         dataEntry%long_name = "vic pet tall"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_natveg") then 
      if(LVT_MOC_VIC_PET_NATVEG(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_NATVEG(source) = var_count
         dataEntry%standard_name ="vic_pet_natveg"
         dataEntry%long_name = "vic pet natveg"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."vic_pet_vegnocr") then 
      if(LVT_MOC_VIC_PET_VEGNOCR(source).eq.LVT_rc%udef) then 
         LVT_MOC_VIC_PET_VEGNOCR(source) = var_count
         dataEntry%standard_name ="vic_pet_vegnocr"
         dataEntry%long_name = "vic pet vegnocr"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2 ","kg/m2s"/)
         dataEntry%valid_min = (/-9999.0,-9999.0/)
         dataEntry%valid_max = (/-9999.0,-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Tsnow") then 
      if(LVT_MOC_LAKE_T_SNOW(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_T_SNOW(source) = var_count
         dataEntry%standard_name ="lake_temperature_at_air_snow_interface"
         dataEntry%long_name = "lake temperature at air snow interface"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Tmnw") then 
      if(LVT_MOC_LAKE_T_MNW(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_T_MNW(source) = var_count
         dataEntry%standard_name ="mean_temperature_of_the_water_column"
         dataEntry%long_name = "mean temperature of the water column"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Twml") then
      if(LVT_MOC_LAKE_T_WML(source).eq.LVT_rc%udef) then  
         LVT_MOC_LAKE_T_WML(source) = var_count
         dataEntry%standard_name ="Lake_temperature_of_the_mixed_layer"
         dataEntry%long_name = "Lake temperature of the mixed layer"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Tbot") then 
      if(LVT_MOC_LAKE_T_BOT(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_T_BOT(source) = var_count
         dataEntry%standard_name ="Lake_temperature_at_the_water_bottom"
         dataEntry%long_name = "Lake temperature at the water bottom"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Tb1") then 
      if(LVT_MOC_LAKE_T_B1(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_T_B1(source) = var_count
         dataEntry%standard_name ="Temperature_at_the_bottom_of_upper_layer_of_sediments"
         dataEntry%long_name = "Temperature at the bottom of upper layer of sediments"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_CT") then 
      if(LVT_MOC_LAKE_C_T(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_C_T(source) = var_count
         dataEntry%standard_name ="Thermocline_shape_factor_of_lake"
         dataEntry%long_name = "Thermocline shape factor of lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Hice") then 
      if(LVT_MOC_LAKE_H_ICE(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_H_ICE(source) = var_count
         dataEntry%standard_name ="Ice_thickness_above_lake"
         dataEntry%long_name = "Ice thickness above lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Hml") then 
      if(LVT_MOC_LAKE_H_ML(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_H_ML(source) = var_count
         dataEntry%standard_name ="Thickness_of_mixed_layer_of_lake"
         dataEntry%long_name = "Thickness of mixed layer of lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Hb1") then 
      if(LVT_MOC_LAKE_H_B1(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_H_B1(source) = var_count
         dataEntry%standard_name ="Thickness_of_upper_layer_of_bottom_sediments"
         dataEntry%long_name = "Thickness of upper layer of bottom sediments"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Walbedo") then 
      if(LVT_MOC_LAKE_ALBEDO_WATER(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_ALBEDO_WATER(source) = var_count
         dataEntry%standard_name ="Water_surface_albedo_over_lake"
         dataEntry%long_name = "Water surface albedo over lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_IceAlbedo") then 
      if(LVT_MOC_LAKE_ALBEDO_ICE(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_ALBEDO_ICE(source) = var_count
         dataEntry%standard_name ="Ice_surface_albedo_over_lake"
         dataEntry%long_name = "Ice surface albedo over lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_SnowAlbedo") then 
      if(LVT_MOC_LAKE_ALBEDO_SNOW(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_ALBEDO_SNOW(source) = var_count
         dataEntry%standard_name ="Snow_surface_albedo_over_lake"
         dataEntry%long_name = "Water surface albedo over lake"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_UFRa") then 
      if(LVT_MOC_LAKE_UFR_A(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_UFR_A(source) = var_count
         dataEntry%standard_name ="Lake_friction_velocity_in_air"
         dataEntry%long_name = "Lake friction velocity in air"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_UFRw") then 
      if(LVT_MOC_LAKE_UFR_W(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_UFR_W(source) = var_count
         dataEntry%standard_name ="Lake_friction_velocity_in_surface_water"
         dataEntry%long_name = "Lake friction velocity in surface water"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_WConv") then 
      if(LVT_MOC_LAKE_WCONV(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_WCONV(source) = var_count
         dataEntry%standard_name ="Lake_convective_velocity_scale"
         dataEntry%long_name = "Lake convective velocity scale"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_IW") then 
      if(LVT_MOC_LAKE_I_W(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_I_W(source) = var_count
         dataEntry%standard_name ="Lake_radiation_flux_at_the_interface"
         dataEntry%long_name = "Lake radiation flux at the interface"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Lake_Qbot") then 
      if(LVT_MOC_LAKE_Q_BOT(source).eq.LVT_rc%udef) then 
         LVT_MOC_LAKE_Q_BOT(source) = var_count
         dataEntry%standard_name ="Lake_heat_flux_across_water_sediment_boundary"
         dataEntry%long_name = "Lake heat flux across water sediment boundary"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RelSMC") then 
      if(LVT_MOC_RELSMC(source).eq.LVT_rc%udef) then 
         LVT_MOC_RELSMC(source) = var_count
         dataEntry%standard_name ="relative_soil_moisture"
         dataEntry%long_name = "relative soil moisture"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/m3","%    ","-    "/)
         dataEntry%valid_min = (/0.0, 0.0, 0.0/)
         dataEntry%valid_max = (/1.0, 100.0, 1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RHMin") then 
      if(LVT_MOC_RHMIN(source).eq.LVT_rc%udef) then 
         LVT_MOC_RHMIN(source) = var_count
         dataEntry%standard_name ="min_relative_humidity"
         dataEntry%long_name = "minimum relative humidity"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-","%"/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TotalPrecip") then 
      if(LVT_MOC_TOTALPRECIP(source).eq.LVT_rc%udef) then 
         LVT_MOC_TOTALPRECIP(source) = var_count
         dataEntry%standard_name ="total_precipitation_amount"
         dataEntry%long_name = "total precipitation amount"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s","kg/m2 "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"UP","DN"/)
      endif
   elseif(name.eq."Emiss_f") then 
      if(LVT_MOC_EMISSFORC(source).eq.LVT_rc%udef) then 
         LVT_MOC_EMISSFORC(source) = var_count
         dataEntry%standard_name ="emissivity"
         dataEntry%long_name = "emissivity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Roughness") then 
      if(LVT_MOC_ROUGHNESS(source).eq.LVT_rc%udef) then 
         LVT_MOC_ROUGHNESS(source) = var_count
         dataEntry%standard_name ="surface_roughness_length"
         dataEntry%long_name = "surface roughness length"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Runoff") then 
      if(LVT_MOC_RUNOFF(source).eq.LVT_rc%udef) then 
         LVT_MOC_RUNOFF(source) = var_count
         dataEntry%standard_name = "total_runoff_amount"
         dataEntry%long_name = "total runoff amount"
         dataEntry%nunits = 3
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s  ","kg/m2   ","mm/month"/)
         dataEntry%valid_min = (/0.0, 0.0, 0.0/)
         dataEntry%valid_max = (/5.0, 43200.0, 43200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."dS") then 
      if(LVT_MOC_dS(source).eq.LVT_rc%udef) then 
         LVT_MOC_dS(source) = var_count
         dataEntry%standard_name = "precip_minus_et_minus_runoff"
         dataEntry%long_name = "precip minus evapotranspiration minus runoff"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s  ","kg/m2   "/)
         dataEntry%valid_min = (/0.0, 0.0/)
         dataEntry%valid_max = (/5.0, 43200.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."EBAL") then 
      if(LVT_MOC_EBAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_EBAL(source) = var_count
         dataEntry%standard_name = "energy_balance"
         dataEntry%long_name = "energy balance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"W/m2"/)
         dataEntry%valid_min = (/-5.0/)
         dataEntry%valid_max = (/5.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."WBAL") then 
      if(LVT_MOC_WBAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_WBAL(source) = var_count
         dataEntry%standard_name = "water_balance"
         dataEntry%long_name = "water balance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."EVAPBAL") then 
      if(LVT_MOC_EVAPBAL(source).eq.LVT_rc%udef) then 
         LVT_MOC_EVAPBAL(source) = var_count
         dataEntry%standard_name = "evaporation_balance"
         dataEntry%long_name = "evaporation balance"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 2
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"IN ", "OUT"/)
      endif
   elseif(name.eq."SWE/P") then 
      if(LVT_MOC_SWEOVERP(source).eq.LVT_rc%udef) then 
         LVT_MOC_SWEOVERP(source) = var_count
         dataEntry%standard_name = "SWE_normalized_by_precipitation"
         dataEntry%long_name = "SWE normalized by precipitation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ET/P") then 
      if(LVT_MOC_ETOVERP(source).eq.LVT_rc%udef) then 
         LVT_MOC_ETOVERP(source) = var_count
         dataEntry%standard_name = "ET_normalized_by_precipitation"
         dataEntry%long_name = "ET normalized by precipitation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Qs/P") then 
      name = "Qs_Over_P"
      if(LVT_MOC_QSOVERP(source).eq.LVT_rc%udef) then 
         LVT_MOC_QSOVERP(source) = var_count
         dataEntry%standard_name = "surface_runoff_normalized_by_precipitation"
         dataEntry%long_name = "surface runoff normalized by precipitation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Qsb/P") then 
      name = "Qsb_Over_P"
      if(LVT_MOC_QSBOVERP(source).eq.LVT_rc%udef) then 
         LVT_MOC_QSBOVERP(source) = var_count
         dataEntry%standard_name = "subsurface_runoff_normalized_by_precipitation"
         dataEntry%long_name = "subsurface runoff normalized by precipitation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ECanop/Qle") then 
      name = "ECanop_Over_Qle"
      if(LVT_MOC_ECANOPOVERQLE(source).eq.LVT_rc%udef) then 
         LVT_MOC_ECANOPOVERQLE(source) = var_count
         dataEntry%short_name = "ECanop_Over_Qle"
         dataEntry%standard_name = "Interception_evaporation_normalized_by_latentheatflux"
         dataEntry%long_name = "Interception evaporation normalized by latentheatflux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."TVeg/Qle") then 
      name = "TVeg_Over_Qle"
      if(LVT_MOC_TVEGOVERQLE(source).eq.LVT_rc%udef) then 
         LVT_MOC_TVEGOVERQLE(source) = var_count
         dataEntry%short_name = "TVeg_Over_Qle"
         dataEntry%standard_name = "Vegetation_transpiration_normalized_by_latentheatflux"
         dataEntry%long_name = "Vegetation transpiration normalized by latentheatflux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ESoil/Qle") then 
      name = "ESoil_Over_Qle"
      if(LVT_MOC_ESOILOVERQLE(source).eq.LVT_rc%udef) then 
         LVT_MOC_ESOILOVERQLE(source) = var_count
         dataEntry%short_name = "ESoil_Over_Qle"
         dataEntry%standard_name = "Bare_soil_evaporation_normalized_by_latentheatflux"
         dataEntry%long_name = "Bare soil evaporation normalized by latentheatflux"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Streamflow") then 
      if(LVT_MOC_STREAMFLOW(source).eq.LVT_rc%udef) then 
         LVT_MOC_STREAMFLOW(source) = var_count
         dataEntry%standard_name = "streamflow"
         dataEntry%long_name = "streamflow"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/s"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RiverStor") then 
      if(LVT_MOC_RIVSTO(source).eq.LVT_rc%udef) then 
         LVT_MOC_RIVSTO(source) = var_count
         dataEntry%standard_name = "river_water_storage"
         dataEntry%long_name = "river water storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RiverDepth") then 
      if(LVT_MOC_RIVDPH(source).eq.LVT_rc%udef) then 
         LVT_MOC_RIVDPH(source) = var_count
         dataEntry%standard_name = "river_depth"
         dataEntry%long_name = "river depth"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RiverVelocity") then 
      if(LVT_MOC_RIVVEL(source).eq.LVT_rc%udef) then 
         LVT_MOC_RIVVEL(source) = var_count
         dataEntry%standard_name ="river_flow_velocity"
         dataEntry%long_name = "river_flow_velocity"
         dataEntry%short_name = "RiverFlowVelocity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodQ") then 
      if(LVT_MOC_FLDOUT(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDOUT(source) = var_count
         dataEntry%standard_name = "floodplain_water_discharge"
         dataEntry%long_name = "floodplain water discharge"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodEvap") then 
      if(LVT_MOC_FLDEVAP(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDEVAP(source) = var_count
         dataEntry%standard_name = "floodplain_evaporation"
         dataEntry%long_name = "floodplain evaporation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodStor") then 
      if(LVT_MOC_FLDSTO(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDSTO(source) = var_count
         dataEntry%standard_name = "floodplain_water_storage"
         dataEntry%long_name = "floodplain water storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodDepth") then 
      if(LVT_MOC_FLDDPH(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDDPH(source) = var_count
         dataEntry%standard_name = "floodplain_depth"
         dataEntry%long_name = "floodplain depth"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodVelocity") then 
      if(LVT_MOC_FLDVEL(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDVEL(source) = var_count
         dataEntry%standard_name = "floodplain_flow_velocity"
         dataEntry%long_name = "floodplain flow velocity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m/s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodedFrac") then 
      if(LVT_MOC_FLDFRC(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDFRC(source) = var_count
         dataEntry%standard_name = "flooded fraction"
         dataEntry%long_name = "flooded fraction"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."FloodedArea") then 
      if(LVT_MOC_FLDARE(source).eq.LVT_rc%udef) then 
         LVT_MOC_FLDARE(source) = var_count
         dataEntry%standard_name = "flooded_area"
         dataEntry%long_name = "flooded area"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m2"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."SurfElev") then 
      if(LVT_MOC_SFCELV(source).eq.LVT_rc%udef) then 
         LVT_MOC_SFCELV(source) = var_count
         dataEntry%standard_name = "surface_water_elevation"
         dataEntry%long_name = "surface water elevation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RunoffStor") then 
      if(LVT_MOC_RNFSTO(source).eq.LVT_rc%udef) then 
         LVT_MOC_RNFSTO(source) = var_count
         dataEntry%standard_name = "runoff_reservoir_storage"
         dataEntry%long_name = "runoff reservoir storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."BaseflowStor") then 
      if(LVT_MOC_BSFSTO(source).eq.LVT_rc%udef) then 
         LVT_MOC_BSFSTO(source) = var_count
         dataEntry%standard_name = "baseflow_reservoir_storage"
         dataEntry%long_name = "baseflow reservoir storage"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/500000.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RTM_emissivity") then 
      if(LVT_MOC_RTM_EMISSIVITY(source).eq.LVT_rc%udef) then 
         LVT_MOC_RTM_EMISSIVITY(source) = var_count
         dataEntry%standard_name = "rtm_emissivity"
         dataEntry%long_name = "rtm emissivity"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/0.7/)
         dataEntry%valid_max = (/1.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RTM_Tb") then 
      if(LVT_MOC_RTM_TB(source).eq.LVT_rc%udef) then 
         LVT_MOC_RTM_TB(source) = var_count
         dataEntry%standard_name = "rtm_brightness_temperature"
         dataEntry%long_name = "rtm brightness temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/213.0/)
         dataEntry%valid_max = (/353.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RTM_SoilMoist") then 
      if(LVT_MOC_RTM_SM(source).eq.LVT_rc%udef) then 
         LVT_MOC_RTM_SM(source) = var_count
         dataEntry%standard_name = "rtm_soil_moisture"
         dataEntry%long_name = "rtm soil moisture"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"m3/m3"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/0.5/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."IrrigatedWater") then 
      if(LVT_MOC_IRRIGATEDWATER(source).eq.LVT_rc%udef) then 
         LVT_MOC_IRRIGATEDWATER(source) = var_count
         dataEntry%standard_name = "irrigated_water_amount"
         dataEntry%long_name = "irrigated water amount"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"kg/m2s"/)
         dataEntry%valid_min = (/0.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WT") then 
      if(LVT_MOC_WT(source).eq.LVT_rc%udef) then 
         LVT_MOC_WT(source) = var_count
         dataEntry%standard_name = "water_in_aquifer_and_saturated_soil" 
         dataEntry%long_name = "water in aquifer and saturated soil"
         dataEntry%nunits = 2
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"mm", "m "/)
         dataEntry%valid_min = (/-9999.0, -9999.0/)
         dataEntry%valid_max = (/-9999.0, -9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."LeafMass") then 
      if(LVT_MOC_LEAFMASS(source).eq.LVT_rc%udef) then 
         LVT_MOC_LEAFMASS(source) = var_count
         dataEntry%standard_name = "leaf_mass"
         dataEntry%long_name = "leaf mass"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."RootMass") then 
      if(LVT_MOC_ROOTMASS(source).eq.LVT_rc%udef) then 
         LVT_MOC_ROOTMASS(source) = var_count
         dataEntry%standard_name = "mass_of_fine_roots"
         dataEntry%long_name = "mass of fine roots"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."StemMass") then 
      if(LVT_MOC_STEMMASS(source).eq.LVT_rc%udef) then 
         LVT_MOC_STEMMASS(source) = var_count
         dataEntry%standard_name = "mass_of_wood_stem"
         dataEntry%long_name = "mass of wood stem"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."WoodMass") then 
      if(LVT_MOC_WOODMASS(source).eq.LVT_rc%udef) then 
         LVT_MOC_WOODMASS(source) = var_count
         dataEntry%standard_name = "mass_of_wood_including_woody_roots"
         dataEntry%long_name = "mass of wood including woody roots"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DeepSoilCarbon") then 
      if(LVT_MOC_CARBON_DEEPSOIL(source).eq.LVT_rc%udef) then 
         LVT_MOC_CARBON_DEEPSOIL(source) = var_count
         dataEntry%standard_name = "stable_carbon_in_deep_soil"
         dataEntry%long_name = "stable carbon in deep soil"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ShallowSoilCarbon") then 
      if(LVT_MOC_CARBON_SHALLOWSOIL(source).eq.LVT_rc%udef) then 
         LVT_MOC_CARBON_SHALLOWSOIL(source) = var_count
         dataEntry%standard_name = "short-lived_carbon_in_shallow_soil"
         dataEntry%long_name = "short-lived carbon in shallow soil"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"g/m2"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DAensspread") then 
      if(LVT_MOC_DA_ENSSPREAD(source).eq.LVT_rc%udef) then 
         LVT_MOC_DA_ENSSPREAD(source) = var_count
         dataEntry%standard_name = "DA ensemble spread"
         dataEntry%long_name = "DA ensemble spread"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DAincr") then 
      if(LVT_MOC_DA_INCR(source).eq.LVT_rc%udef) then 
         LVT_MOC_DA_INCR(source) = var_count
         dataEntry%standard_name = "DA analysis increments"
         dataEntry%long_name = "DA analysis increments"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DAninnov") then 
      if(LVT_MOC_DA_NINNOV(source).eq.LVT_rc%udef) then 
         LVT_MOC_DA_NINNOV(source) = var_count
         dataEntry%standard_name = "DA ninnov"
         dataEntry%long_name = "DA normalized innovation"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."DAobscount") then 
      if(LVT_MOC_DA_OBSCOUNT(source).eq.LVT_rc%udef) then 
         LVT_MOC_DA_OBSCOUNT(source) = var_count
         dataEntry%standard_name = "DA obscount"
         dataEntry%long_name = "DA observation count"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Tair_f_min") then 
      if(LVT_MOC_TAIRFORC_MIN(source).eq.LVT_rc%udef) then 
         LVT_MOC_TAIRFORC_MIN(source) = var_count
         dataEntry%standard_name = "min air temperature"
         dataEntry%long_name = "minimum air temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."Tair_f_max") then 
      if(LVT_MOC_TAIRFORC_MAX(source).eq.LVT_rc%udef) then 
         LVT_MOC_TAIRFORC_MAX(source) = var_count
         dataEntry%standard_name = "max air temperature"
         dataEntry%long_name = "maximum air temperature"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"K"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."ESI") then 
      if(LVT_MOC_ESI(source).eq.LVT_rc%udef) then 
         LVT_MOC_ESI(source) = var_count
         dataEntry%standard_name = "ESI"
         dataEntry%long_name = "ESI"
         dataEntry%nunits = 1
         allocate(dataEntry%unittypes(dataEntry%nunits))
         allocate(dataEntry%valid_min(dataEntry%nunits))
         allocate(dataEntry%valid_max(dataEntry%nunits))
         dataEntry%unittypes = (/"-"/)
         dataEntry%valid_min = (/-9999.0/)
         dataEntry%valid_max = (/-9999.0/)
         dataEntry%ndirs = 1
         allocate(dataEntry%dirtypes(dataEntry%ndirs))
         dataEntry%dirtypes = (/"-"/)
      endif
   elseif(name.eq."none") then 
      dataEntry%standard_name = "none"
      dataEntry%long_name = "none"
      dataEntry%nunits = 1
      allocate(dataEntry%unittypes(dataEntry%nunits))
      allocate(dataEntry%valid_min(dataEntry%nunits))
      allocate(dataEntry%valid_max(dataEntry%nunits))
      dataEntry%unittypes = (/"none"/)
      dataEntry%valid_min = (/0.0/)
      dataEntry%valid_max = (/-9999.0/)
      dataEntry%ndirs = 1
      allocate(dataEntry%dirtypes(dataEntry%ndirs))
      dataEntry%dirtypes = (/"-"/)
   endif
      
 end subroutine updateDataEntryMetaData
  

!BOP
! !ROUTINE: set_ptr_into_list
! \label{set_ptr_into_list}
! 
! !INTERFACE: 
  subroutine set_ptr_into_list(var_count, &
       ds1head, ds1array, &
       ds2Head, ds2Array, &
       ds3Head, ds3Array)
     implicit none
! !ARGUMENTS:
     integer                                 :: var_count(3)
     type(LVT_metadataEntry), pointer        :: ds1head
     type(dep), dimension(var_count(1))      :: ds1array
     type(LVT_metadataEntry), pointer        :: ds2Head
     type(dep), dimension(var_count(2))      :: ds2Array
     type(LVT_metadataEntry), pointer        :: ds3Head
     type(dep), dimension(var_count(3))      :: ds3Array
! 
! !DESCRIPTION: 
!  This routine takes an array of pointers and sets each element of the
!  array to point directly to its corresponding element in the given
!  datastream linked list.  This allows for direct access to the
!  elements in the history output linked list.
!
!   The arguments are: 
!   \begin{description}
!   \item[var_count] number of elements in the given history output linked list
!   \item[ds1head] head of the datastream 1 linked list
!   \item[ds1array] array of pointers to point directly to the elements of \newline
!                the datastream 1 linked list \newline
!   \item[ds2head] head of the datastream 2 linked list
!   \item[ds2array] array of pointers to point directly to the elements of \newline
!                the datastream 2 linked list \newline
!   \item[ds3head] head of the datastream 3 linked list
!   \item[ds3array] array of pointers to point directly to the elements of \newline
!                the datastream 3 linked list \newline
!   \end{description}
!EOP

     type(LVT_metadataEntry), pointer :: ds1dataEntry
     type(LVT_metadataEntry), pointer :: ds2dataEntry
     type(LVT_metadataEntry), pointer :: ds3dataEntry
     integer :: i

     ds1dataEntry    => ds1Head
     ds2dataEntry    => ds2Head
     ds3dataEntry    => ds3Head

     do i = 1, var_count(1)
        ds1array(i)%dataEntryPtr => ds1dataEntry
        ds1dataEntry => ds1dataEntry%next
     enddo

     do i = 1, var_count(2)
        ds2Array(i)%dataEntryPtr =>ds2dataEntry
        ds2dataEntry => ds2dataEntry%next
     enddo

     do i = 1, var_count(3)
        ds3Array(i)%dataEntryPtr =>ds3dataEntry
        ds3dataEntry => ds3dataEntry%next
     enddo

  end subroutine set_ptr_into_list


!BOP
! 
! !ROUTINE: LVT_logSingleDataStreamVar
! \label{LVT_logSingleDataStreamVar}
!
! !INTERFACE:
  subroutine LVT_logSingleDataStreamVar(index, source, value_inp,vlevel,units,&
       ens_index,dir,stdev)
! 
! !USES: 

    implicit none
! !ARGUMENTS: 
    integer                    :: index(3)
    integer                    :: source
    real                       :: value_inp(LVT_rc%lnc, LVT_rc%lnr)
    integer                    :: vlevel
    character(len=*)           :: units
    integer,          optional :: ens_index
    character(len=*), optional :: dir
    real,             optional :: stdev(LVT_rc%lnc,LVT_rc%lnr)
!
! !DESCRIPTION: 
!  This subroutine maps the processed LSM related data values 
!  onto the LVT data structures for future temporal averaging 
!  and comparisons. 
! 
!   The arguments are: 
!   \begin{description}
!   \item[index]
!       index of the variable to be mapped to the LVT data structures
!   \item[source] 
!       index of the data source (1 or 2)
!   \item[value\_inp]
!       array of input data values to be mapped into the LVT data structures
!   \item[vlevel]
!       vertical level at which the input data is mapped
!   \item[units]
!       units of the input data
!   \item[dir]
!       direction type of the input data
!   \item[stdev]
!       standard deviation values of the input data
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[logSingleVar] (\ref{logSingleVar}) \newline
!      Generic routine that maps the input data to LVT data structures. 
!      Application of external masks to the mapped data is also performed
!      in this routine. 
!   \end{description}
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer      :: ensem
    logical      :: stdev_flag
    logical      :: dir_flag
    character*20 :: dir_tmp
    real         :: stdev_tmp(LVT_rc%lnc,LVT_rc%lnr)

    ! If index is less than zero, the variable was not present in the
    ! LIS output or it was not specified for computing within LVT. 

    if(index(source).gt.0) then 

       stdev_flag = .false.
       stdev_tmp  = LVT_rc%udef
       if(present(stdev)) then 
          stdev_flag = .true.
          stdev_tmp = stdev
       endif
    
       if(present(dir)) then 
          dir_flag = .true.
          dir_tmp = dir
       else
          dir_flag = .false.
       endif

       if(present(ens_index)) then 
          ensem = ens_index
       else
          ensem = 1
       endif
       
       if(source.eq.1) then 
          call logSingleVar(LVT_MOC_COUNT(1),&
               LVT_histData%ptr_into_ds1_list,&
               index(source), value_inp, vlevel, ensem, units, dir_flag, dir,&
               stdev_flag,stdev_tmp)
       elseif(source.eq.2) then 
          call logSingleVar(LVT_MOC_COUNT(2),&
               LVT_histData%ptr_into_ds2_list,&
               index(source), value_inp, vlevel, ensem, units, dir_flag, dir,&
               stdev_flag,stdev_tmp)
       elseif(source.eq.3) then 
          call logSingleVar(LVT_MOC_COUNT(3),&
               LVT_histData%ptr_into_ds3_list,&
               index(source), value_inp, vlevel, ensem, units, dir_flag, dir,&
               stdev_flag,stdev_tmp)
       endif
    endif

  end subroutine LVT_logSingleDataStreamVar


!BOP
! 
! !ROUTINE: logSingleVar
! \label{logSingleVar}
!
! !INTERFACE:
  subroutine logSingleVar(var_count,ptr_into_list, index, value_inp,&
       vlevel, ensem,units,dir_flag, dir,stdev_flag,stdev)
! 
! !USES: 
    use LVT_statsDataMod,only : LVT_stats

    implicit none
! !ARGUMENTS: 
    integer                         :: var_count
    type(dep), dimension(var_count) :: ptr_into_list
    integer                         :: index
    real                            :: value_inp(LVT_rc%lnc, LVT_rc%lnr)
    integer                         :: vlevel
    integer                         :: ensem
    character(len=*)                :: units
    logical                         :: dir_flag
    character(len=*)                :: dir
    logical                         :: stdev_flag
    real                            :: stdev(LVT_rc%lnc,LVT_rc%lnr)

!
! !DESCRIPTION: 
!  This subroutine maps the processed observations onto the LVT data
!  structures for future temporal averaging and comparisons. The 
!  data are also filtered using the specified external mask. 
! 
!   The arguments are: 
!   \begin{description}
!   \item[var\_count]
!       total variable count for the data stream type 
!   \item[ptr\_into\_list] 
!       pointer to the datastream linked list
!   \item[index]
!       index of the variable to be mapped to the LVT data structures
!   \item[value\_inp]
!       array of input data values to be mapped into the LVT data structures
!   \item[vlevel]
!       vertical level at which the input data is mapped
!   \item[units]
!       units of the input data
!   \item[dir\_flag]
!       flag indicating if direction types are specified or not
!   \item[dir]
!       direction type of the input data
!   \item[stdev\_flag]
!       flag indicating if standard deviation values are specified or not
!   \item[stdev]
!       standard deviation values of the input data
!   \end{description}
!
!EOP
    integer                     :: k,i,m,c,r,gid
    integer                     :: nsize
    real                        :: mfactor
    logical                     :: unit_status
    logical                     :: dir_status
    logical                     :: log_var

    type(LVT_metadataEntry), pointer :: dataEntry

    dataEntry => ptr_into_list(index)%dataEntryPtr

    if(dataEntry%startNlevs.gt.0) then 

       if(stdev_flag) then 
          dataEntry%stdev_flag = .true. 
    
          if(.not.allocated(dataEntry%stdev)) then 
             if(LVT_rc%computeEnsMetrics.eq.1) then 
                if(LVT_rc%obssource(1).eq."LIS output") then 
                   nsize = LVT_LIS_rc(1)%ntiles
                elseif(LVT_rc%obssource(2).eq."LIS output") then 
                   nsize = LVT_LIS_rc(2)%ntiles
                elseif(LVT_rc%obssource(3).eq."LIS output") then 
                   nsize = LVT_LIS_rc(3)%ntiles
                endif
             else
                nsize = LVT_rc%ngrid
             endif

             allocate(dataEntry%stdev(LVT_rc%ngrid,dataEntry%vlevels))
             allocate(dataEntry%count_stdev(LVT_rc%ngrid,dataEntry%vlevels))
             
             dataEntry%stdev = 0 
             dataEntry%count_stdev = 0 
          endif
       endif

       k = vlevel
       log_var = .true. 
       
       unit_status = .false. 
       do i=1, dataEntry%nunits
          if(trim(units).eq.trim(dataEntry%unittypes(i))) then 
             unit_status = .true. 
          endif
       enddo
       
       if(unit_status) then 
          if(trim(dataEntry%units).eq.trim(units)) then 
             log_var = .true. 
          else
             log_var = .false. 
          endif
       else
          log_var = .false. 
       endif

!first check if the direction types are consistent
       mfactor = 1.0 
       
       if(dir_flag) then 
          dir_status = .false. 
          do i=1,dataEntry%ndirs
             if(trim(dir).eq.dataEntry%dirtypes(i)) then 
                dir_status = .true. 
             endif
          enddo
          
          if(dir_status) then 
             if(trim(dataEntry%dir).eq.trim(dir)) then 
                mfactor = 1.0
             else
                mfactor = -1.0
             endif
          else
             write(LVT_logunit,*) 'The observation direction type ',trim(dir),&
                  'does not map to the model output types'
             write(LVT_logunit,*) 'program stopping...'
             call LVT_endrun
          endif
       endif

       if(dataEntry%startNlevs.gt.0.and.log_var) then           
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_stats%datamask(c,r).eq.1) then
                   if(value_inp(c,r).ne.LVT_rc%udef) then 
                      if(LVT_domain%gindex(c,r).ne.-1) then 
                         gid = LVT_domain%gindex(c,r)
! EMK...Make sure we add accumulations
                         dataEntry%value(gid,:,k) = dataEntry%value(gid,:,k) + &
                              value_inp(c,r)*mfactor
                         dataEntry%count_status(gid,:,k) = 1
! EMK...For accumulations, only update the count at the compute time.
                         if (dataEntry%timeAvgOpt .ne. 3) then 
                            dataEntry%count(gid,:,k) = &
                                 dataEntry%count(gid,:,k) + 1
                            
                         end if
                      endif
                   endif
                   if(stdev_flag) then 
                      if(stdev(c,r).ne.LVT_rc%udef) then 
                         if(LVT_domain%gindex(c,r).ne.-1) then 
                            gid = LVT_domain%gindex(c,r)
                            dataEntry%stdev(gid,k) = dataEntry%stdev(gid,k) + & 
                                 stdev(c,r)
                            dataEntry%count_stdev(gid,k) = dataEntry%count_stdev(gid,k)+1
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
! EMK...For accumulations, only update the count at the compute time.
          if(LVT_rc%computeFlag.and.dataEntry%timeAvgOpt.eq.3) then 
             do r=1,LVT_rc%lnr
                do c=1,LVT_rc%lnc
                   if(LVT_stats%datamask(c,r).eq.1) then
!                      if(value_inp(c,r).ne.LVT_rc%udef) then 
                      if(LVT_domain%gindex(c,r).ne.-1) then 
                         gid = LVT_domain%gindex(c,r)
                         do m=1,LVT_rc%nensem
                            if(dataEntry%count_status(gid,m,k).eq.1) then 
                               if(dataEntry%value(gid,m,k).ne.&
                                    LVT_rc%udef) then
                                  dataEntry%count(gid,m,k) = &
                                       dataEntry%count(gid,m,k) + 1
                               end if
                            endif
                         enddo
                      endif
                   endif
                enddo
             enddo
          endif
       endif
    endif
  end subroutine logSingleVar

!BOP
!
! !ROUTINE: LVT_getDataStream1Ptr
! \label{LVT_getDataStream1Ptr}
! 
! !INTERFACE: 
  subroutine LVT_getDataStream1Ptr(dataEntry)
! !ARGUMENTS:     
    type(LVT_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
! 
!  This routine returns the data object at the head of the
!  datastream 1 linked list
!
!   The arguments are: 
!   \begin{description}
!   \item[dataEntry]
!      data entry object at the head of the datastream 1 linked list
!   \end{description}
!EOP
    dataEntry => LVT_histData%head_ds1_list
  end subroutine LVT_getDataStream1Ptr

!BOP
!
! !ROUTINE: LVT_getDataStream2Ptr
! \label{LVT_getDataStream2Ptr}
! 
! !INTERFACE: 
  subroutine LVT_getDataStream2Ptr(dataEntry)
! !ARGUMENTS:     
    type(LVT_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
!
!  This routine returns the data object at the head of the
!  datastream 2 linked list
!
!   The arguments are: 
!   \begin{description}
!   \item[dataEntry]
!      data entry object at the head of the datastream 1 linked list
!   \end{description}
!EOP
    dataEntry => LVT_histData%head_ds2_list
  end subroutine LVT_getDataStream2Ptr

!BOP
!
! !ROUTINE: LVT_getDataStream3Ptr
! \label{LVT_getDataStream3Ptr}
! 
! !INTERFACE: 
  subroutine LVT_getDataStream3Ptr(dataEntry)
! !ARGUMENTS:     
    type(LVT_metadataEntry), pointer :: dataEntry
! 
! !DESCRIPTION: 
! 
!  This routine returns the data object at the head of the
!  datastream 1 linked list
!
!   The arguments are: 
!   \begin{description}
!   \item[dataEntry]
!      data entry object at the head of the datastream 1 linked list
!   \end{description}
!EOP
    dataEntry => LVT_histData%head_ds3_list
  end subroutine LVT_getDataStream3Ptr
!BOP
! 
! !ROUTINE: LVT_getstatsEntryPtr
! \label{LVT_getstatsEntryPtr}
! 
! !INTERFACE: 
  subroutine LVT_getstatsEntryPtr(statsEntry)
! !ARGUMENTS:     
    type(LVT_statsEntry), pointer :: statsEntry
! !DESCRIPTION: 
!
!  This routine returns the data object at the head of the
!  stats linked list
!
!   The arguments are: 
!   \begin{description}
!   \item[statsEntry]
!      object at the head of the stats linked list
!   \end{description}
!EOP
    statsEntry => LVT_histData%head_stats_list
    
  end subroutine LVT_getstatsEntryPtr

!BOP
! 
! !ROUTINE: LVT_getdataEntryUnits
! \label{LVT_getdataEntryUnits}
! 
! !INTERFACE:   
  function LVT_getdataEntryUnits(index)
! !ARGUMENTS:     
    integer      :: index
    character*20 :: LVT_getdataEntryUnits
!
! !DESCRIPTION: 
! This routine returns the units of the data object at the 
! specified index in the datastream 1 linked list
!
!   The arguments are: 
!   \begin{description}
!   \item[index]
!      index of the data object 
!   \end{description}
!EOP

    if(index.gt.0) then 
       LVT_getdataEntryUnits = &
            LVT_histData%ptr_into_ds1_list(index)%dataEntryPtr%units
    else
       LVT_getdataEntryUnits = ""
    endif

  end function LVT_getdataEntryUnits

!BOP
! !ROUTINE: LVT_checkDataStreamSetup
! \label{LVT_checkDataStreamSetup}
!
! !INTERFACE: 
  subroutine LVT_checkDatastreamSetup
!
! !DESCRIPTION: 
!  This check examines if the model and obs linked
!  lists are in sync. For data intercomparison, they must 
!  be constructed such that all selected variables from the
!  model list must match up with all the sected variables
!  from the obs list. 
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    integer                          :: checkDiff
!-------------------------------------------------------------------
! If the counts for the model and obs lists are different, then 
! the datastream entries are setup incorrectly. 
!-------------------------------------------------------------------
    checkDiff = LVT_MOC_COUNT(1)-LVT_MOC_COUNT(2)

    if(checkDiff.ne.0) then 
       write(LVT_logunit,*) '[ERR] '
       write(LVT_logunit,*) '[ERR] The data stream attributes are specified incorrectly.'
       write(LVT_logunit,*) '[ERR] Please correct the LVT datastream attributes table'
       write(LVT_logunit,*) '[ERR] so that the model and obs entries are both enabled'
       write(LVT_logunit,*) '[ERR] simultaneously for the selected variables.'
       call LVT_endrun()
    endif

  end subroutine LVT_checkDatastreamSetup

end module LVT_histDataMod
