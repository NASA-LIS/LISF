!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define FILENAME "LIS_NUOPC_Gluecode.F90"
#define MODNAME "LIS_NUOPC_Gluecode"
#include "LIS_NUOPC_Macros.h"

#ifndef GSM_EXTLND
#ifdef ESMF_TRACE
#define T_ENTER(region) call ESMF_TraceRegionEnter(region)
#define T_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define T_ENTER(region)
#define T_EXIT(region)
#endif
#else
#define T_ENTER(region)
#define T_EXIT(region)
#endif

!> @file LIS_NUOPC_Gluecode.F90 LIS NUOPC Gluecode interfaces
module LIS_NUOPC_Gluecode
!BOP
!
! !MODULE: LIS_NUOPC_Gluecode
!
! !DESCRIPTION:
!   This module connects NUOPC initialize, advance,
!   and finalize to LIS.
!
! !REVISION HISTORY:
!  13Oct15    Dan Rosen  Initial Specification
!
! !USES:
  use ESMF
  use NUOPC
  use LIS_coreMod, only: &
    LIS_config_init, &
    LIS_core_init, &
    LIS_finalize, &
    LIS_rc, &
    LIS_domain, &
    LIS_initialized, &
    LIS_vm, &
    LIS_masterproc, &
    LIS_localPet, &
    LIS_npes, &
    LIS_ews_ind, &
    LIS_ewe_ind, &
    LIS_nss_ind, &
    LIS_nse_ind, &
    LIS_ews_halo_ind, &
    LIS_ewe_halo_ind, &
    LIS_nss_halo_ind, &
    LIS_nse_halo_ind
  use LIS_domainMod, only: &
    LIS_domain_init         !initialize specified domains
  use LIS_surfaceModelMod, only: &
    LIS_surfaceModel_setexport
  use LIS_metforcingMod, only: &
    LIS_FORC_State
  use LISWRFGridCompMod, only: &
    LISWRF_alloc_states, &
    LISWRF_reset_states, &
    LISWRF_export
  use LIS_timeMgrMod, only: &
    LIS_timemgr_set      ! set the current time
  use LIS_FORC_AttributesMod, only: &
    forc_attrib_type, &     ! forcing arrribute type
    LIS_FORC_Tair, &
    LIS_FORC_Qair, &
    LIS_FORC_SWdown, &
    LIS_FORC_SWdirect, &
    LIS_FORC_SWdiffuse, &
    LIS_FORC_LWdown, &
    LIS_FORC_Wind_E, &
    LIS_FORC_Wind_N, &
    LIS_FORC_Psurf, &
    LIS_FORC_Rainf, &
    LIS_FORC_Snowf, &
    LIS_FORC_CRainf, &
    LIS_FORC_Forc_Hgt, &
    LIS_FORC_Ch, &
    LIS_FORC_Cm, &
    LIS_FORC_Emiss, &
    LIS_FORC_Q2sat, &
    LIS_FORC_Cosz, &
    LIS_FORC_Alb, &
    LIS_FORC_XICE, &
    LIS_FORC_QSFC, &
    LIS_FORC_CHS2, &
    LIS_FORC_CQS2, &
    LIS_FORC_T2, &
    LIS_FORC_Q2, &
    LIS_FORC_TH2, &
    LIS_FORC_TMN, &
    LIS_FORC_LPressure, &
    LIS_FORC_O3, &    !absorber
    LIS_FORC_PET, & ! SY for FEWSNET
    LIS_FORC_RefET, & ! SY for FEWSNET
    LIS_FORC_CMFORC, & ! dmm
    LIS_FORC_CHFORC, & ! dmm for NLDAS-2
    LIS_FORC_CAPE, & ! dmm for NLDAS-2
    LIS_FORC_PARDR, &
    LIS_FORC_PARDF, &
    LIS_FORC_SWNET, &
    LIS_FORC_SNOWFLAG, &
    LIS_FORC_DENSITY, &
    LIS_FORC_VAPORPRESS, &
    LIS_FORC_VAPORPRESSDEFICIT, &
    LIS_FORC_WIND, &
    LIS_FORC_Z0, &
    LIS_FORC_GVF, &
    LIS_FORC_CO2
  use LIS_ESMF_Extensions
  use LIS_NUOPC_DataCopy
  use LIS_NUOPC_Flags

  IMPLICIT NONE

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_NUOPC_Init       ! init method for nuopc cpl mode
  public :: LIS_NUOPC_DataInit   ! Copy data from internal state to export
  public :: LIS_ImportFieldsCopy ! Copy data from import to internal state
  public :: LIS_NUOPC_Run        ! run method for nuopc cpl mode
  public :: LIS_NUOPC_Final      ! finalize method for nuopc cpl mode
  public :: LIS_GridCreate
  public :: LIS_TimestepGet
  public :: LIS_NestCntGet
  public :: LIS_EnsMemberCntGet
  public :: LIS_RunModeGet
  public :: LIS_Unknown
  public :: LIS_Offline
  public :: LIS_Coupled
  public :: LIS_Hybrid
  public :: LIS_FieldDictionaryAdd
  public :: LIS_Log
  public :: LIS_FieldList
  public :: LIS_FieldListLog
  public :: LIS_TestFillImport
  public :: LIS_TestFillExport

!-----------------------------------------------------------------------------
! !LOCAL VARIABLES:
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: LIS_Unknown = -1
  INTEGER, PARAMETER :: LIS_Offline =  0
  INTEGER, PARAMETER :: LIS_Coupled =  1
  INTEGER, PARAMETER :: LIS_Hybrid  =  2

  type LIS_FieldHookup
    real,pointer,dimension(:,:) :: exportArray   => null()
    real,pointer,dimension(:)   :: exportArray_t => null()
  end type

  type LIS_Field
    character(len=64)                  :: stdName        = ""
    character(len=16)                  :: stateName      = ""
    real                               :: ampValue       = 1.d0
    real                               :: meanValue      = 0.d0
    character(len=10)                  :: units          = ""
    character(len=20)                  :: transferOffer  = ""
    logical                            :: adImport       = .FALSE.
    logical                            :: reqImport      = .FALSE.
    logical                            :: realizedImport = .FALSE.
    logical                            :: directConn     = .FALSE.
    logical                            :: sharedMem      = .FALSE.
    logical                            :: adExport       = .FALSE.
    logical                            :: realizedExport = .FALSE.
    character(len=100)                 :: lisForcVarname = ""
    type(LIS_FieldHookup), allocatable :: hookup(:)
  end type

  type(LIS_Field),dimension(83)  :: LIS_FieldList = (/ &
    LIS_Field(stdName='2m_air_temperature', &
      stateName='t2_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='2m_heat_exchange_coefficient', &
      stateName='chs2_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='2m_moisture_exchange_coefficient', &
      stateName='cqs2_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='2m_potential_temperature', &
      stateName='th2_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='2m_specific_humidity', &
      stateName='q2_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='air_temperature', &
      stateName='tair_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='albedo', &
      stateName='albedo_f', &
      ampValue=0.02d0, meanValue=0.14d0, &
      units='1',transferOffer='will provide'), &
    LIS_Field(stdName='albedo_w_snow_effect', &
      stateName='albedo_snwff', &
      ampValue=0.02d0, meanValue=0.14d0, &
      units='1',transferOffer='will provide'), &
    LIS_Field(stdName='atmospheric_density', &
      stateName='density_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-3',transferOffer='will provide'), &
    LIS_Field(stdName='canopy_moisture', &
      stateName='canopmoist', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdName='convective_available_potential_energy', &
      stateName='cape_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='J kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='convective_rainfall_flux', &
      stateName='crainf_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='cosine_solar_zenith_angle', &
      stateName='coszenith_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='downward_heat_flux_in_soil', &
      stateName='qg', &
      ampValue=90.d0, meanValue=-100.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='eastward_wind', &
      stateName='ewind_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='effective_mixing_ratio', &
      stateName='effmixratio', &
      ampValue=0.2d0, meanValue=-0.35d0, &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='emissivity', &
      stateName='emiss_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='forcing_height', &
      stateName='fheight_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdName='green_vegetation_fraction', &
      stateName='greenness_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='1',transferOffer='will provide'), &
    LIS_Field(stdName='heat_exchange_coefficient_in_air', &
      stateName='ch_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='latent_heat_flux_kinematic', &
      stateName='qlekinematic', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='level_pressure', &
      stateName='lpressure_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='Pa',transferOffer='will provide'), &
    LIS_Field(stdName='liquid_fraction_of_soil_moisture_layer_1', &
      stateName='smliqfracl1', &
      ampValue=0.02d0, meanValue=0.48d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='liquid_fraction_of_soil_moisture_layer_2', &
      stateName='smliqfracl2', &
      ampValue=0.02d0, meanValue=0.19d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='liquid_fraction_of_soil_moisture_layer_3', &
      stateName='smliqfracl3', &
      ampValue=1.d0, meanValue=0.17d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='liquid_fraction_of_soil_moisture_layer_4', &
      stateName='smliqfracl4', &
      ampValue=0.02d0, meanValue=0.22d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='liquid_water_content_of_surface_snow', &
      stateName='swe', &
      ampValue=0.02d0, meanValue=0.d0, &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdName='mixing_ratio', &
      stateName='mixratio_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='momentum_exchange_coefficient_in_air', &
      stateName='cm_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='northward_wind', &
      stateName='nwind_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='ozone_concentration', &
      stateName='o3_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='porosity', &
      stateName='porosity', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='potential_evaporation', &
      stateName='pet_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='rainfall_flux', &
      stateName='rainf_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='reference_et', &
      stateName='refet_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdName='relative_soil_moisture', &
      stateName='relsmc', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='root_zone_soil_moisture', &
      stateName='rootmoist', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='seaicemask', &
      stateName='xice_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='snow_depth', &
      stateName='snowdepth', &
      ampValue=0.d0, meanValue=0.d0, &
      units='m ',transferOffer='will provide'), &
    LIS_Field(stdName='snowfall_flux', &
      stateName='snowf_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='snowflag', &
      stateName='snowflag_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='snowmelt', &
      stateName='qsm', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
    LIS_Field(stdName='soil_moisture_content', &
      stateName='soilmoist', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdName='soil_moisture_fraction_layer_1', &
      stateName='smfracl1', &
      ampValue=0.1d0, meanValue=0.20d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='soil_moisture_fraction_layer_2', &
      stateName='smfracl2', &
      ampValue=0.1d0, meanValue=0.19d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='soil_moisture_fraction_layer_3', &
      stateName='smfracl3', &
      ampValue=0.1d0, meanValue=0.17d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='soil_moisture_fraction_layer_4', &
      stateName='smfracl4', &
      ampValue=0.1d0, meanValue=0.22d0, &
      units='m3 m-3',transferOffer='will provide'), &
    LIS_Field(stdName='soil_temperature_layer_1', &
      stateName='soiltempl1', &
      ampValue=5.d0, meanValue=300.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='soil_temperature_layer_2', &
      stateName='soiltempl2', &
      ampValue=5.d0, meanValue=295.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='soil_temperature_layer_3', &
      stateName='soiltempl3', &
      ampValue=5.d0, meanValue=293.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='soil_temperature_layer_4', &
      stateName='soiltempl4', &
      ampValue=5.d0, meanValue=290.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='soil_temperature_lower_boundary', &
      stateName='tmn_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdName='specific_humidity', &
      stateName='qair_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdName='subsurface_runoff_amount', &
      stateName='qsb', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_air_pressure', &
      stateName='psurf_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='Pa',transferOffer='will provide'), &
    LIS_Field(stdName='surface_diffuse_downwelling_shortwave_flux_in_air', &
      stateName='diffusesw_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_direct_downwelling_shortwave_flux_in_air', &
      stateName='directsw_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_downward_par_diffuse', &
      stateName='pardf_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_downward_par_direct', &
      stateName='pardr_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_downwelling_longwave_flux_in_air', &
      stateName='lwdown_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_downwelling_shortwave_flux_in_air', &
      stateName='swdown_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_net_downward_shortwave_flux', &
      stateName='swnet_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_roughness_length', &
      stateName='roughness_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdName='surface_runoff_amount', &
      stateName='qs', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_snow_area_fraction', &
      stateName='snowcover', &
      ampValue=0.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='surface_specific_humidity', &
      stateName='qsfc_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg kg-1',transferOffer='will provide'), &
#ifdef GSM_EXTLND
    LIS_Field(stdName='surface_temperature_land', &
      stateName='avgsurft', &
      ampValue=10.d0, meanValue=295.d0, &
      units='K',transferOffer='will provide'), &
#else
    LIS_Field(stdName='surface_temperature', &
      stateName='avgsurft', &
      ampValue=10.d0, meanValue=295.d0, &
      units='K',transferOffer='will provide'), &
#endif
    LIS_Field(stdName='surface_upward_latent_heat_flux', &
      stateName='qle', &
      ampValue=500.d0, meanValue=-100.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='surface_upward_sensible_heat_flux', &
      stateName='qh', &
      ampValue=500.d0, meanValue=450.d0, &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdName='vapor_pressure', &
      stateName='vaporpress_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='vapor_pressure_deficit', &
      stateName='vaporpressdef_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='vegetationtype', &
      stateName='vegtype', &
      ampValue=1.d0, meanValue=0.d0, &
      units='-',transferOffer='will provide'), &
    LIS_Field(stdName='wind_speed', &
      stateName='wind_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdName='surface_water_depth', &
      stateName='sfcheadrt_f', &
      ampValue=1.d0, meanValue=0.d0, &
      units='mm',transferOffer='will provide'), &
   LIS_Field(stdName='time_step_infiltration_excess', &
      stateName='infxsrt', &
      ampValue=1.d0, meanValue=0.d0, &
      units='mm',transferOffer='will provide'), &
   LIS_Field(stdName='soil_column_drainage', &
      stateName='soldrain', &
      ampValue=1.d0, meanValue=0.d0, &
      units='mm',transferOffer='will provide'), &
   LIS_Field(stdName='final_potential_evaporation', &
      stateName='etp', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
   LIS_Field(stdName='accum_plant_transpiration', &
      stateName='ett', &
      ampValue=1.d0, meanValue=0.d0, &
      units='W m-2',transferOffer='will provide'), &
   LIS_Field(stdName='total_water_flux_layer_1', &
      stateName='wtrflx1', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
   LIS_Field(stdName='total_water_flux_layer_2', &
      stateName='wtrflx2', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
   LIS_Field(stdName='total_water_flux_layer_3', &
      stateName='wtrflx3', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
   LIS_Field(stdName='total_water_flux_layer_4', &
      stateName='wtrflx4', &
      ampValue=1.d0, meanValue=0.d0, &
      units='kg m-2 s-1',transferOffer='will provide'), &
   LIS_Field(stdName='ground_water_storage', &
      stateName='wa', &
      ampValue=1.d0, meanValue=0.d0, &
      units='mm',transferOffer='will provide')/)

!EOP

contains

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Init: Allocate memory and initialize domain and start values
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_NUOPC_Init"

!BOP
! !ROUTINE: LIS_NUOPC_Init
!
! !INTERFACE:
  subroutine LIS_NUOPC_Init(vm,configFile,rc)
! !ARGUMENTS:
    type(ESMF_VM),intent(in)                :: vm
    integer,intent(out)                     :: rc
    character(len=*),intent(in),optional    :: configFile
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS during the
!  initialization of the LIS-NUOPC system.
! \begin{description}
!  \item[LIS\_config\_init] (\ref{LIS_config_init}) \\
!    initialize the LIS configuration
!  \item[LIS\_domain\_init] (\ref{LIS_domain_init}) \\
!    initialize the LIS domains
!  \item[LIS\_core\_init] (\ref{LIS_core_init}) \\
!    no description
! \end{description}
!
!EOP
!
! ! LOCAL VARIABLES

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(.not.LIS_initialized) then
       if (present(configFile)) LIS_rc%lis_config_file = trim(configFile)
       call LIS_config_init(vm=vm)
       call lisinit(trim(LIS_rc%runmode)//char(0))

       call LISWRF_alloc_states
       call LISWRF_reset_states

       call LIS_HookupInit(rc)
       if(ESMF_STDERRORCHECK(rc)) return

       LIS_initialized = .true.
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_DataInit: Copy Data from internal state to nuopc state
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_NUOPC_DataInit"

!BOP
! !ROUTINE: LIS_NUOPC_DataInit
!
! !INTERFACE:
  subroutine LIS_NUOPC_DataInit(nest,exportState,rc)
! !ARGUMENTS:
    integer,intent(in)                     :: nest
    type(ESMF_State),intent(inout)         :: exportState
    integer,intent(out)                    :: rc

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    ! See LIS_lsmcpl_pluginMod for subroutines registered to:
    ! "retrospective", "WRF coupling", "NUOPC coupling"
    call LIS_surfaceModel_setexport(nest)

    call LIS_ExportFieldsCopy(nest,exportState,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Run: Advance LIS Model
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_NUOPC_Run"

!BOP
! !ROUTINE: LIS_NUOPC_Run
!
! !INTERFACE:
  subroutine LIS_NUOPC_Run(impStateList,expStateList,clock, &
  misg_import,rc)
! !ARGUMENTS:
    type(ESMF_State),intent(inout)         :: impStateList(:)
    type(ESMF_State),intent(inout)         :: expStateList(:)
    type(ESMF_Clock),intent(in)            :: clock
    type(missingval_flag),intent(in)       :: misg_import
    integer,intent(out)                    :: rc
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS during a coupled
!  NUOPC run. The routine translates the NUOPC import fields to NUOPC forcing
!  variables, and then calls the LIS LSM physics routines. Finally the
!  export state from LIS to NUOPC is generated.
!
!EOP
!
! !LOCAL VARIABLES:
    type(ESMF_TimeInterval)        :: timeStep
    type(ESMF_StateItem_Flag)      :: itemType
    type(ESMF_Grid)                :: grid
    type(ESMF_Field)               :: exportField
    real(ESMF_KIND_R4),pointer     :: exportFarray(:,:)
    integer                        :: elb(2), eub(2)
    integer                        :: i, j
    integer                        :: timeStepSecs
    integer                        :: sIndex

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    T_ENTER("datacopy")
    do sIndex=1, size(impStateList)
      call LIS_ImportFieldsCopy(sIndex,impStateList(sIndex),misg_import,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return
    enddo
    T_EXIT("datacopy")

    call lisstep(trim(LIS_rc%runmode)//char(0))

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return
    call ESMF_TimeIntervalGet(timeStep, s=timeStepSecs, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    do sIndex=1, size(expStateList)
      ! Write LIS output data to export state
      ! See LIS_lsmcpl_pluginMod for subroutines registered to:
      ! "retrospective", "WRF coupling", "NUOPC coupling"
      call LIS_surfaceModel_setexport(sIndex)
      call LIS_ExportFieldsCopy(sIndex,expStateList(sIndex),rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return

      ! convert runoff fields qs qsb from (kg m-2 s-1) to (kg m-2)
      call ESMF_StateGet(expStateList(sIndex), &
        itemName="qs", itemType=itemType, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (itemType == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(expStateList(sIndex), &
          itemName="qs", field=exportField, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_FieldGet(exportField, farrayPtr=exportFarray, &
          exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        do j=elb(2), eub(2)
        do i=elb(1), eub(1)
          if (exportFarray(i,j) /= MISSINGVALUE) then
            exportFarray(i,j) = exportFarray(i,j) * timeStepSecs
          endif
        enddo
        enddo
      endif

      call ESMF_StateGet(expStateList(sIndex), &
        itemName="qsb", itemType=itemType, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      if (itemType == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(expStateList(sIndex), &
          itemName="qsb", field=exportField, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_FieldGet(exportField, farrayPtr=exportFarray, &
          exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        do j=elb(2), eub(2)
        do i=elb(1), eub(1)
          if (exportFarray(i,j) /= MISSINGVALUE) then
            exportFarray(i,j) = exportFarray(i,j) * timeStepSecs
          endif
        enddo
        enddo
      endif

    enddo ! end nest loop

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Final: Cleanup allocated memory and set endtime.
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_NUOPC_Final"

!BOP
! !ROUTINE: LIS_NUOPC_Final
!
! !INTERFACE:
  subroutine LIS_NUOPC_Final(nest,clock,rc)
! !ARGUMENTS:
    integer,intent(in)                      :: nest
    type(ESMF_Clock), intent(in)            :: clock
    integer,intent(out)                     :: rc
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS to finalize
!  a coupled LIS-NUOPC run. The routine writes LIS output and restart files
!  valid at the end of the run.
!
!EOP
!
! !LOCAL VARIABLES:
    integer                    :: stat
    type(ESMF_Time)            :: currTime
    integer                    :: yy, mm, dd, h, m, s

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    ! use incoming clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    ! Confirm if the timemgr should receive current time or stop time
    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    call LIS_timemgr_set(LIS_rc, yy, mm, dd, h, m, s, 0, 0.0)

    LIS_rc%endtime = 1

    call lisfinalize(trim(LIS_rc%runmode)//char(0))

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_HookupInit: Initialize the field list hookups
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_HookupInit"

  subroutine LIS_HookupInit(rc)
    ! ARGUMENTS
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    integer :: stat
    integer :: fIndex
    integer :: nIndex
    character(len=10)  :: nStr

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if (LIS_rc%nnest > 999999999) then
      nStr = '999999999+'
    else
      write (nStr,"(I0)") LIS_rc%nnest
    endif

    if (.NOT.allocated(LISWRF_export)) then
      call ESMF_LogSetError(ESMF_RC_OBJ_INIT, &
        msg="LISWRF_export is not initialized.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (size(LISWRF_export) < LIS_rc%nnest) then
      call ESMF_LogSetError(ESMF_RC_OBJ_INIT, &
        msg="LISWRF_export is not initialized for each nest.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    do fIndex=1, size(LIS_FieldList)
      allocate(LIS_FieldList(fIndex)%hookup(LIS_rc%nnest), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of field hookup memory failed.", &
        line=__LINE__,file=__FILE__)) &
      return
    enddo

    do nIndex=1, LIS_rc%nnest
    do fIndex=1, size(LIS_FieldList)
      select case (trim(LIS_FieldList(fIndex)%stdName))
        case ('2m_air_temperature')              ! (01)
          if (allocated(LIS_FORC_T2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_T2%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_T2%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('2m_heat_exchange_coefficient')              ! (02)
          if (allocated(LIS_FORC_CHS2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_CHS2%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_CHS2%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%chs2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%chs2_t
        case ('2m_moisture_exchange_coefficient')              ! (03)
          if (allocated(LIS_FORC_CQS2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_CQS2%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_CQS2%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%cqs2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%cqs2_t
        case ('2m_potential_temperature')              ! (04)
          if (allocated(LIS_FORC_TH2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_TH2%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_TH2%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('2m_specific_humidity')              ! (05)
          if (allocated(LIS_FORC_Q2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Q2%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Q2%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('air_temperature')              ! (06)
          if (allocated(LIS_FORC_Tair%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Tair%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Tair%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('albedo')              ! (07)
          if (allocated(LIS_FORC_Alb%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Alb%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Alb%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%albedo
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%albedo_t
        case ('albedo_w_snow_effect')              ! (08)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%albedo
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%albedo_t
        case ('atmospheric_density')              ! (09)
          if (allocated(LIS_FORC_DENSITY%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_DENSITY%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_DENSITY%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('canopy_moisture')              ! (10)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%cmc
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%cmc_t
        case ('convective_available_potential_energy')              ! (11)
          if (allocated(LIS_FORC_CAPE%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_CAPE%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_CAPE%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('convective_rainfall_flux')              ! (12)
          if (allocated(LIS_FORC_CRainf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_CRainf%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_CRainf%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('cosine_solar_zenith_angle')              ! (13)
          if (allocated(LIS_FORC_Cosz%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Cosz%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Cosz%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('downward_heat_flux_in_soil')              ! (14)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qg
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qg_t
        case ('eastward_wind')              ! (15)
          if (allocated(LIS_FORC_Wind_E%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Wind_E%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Wind_E%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('effective_mixing_ratio')              ! (16)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%q1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%q1_t
        case ('emissivity')              ! (17)
          if (allocated(LIS_FORC_Emiss%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Emiss%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Emiss%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%emiss
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%emiss_t
        case ('forcing_height')              ! (18)
          if (allocated(LIS_FORC_Forc_Hgt%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Forc_Hgt%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Forc_Hgt%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('green_vegetation_fraction')              ! (19)
          if (allocated(LIS_FORC_GVF%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_GVF%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_GVF%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('heat_exchange_coefficient_in_air')              ! (20)
          if (allocated(LIS_FORC_Ch%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Ch%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Ch%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('latent_heat_flux_kinematic')              ! (21)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%eta_kinematic
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%eta_kinematic_t
        case ('level_pressure')              ! (22)
          if (allocated(LIS_FORC_lpressure%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_lpressure%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_lpressure%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('liquid_fraction_of_soil_moisture_layer_1')              ! (23)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o1_t
        case ('liquid_fraction_of_soil_moisture_layer_2')              ! (24)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o2_t
        case ('liquid_fraction_of_soil_moisture_layer_3')              ! (25)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o3_t
        case ('liquid_fraction_of_soil_moisture_layer_4')              ! (26)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o4_t
        case ('liquid_water_content_of_surface_snow')              ! (27)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snow
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snow_t
        case ('mixing_ratio')              ! (28)
          if (allocated(LIS_FORC_Q2sat%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Q2sat%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Q2sat%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('momentum_exchange_coefficient_in_air')              ! (29)
          if (allocated(LIS_FORC_Cm%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Cm%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Cm%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('northward_wind')              ! (30)
          if (allocated(LIS_FORC_Wind_N%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Wind_N%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Wind_N%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('ozone_concentration')              ! (31)
          if (allocated(LIS_FORC_o3%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_o3%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_o3%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('porosity')              ! (32)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.FALSE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%lispor
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%lispor_t
        case ('potential_evaporation')              ! (33)
          if (allocated(LIS_FORC_PET%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_PET%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_PET%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('rainfall_flux')              ! (34)
          if (allocated(LIS_FORC_Rainf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Rainf%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Rainf%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('reference_et')              ! (35)
          if (allocated(LIS_FORC_RefET%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_RefET%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_RefET%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('relative_soil_moisture')              ! (36)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%relsmc1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%relsmc1_t
        case ('root_zone_soil_moisture')              ! (37)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%rootmoist
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%rootmoist_t
        case ('seaicemask')              ! (38)
          if (allocated(LIS_FORC_XICE%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_XICE%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_XICE%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%xice
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%xice_t
        case ('snow_depth')              ! (39)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snowh
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snowh_t
        case ('snowfall_flux')              ! (40)
          if (allocated(LIS_FORC_Snowf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Snowf%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Snowf%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('snowflag')              ! (41)
          if (allocated(LIS_FORC_SNOWFLAG%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SNOWFLAG%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_SNOWFLAG%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('snowmelt')              ! (42)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qsm
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qsm_t
        case ('soil_moisture_content')              ! (43)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%soilm
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%soilm_t
        case ('soil_moisture_fraction_layer_1')              ! (44)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc1_t
        case ('soil_moisture_fraction_layer_2')              ! (45)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc2_t
        case ('soil_moisture_fraction_layer_3')              ! (46)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc3_t
        case ('soil_moisture_fraction_layer_4')              ! (47)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc4_t
        case ('soil_temperature_layer_1')              ! (48)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc1_t
        case ('soil_temperature_layer_2')              ! (49)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc2_t
        case ('soil_temperature_layer_3')              ! (50)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc3_t
        case ('soil_temperature_layer_4')              ! (51)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc4_t
        case ('soil_temperature_lower_boundary')              ! (52)
          if (allocated(LIS_FORC_TMN%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_TMN%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_TMN%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('specific_humidity')              ! (53)
          if (allocated(LIS_FORC_Qair%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Qair%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Qair%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('subsurface_runoff_amount')              ! (54)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qsb
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qsb_t
        case ('surface_air_pressure')              ! (55)
          if (allocated(LIS_FORC_Psurf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Psurf%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Psurf%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_diffuse_downwelling_shortwave_flux_in_air')              ! (56)
          if (allocated(LIS_FORC_SWdiffuse%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SWdiffuse%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_SWdiffuse%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_direct_downwelling_shortwave_flux_in_air')              ! (57)
          if (allocated(LIS_FORC_SWdirect%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SWdirect%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_SWdirect%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_downward_par_diffuse')              ! (58)
          if (allocated(LIS_FORC_Pardf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Pardf%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Pardf%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_downward_par_direct')              ! (59)
          if (allocated(LIS_FORC_Pardr%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Pardr%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Pardr%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_downwelling_longwave_flux_in_air')              ! (60)
          if (allocated(LIS_FORC_LWdown%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_LWdown%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_LWdown%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_downwelling_shortwave_flux_in_air')              ! (61)
          if (allocated(LIS_FORC_SWdown%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SWdown%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_SWdown%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_net_downward_shortwave_flux')              ! (62)
          if (allocated(LIS_FORC_SWnet%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SWnet%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_SWnet%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('surface_roughness_length')              ! (63)
          if (allocated(LIS_FORC_Z0%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Z0%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_Z0%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%znt
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%znt_t
        case ('surface_runoff_amount')              ! (64)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qs
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qs_t
        case ('surface_snow_area_fraction')              ! (65)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snocvr
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snocvr_t
        case ('surface_specific_humidity')              ! (66)
          if (allocated(LIS_FORC_QSFC%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_QSFC%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_QSFC%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
#ifdef GSM_EXTLND
        case ('surface_temperature_land')              ! (67)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%avgsurft
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%avgsurft_t
#else
        case ('surface_temperature')              ! (67)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%avgsurft
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%avgsurft_t
#endif
        case ('surface_upward_latent_heat_flux')              ! (68)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qle
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qle_t
        case ('surface_upward_sensible_heat_flux')              ! (69)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qh
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qh_t
        case ('vapor_pressure')              ! (70)
          if (allocated(LIS_FORC_VAPORPRESS%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_VAPORPRESS%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_VAPORPRESS%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('vapor_pressure_deficit')              ! (71)
          if (allocated(LIS_FORC_VAPORPRESSDEFICIT%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_VAPORPRESSDEFICIT%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_VAPORPRESSDEFICIT%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('vegetationtype')              ! (72)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%xland
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%xland_t
        case ('wind_speed')              ! (73)
          if (allocated(LIS_FORC_WIND%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_WIND%varname(1)
          LIS_FieldList(fIndex)%reqImport=(LIS_FORC_WIND%selectOpt == 1)
          LIS_FieldList(fIndex)%adImport=.TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
#ifdef WRF_HYDRO
        case ('surface_water_depth')              ! (74)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%#NOTAVAILABLE#_t
        case ('time_step_infiltration_excess')              ! (75)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%infxsrt
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%infxsrt_t
        case ('soil_column_drainage')              ! (76)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%soldrain
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%soldrain_t
#endif
#ifdef PARFLOW
        case ('total_water_flux_layer_1')              ! (?)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%wtrflx1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%wtrflx1_t
        case ('total_water_flux_layer_2')              ! (?)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%wtrflx2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%wtrflx2_t
        case ('total_water_flux_layer_3')              ! (?)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%wtrflx3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%wtrflx3_t
        case ('total_water_flux_layer_4')              ! (?)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%reqImport=(#NOTAVAILABLE#%selectOpt == 1)
!          LIS_FieldList(fIndex)%adImport=.FALSE.
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%wtrflx4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%wtrflx4_t
        case ('ground_water_storage')              ! (?)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
          LIS_FieldList(fIndex)%adImport = .TRUE.
          LIS_FieldList(fIndex)%directConn = .TRUE.
          LIS_FieldList(fIndex)%sharedMem = .TRUE.
!          LIS_FieldList(fIndex)%adExport=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%%#NOTAVAILABLE#
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%%#NOTAVAILABLE#_t
#endif
#ifdef GSM_EXTLND
        case ('accum_plant_transpiration')              ! (77)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%adImport=(#NOTAVAILABLE#%selectOpt == 1)
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%ett
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%ett_t
        case ('final_potential_evaporation')              ! (78)
!          if (allocated(#NOTAVAILABLE#%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=#NOTAVAILABLE#%varname(1)
!          LIS_FieldList(fIndex)%adImport=(#NOTAVAILABLE#%selectOpt == 1)
          LIS_FieldList(fIndex)%adExport=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%etp
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%etp_t
#endif
        case default
          call ESMF_LogWrite("LIS: Field hookup information missing. " //&
            "Skipping hookup: "//trim(LIS_FieldList(fIndex)%stdName), &
            ESMF_LOGMSG_WARNING)
      end select
    enddo
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy NUOPC Import Fields to LIS Forcing State
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_ImportFieldsCopy"

  subroutine LIS_ImportFieldsCopy(nest,importState,missing,label,rc)
    ! ARGUMENTS
    integer,intent(in)                :: nest
    type(ESMF_State),intent(inout)    :: importState
    type(missingval_flag),intent(in)  :: missing
    character(*),intent(in),optional  :: label
    integer,intent(out)               :: rc
    ! LOCAL VARIABLES
    character(len=64)          :: l_label
    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_StateItem_Flag)  :: lisItemType
    type(ESMF_Field)           :: importField
    type(ESMF_Field)           :: lisImportField
    integer                    :: fIndex

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_ImportFieldsCopy:'
    endif

    do fIndex = 1,size(LIS_FieldList)
      ! Skip over field if it is not realized in any import state
      if (LIS_FieldList(fIndex)%realizedImport) then
        ! Check itemType to see if field exists in import state
        call ESMF_StateGet(importState, &
          itemName=trim(LIS_FieldList(fIndex)%stateName), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (itemType == ESMF_STATEITEM_FIELD) then

          call ESMF_StateGet(importState, &
            itemName=trim(LIS_FieldList(fIndex)%stateName), &
            field=importField,rc=rc)

          if (LIS_FieldList(fIndex)%directConn) then
            if (LIS_rc%lsm.eq."Noah.3.3") then
              call LIS_CopyToNoah_3_3(field=importField, &
                stdName=LIS_FieldList(fIndex)%stdName, &
                nest=nest,missing=missing,rc=rc)
              if(ESMF_STDERRORCHECK(rc)) return
            else if (LIS_rc%lsm.eq."NoahMP.3.6") then
              call LIS_CopyToNoahMP_3_6(field=importField, &
                stdName=LIS_FieldList(fIndex)%stdName, &
                nest=nest,missing=missing,rc=rc)
              if(ESMF_STDERRORCHECK(rc)) return
            else if (LIS_rc%lsm.eq."Noah-MP.4.0.1") then
              call LIS_CopyToNoahMP_4_0_1(field=importField, &
                stdName=LIS_FieldList(fIndex)%stdName, &
                nest=nest,missing=missing,rc=rc)
              if(ESMF_STDERRORCHECK(rc)) return
            else
              call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                msg="Direct coupling is not implemented for "//trim(LIS_rc%lsm), &
                line=__LINE__, file=FILENAME, rcToReturn=rc)
              return
            endif
          else
            call LIS_ForcFieldGet(LIS_FieldList(fIndex)%lisForcVarname, &
              nest=nest,itemType=lisItemType,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return

            if (lisItemType == ESMF_STATEITEM_FIELD ) then
              call LIS_ForcFieldGet(LIS_FieldList(fIndex)%lisForcVarname, &
                nest=nest,field=lisImportField,rc=rc)
              if(ESMF_STDERRORCHECK(rc)) return
              call LIS_CopyToLIS(field=importField, &
                fieldLIS=lisImportField, &
                nest=nest,rc=rc)
              if(ESMF_STDERRORCHECK(rc)) return
            else ! not present in LIS_Forc
              call ESMF_LogWrite( trim(l_label)// &
                " field is not present in LIS_Forc state="// &
                trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
            endif ! check LIS_Forc
          endif ! direct connection vs LIS_FORC
        else ! not present in NOUPC import
          call ESMF_LogWrite( trim(l_label)// &
           " field is not present in NUOPC import state="// &
            trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
          cycle
        endif ! Not present in import state
      endif ! realizedImport
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy LIS Export Fields to NUOPC Export State
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_ExportFieldsCopy"

  subroutine LIS_ExportFieldsCopy(nest,exportState,label,rc)
! !ARGUMENTS:
    integer,intent(in)                :: nest
    type(ESMF_State),intent(inout)    :: exportState
    character(*),intent(in),optional  :: label
    integer,intent(out)               :: rc
!EOP
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    type(ESMF_Field)           :: exportField
    integer                    :: fIndex
    type(ESMF_StateItem_Flag)  :: itemType

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_ExportFieldsCopy:'
    endif

    do fIndex = 1,size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%realizedExport) then
        ! Check itemType to see if field exists in import state
        call ESMF_StateGet(exportState, &
          itemName=trim(LIS_FieldList(fIndex)%stateName), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        if (itemType == ESMF_STATEITEM_FIELD) then
          if (.NOT. allocated(LIS_FieldList(fIndex)%hookup)) then
            call ESMF_LogWrite( trim(l_label)// &
              " export field hookups have not been allocated="// &
              trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
            cycle
          endif
          if (.NOT.associated(LIS_FieldList(fIndex)%hookup(nest)%exportArray_t)) then
            call ESMF_LogWrite( trim(l_label)// &
              " export array_t is missing="// &
              trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
            cycle
          endif
          call ESMF_StateGet(exportState, &
            itemName=trim(LIS_FieldList(fIndex)%stateName), &
            field=exportField,rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return
          if (LIS_FieldList(fIndex)%stateName .eq. "smliqfracl1") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smliqfracl2") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smliqfracl3") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smliqfracl4") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smfracl1") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smfracl2") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smfracl3") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "smfracl4") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=1.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "infxsrt") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=0.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          elseif (LIS_FieldList(fIndex)%stateName .eq. "soldrain") then
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,fillVal=0.0,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          else
            call LIS_CopyFromLIS( &
              farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
              field=exportField,nest=nest,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          endif
        else
          call ESMF_LogWrite( trim(l_label)// &
            " field is not present in NUOPC export state="// &
            trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
          cycle
        endif
      endif
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Fill LIS Import Fields
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_TestFillImport"

  subroutine LIS_TestFillImport(nest,importState,step,label,rc)
! !ARGUMENTS:
    integer,intent(in)                :: nest
    type(ESMF_State),intent(inout)    :: importState
    integer,intent(in)                :: step
    character(*),intent(in),optional  :: label
    integer,intent(out)               :: rc
!EOP
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    type(ESMF_Field)           :: importField
    integer                    :: fIndex
    type(ESMF_StateItem_Flag)  :: itemType

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_TestFillImport:'
    endif

    do fIndex = 1,size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%realizedImport) then
        ! Check itemType to see if field exists in import state
        call ESMF_StateGet(importState, &
          itemName=trim(LIS_FieldList(fIndex)%stateName), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        if (itemType == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(importState, &
            itemName=trim(LIS_FieldList(fIndex)%stateName), &
            field=importField,rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return
          call LIS_ESMF_FieldFill(importField, dataFillScheme="sincos", &
            step=step, amplitude=LIS_FieldList(fIndex)%ampValue, &
            meanValue=LIS_FieldList(fIndex)%meanValue, rc=rc)
        else
          call ESMF_LogWrite( trim(l_label)// &
            " field is not present in NUOPC import state="// &
            trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
          cycle
        endif
      endif
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Fill LIS Export Fields
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_TestFillExport"

  subroutine LIS_TestFillExport(nest,exportState,step,label,rc)
! !ARGUMENTS:
    integer,intent(in)                :: nest
    type(ESMF_State),intent(inout)    :: exportState
    integer,intent(in)                :: step
    character(*),intent(in),optional  :: label
    integer,intent(out)               :: rc
!EOP
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    type(ESMF_Field)           :: exportField
    integer                    :: fIndex
    type(ESMF_StateItem_Flag)  :: itemType

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_TestFillExport:'
    endif

    do fIndex = 1,size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%realizedExport) then
        ! Check itemType to see if field exists in export state
        call ESMF_StateGet(exportState, &
          itemName=trim(LIS_FieldList(fIndex)%stateName), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        if (itemType == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(exportState, &
            itemName=trim(LIS_FieldList(fIndex)%stateName), &
            field=exportField,rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return
          call LIS_ESMF_FieldFill(exportField, dataFillScheme="sincos", &
            step=step, amplitude=LIS_FieldList(fIndex)%ampValue, &
            meanValue=LIS_FieldList(fIndex)%meanValue, rc=rc)
        else
          call ESMF_LogWrite( trim(l_label)// &
            " field is not present in NUOPC export state="// &
            trim(LIS_FieldList(fIndex)%stateName),ESMF_LOGMSG_WARNING)
          cycle
        endif
      endif
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_GridCreate"

  function LIS_GridCreate(nest,rc)
    !RETURN VALUE
    type(ESMF_Grid) :: LIS_GridCreate
    !ARGUEMENTS
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc

    ! LOCAL VARIABLES
    integer                                    :: stat
    integer,allocatable                        :: deBlockList(:,:,:)
    integer,allocatable                        :: petMap(:)
    type(ESMF_DELayout)                        :: deLayout
    type(ESMF_DistGrid)                        :: distGrid
    character(len=10)   :: did
    integer             :: petID, deID
    integer             :: start(2)
    integer             :: id, col, row
    real                :: lat_cor, lon_cor
    real                :: lat_cen, lon_cen
    integer             :: lbnd(2),ubnd(2)
    real(ESMF_KIND_COORD), pointer    :: coordXcenter(:,:)
    real(ESMF_KIND_COORD), pointer    :: coordYcenter(:,:)
    real(ESMF_KIND_COORD), pointer    :: coordXcorner(:,:)
    real(ESMF_KIND_COORD), pointer    :: coordYcorner(:,:)
    integer(ESMF_KIND_I4), pointer :: gridmask(:,:)

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    ! Setup list of 2-D decompositions
    ! One block per PET
    allocate(deBlockList(2,2,LIS_npes),petMap(LIS_npes),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of deBlockList and petMap memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return
    do petID=0, LIS_npes-1
      deID = petID + 1
      deBlockList(1,1,deID) = LIS_ews_ind(nest,deID) ! East-West Start
      deBlockList(1,2,deID) = LIS_ewe_ind(nest,deID) ! East-West End
      deBlockList(2,1,deID) = LIS_nss_ind(nest,deID) ! North-South Start
      deBlockList(2,2,deID) = LIS_nse_ind(nest,deID) ! North-South End
      petMap(deID)=petID
    enddo

    deLayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    distGrid = ESMF_DistGridCreate( &
      minIndex=(/1,1/), maxIndex=(/LIS_rc%gnc(nest),LIS_rc%gnr(nest)/), &
      indexflag = ESMF_INDEX_DELOCAL, &
      deBlockList=deBlockList, &
      delayout=deLayout, &
      rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call LIS_DecompGet(distgrid,istart=start(1),jstart=start(2),rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    deallocate(deBlockList,petMap,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of deBlockList and petMap memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return

    write(did,"(I0)") nest
    LIS_GridCreate = ESMF_GridCreate(name='LIS_Grid_'//trim(did), &
      distgrid=distgrid, &
      coordSys=ESMF_COORDSYS_SPH_DEG, &
      coordTypeKind=ESMF_TYPEKIND_COORD, &
      gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
      rc = rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Add Center Latitude & Longitude Coordinates
    call ESMF_GridAddCoord(LIS_GridCreate, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Get pointer to Center Latitude array
    call ESMF_GridGetCoord(LIS_GridCreate, coordDim=1, &
      staggerloc=ESMF_STAGGERLOC_CENTER, localDE=0, &
      computationalLBound=lbnd, computationalUBound=ubnd, &
      farrayPtr=coordXcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    ! Get pointer to Center Longitude array
    call ESMF_GridGetCoord(LIS_GridCreate, coordDim=2, &
      staggerloc=ESMF_STAGGERLOC_CENTER, localDE=0, &
      farrayPtr=coordYcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! error checking
    ! lnc_red and lnr_red are "reduced" to not include the halo region
    if ((ubnd(1)-lbnd(1)+1 /= LIS_rc%lnc_red(nest)) .OR. &
    (ubnd(2)-lbnd(2)+1 /= LIS_rc%lnr_red(nest))) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="Size of local coordinate arrays not equal to size of local DE", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call LIS_ESMF_NetcdfReadIXJX("lon",trim(LIS_rc%paramfile(nest)), &
      start,coordXcenter,rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call LIS_ESMF_NetcdfReadIXJX("lat",trim(LIS_rc%paramfile(nest)), &
      start,coordYcenter,rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Add Grid Mask
    call ESMF_GridAddItem(LIS_GridCreate, itemFlag=ESMF_GRIDITEM_MASK, &
      itemTypeKind=ESMF_TYPEKIND_I4, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    ! Get pointer to Grid Mask array
    call ESMF_GridGetItem(LIS_GridCreate, itemflag=ESMF_GRIDITEM_MASK, &
      localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=gridmask, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    call LIS_ESMF_NetcdfReadIXJX("LANDMASK",trim(LIS_rc%paramfile(nest)), &
      start,gridmask,rc)
    if (ESMF_STDERRORCHECK(rc)) return

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_DecompGet"

  subroutine LIS_DecompGet(distgrid,istart,iend,jstart,jend,rc)
    type(ESMF_DistGrid),intent(in)          :: distGrid
    integer,intent(out),optional            :: istart,iend,jstart,jend
    integer,intent(out)                     :: rc

    ! LOCAL VARIABLES
    integer                :: stat
    integer,allocatable    :: indexCountPDe(:,:)
    integer,allocatable    :: iIndexList(:), jIndexList(:)

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    !! Get the grid distribution for this pet
    allocate(indexCountPDe(2, 0:(LIS_npes - 1)),stat=stat) ! (dimCount, deCount)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of indexCountPDE memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return
    call ESMF_DistGridGet(distgrid, indexCountPDe=indexCountPDe, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(iIndexList(indexCountPDe(1, LIS_localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of iIndexList memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return
    call ESMF_DistGridGet(distgrid, localDe=0, dim=1, &
      indexList=iIndexList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(jIndexList(indexCountPDe(2, LIS_localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of jIndexList memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return
    call ESMF_DistGridGet(distgrid, localDe=0, dim=2, &
      indexList=jIndexList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (present(istart)) istart = minVal(iIndexList)
    if (present(iend)) iend = maxVal(iIndexList)
    if (present(jstart)) jstart = minVal(jIndexList)
    if (present(jend)) jend = maxVal(jIndexList)

    deallocate(indexCountPDe,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of indexCountPDe memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return
    deallocate(iIndexList,jIndexList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of iIndexList and jIndexList memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Retrieve ESMF Field from LIS_FORC_State
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_ForcFieldGet"

  !BOP
! !ROUTINE: LIS_ForcFieldGet(varname,nest,field,itemType,rc)
!
! !INTERFACE:
  subroutine LIS_ForcFieldGet(varname,nest,field,itemType,rc)
! !ARGUMENTS:
    character(len=*),intent(in)                     :: varname
    integer,intent(in)                              :: nest
    type(ESMF_Field),intent(out),optional           :: field
    type(ESMF_StateItem_Flag),intent(out),optional  :: itemType
    integer,intent(out)                             :: rc
! !LOCAL VARIABLES:
    type(ESMF_StateItem_Flag)                       :: l_itemType
! !DESCRIPTION:
!   Return ESMF_Field for forcing variable
!
!EOP

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if (nest > size(LIS_FORC_State)) then
      if (present(itemType)) itemType = ESMF_STATEITEM_NOTFOUND
      return
    endif

    call ESMF_StateGet(LIS_FORC_State(nest), &
      itemName=trim(varname), &
      itemType=l_itemType, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (present(itemType)) itemType = l_itemType

    if (l_itemType == ESMF_STATEITEM_FIELD) then
      if (present(field)) then
        call ESMF_StateGet(LIS_FORC_State(nest), &
          itemName=trim(varname), &
          field=field,rc=rc)
        if(ESMF_STDERRORCHECK(rc)) return
      endif
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Retrieve timestep of nest from LIS_rc
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_TimestepGet"

!BOP
! !FUNCTION: LIS_TimestepGet(rc)
! !INTERFACE:
  function LIS_TimestepGet(rc)
! !RETURN VALUE:
    real :: LIS_TimestepGet
! !ARGUMENTS:
    integer,intent(out)                    :: rc
! !DESCRIPTION:
!   Return timestep for LIS
!
!EOP
!
! !LOCAL VARIABLES:

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    LIS_TimestepGet = LIS_rc%ts

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------
  ! Retrieve nest count from LIS_rc
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_NestCntGet"

!BOP
! !FUNCTION: LIS_NestCntGet(rc)
! !INTERFACE:
  function LIS_NestCntGet(rc)
! !RETURN VALUE:
    integer               :: LIS_NestCntGet
! !ARGUMENTS:
    integer,intent(out)           :: rc
! !DESCRIPTION:
!   Return timestep for LIS nest
!
!EOP
!
! !LOCAL VARIABLES:

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    LIS_NestCntGet = LIS_rc%nnest

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------
  ! Retrieve ensemble information from LIS_rc
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_EnsMemberCntGet"

!BOP
! !FUNCTION: LIS_EnsMemberCntGet(nest,ensemble,ensMemberCnt,rc)
! !INTERFACE:
  subroutine LIS_EnsMemberCntGet(nest,ensMemberCnt,rc)
! !ARGUMENTS:
    integer,intent(in)            :: nest
    integer,intent(out)           :: ensMemberCnt
    integer,intent(out)           :: rc
! !DESCRIPTION:
!   Return  ensemble information
!
!EOP
!
! !LOCAL VARIABLES:

    rc = ESMF_SUCCESS

    ensMemberCnt = LIS_rc%nensem(nest)

  end subroutine

  !-----------------------------------------------------------------------------
  ! Retrieve NUOPC coupling mode
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_RunModeGet"

!BOP
! !FUNCTION: LIS_RunModeGet(fieldList,rc)
! !INTERFACE:
  function LIS_RunModeGet(fieldList,stateList,rc)
! !RETURN VALUE:
    integer :: LIS_RunModeGet
! !ARGUMENTS:
    type(LIS_Field),intent(in)      :: fieldList(:)
    type(ESMF_State),intent(in)     :: stateList(:)
    integer,intent(out)             :: rc
! !DESCRIPTION:
!   Return NUOPC coupling for LIS nest
!
!EOP
!
! !LOCAL VARIABLES:
    integer                   :: sIndex
    integer                   :: fIndex
    logical                   :: allConnected
    logical                   :: isConnected
    integer                   :: connectedCount
    type(ESMF_StateItem_Flag) :: itemType

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    LIS_RunModeGet = LIS_Unknown
    connectedCount = 0
    allConnected = .true.

    do fIndex=1, size(fieldList)
      if(fieldList(fIndex)%reqImport) then
        ! Check to see if field is connected in each state
        do sIndex=1, size(stateList)
          call ESMF_StateGet(stateList(sIndex), &
            itemName=trim(fieldList(fIndex)%stateName), &
            itemType=itemType, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return

          if (itemType == ESMF_STATEITEM_FIELD) then
            isConnected = NUOPC_IsConnected(stateList(sIndex), &
              fieldName=trim(fieldList(fIndex)%stateName), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            isConnected = .false.
          endif

          if (isConnected) then
            connectedCount = connectedCount + 1
          else
            allConnected = .false.
          endif
        enddo
      endif
    enddo

    if ( allConnected ) then
      LIS_RunModeGet = LIS_Coupled
      LIS_rc%metforc_blend_alg = "overlay"
      LIS_rc%metforc(:) = "none"
    elseif( connectedCount > 0 ) then
      LIS_RunModeGet = LIS_Hybrid
    elseif( connectedCount == 0 ) then
      LIS_RunModeGet = LIS_Offline
    else
      LIS_RunModeGet = LIS_Unknown
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------
  ! Dictionary Utility
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_FieldDictionaryAdd"

  subroutine LIS_FieldDictionaryAdd(rc)
    ! ARGUMENTS
    integer,intent(out)                     :: rc
    ! LOCAL VARIABLES
    integer                    :: fIndex
    logical                    :: isPresent

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    do fIndex=1,size(LIS_FieldList)
      isPresent = NUOPC_FieldDictionaryHasEntry( &
        trim(LIS_FieldList(fIndex)%stdName), &
        rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (.not.isPresent) then
        call NUOPC_FieldDictionaryAddEntry( &
          trim(LIS_FieldList(fIndex)%stdName), &
          trim(LIS_FieldList(fIndex)%units), &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Log Utilities
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_Log"

  subroutine LIS_Log(label,rc)
! !INTERFACE:
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    integer                    :: nIndex
    integer                    :: mIndex
    character(len=4)           :: nestStr
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_Log'
    endif

    ! LIS State
    write (logMsg,"(A,A,L1)") trim(l_label), &
      " Initialized: ",LIS_initialized
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,L1)") trim(l_label), &
      " Master process: ",LIS_masterproc
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Local PET: ",LIS_localPet
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of PETs: ",LIS_npes
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    ! Run Config
    write (logMsg,"(A,A,A)") trim(l_label), &
      " Running mode: ",LIS_rc%runmode
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,A)") trim(l_label), &
      " Start code: ",LIS_rc%startcode
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of nests: ",LIS_rc%nnest
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Number of processors (East-West,North-South): (", &
      LIS_rc%npesx,",",LIS_rc%npesy,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    ! Run time
    write (logMsg,"(A,A,6(I0,A))") trim(l_label), &
      " Current time LIS (yr-mo-da-hr:mn:ss): (", &
      LIS_rc%yr,"-",LIS_rc%mo,"-",LIS_rc%da,"_", &
      LIS_rc%hr,":",LIS_rc%mn,":",LIS_rc%ss,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of forcing variables: ", &
      LIS_rc%nf
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of meteorological forcing datasets: ", &
      LIS_rc%nmetforc
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

    if (allocated(LIS_rc%met_nf)) then
      do mIndex=1,size(LIS_rc%met_nf)
        write (logMsg,"(A,A,I0,A,I0)") trim(l_label), &
          " Number of forcing variables in met forcing dataset (", &
          mIndex,"): ",LIS_rc%met_nf(mIndex)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      enddo
    endif

    if (allocated(LIS_rc%metforc)) then
      do mIndex=1,size(LIS_rc%metforc)
        write (logMsg,"(A,A,I0,A,A)") trim(l_label), &
          " Met forcing (", &
          mIndex,"): ",trim(LIS_rc%metforc(mIndex))
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      enddo
      write (logMsg,"(A,A,A)") trim(l_label), &
        " Met forcing blending: ", trim(LIS_rc%metforc_blend_alg)
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    endif

    do nIndex=1,LIS_rc%nnest
      if (nIndex > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nIndex
      endif
      call LIS_LogNest(nIndex,trim(l_label)//" Nest: "//trim(nestStr),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_LogNest"

  subroutine LIS_LogNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    integer                    :: lbounds(3), ubounds(3)
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogNest_'//trim(nestStr)
    endif

    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Grid dimension (East-West,North-South)=(", &
      LIS_rc%gnc(nest),",",LIS_rc%gnr(nest),")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Size of gridspace: ",LIS_rc%glbngrid(nest)
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.5)") trim(l_label), &
      " Timestep for nest: ",LIS_rc%nts(nest)
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_LogDomainForNest"

  subroutine LIS_LogDomainForNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    integer                    :: lbounds(3), ubounds(3)
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogDomainForNest_'//trim(nestStr)
    endif

    ! LIS Domain
    lbounds(1) = lbound(LIS_domain(nest)%lat,1)
    ubounds(1) = ubound(LIS_domain(nest)%lat,1)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Latitude array values=(", &
      LIS_domain(nest)%lat(lbounds(1)), ":", &
      LIS_domain(nest)%lat(ubounds(1)), ")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    lbounds(1) = lbound(LIS_domain(nest)%lon,1)
    ubounds(1) = ubound(LIS_domain(nest)%lon,1)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Longitude array values=(", &
      LIS_domain(nest)%lon(lbounds(1)), ":", &
      LIS_domain(nest)%lon(ubounds(1)), ")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Grid index (min,max)=(", &
      MINVAL(LIS_domain(nest)%gindex),":",MAXVAL(LIS_domain(nest)%gindex),")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Latitude (min,max)=(", &
      LIS_domain(nest)%minLat,",",LIS_domain(nest)%maxLat,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Longitude (min,max)=(", &
      LIS_domain(nest)%minLon,",",LIS_domain(nest)%maxLon,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.3,A))") trim(l_label), &
      " Domain Resolution (DX,DY)=(", &
      LIS_domain(nest)%dx,",",LIS_domain(nest)%dy,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.3)") trim(l_label), &
      " Lat increment for lat/lon: ",LIS_domain(nest)%lisproj%dlat
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.3)") trim(l_label), &
      " Lon increment for lat/lon: ",LIS_domain(nest)%lisproj%dlon
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_LogGridDescForNest"

  subroutine LIS_LogGridDescForNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogGridDescForNest'//trim(nestStr)
    endif

    ! Grid Description
    select case ( INT(LIS_rc%gridDesc(nest,1)) )
      case (0) ! 0 = indicateds lat/lon projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(0) lat/lon", ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of columns in the domain: ",LIS_rc%gridDesc(nest,2)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
           " Number of rows in the domain: ",LIS_rc%gridDesc(nest,3)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner latitude: ",LIS_rc%gridDesc(nest,4)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner longitude: ",LIS_rc%gridDesc(nest,5)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner latitude: ",LIS_rc%gridDesc(nest,7)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner longitude: ",LIS_rc%gridDesc(nest,8)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " East-West spatial resolution: ",LIS_rc%gridDesc(nest,9)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " North-South spatial resolution: ",LIS_rc%gridDesc(nest,10)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      case (1) ! 1 = indicates mercator projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(1) mercator", ESMF_LOGMSG_INFO)
      case (3) ! 3 = indicates lambert projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(3) lambert", ESMF_LOGMSG_INFO)
      case (4) ! 4 = indicates gaussian projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(4) gaussian", ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of columns in the domain: ",LIS_rc%gridDesc(nest,2)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of rows in the domain: ",LIS_rc%gridDesc(nest,3)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner latitude: ",LIS_rc%gridDesc(nest,4)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner longitude: ",LIS_rc%gridDesc(nest,5)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner latitude: ",LIS_rc%gridDesc(nest,7)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner longitude: ",LIS_rc%gridDesc(nest,8)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " East-West spatial resolution: ",LIS_rc%gridDesc(nest,9)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " Number of latitude circles: ",LIS_rc%gridDesc(nest,10)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      case (5) ! 5 = indicates polar stereographic projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(5) polar stereographic", ESMF_LOGMSG_INFO)
      case (7) ! 7 = indicates UTM projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(7) UTM", ESMF_LOGMSG_INFO)
      case default
        write (logMsg,"(A,A,I0,A)") trim(l_label), &
          " Grid projection=(",LIS_rc%gridDesc(nest,1),") UNKNOWN"
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    end select

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_FieldListLog"

  subroutine LIS_FieldListLog(label)
    ! ARGUMENTS
    character(len=*),intent(in) :: label
    ! LOCAL VARIABLES
    integer                     :: fIndex

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    do fIndex=1, size(LIS_FieldList)
      call LIS_FieldLog(LIS_FieldList(fIndex),label=label)
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "LIS_FieldLog"

  subroutine LIS_FieldLog(field,label)
    ! ARGUMENTS
    type(LIS_Field),intent(in)  :: field
    character(len=*),intent(in) :: label
    ! LOCAL VARIABLES
    integer                     :: nIndex
    character(len=10)           :: nStr
    character(ESMF_MAXSTR)      :: logMsg
    logical                     :: iAssoc
    logical                     :: eAssoc
    logical                     :: tAssoc
    type(ESMF_StateItem_Flag)   :: itemType
    type(ESMF_Field)            :: importField
    integer                     :: rc

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    write (logMsg, "(A,A,(A,A))") trim(label)//': ', &
      trim(field%stateName), &
      ' Transfer offer: ',trim(field%transferOffer)
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,A,(A,L1))") trim(label)//': ', &
      trim(field%stateName), &
      ' Import advertised: ',field%adImport
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,A,(A,L1))") trim(label)//': ', &
      trim(field%stateName), &
      ' Import realized: ',field%realizedImport
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,A,(A,L1))") trim(label)//': ', &
      trim(field%stateName), &
      ' Export advertised: ',field%adExport
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,A,(A,L1))") trim(label)//': ', &
      trim(field%stateName), &
      ' Export realized: ',field%realizedExport
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

    if (allocated(field%hookup)) then
      do nIndex=1, size(field%hookup)
        if (nIndex > 999999999) then
          write (nStr,"(A)") '999999999+'
        else
          write (nStr,"(I0)") nIndex
        endif
        call LIS_ForcFieldGet(field%lisForcVarname, &
          nest=nIndex,itemType=itemType,rc=rc)
        if(ESMF_STDERRORCHECK(rc)) return
        if (itemType == ESMF_STATEITEM_FIELD) then
          iAssoc = .TRUE.
        else
          iAssoc = .FALSE.
        endif
        eAssoc = associated(field%hookup(nIndex)%exportArray)
        tAssoc = associated(field%hookup(nIndex)%exportArray_t)
        write (logMsg,"(A,A,(A,A),(A,L1))") trim(label)//': ', &
          trim(field%stateName), &
          " Nest=",trim(nStr), &
          " Import hookup=",iAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,(A,A),(A,L1))") trim(label)//': ', &
          trim(field%stateName), &
          " Nest=",trim(nStr), &
          " Export hookup_2D=",eAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,(A,A),(A,L1))") trim(label)//': ', &
          trim(field%stateName), &
          " Nest=",trim(nStr), &
          " Export hookup_1D=",tAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

        if (iAssoc) then
          call LIS_ForcFieldGet(field%lisForcVarname, &
            nest=nIndex,field=importField,rc=rc)
          if (.NOT.(rc.eq.ESMF_SUCCESS)) then
            call ESMF_LogWrite(trim(label)// &
            ': LIS_ForcFieldGet failed for '//trim(field%lisForcVarname), &
            ESMF_LOGMSG_ERROR)
          else
            call LIS_ESMF_LogFieldLclVal(importField, &
              label=trim(label)//": Nest="//trim(nStr)//" import")
          endif
        endif
        if (eAssoc) then
          call LIS_ESMF_LogFarrayLclVal(field%hookup(nIndex)%exportArray, &
            label=trim(label)//": Nest="//trim(nStr)//" export2D")
        endif
        if (tAssoc) then
          call LIS_ESMF_LogFarrayLclVal(field%hookup(nIndex)%exportArray_t, &
            label=trim(label)//": Nest="//trim(nStr)//" export1D")
        endif
      enddo
    else
      write (logMsg,"(A,A)") trim(label)//": ", &
        field%stateName//" does not have field hookups."
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO, &
          line=__LINE__, file=FILENAME)
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

end module
