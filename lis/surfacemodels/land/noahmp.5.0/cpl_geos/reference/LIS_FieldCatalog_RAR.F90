module LIS_FieldCatalog

   !DESCRIPTION:
   ! Self-contained copy of the stdName/stateName/units metadata and role
   ! classification (adImport/adExport/directConn) that LIS_Fields
   ! needs from LIS_NUOPC_Gluecode's LIS_FieldList/LIS_HookupInit -- kept here
   ! instead so the Plug does not have to compile/link LISF's NUOPC cap
   ! (LIS_NUOPC_Cap/LIS_NUOPC_Gluecode/LIS_NUOPC_DataCopy/LIS_ESMF_Extensions/
   ! LIS_NUOPC_Flags, and transitively LISWRFGridCompMod and NUOPC itself) or
   ! LIS_FORC_AttributesMod just to get this static data.
   !
   ! Roles mirror LIS_HookupInit's select case in
   ! @LISF/lis/runmodes/nuopc_cpl_mode/LIS_NUOPC_Gluecode.F90 for the default
   ! build only (no WRF_HYDRO/PARFLOW/GSM_EXTLND); keep in sync if that
   ! select case changes upstream.
   !
   ! Reference only: this file is intentionally kept out of the cpl_geos/
   ! build directory (moved to this reference/ subdirectory) so it is never
   ! compiled. Implementation status: every adExport=.true. field below is
   ! now wired into noahMP50_setgeosexport.F90/NoahMP50_geosCplMod.F90,
   ! except Landcover, which is exported as a raw NoahMP vegetype passthrough
   ! pending a NoahMP->CLSM land-cover category remap table.

   implicit none
   private

   public :: LIS_Field
   public :: LIS_FieldList

   type :: LIS_Field
      character(len=:), allocatable :: stdName
      character(len=:), allocatable :: stateName
      character(len=:), allocatable :: units
      logical :: adImport = .false.
      logical :: adExport = .false.
      logical :: directConn = .false.
   end type LIS_Field

   type :: LIS_FieldList
      type(LIS_Field), allocatable :: fields(:)
   end type LIS_FieldList

   interface LIS_FieldList
      module procedure new_LIS_FieldList
   end interface LIS_FieldList

contains

   ! Constructs the catalog. Roles (adImport/adExport/directConn) are baked
   ! in directly, mirroring LIS_HookupInit's select case (default build only,
   ! no WRF_HYDRO/PARFLOW/GSM_EXTLND -- those fields keep all-.FALSE. defaults).
   function new_LIS_FieldList() result(this)
      type(LIS_FieldList) :: this

!!=========================================================================================
!! LSM exports to routing 
!!=========================================================================================
LIS_Field( &
                stdName='surface_runoff_amount', &    ! Noah-MP: RunoffSurface
                stateName='Qs', & 					! NoahMP50_struc(n)%noahmp50(t)%sfcrunoff
                units='mm/soil_timestep', &			! CLSM exports runoff_total_flux (RUNOFF )
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='subsurface_runoff_amount', &         ! Noah-MP: RunoffSubsurface
                stateName='Qsb', &						   ! NoahMP50_struc(n)%noahmp50(t)%udrrunoff
                units='mm/soil_timestep', &					! CLSM exports runoff_total_flux (RUNOFF )
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='totalrunoff', &         ! RunoffSubsurface + RunoffSurface
                stateName='Qsum', &						   
                units='mm/soil_timestep', &					! CLSM exports runoff_total_flux (RUNOFF )
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!!=========================================================================================
!! exports to Moist / Turbulence / Radiation
!!=========================================================================================
!! --- Surface Skin Temperature  ---
LIS_Field( &
                stdName='surface_temperature', & ! Noah-MP: TemperatureSfc
                stateName='AvgSurfT', &			! NoahMP50_struc(n)%noahmp50(t)%tg
                units='K', &					! CLSM exports: TH (turbulence_surface_skin_temperature)
                adImport=.false., &			     ! This is used to calculate THAT and SHAT in SurfaceGridComp
                adExport=.true., &
                directConn=.false.), &


!! --- Surface Specific Humidity  ---
LIS_Field( &
                stdName='mean_surface_specific_humidity', & ! Noah-MP: SpecHumiditySfcMean
                stateName='AvgSurfQ', &			! NoahMP50_struc(n)%noahmp50(t)%qsfcmean
                units='kg/kg', &					! CLSM exports: QH (turbulence_surface_skin_specific_hum)
                adImport=.false., &			! This is used to calculate QHAT in SurfaceGridComp
                adExport=.true., &
                directConn=.false.), &
! This is used to calculate QHAT in SurfaceGridComp

!! --- Snow Mass ---
LIS_Field( &
                stdName='liquid_water_content_of_surface_snow', &     ! Noah-MP: SnowWaterEquiv
                stateName='SWE', &								  ! NoahMP50_struc(n)%noahmp50(t)%sneqv
                units='kg/m2', &							! CLSM exports: snow_mass  (SNOWMASS)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!! --- 2-Meter Air Temperature  ---
LIS_Field( &
                stdName='2m_air_temperature', &  ! Noah-MP: TemperatureAir2m
                stateName='t2_f', &			! NoahMP50_struc(n)%noahmp50(t)% t2mean
                units='K', &					! CLSM MOT2M (temperature 2m wind from MO sfc)
                adImport=.false., &.               
                adExport=.true., &
                directConn=.false.), &

!! --- Sensible Heat Flux  ---
LIS_Field( &
                stdName='surface_upward_sensible_heat_flux', &        ! Noah-MP: HeatSensibleSfc
                stateName='Qh', &								  
                units='W/m2', &								 ! CLSM: SHOUT (upward_sensible_heat_flux)
                adImport=.false., &								
                adExport=.true., &
                directConn=.false.), &

!! --- Surface Roughness, ZPD, Exchange Coefficients (Partitioned) ---
LIS_Field( &
                stdName='surface_roughness_length', &          ! Noah-MP: RoughLenMomSfcToAtm (maps to Z0 and Z0H)
                stateName='roughness_f', &					! NoahMP50_struc(n)%noahmp50(t)%z0 
                units='m', &							!one to many mapping: CLSM Z0 (surface_roughness), Z0H (surface_roughness_for_heat)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='heat_exchange_coefficient_in_air', & ! Noah-MP: ExchCoeffShSfc 
                stateName='ch_f', &						! NoahMP50_struc(n)%noahmp50(t)%ch
                units='m/s', &						! (maps to GEOS heat and momentum exchange coefficients: CT(surface_exchange_coefficient_for_heat), CQ(surface_exchange_coefficient_for_moisture)
                adImport=.false., &					
                adExport=.true., &					
                directConn=.false.), &
LIS_Field( &
                stdName='momentum_exchange_coefficient_in_air', & ! Noah-MP: ExchCoeffMomSfc      
                stateName='cm_f', &						! NoahMP50_struc(n)%noahmp50(t)%cm
                units='m/s', &						
                adImport=.false., &					
                adExport=.true., &					
                directConn=.false.), &

!! --- 2-Meter Specific Humidity ---
LIS_Field( &
                stdName='2m_specific_humidity', &      ! Noah-MP: SpecificHumidity2m (Derived from CLSM Q2M)
                stateName='QairT2m', &				  ! NoahMP50_struc(n)%noahmp50(t)%q2mv & %q2mb
                units='kg/kg', &					!CLSM: Q2M (2-meter_specific_humidity)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &
!! --- Evaporation ---
LIS_Field( &
                stdName='total_evapotranspiration', &  ! Noah-MP: WaterToAtmosTotal
                stateName='Evap', &					! NoahMP50_struc(n)%noahmp50(t)%qfx
                units='kg/m2s', &					! CLSM: EVLAND, 'total_evapotranspiration_land'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!! --- displacement height (does not look like this needs to be passed to other components)---
LIS_Field( &
                stdName='displacement_height', & ! Noah-MP: ZeroPlaneDispSfc 
                stateName='zpd', &				! NoahMP50_struc(n)%noahmp50(t)%zpd
                units='m', &					! This is considered a MAPL_AddInternalSpec in CLSM, D0 (zero_plane_displacement_height)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!! --- surface emissivity ---
LIS_Field( &
                stdName='emissivity', &                         ! Noah-MP: EmissivitySfc
                stateName='Emiss_f', &						! NoahMP50_struc(n)%noahmp50(t)%emiss
                units='-', &						! CLSM: EMIS, 'surface_emissivity'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!!--- surface albedos ---
LIS_Field( &
                stdName='surface_reflectivity_visible_beam', &          ! gx%albedo_dir(1) = nm%albsfcdir(1)   ! AlbedoVisDir
                stateName='AlbedoVisDir', &                             
                units='-', &								  ! CLSM: ALBVR (surface_reflectivity_visible_beam) 
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_reflectivity_visible_diffuse', &       ! gx%albedo_dif(1) = nm%albsfcdif(1)   ! AlbedoVisDif
                stateName='AlbedoVisDif', &                             
                units='-', &								  ! CLSM: ALBVF (surface_reflectivity_visible_diffuse)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_reflectivity_near_infrared_beam', &    ! gx%albedo_dir(2) = nm%albsfcdir(2)   ! AlbedoNIRDir
                stateName='AlbedoNIRDir', &                             
                units='-', &								  ! CLSM: ALBNR, 'surface_reflectivity_near_infrared_beam' 
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_reflectivity_near_infrared_diffuse', & ! gx%albedo_dif(2) = nm%albsfcdif(2)   ! AlbedoNIRDif
                stateName='AlbedoNIRDif', &                             
                units='-', &								! CLSM: ALBNF, 'surface_reflectivity_near_infrared_diffuse' 								
                adImport=.false., &	
                adExport=.true., &
                directConn=.false.), &

!--- surfacelayer variable ---
LIS_Field( &
                stdName='friction_velocity_vegetated', &        ! Noah-MP: FrictionVelVeg, CLSM exports bulk USTAR=FrictionVelVeg*VegFrac + FrictionVelBare*(1-VegFrac)
                stateName='UstarVeg', &                         ! NoahMP50_struc(n)%noahmp50(t)%ustarveg
                units='m/s', &						! GEOS/CLSM Variable: USTAR (friction_velocity)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='friction_velocity_bare', &             ! Noah-MP: FrictionVelBare, CLSM exports bulk USTAR=FrictionVelVeg*VegFrac + FrictionVelBare*(1-VegFrac)
                stateName='UstarBare', &                        ! NoahMP50_struc(n)%noahmp50(t)%ustarbare 
                units='m/s', &						! GEOS/CLSM Variable: USTAR (friction_velocity)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='friction_velocity', &             
                stateName='Ustar', &                        ! Calculated from conversion script: gx%ustar = fveg*nm%ustarveg + (1.0-fveg)*nm%ustarbare   
                units='m/s', &						! GEOS/CLSM Variable: USTAR (friction_velocity)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='eastward_surface_stress_on_air', &     ! Noah-MP: WindStressEwSfc 
                stateName='Taux', &                             ! NoahMP50_struc(n)%noahmp50(t)%taux
                units='N/m2', &                                 ! CLSM: TAUX (eastward_surface_stress_on_air)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='northward_surface_stress_on_air', &    ! Noah-MP: WindStressNsSfc 
                stateName='Tauy', &                             !  NoahMP50_struc(n)%noahmp50(t)%tauy
                units='N/m2', &                                 ! CLSM: northward_surface_stress_on_air (Tauy)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

!!=========================================================================================
!! Met Forcings Imports 
!!=========================================================================================
LIS_Field( &
                stdName='surface_air_pressure', &                  ! Noah-MP: PressureAirRefHeight
                stateName='Psurf_f', &						! NoahMP50_struc(n)%noahmp50(t)%psurf
                units='Pa', &						! CLSM: surface_air_pressure (Psurf_f)
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='air_temperature', &           ! Noah-MP: TemperatureAirRefHeight
                stateName='Tair_f', &				! NoahMP50_struc(n)%noahmp50(t)%tair
                units='K', &						! CLSM: TA (surface_air_temperature)
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='specific_humidity', &     ! Noah-MP: SpecHumidityRefHeight
                stateName='Qair_f', &				! NoahMP50_struc(n)%noahmp50(t)%qair
                units='kg/kg', &				! CLSM: QA (surface_air_specific_humidity)
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_downwelling_longwave_flux_in_air', &    ! Noah-MP: RadLwDownRefHeight
                stateName='LWdown_f', &								! NoahMP50_struc(n)%noahmp50(t)%lwdown
                units='W/m2', &								! CLSM: LWDNSRF, 'surface_absorbed_longwave_flux' (imported by surface)
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_downwelling_shortwave_flux_in_air', & ! Noah-MP: RadSwDownRefHeight
                stateName='SWdown_f', &								! NoahMP50_struc(n)%noahmp50(t)%swdown 
                units='W/m2', &							
                adImport=.true., &
                adExport=. true., &						! CLSM: (DRPAR + DFPAR + DRNIR + DFNIR + DRUVR + DFUVR) (from Rad), in CLSM: SWDOWNLAND = DRPAR + DFPAR + DRUVR + DFUVR + DRNIR + DFNIR, conversion needed to bring GEOS SWdown partitions to bulk SWdown required by NoahMP 
                directConn=.false.), &

LIS_Field( &
                stdName='eastward_wind', &        ! Noah-MP: WindEastwardRefHeight
                stateName='EWind_f', &			! NoahMP50_struc(n)%noahmp50(t)%wind_e
                units='m/s', &				! CLSM: UL (eastward_wind_bottom_level)
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='northward_wind', &       ! Noah-MP: WindNorthwardRefHeight
                stateName='NWind_f', &				! NoahMP50_struc(n)%noahmp50(t)%wind_n
                units='m/s', &					! CLSM: VL (northward_wind_bottom_level)
                adImport=.true., &	
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='snowfall_flux', &                          ! Noah-MP: PrecipSnowRefHeight
                stateName='snowf_f', &							! NoahMP50_struc(n)%noahmp50(t)%snowf
                units='kg/m2s', &						! CLSM exports: snowfall_land
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

!! --- Precipitation Forcing (All other physics options other than 4) ---
LIS_Field( &
                stdName='total_precipitation', &               ! Noah-MP: PrecipTotRefHeight
                stateName='TotalPrecip', &					! NoahMP50_struc(n)%noahmp50(t)%prcp 
                units='kg/m2s', &						! CLSM: TPREC, total_precipitation
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &
!!==============================================================================================================================
!!-------------did not find these variable names from the LIS side (talk to Jim) start--------------------------------------------
!!==============================================================================================================================
!! --- Precipitation Forcings (Microphysics Option 4) ---

LIS_Field( &
                stdName='convective_precipitation', &          ! Noah-MP: PrecipConvRefHeight
                stateName='RainfConv_f', &					! Not mapped in NoahMP50_main.F90
                units='kg/m2s', &
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='shallow_convective_precipitation', &  ! Noah-MP: PrecipShConvRefHeight
                stateName='RainfShConv_f', &				! Not mapped in NoahMP50_main.F90
                units='kg/m2s', &
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='liquid_water_large_scale_precipitation', & ! Noah-MP: PrecipNonConvRefHeight
                stateName='RainfLS_f', &						! Not mapped in NoahMP50_main.F90
                units='kg/m2s', &
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='freezing_rain_fall', &                ! Noah-MP: PrecipGraupelRefHeight
                stateName='Graupelf', &						! Not mapped in NoahMP50_main.F90
                units='kg/m2s', &
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &

LIS_Field( &
                stdName='icefall', &                           ! Noah-MP: PrecipHailRefHeight
                stateName='Hailf', &						! Not mapped in NoahMP50_main.F90
                units='kg/m2s', &
                adImport=.true., &
                adExport=.false., &
                directConn=.false.), &
!!==============================================================================================================================
!!-------------did not find these variable names from the LIS side (talk to Jim) end--------------------------------------------
!!==============================================================================================================================

!!=========================================================================================
!! Surface exports to Chem
!!=========================================================================================
LIS_Field( &
                stdName='leaf_area_index', &                   ! Noah-MP: LeafAreaIndex
                stateName='LAI', &							! NoahMP50_struc(n)%noahmp50(t)%lai
                units='-', &				
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='stem_area_index', &                   ! Noah-MP: StemAreaIndex
                stateName='SAI', &							! NoahMP50_struc(n)%noahmp50(t)%sai
                units='-', &					
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &


LIS_Field( &
                stdName='vegetation_area_index', &                   ! Conversion script: gx%lai = nm%LAI + nm%SAI
                stateName='VAI', &							
                units='-', &					! CLSM: LAI, 'leaf_area_index'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='greenness_fraction', &                ! Noah-MP: VegFrac
                stateName='Greenness', &					! NoahMP50_struc(n)%noahmp50(t)%fveg
                units='-', &						! CLSM: GRN, 'vegetation_greenness_fraction'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='vegetation_type', &                   ! Noah-MP: VegType - we need to determine how NoahMP LCs map to GEOS LCs
                stateName='Landcover', &					! 	NoahMP50_struc(n)%noahmp50(t)%vegetype
                units='-', &							! CLSM: ITY (vegetation_type) or NITY (NCEP_vegetation_type) - conversion needed to map NoahMP categories to CLSM categories (Rolf may be able to inform this?)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='latent_heat_flux', &                  ! Noah-MP: HeatLatentCanopy + HeatLatentGrd + HeatLatentTransp + HeatLatentIrriEvap
                stateName='Qle', &							! NoahMP50_struc(n)%noahmp50(t)%lh
                units='W/m2', &						! CLSM: GC, 'total_latent_heat_flux_consistent_with_evaporation_from_turbulence' 
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='canopy_water_content', &              ! Noah-MP: CanopyTotalWater
                stateName='CanopInt', &						! Conversion calculated as: gx%canopint = nm%canliq + nm%canice
                units='kg/m2', &							
                adImport=.false., &							! CLSM: CAPAC, 'vegetation_interception_water_storage'  
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='net_ecosystem_exchange', &            ! Noah-MP: NetEcoExchange
                stateName='NEE', &							! NoahMP50_struc(n)%noahmp50(t)%nee
                units='g/m2s', &						!CLSM, CNNEE, 'CN_net_ecosystem_exchange'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='soil_temperature_layer_1', &          ! Noah-MP: TemperatureSoilSnow(1)
                stateName='SoilTemp1', &					! NoahMP50_struc(n)%noahmp50(t)%tslb(1)
                units='K', &							! CLSM: TSOIL1, soil_temperature_layer_1
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='soil_moisture_layer_1', &             ! Noah-MP: SoilMoisture(1)
                stateName='SoilMoist1', &					! NoahMP50_struc(n)%noahmp50(t)%smc(1)
                units='m3/m3', &						! CLSM, WCSF , 'soil_moisture_surface'
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='soil_moisture_rootzone', &             ! Noah-MP: SoilMoisture
                stateName='SoilMoistRZ', &				! NoahMP50_struc(n)%noahmp50(t)%smc(1)*0.1 + NoahMP50_struc(n)%noahmp50(t)%smc(2)*0.3 + NoahMP50_struc(n)%noahmp50(t)%smc*0.6 =
                units='m3/m3', &						!conversion needed for CLSM: WCRZ (soil_moisture_rootzone)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='snow_depth', &                        ! Noah-MP: SnowDepth
                stateName='SnowDepth', &					! NoahMP50_struc(n)%noahmp50(t)%snowh
                units='m', &							! CLSM equivalent: snow_depth_within_snow_covered_area_fraction_on_land, conversion needed using SnowCoverFrac b/c noah-mp represents pixel-mean snow depth.
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='snow_cover_fraction', &               ! Noah-MP: SnowCoverFrac
                stateName='Snowcover', &					! NoahMP50_struc(n)%noahmp50(t)%snowc
                units='-', &							! CLSM exports: fractional_area_of_snow_on_land
                adImport=.false., &
                adExport=.true., &
                directConn=.false.), &

LIS_Field( &
                stdName='surface_air_density', &               ! Noah-MP: DensityAirRefHeight
                stateName='rhos', &						! NoahMP50_struc(n)%noahmp50(t)%rhoair
                units='-', &						! CLSM: RHOS (Surface air density)
                adImport=.false., &
                adExport=.true., &
                directConn=.false.)

!!==============================================================================================================================
!! Variables that GEOS-surface exports that Noah-MP does not solve for - all of which are for Chem component:
!!==============================================================================================================================
! BSTAR (surface buoyancy scale) | This is calculated in Surface GridComp w/ variables that NoahMP does export.
! Land/ocean fractions (FRLAND, FROCEAN)
! U10M / V10M (10-meter wind speed components) - solved with surfacelayer.F90
! U10N / V10N (10-meter neutral wind speed components) - solved with GEOS_SurfaceGridComp.F90
! T10M (10-meter air temperature) - solved with surfacelayer.F90
! Q10M (10-meter specific humidity) - solved with surfacelayer.F90
! Surface: CN_PRCP (Convective precipitation) - solved with moisture component


!!==============================================================================================================================
!! Variables that GEOS-surface imports that Noah-MP solves for or does not require - do not import these:
!!==============================================================================================================================
! CLSM import DEWL, whereas it looks like Noah-MP prognostically solves for the components of this variable: CondenseVapGrd, DewCanopyLiq, DewSoilSfcLiq
! CLSM imports SH and EVAP whereas noah-mp prognostically solves for these and will provide them to turbulence via export (see above)
! CLSM import FRSL, whereas noah-mp prognostically solves for its partitions: FrostCanopyIce & FrostSnowSfcIce
! CLSM imports a suite of derivatives used in the implicit solver which noah-mp will not need: DSH, DFU, DFV, DEVAP, DDEWL, DFRSL
           
   end function new_LIS_FieldList

end module LIS_FieldCatalog
