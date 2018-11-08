! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE var_lookup
 ! defines named variables used to index array elements
 USE nrtype, integerMissing=>nr_integerMissing
 implicit none
 private
 ! local variables
 integer(i4b),parameter,public     :: numStats = 7                 ! number of output stats
 integer(i4b),parameter,public     :: maxFreq = 10                 ! maximum number of output streams
 integer(i4b),parameter            :: ixVal=1                      ! an example integer
 integer(i4b),parameter            :: iLength=storage_size(ixVal)  ! size of the example integer

 ! ***************************************************************************************
 ! (0) define model decisions
 ! ***************************************************************************************
 type, public  ::  iLook_decision
  integer(i4b)    :: simulStart = integerMissing     ! simulation start time
  integer(i4b)    :: simulFinsh = integerMissing     ! simulation end time
  integer(i4b)    :: soilCatTbl = integerMissing     ! soil-category dateset
  integer(i4b)    :: vegeParTbl = integerMissing     ! vegetation category dataset
  integer(i4b)    :: soilStress = integerMissing     ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist = integerMissing     ! choice of function for stomatal resistance
  integer(i4b)    :: bbTempFunc = integerMissing     ! Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  integer(i4b)    :: bbHumdFunc = integerMissing     ! Ball-Berry: humidity controls on stomatal resistance
  integer(i4b)    :: bbElecFunc = integerMissing     ! Ball-Berry: dependence of photosynthesis on PAR
  integer(i4b)    :: bbCO2point = integerMissing     ! Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  integer(i4b)    :: bbNumerics = integerMissing     ! Ball-Berry: iterative numerical solution method
  integer(i4b)    :: bbAssimFnc = integerMissing     ! Ball-Berry: controls on carbon assimilation
  integer(i4b)    :: bbCanIntg8 = integerMissing     ! Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  integer(i4b)    :: num_method = integerMissing     ! choice of numerical method
  integer(i4b)    :: fDerivMeth = integerMissing     ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method = integerMissing     ! method used to determine LAI and SAI
  integer(i4b)    :: cIntercept = integerMissing     ! choice of parameterization for canopy interception
  integer(i4b)    :: f_Richards = integerMissing     ! form of richards' equation
  integer(i4b)    :: groundwatr = integerMissing     ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile = integerMissing     ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn = integerMissing     ! type of upper boundary condition for thermodynamics
  integer(i4b)    :: bcLowrTdyn = integerMissing     ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH = integerMissing     ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH = integerMissing     ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits = integerMissing     ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: rootProfil = integerMissing     ! choice of parameterization for the rooting profile
  integer(i4b)    :: canopyEmis = integerMissing     ! choice of parameterization for canopy emissivity
  integer(i4b)    :: snowIncept = integerMissing     ! choice of parameterization for snow interception
  integer(i4b)    :: windPrfile = integerMissing     ! choice of canopy wind profile
  integer(i4b)    :: astability = integerMissing     ! choice of stability function
  integer(i4b)    :: canopySrad = integerMissing     ! choice of method for canopy shortwave radiation
  integer(i4b)    :: alb_method = integerMissing     ! choice of albedo representation
  integer(i4b)    :: snowLayers = integerMissing     ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction = integerMissing     ! choice of compaction routine
  integer(i4b)    :: thCondSnow = integerMissing     ! choice of thermal conductivity representation for snow
  integer(i4b)    :: thCondSoil = integerMissing     ! choice of thermal conductivity representation for soil
  integer(i4b)    :: spatial_gw = integerMissing     ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting = integerMissing     ! choice of method for sub-grid routing
  integer(i4b)    :: snowDenNew = integerMissing     ! choice of method for new snow density
 endtype iLook_decision

 ! ***********************************************************************************************************
 ! (1) define model time
 ! ***********************************************************************************************************
 type, public  ::  iLook_time
  integer(i4b)    :: iyyy       = integerMissing     ! year
  integer(i4b)    :: im         = integerMissing     ! month
  integer(i4b)    :: id         = integerMissing     ! day
  integer(i4b)    :: ih         = integerMissing     ! hour
  integer(i4b)    :: imin       = integerMissing     ! minute
 endtype iLook_time

 ! ***********************************************************************************************************
 ! (2) define model forcing data
 ! ***********************************************************************************************************
 type, public  ::  iLook_force
  integer(i4b)    :: time       = integerMissing     ! time since time reference       (s)
  integer(i4b)    :: pptrate    = integerMissing     ! precipitation rate              (kg m-2 s-1)
  integer(i4b)    :: airtemp    = integerMissing     ! air temperature                 (K)
  integer(i4b)    :: spechum    = integerMissing     ! specific humidity               (g/g)
  integer(i4b)    :: windspd    = integerMissing     ! windspeed                       (m/s)
  integer(i4b)    :: SWRadAtm   = integerMissing     ! downwelling shortwave radiaiton (W m-2)
  integer(i4b)    :: LWRadAtm   = integerMissing     ! downwelling longwave radiation  (W m-2)
  integer(i4b)    :: airpres    = integerMissing     ! pressure                        (Pa)
 endtype iLook_force

 ! ***********************************************************************************************************
 ! (3) define local attributes
 ! ***********************************************************************************************************
 type, public  ::  iLook_attr
  integer(i4b)    :: latitude      = integerMissing  ! latitude (degrees north)
  integer(i4b)    :: longitude     = integerMissing  ! longitude (degrees east)
  integer(i4b)    :: elevation     = integerMissing  ! elevation (m)
  integer(i4b)    :: tan_slope     = integerMissing  ! tan water table slope, taken as tan local ground surface slope (-)
  integer(i4b)    :: contourLength = integerMissing  ! length of contour at downslope edge of HRU (m)
  integer(i4b)    :: HRUarea       = integerMissing  ! area of each HRU  (m2)
  integer(i4b)    :: mHeight       = integerMissing  ! measurement height above bare ground (m)
 end type iLook_attr

 ! ***********************************************************************************************************
 ! (4) define local classification of veg, soil, etc.
 ! ***********************************************************************************************************
 type, public  ::  iLook_type
  integer(i4b)    :: hruIndex      = integerMissing  ! index defining hydrologic response unit (-)
  integer(i4b)    :: vegTypeIndex  = integerMissing  ! index defining vegetation type (-)
  integer(i4b)    :: soilTypeIndex = integerMissing  ! index defining soil type (-)
  integer(i4b)    :: slopeTypeIndex= integerMissing  ! index defining slope (-)
  integer(i4b)    :: downHRUindex  = integerMissing  ! index of downslope HRU (0 = basin outlet)
 end type iLook_type

 ! ***********************************************************************************************************
 ! (5) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead        = integerMissing    ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead        = integerMissing    ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta       = integerMissing    ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta       = integerMissing    ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp        = integerMissing    ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp        = integerMissing    ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain          = integerMissing    ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep     = integerMissing    ! temperature range over the time step (K)
  integer(i4b)    :: frozenPrecipMultip    = integerMissing    ! frozen precipitation multiplier (-)
  ! snow properties
  integer(i4b)    :: snowfrz_scale         = integerMissing    ! scaling parameter for the freezing curve for snow (K-1)
  integer(i4b)    :: fixedThermalCond_snow = integerMissing    ! fixed thermal conductivity for snow (W m-1 K-1)
  ! snow albedo
  integer(i4b)    :: albedoMax             = integerMissing    ! maximum snow albedo for a single spectral band (-)
  integer(i4b)    :: albedoMinWinter       = integerMissing    ! minimum snow albedo during winter for a single spectral band (-)
  integer(i4b)    :: albedoMinSpring       = integerMissing    ! minimum snow albedo during spring for a single spectral band (-)
  integer(i4b)    :: albedoMaxVisible      = integerMissing    ! maximum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMinVisible      = integerMissing    ! minimum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMaxNearIR       = integerMissing    ! maximum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoMinNearIR       = integerMissing    ! minimum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoDecayRate       = integerMissing    ! albedo decay rate (s)
  integer(i4b)    :: albedoSootLoad        = integerMissing    ! soot load factor (-)
  integer(i4b)    :: albedoRefresh         = integerMissing    ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow           = integerMissing    ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: directScale           = integerMissing    ! scaling factor for fractional driect radiaion parameterization (-)
  integer(i4b)    :: Frad_direct           = integerMissing    ! maximum fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis              = integerMissing    ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin         = integerMissing    ! minimum new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMult        = integerMissing    ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal        = integerMissing    ! scaling factor for new snow density (K)
  integer(i4b)    :: constSnowDen          = integerMissing    ! constDens, Constant new snow density (kg m-3)
  integer(i4b)    :: newSnowDenAdd         = integerMissing    ! Pahaut 1976, additive factor for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMultTemp    = integerMissing    ! Pahaut 1976, multiplier for new snow density applied to air temperature (kg m-3 K-1)
  integer(i4b)    :: newSnowDenMultWind    = integerMissing    ! Pahaut 1976, multiplier for new snow density applied to wind speed (kg m-7/2 s-1/2)
  integer(i4b)    :: newSnowDenMultAnd     = integerMissing    ! Anderson 1976, multiplier for new snow density for Anderson function (K-1)
  integer(i4b)    :: newSnowDenBase        = integerMissing    ! Anderson 1976, base value that is rasied to the (3/2) power (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth        = integerMissing    ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth        = integerMissing    ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate       = integerMissing    ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn        = integerMissing    ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn        = integerMissing    ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: baseViscosity         = integerMissing    ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil                = integerMissing    ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow                = integerMissing    ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp                = integerMissing    ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow                = integerMissing    ! roughness length of snow (m)
  integer(i4b)    :: z0Soil                = integerMissing    ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy              = integerMissing    ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction           = integerMissing    ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber        = integerMissing    ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam        = integerMissing    ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar         = integerMissing    ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale        = integerMissing    ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff     = integerMissing    ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam    = integerMissing    ! canopy wind reduction parameter (-)
  ! stomatal conductance
  integer(i4b)    :: Kc25                  = integerMissing    ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  integer(i4b)    :: Ko25                  = integerMissing    ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  integer(i4b)    :: Kc_qFac               = integerMissing    ! factor in the q10 function defining temperature controls on Kc (-)
  integer(i4b)    :: Ko_qFac               = integerMissing    ! factor in the q10 function defining temperature controls on Ko (-)
  integer(i4b)    :: kc_Ha                 = integerMissing    ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  integer(i4b)    :: ko_Ha                 = integerMissing    ! activation energy for the Michaelis-Menten constant for O2 (J mol-1)
  integer(i4b)    :: vcmax25_canopyTop     = integerMissing    ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  integer(i4b)    :: vcmax_qFac            = integerMissing    ! factor in the q10 function defining temperature controls on vcmax (-)
  integer(i4b)    :: vcmax_Ha              = integerMissing    ! activation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Hd              = integerMissing    ! deactivation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Sv              = integerMissing    ! entropy term in the vcmax function (J mol-1 K-1)
  integer(i4b)    :: vcmax_Kn              = integerMissing    ! foliage nitrogen decay coefficient (-)
  integer(i4b)    :: jmax25_scale          = integerMissing    ! scaling factor to relate jmax25 to vcmax25 (-)
  integer(i4b)    :: jmax_Ha               = integerMissing    ! activation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Hd               = integerMissing    ! deactivation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Sv               = integerMissing    ! entropy term in the jmax function (J mol-1 K-1)
  integer(i4b)    :: fractionJ             = integerMissing    ! fraction of light lost by other than the chloroplast lamellae (-)
  integer(i4b)    :: quantamYield          = integerMissing    ! quantam yield (mol e mol-1 q)
  integer(i4b)    :: vpScaleFactor         = integerMissing    ! vapor pressure scaling factor in stomatal conductance function (Pa)
  integer(i4b)    :: cond2photo_slope      = integerMissing    ! slope of conductance-photosynthesis relationship (-)
  integer(i4b)    :: minStomatalConductance= integerMissing    ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  integer(i4b)    :: winterSAI             = integerMissing    ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI             = integerMissing    ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootScaleFactor1      = integerMissing    ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootScaleFactor2      = integerMissing    ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootingDepth          = integerMissing    ! rooting depth (m)
  integer(i4b)    :: rootDistExp           = integerMissing    ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi          = integerMissing    ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam       = integerMissing    ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting       = integerMissing    ! critical vol. liq. water content when plants are wilting (-)
  integer(i4b)    :: critSoilTranspire     = integerMissing    ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire  = integerMissing    ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: minStomatalResistance = integerMissing    ! minimum canopy resistance (s m-1)
  integer(i4b)    :: leafDimension         = integerMissing    ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop       = integerMissing    ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom    = integerMissing    ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: specificHeatVeg       = integerMissing    ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation     = integerMissing    ! maximum mass of vegetation (full foliage) (kg m-2)
  integer(i4b)    :: throughfallScaleSnow  = integerMissing    ! scaling factor for throughfall (snow) (-)
  integer(i4b)    :: throughfallScaleRain  = integerMissing    ! scaling factor for throughfall (rain) (-)
  integer(i4b)    :: refInterceptCapSnow   = integerMissing    ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  integer(i4b)    :: refInterceptCapRain   = integerMissing    ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  integer(i4b)    :: snowUnloadingCoeff    = integerMissing    ! time constant for unloading of snow from the forest canopy (s-1)
  integer(i4b)    :: canopyDrainageCoeff   = integerMissing    ! time constant for drainage of liquid water from the forest canopy (s-1)
  integer(i4b)    :: ratioDrip2Unloading   = integerMissing    ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  integer(i4b)    :: canopyWettingFactor   = integerMissing    ! maximum wetted fraction of the canopy (-)
  integer(i4b)    :: canopyWettingExp      = integerMissing    ! exponent in canopy wetting function (-)
  ! soil properties
  integer(i4b)    :: soil_dens_intr        = integerMissing    ! intrinsic soil density (kg m-3)
  integer(i4b)    :: thCond_soil           = integerMissing    ! thermal conductivity of soil (W m-1 K-1)
  integer(i4b)    :: frac_sand             = integerMissing    ! fraction of sand (-)
  integer(i4b)    :: frac_silt             = integerMissing    ! fraction of silt (-)
  integer(i4b)    :: frac_clay             = integerMissing    ! fraction of clay (-)
  integer(i4b)    :: fieldCapacity         = integerMissing    ! field capacity (-)
  integer(i4b)    :: wettingFrontSuction   = integerMissing    ! Green-Ampt wetting front suction (m)
  integer(i4b)    :: theta_mp              = integerMissing    ! volumetric liquid water content when macropore flow begins (-)
  integer(i4b)    :: theta_sat             = integerMissing    ! porosity (-)
  integer(i4b)    :: theta_res             = integerMissing    ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha             = integerMissing    ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                 = integerMissing    ! van Genuchten "n" parameter (-)
  integer(i4b)    :: mpExp                 = integerMissing    ! empirical exponent in macropore flow equation (-)
  integer(i4b)    :: k_soil                = integerMissing    ! hydraulic conductivity of soil (m s-1)
  integer(i4b)    :: k_macropore           = integerMissing    ! saturated hydraulic conductivity for macropores (m s-1)
  integer(i4b)    :: kAnisotropic          = integerMissing    ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL       = integerMissing    ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth        = integerMissing    ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor    = integerMissing    ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp    = integerMissing    ! baseflow exponent (-)
  integer(i4b)    :: qSurfScale            = integerMissing    ! scaling factor in the surface runoff parameterization (-)
  integer(i4b)    :: specificYield         = integerMissing    ! specific yield (-)
  integer(i4b)    :: specificStorage       = integerMissing    ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede              = integerMissing    ! ice impedence factor (-)
  integer(i4b)    :: soilIceScale          = integerMissing    ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  integer(i4b)    :: soilIceCV             = integerMissing    ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind               = integerMissing    ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep               = integerMissing    ! minimum length of the time step
  integer(i4b)    :: maxstep               = integerMissing    ! maximum length of the time step
  integer(i4b)    :: wimplicit             = integerMissing    ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter               = integerMissing    ! maximum number of iteration
  integer(i4b)    :: relConvTol_liquid     = integerMissing    ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid     = integerMissing    ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric     = integerMissing    ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric     = integerMissing    ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy     = integerMissing    ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy     = integerMissing    ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr     = integerMissing    ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr     = integerMissing    ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                  = integerMissing    ! minimum layer depth (m)
  integer(i4b)    :: zmax                  = integerMissing    ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1            = integerMissing    ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2            = integerMissing    ! minimum layer depth for the 2nd layer (m)
  integer(i4b)    :: zminLayer3            = integerMissing    ! minimum layer depth for the 3rd layer (m)
  integer(i4b)    :: zminLayer4            = integerMissing    ! minimum layer depth for the 4th layer (m)
  integer(i4b)    :: zminLayer5            = integerMissing    ! minimum layer depth for the 5th (bottom) layer (m)
  integer(i4b)    :: zmaxLayer1_lower      = integerMissing    ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  integer(i4b)    :: zmaxLayer2_lower      = integerMissing    ! maximum layer depth for the 2nd layer when only 2 layers (m)
  integer(i4b)    :: zmaxLayer3_lower      = integerMissing    ! maximum layer depth for the 3rd layer when only 3 layers (m)
  integer(i4b)    :: zmaxLayer4_lower      = integerMissing    ! maximum layer depth for the 4th layer when only 4 layers (m)
  integer(i4b)    :: zmaxLayer1_upper      = integerMissing    ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  integer(i4b)    :: zmaxLayer2_upper      = integerMissing    ! maximum layer depth for the 2nd layer when > 2 layers (m)
  integer(i4b)    :: zmaxLayer3_upper      = integerMissing    ! maximum layer depth for the 3rd layer when > 3 layers (m)
  integer(i4b)    :: zmaxLayer4_upper      = integerMissing    ! maximum layer depth for the 4th layer when > 4 layers (m)
 endtype ilook_param


 ! ***********************************************************************************************************
 ! (6) define model prognostic (state) variables
 ! ***********************************************************************************************************
 type, public :: iLook_prog
  ! variables for time stepping
  integer(i4b)    :: dt_init                     = integerMissing    ! length of initial time step at start of next data interval (s)
  ! state variables for vegetation
  integer(i4b)    :: scalarCanopyIce             = integerMissing    ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq             = integerMissing    ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyWat             = integerMissing    ! mass of total water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanairTemp            = integerMissing    ! temperature of the canopy air space (Pa)
  integer(i4b)    :: scalarCanopyTemp            = integerMissing    ! temperature of the vegetation canopy (K)
  ! state variables for snow
  integer(i4b)    :: spectralSnowAlbedoDiffuse   = integerMissing    ! diffuse snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarSnowAlbedo            = integerMissing    ! snow albedo for the entire spectral band (-)
  integer(i4b)    :: scalarSnowDepth             = integerMissing    ! total snow depth (m)
  integer(i4b)    :: scalarSWE                   = integerMissing    ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond           = integerMissing    ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  ! state variables for the snow+soil domain
  integer(i4b)    :: mLayerTemp                  = integerMissing    ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracIce            = integerMissing    ! volumetric fraction of ice in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq            = integerMissing    ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolFracWat            = integerMissing    ! volumetric fraction of total water in each layer (-)
  integer(i4b)    :: mLayerMatricHead            = integerMissing    ! matric head of water in the soil (m)
  ! other state variables
  integer(i4b)    :: scalarAquiferStorage        = integerMissing    ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarSurfaceTemp           = integerMissing    ! surface temperature (K)
  ! coordinate variables
  integer(i4b)    :: mLayerDepth                 = integerMissing    ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                = integerMissing    ! height at the mid-point of each layer (m)
  integer(i4b)    :: iLayerHeight                = integerMissing    ! height of the layer interface; top of soil = 0 (m)
 endtype iLook_prog

 ! ***********************************************************************************************************
 ! (7) define diagnostic variables
 ! ***********************************************************************************************************
 type, public :: iLook_diag
  ! local properties
  integer(i4b)    :: scalarCanopyDepth               = integerMissing ! canopy depth (m)
  integer(i4b)    :: scalarGreenVegFraction          = integerMissing ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg         = integerMissing ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarCanopyEmissivity          = integerMissing ! effective canopy emissivity (-)
  integer(i4b)    :: scalarRootZoneTemp              = integerMissing ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                       = integerMissing ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                       = integerMissing ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                = integerMissing ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                = integerMissing ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarCanopyIceMax              = integerMissing ! maximum interception storage capacity for ice (kg m-2)
  integer(i4b)    :: scalarCanopyLiqMax              = integerMissing ! maximum interception storage capacity for liquid water (kg m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex        = integerMissing ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVolHtCap_air              = integerMissing ! volumetric heat capacity air (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = integerMissing ! volumetric heat capacity ice (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = integerMissing ! volumetric heat capacity dry soil (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = integerMissing ! volumetric heat capacity liquid wat (J m-3 K-1)
  integer(i4b)    :: mLayerVolHtCapBulk              = integerMissing ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = integerMissing ! thermal conductivity of dry soil     (W m-1 K-1)
  integer(i4b)    :: scalarLambda_wetsoil            = integerMissing ! thermal conductivity of wet soil     (W m-1 K-1)
  integer(i4b)    :: mLayerThermalC                  = integerMissing ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerThermalC                  = integerMissing ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! forcing
  integer(i4b)    :: scalarVPair                     = integerMissing ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarVP_CanopyAir              = integerMissing ! vapor pressure of the canopy air space (Pa)
  integer(i4b)    :: scalarTwetbulb                  = integerMissing ! wet bulb temperature (K)
  integer(i4b)    :: scalarSnowfallTemp              = integerMissing ! temperature of fresh snow (K)
  integer(i4b)    :: scalarNewSnowDensity            = integerMissing ! density of fresh snow (kg m-3)
  integer(i4b)    :: scalarO2air                     = integerMissing ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                    = integerMissing ! atmospheric co2 concentration (Pa)
  ! shortwave radiation
  integer(i4b)    :: scalarCosZenith                 = integerMissing ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: scalarFractionDirect            = integerMissing ! fraction of direct radiation (0-1)
  integer(i4b)    :: scalarCanopySunlitFraction      = integerMissing ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI           = integerMissing ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI           = integerMissing ! shaded leaf area (-)
  integer(i4b)    :: spectralAlbGndDirect            = integerMissing ! direct  albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: spectralAlbGndDiffuse           = integerMissing ! diffuse albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: scalarGroundAlbedo              = integerMissing ! albedo of the ground surface (-)
  ! turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy       = integerMissing ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround       = integerMissing ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp          = integerMissing ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp          = integerMissing ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                  = integerMissing ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor       = integerMissing ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement     = integerMissing ! zero plane displacement (m)
  integer(i4b)    :: scalarRiBulkCanopy              = integerMissing ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround              = integerMissing ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection = integerMissing ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection = integerMissing ! stability correction for the ground surface (-)
  ! evapotranspiration
  integer(i4b)    :: scalarIntercellularCO2Sunlit    = integerMissing ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  integer(i4b)    :: scalarIntercellularCO2Shaded    = integerMissing ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  integer(i4b)    :: scalarTranspireLim              = integerMissing ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr          = integerMissing ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor     = integerMissing ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarSoilRelHumidity           = integerMissing ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: mLayerTranspireLim              = integerMissing ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerRootDensity               = integerMissing ! fraction of roots in each soil layer (-)
  integer(i4b)    :: scalarAquiferRootFrac           = integerMissing ! fraction of roots below the soil profile (-)
  ! canopy hydrology
  integer(i4b)    :: scalarFracLiqVeg                = integerMissing ! fraction of liquid water on vegetation (-)
  integer(i4b)    :: scalarCanopyWetFraction         = integerMissing ! fraction of canopy that is wet
  ! snow hydrology
  integer(i4b)    :: scalarSnowAge                   = integerMissing ! non-dimensional snow age (-)
  integer(i4b)    :: scalarGroundSnowFraction        = integerMissing ! fraction of ground that is covered with snow (-)
  integer(i4b)    :: spectralSnowAlbedoDirect        = integerMissing ! direct snow albedo for individual spectral bands (-)
  integer(i4b)    :: mLayerFracLiqSnow               = integerMissing ! fraction of liquid water in each snow layer (-)
  integer(i4b)    :: mLayerThetaResid                = integerMissing ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = integerMissing ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerMeltFreeze                = integerMissing ! change in ice content due to melt/freeze in each layer (kg m-3)
  ! soil hydrology
  integer(i4b)    :: scalarInfilArea                 = integerMissing ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea                = integerMissing ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: scalarSoilControl               = integerMissing ! soil control on infiltration: 1=controlling; 0=not (-)
  integer(i4b)    :: mLayerVolFracAir                = integerMissing ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerTcrit                     = integerMissing ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerCompress                  = integerMissing ! change in volumetric water content due to compression of soil (-)
  integer(i4b)    :: scalarSoilCompress              = integerMissing ! change in total soil storage due to compression of the soil matrix (kg m-2)
  integer(i4b)    :: mLayerMatricHeadLiq             = integerMissing ! matric potential of liquid water (m)
  ! mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = integerMissing ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = integerMissing ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = integerMissing ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = integerMissing ! total mass of ice in the soil (kg m-2)
  ! variable shortcuts
  integer(i4b)    :: scalarVGn_m                     = integerMissing ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = integerMissing ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolLatHt_fus              = integerMissing ! volumetric latent heat of fusion     (J m-3)
  ! number of function evaluations
  integer(i4b)    :: numFluxCalls                    = integerMissing ! number of flux calls (-)
 endtype iLook_diag

 ! ***********************************************************************************************************
 ! (8) define model fluxes
 ! ***********************************************************************************************************
 type, public :: iLook_flux
  ! net energy and mass fluxes for the vegetation domain
  integer(i4b)    :: scalarCanairNetNrgFlux          = integerMissing ! net energy flux for the canopy air space (W m-2)
  integer(i4b)    :: scalarCanopyNetNrgFlux          = integerMissing ! net energy flux for the vegetation canopy (W m-2)
  integer(i4b)    :: scalarGroundNetNrgFlux          = integerMissing ! net energy flux for the ground surface (W m-2)
  integer(i4b)    :: scalarCanopyNetLiqFlux          = integerMissing ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
  ! forcing
  integer(i4b)    :: scalarRainfall                  = integerMissing ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                  = integerMissing ! computed snowfall rate (kg m-2 s-1)
  ! shortwave radiation
  integer(i4b)    :: spectralIncomingDirect          = integerMissing ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse         = integerMissing ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarCanopySunlitPAR           = integerMissing ! average absorbed par for sunlit leaves (W m-2)
  integer(i4b)    :: scalarCanopyShadedPAR           = integerMissing ! average absorbed par for shaded leaves (W m-2)
  integer(i4b)    :: spectralBelowCanopyDirect       = integerMissing ! downward direct flux below veg layer for each spectral band  (W m-2)
  integer(i4b)    :: spectralBelowCanopyDiffuse      = integerMissing ! downward diffuse flux below veg layer for each spectral band (W m-2)
  integer(i4b)    :: scalarBelowCanopySolar          = integerMissing ! solar radiation transmitted below the canopy (W m-2)
  integer(i4b)    :: scalarCanopyAbsorbedSolar       = integerMissing ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar       = integerMissing ! solar radiation absorbed by ground (W m-2)
  ! longwave radiation
  integer(i4b)    :: scalarLWRadCanopy               = integerMissing ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround               = integerMissing ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy        = integerMissing ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground        = integerMissing ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound        = integerMissing ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound        = integerMissing ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground        = integerMissing ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy        = integerMissing ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound        = integerMissing ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy        = integerMissing ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy               = integerMissing ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround               = integerMissing ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound               = integerMissing ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! turbulent heat transfer
  integer(i4b)    :: scalarEddyDiffusCanopyTop       = integerMissing ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity          = integerMissing ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop          = integerMissing ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom       = integerMissing ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance          = integerMissing ! below canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarCanopyResistance          = integerMissing ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance            = integerMissing ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance            = integerMissing ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSenHeatTotal              = integerMissing ! sensible heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarSenHeatCanopy             = integerMissing ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround             = integerMissing ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatTotal              = integerMissing ! latent heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap         = integerMissing ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyTrans        = integerMissing ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround             = integerMissing ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarCanopyAdvectiveHeatFlux   = integerMissing ! heat advected to the canopy surface with rain + snow (W m-2)
  integer(i4b)    :: scalarGroundAdvectiveHeatFlux   = integerMissing ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  integer(i4b)    :: scalarCanopySublimation         = integerMissing ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation           = integerMissing ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b)    :: scalarStomResistSunlit          = integerMissing ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded          = integerMissing ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit      = integerMissing ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded      = integerMissing ! shaded photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarCanopyTranspiration       = integerMissing ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation         = integerMissing ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation         = integerMissing ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = integerMissing ! transpiration loss from each soil layer (kg m-2 s-1)
  ! liquid and solid water fluxes through the canopy
  integer(i4b)    :: scalarThroughfallSnow           = integerMissing ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarThroughfallRain           = integerMissing ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySnowUnloading       = integerMissing ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyLiqDrainage         = integerMissing ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyMeltFreeze          = integerMissing ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! energy fluxes and for the snow and soil domains
  integer(i4b)    :: iLayerConductiveFlux            = integerMissing ! conductive energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux             = integerMissing ! advective energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerNrgFlux                   = integerMissing ! energy flux at layer interfaces (W m-2)
  integer(i4b)    :: mLayerNrgFlux                   = integerMissing ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes for the snow domain
  integer(i4b)    :: scalarSnowDrainage              = integerMissing ! drainage from the bottom of the snow profile (m s-1)
  integer(i4b)    :: iLayerLiqFluxSnow               = integerMissing ! liquid flux at snow layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSnow               = integerMissing ! net liquid water flux for each snow layer (s-1)
  ! liquid water fluxes for the soil domain
  integer(i4b)    :: scalarRainPlusMelt              = integerMissing ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarMaxInfilRate              = integerMissing ! maximum infiltration rate (m s-1)
  integer(i4b)    :: scalarInfiltration              = integerMissing ! infiltration of water into the soil profile (m s-1)
  integer(i4b)    :: scalarExfiltration              = integerMissing ! exfiltration of water from the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = integerMissing ! surface runoff (m s-1)
  integer(i4b)    :: mLayerSatHydCondMP              = integerMissing ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  integer(i4b)    :: mLayerSatHydCond                = integerMissing ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: iLayerSatHydCond                = integerMissing ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: mLayerHydCond                   = integerMissing ! hydraulic conductivity in each soil layer (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = integerMissing ! liquid flux at soil layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSoil               = integerMissing ! net liquid water flux for each soil layer (s-1)
  integer(i4b)    :: mLayerBaseflow                  = integerMissing ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow              = integerMissing ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow             = integerMissing ! total outflow from each layer in a given soil column (m3 s-1)
  integer(i4b)    :: scalarSoilBaseflow              = integerMissing ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = integerMissing ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = integerMissing ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = integerMissing ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = integerMissing ! baseflow from the aquifer (m s-1)
 endtype iLook_flux

 ! ***********************************************************************************************************
 ! (9) define derivatives
 ! ***********************************************************************************************************
 type, public :: iLook_deriv
  ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanairNetFlux_dCanairTemp      = integerMissing ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dCanopyTemp      = integerMissing ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dGroundTemp      = integerMissing ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanairTemp      = integerMissing ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanopyTemp      = integerMissing ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dGroundTemp      = integerMissing ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanLiq          = integerMissing ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  integer(i4b)    :: dGroundNetFlux_dCanairTemp      = integerMissing ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanopyTemp      = integerMissing ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dGroundTemp      = integerMissing ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanLiq          = integerMissing ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanopyEvaporation_dTCanair     = integerMissing ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTCanopy     = integerMissing ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTGround     = integerMissing ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dCanLiq      = integerMissing ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
  integer(i4b)    :: dGroundEvaporation_dTCanair     = integerMissing ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTCanopy     = integerMissing ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTGround     = integerMissing ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dCanLiq      = integerMissing ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
  ! derivatives in canopy water w.r.t canopy temperature
  integer(i4b)    :: dTheta_dTkCanopy                = integerMissing ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dCanLiq_dTcanopy                = integerMissing ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
  ! derivatives in canopy liquid fluxes w.r.t. canopy water
  integer(i4b)    :: scalarCanopyLiqDeriv            = integerMissing ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarThroughfallRainDeriv      = integerMissing ! derivative in throughfall w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarCanopyLiqDrainageDeriv    = integerMissing ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
  integer(i4b)    :: dNrgFlux_dTempAbove             = integerMissing ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  integer(i4b)    :: dNrgFlux_dTempBelow             = integerMissing ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
  integer(i4b)    :: iLayerLiqFluxSnowDeriv          = integerMissing ! derivative in vertical liquid water flux at layer interfaces (m s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  integer(i4b)    :: dVolTot_dPsi0                   = integerMissing ! derivative in total water content w.r.t. total water matric potential (m-1)
  integer(i4b)    :: dq_dHydStateAbove               = integerMissing ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dHydStateBelow               = integerMissing ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dPsi               = integerMissing ! derivative in the soil water characteristic w.r.t. psi (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = integerMissing ! derivative in the soil water characteristic w.r.t. theta (m)
  integer(i4b)    :: dCompress_dPsi                  = integerMissing ! derivative in compressibility w.r.t matric head (m-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  integer(i4b)    :: dq_dNrgStateAbove               = integerMissing ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dNrgStateBelow               = integerMissing ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dTk                = integerMissing ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dPsiLiq_dTemp                   = integerMissing ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
  integer(i4b)    :: dPsiLiq_dPsi0                   = integerMissing ! derivative in liquid water matric potential w.r.t. the total water matric potential (-)
 endtype iLook_deriv

 ! ***********************************************************************************************************
 ! (10) define model indices
 ! ***********************************************************************************************************
 type, public :: iLook_index
 ! number of model layers, and layer indices
 integer(i4b)     :: nSnow              = integerMissing  ! number of snow layers                                                    (-) 
 integer(i4b)     :: nSoil              = integerMissing  ! number of soil layers                                                    (-) 
 integer(i4b)     :: nLayers            = integerMissing  ! total number of layers                                                   (-) 
 integer(i4b)     :: layerType          = integerMissing  ! index defining type of layer (snow or soil)                              (-) 
 ! number of state variables of different type
 integer(i4b)     :: nCasNrg            = integerMissing  ! number of energy state variables for the canopy air space                (-) 
 integer(i4b)     :: nVegNrg            = integerMissing  ! number of energy state variables for the vegetation canopy               (-) 
 integer(i4b)     :: nVegMass           = integerMissing  ! number of hydrology states for vegetation (mass of water)                (-) 
 integer(i4b)     :: nVegState          = integerMissing  ! number of vegetation state variables                                     (-) 
 integer(i4b)     :: nNrgState          = integerMissing  ! number of energy state variables                                         (-) 
 integer(i4b)     :: nWatState          = integerMissing  ! number of "total water" states (vol. total water content)                (-) 
 integer(i4b)     :: nMatState          = integerMissing  ! number of matric head state variables                                    (-) 
 integer(i4b)     :: nMassState         = integerMissing  ! number of hydrology state variables (mass of water)                      (-) 
 integer(i4b)     :: nState             = integerMissing  ! total number of model state variables                                    (-) 
 ! number of state variables within different domains in the snow+soil system
 integer(i4b)     :: nSnowSoilNrg       = integerMissing  ! number of energy states in the snow+soil domain                          (-) 
 integer(i4b)     :: nSnowOnlyNrg       = integerMissing  ! number of energy states in the snow domain                               (-) 
 integer(i4b)     :: nSoilOnlyNrg       = integerMissing  ! number of energy states in the soil domain                               (-) 
 integer(i4b)     :: nSnowSoilHyd       = integerMissing  ! number of hydrology states in the snow+soil domain                       (-) 
 integer(i4b)     :: nSnowOnlyHyd       = integerMissing  ! number of hydrology states in the snow domain                            (-) 
 integer(i4b)     :: nSoilOnlyHyd       = integerMissing  ! number of hydrology states in the soil domain                            (-) 
 ! type of model state variables
 integer(i4b)     :: ixControlVolume    = integerMissing  ! index of the control volume for different domains (veg, snow, soil)      (-) 
 integer(i4b)     :: ixDomainType       = integerMissing  ! index of the type of domain (iname_veg, iname_snow, iname_soil)          (-) 
 integer(i4b)     :: ixStateType        = integerMissing  ! index of the type of every state variable (iname_nrgCanair, ...)         (-) 
 integer(i4b)     :: ixHydType          = integerMissing  ! index of the type of hydrology states in snow+soil domain                (-) 
 ! type of model state variables (state subset)
 integer(i4b)     :: ixDomainType_subset= integerMissing  ! [state subset] id of domain for desired model state variables            (-) 
 integer(i4b)     :: ixStateType_subset = integerMissing  ! [state subset] type of desired model state variables                     (-) 
 ! mapping between state subset and the full state vector
 integer(i4b)     :: ixMapFull2Subset   = integerMissing  ! list of indices of the state subset in the full state vector             (-) 
 integer(i4b)     :: ixMapSubset2Full   = integerMissing  ! list of indices of the full state vector in the state subset             (-) 
 ! indices of model specific state variables
 integer(i4b)     :: ixCasNrg           = integerMissing  ! index IN THE STATE SUBSET of canopy air space energy state variable      (-) 
 integer(i4b)     :: ixVegNrg           = integerMissing  ! index IN THE STATE SUBSET of canopy energy state variable                (-) 
 integer(i4b)     :: ixVegHyd           = integerMissing  ! index IN THE STATE SUBSET of canopy hydrology state variable (mass)      (-) 
 integer(i4b)     :: ixTopNrg           = integerMissing  ! index IN THE STATE SUBSET of upper-most energy state in snow+soil domain (-) 
 integer(i4b)     :: ixTopHyd           = integerMissing  ! index IN THE STATE SUBSET of upper-most hydrol state in snow+soil domain (-) 
 ! vectors of indices for specific state types
 integer(i4b)     :: ixNrgOnly          = integerMissing  ! indices IN THE STATE SUBSET for all energy states                        (-) 
 integer(i4b)     :: ixHydOnly          = integerMissing  ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-) 
 integer(i4b)     :: ixMatOnly          = integerMissing  ! indices IN THE STATE SUBSET for matric head state variables              (-) 
 integer(i4b)     :: ixMassOnly         = integerMissing  ! indices IN THE STATE SUBSET for hydrology states (mass of water)         (-) 
 ! vectors of indices for specific state types within specific sub-domains
 integer(i4b)     :: ixSnowSoilNrg      = integerMissing  ! indices IN THE STATE SUBSET for energy states in the snow+soil domain    (-) 
 integer(i4b)     :: ixSnowOnlyNrg      = integerMissing  ! indices IN THE STATE SUBSET for energy states in the snow domain         (-) 
 integer(i4b)     :: ixSoilOnlyNrg      = integerMissing  ! indices IN THE STATE SUBSET for energy states in the soil domain         (-) 
 integer(i4b)     :: ixSnowSoilHyd      = integerMissing  ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-) 
 integer(i4b)     :: ixSnowOnlyHyd      = integerMissing  ! indices IN THE STATE SUBSET for hydrology states in the snow domain      (-) 
 integer(i4b)     :: ixSoilOnlyHyd      = integerMissing  ! indices IN THE STATE SUBSET for hydrology states in the soil domain      (-) 
 ! vectors of indices for specfic state types within specific sub-domains
 integer(i4b)     :: ixNrgCanair        = integerMissing  ! indices IN THE FULL VECTOR for energy states in canopy air space domain  (-) 
 integer(i4b)     :: ixNrgCanopy        = integerMissing  ! indices IN THE FULL VECTOR for energy states in the canopy domain        (-) 
 integer(i4b)     :: ixHydCanopy        = integerMissing  ! indices IN THE FULL VECTOR for hydrology states in the canopy domain     (-) 
 integer(i4b)     :: ixNrgLayer         = integerMissing  ! indices IN THE FULL VECTOR for energy states in the snow+soil domain     (-) 
 integer(i4b)     :: ixHydLayer         = integerMissing  ! indices IN THE FULL VECTOR for hydrology states in the snow+soil domain  (-) 
 ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
 integer(i4b)     :: ixVolFracWat       = integerMissing  ! indices IN THE SNOW+SOIL VECTOR for hyd states                           (-) 
 integer(i4b)     :: ixMatricHead       = integerMissing  ! indices IN THE SOIL VECTOR for hyd states                                (-) 
 ! indices within state vectors
 integer(i4b)     :: ixAllState         = integerMissing  ! list of indices for all model state variables                            (-) 
 integer(i4b)     :: ixSoilState        = integerMissing  ! list of indices for all soil layers                                      (-) 
 integer(i4b)     :: ixLayerState       = integerMissing  ! list of indices for all model layers                                     (-) 
 ! indices for the model output files
 integer(i4b)     :: midSnowStartIndex  = integerMissing  ! start index of the midSnow vector for a given timestep                   (-) 
 integer(i4b)     :: midSoilStartIndex  = integerMissing  ! start index of the midSoil vector for a given timestep                   (-) 
 integer(i4b)     :: midTotoStartIndex  = integerMissing  ! start index of the midToto vector for a given timestep                   (-) 
 integer(i4b)     :: ifcSnowStartIndex  = integerMissing  ! start index of the ifcSnow vector for a given timestep                   (-) 
 integer(i4b)     :: ifcSoilStartIndex  = integerMissing  ! start index of the ifcSoil vector for a given timestep                   (-) 
 integer(i4b)     :: ifcTotoStartIndex  = integerMissing  ! start index of the ifcToto vector for a given timestep                   (-) 
 endtype iLook_index

 ! ***********************************************************************************************************
 ! (11) define basin-average model parameters
 ! ***********************************************************************************************************
 type, public :: iLook_bpar
  ! baseflow
  integer(i4b)    :: basin__aquiferHydCond      = integerMissing ! hydraulic conductivity for the aquifer (m s-1)
  integer(i4b)    :: basin__aquiferScaleFactor  = integerMissing ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: basin__aquiferBaseflowExp  = integerMissing ! baseflow exponent for the big bucket (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape          = integerMissing ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale          = integerMissing ! scale parameter in Gamma distribution used for sub-grid routing (s)
 endtype iLook_bpar

 ! ***********************************************************************************************************
 ! (12) define basin-average model variables
 ! ***********************************************************************************************************
 type, public :: iLook_bvar
  ! define derived variables
  integer(i4b)    :: basin__totalArea           = integerMissing ! total basin area (m2)
  ! define fluxes
  integer(i4b)    :: basin__SurfaceRunoff       = integerMissing ! surface runoff (m s-1)
  integer(i4b)    :: basin__ColumnOutflow       = integerMissing ! outflow from all "outlet" HRUs (those with no downstream HRU)
  integer(i4b)    :: basin__AquiferStorage      = integerMissing ! aquifer storage (m s-1)
  integer(i4b)    :: basin__AquiferRecharge     = integerMissing ! recharge to the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferBaseflow     = integerMissing ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferTranspire    = integerMissing ! transpiration from the aquifer (m s-1)
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture        = integerMissing ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture      = integerMissing ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff       = integerMissing ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff        = integerMissing ! routed runoff (m s-1)
 endtype iLook_bvar

 ! ***********************************************************************************************************
 ! (10) structure for looking up the type of a model variable (this is only needed for backward 
 ! compatability, and should be removed eventually)
 ! ***********************************************************************************************************
 type, public :: iLook_varType
  integer(i4b)    :: scalarv   = integerMissing ! scalar variables 
  integer(i4b)    :: wLength   = integerMissing ! # spectral bands
  integer(i4b)    :: midSnow   = integerMissing ! mid-layer snow variables
  integer(i4b)    :: midSoil   = integerMissing ! mid-layer soil variables 
  integer(i4b)    :: midToto   = integerMissing ! mid-layer, both snow and soil
  integer(i4b)    :: ifcSnow   = integerMissing ! interface snow variables
  integer(i4b)    :: ifcSoil   = integerMissing ! interface soil variables
  integer(i4b)    :: ifcToto   = integerMissing ! interface, snow and soil
  integer(i4b)    :: parSoil   = integerMissing ! soil depth
  integer(i4b)    :: routing   = integerMissing ! routing variables
  integer(i4b)    :: outstat   = integerMissing ! output statistic
  integer(i4b)    :: unknown   = integerMissing ! cath-cal alternative type
 endtype iLook_varType

 ! ***********************************************************************************************************
 ! (11) structure for looking up statistics 
 ! ***********************************************************************************************************
 type, public :: iLook_stat
  integer(i4b)    :: totl = integerMissing ! summation 
  integer(i4b)    :: inst = integerMissing ! instantaneous 
  integer(i4b)    :: mean = integerMissing ! mean over period
  integer(i4b)    :: vari = integerMissing ! variance over period
  integer(i4b)    :: mini = integerMissing ! minimum over period 
  integer(i4b)    :: maxi = integerMissing ! maximum over period
  integer(i4b)    :: mode = integerMissing ! mode over period
 endtype iLook_stat

 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************

 ! named variables: model decisions
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39)

 ! named variables: model time
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    (  1,  2,  3,  4,  5)

 ! named variables: model forcing data
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   (  1,  2,  3,  4,  5,  6,  7,  8)

 ! named variables: model attributes
 type(iLook_attr),    public,parameter :: iLookATTR     =iLook_attr    (  1,  2,  3,  4,  5,  6,  7)

 ! named variables: soil and vegetation types
 type(iLook_type),    public,parameter :: iLookTYPE     =iLook_type    (  1,  2,  3,  4,  5)

 ! named variables: model parameters
 type(iLook_param),   public,parameter :: iLookPARAM    =iLook_param   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81, 82, 83, 84, 85, 86, 87, 88, 89, 90,&
                                                                         91, 92, 93, 94, 95, 96, 97, 98, 99,100,&
                                                                        101,102,103,104,105,106,107,108,109,110,&
                                                                        111,112,113,114,115,116,117,118,119,120,&
                                                                        121,122,123,124,125,126,127,128,129,130,&
                                                                        131,132,133,134,135,136,137,138,139,140,&
                                                                        141,142,143,144,145,146,147,148,149,150,&
                                                                        151,152,153,154)

 ! named variables: model prognostic (state) variables
 type(iLook_prog),   public,parameter  :: iLookPROG     =iLook_prog    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21)

 ! named variables: model diagnostic variables
 type(iLook_diag),    public,parameter :: iLookDIAG     =iLook_diag    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81)
 ! named variables: model fluxes
 type(iLook_flux),    public,parameter :: iLookFLUX     =iLook_flux    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81, 82, 83, 84, 85, 86)

 ! named variables: derivatives in model fluxes w.r.t. relevant state variables
 type(iLook_deriv),   public,parameter :: iLookDERIV    =iLook_deriv   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38)

 ! named variables: model indices
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58)

 ! named variables: basin-average parameters
 type(iLook_bpar),    public,parameter :: iLookBPAR     =ilook_bpar    (  1,  2,  3,  4,  5)

 ! named variables: basin-average variables
 type(iLook_bvar),    public,parameter :: iLookBVAR     =ilook_bvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11) 

 ! named variables in varibale type structure
 type(iLook_varType), public,parameter :: iLookVarType  =ilook_varType (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12)

 ! number of possible output statistics
 type(iLook_stat),    public,parameter :: iLookStat     =ilook_stat    (  1,  2,  3,  4,  5,  6,  7)

 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions = storage_size(iLookDECISIONS)/iLength
 integer(i4b),parameter,public :: maxvarTime      = storage_size(iLookTIME)/iLength
 integer(i4b),parameter,public :: maxvarForc      = storage_size(iLookFORCE)/iLength
 integer(i4b),parameter,public :: maxvarAttr      = storage_size(iLookATTR)/iLength
 integer(i4b),parameter,public :: maxvarType      = storage_size(iLookTYPE)/iLength
 integer(i4b),parameter,public :: maxvarMpar      = storage_size(iLookPARAM)/iLength
 integer(i4b),parameter,public :: maxvarProg      = storage_size(iLookPROG)/iLength
 integer(i4b),parameter,public :: maxvarDiag      = storage_size(iLookDIAG)/iLength
 integer(i4b),parameter,public :: maxvarFlux      = storage_size(iLookFLUX)/iLength
 integer(i4b),parameter,public :: maxvarDeriv     = storage_size(iLookDERIV)/iLength
 integer(i4b),parameter,public :: maxvarIndx      = storage_size(iLookINDEX)/iLength
 integer(i4b),parameter,public :: maxvarBpar      = storage_size(iLookBPAR)/iLength
 integer(i4b),parameter,public :: maxvarBvar      = storage_size(iLookBVAR)/iLength
 integer(i4b),parameter,public :: maxvarVarType   = storage_size(iLookVarType)/iLength
 integer(i4b),parameter,public :: maxvarStat      = storage_size(iLookStat)/iLength

 ! ***********************************************************************************************************
 ! (Y) define ancillary look-up structures
 ! ***********************************************************************************************************

 integer(i4b),allocatable,save,public   :: childFLUX_MEAN(:)  ! index of the child data structure: mean flux


END MODULE var_lookup
