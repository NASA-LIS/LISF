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

module get_ixname_module
! used to get the index of a named variable
USE nrtype, integerMissing=>nr_integerMissing
implicit none
private
public::get_ixdecisions
public::get_ixtime
public::get_ixattr
public::get_ixtype
public::get_ixforce
public::get_ixparam
public::get_ixprog
public::get_ixdiag
public::get_ixflux
public::get_ixderiv
public::get_ixindex
public::get_ixbpar
public::get_ixbvar
public::get_ixVarType
public::get_varTypeName
public::get_ixUnknown
public::get_statName
contains

 ! *******************************************************************************************************************
 ! public function get_ixdecisions: get the index of the named variables for the model decisions
 ! *******************************************************************************************************************
 function get_ixdecisions(varName)
 USE var_lookup,only:iLookDECISIONS                  ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixdecisions         ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('simulStart'      ); get_ixdecisions=iLookDECISIONS%simulStart  ! ( 1) simulation start time
  case('simulFinsh'      ); get_ixdecisions=iLookDECISIONS%simulFinsh  ! ( 2) simulation end time
  case('soilCatTbl'      ); get_ixdecisions=iLookDECISIONS%soilCatTbl  ! ( 3) soil-category dateset
  case('vegeParTbl'      ); get_ixdecisions=iLookDECISIONS%vegeParTbl  ! ( 4) vegetation category dataset
  case('soilStress'      ); get_ixdecisions=iLookDECISIONS%soilStress  ! ( 5) choice of function for the soil moisture control on stomatal resistance
  case('stomResist'      ); get_ixdecisions=iLookDECISIONS%stomResist  ! ( 6) choice of function for stomatal resistance
  case('bbTempFunc'      ); get_ixdecisions=iLookDECISIONS%bbTempFunc  ! ( 7) Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  case('bbHumdFunc'      ); get_ixdecisions=iLookDECISIONS%bbHumdFunc  ! ( 8) Ball-Berry: humidity controls on stomatal resistance
  case('bbElecFunc'      ); get_ixdecisions=iLookDECISIONS%bbElecFunc  ! ( 9) Ball-Berry: dependence of photosynthesis on PAR
  case('bbCO2point'      ); get_ixdecisions=iLookDECISIONS%bbCO2point  ! (10) Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  case('bbNumerics'      ); get_ixdecisions=iLookDECISIONS%bbNumerics  ! (11) Ball-Berry: iterative numerical solution method
  case('bbAssimFnc'      ); get_ixdecisions=iLookDECISIONS%bbAssimFnc  ! (12) Ball-Berry: controls on carbon assimilation
  case('bbCanIntg8'      ); get_ixdecisions=iLookDECISIONS%bbCanIntg8  ! (13) Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  case('num_method'      ); get_ixdecisions=iLookDECISIONS%num_method  ! (14) choice of numerical method
  case('fDerivMeth'      ); get_ixdecisions=iLookDECISIONS%fDerivMeth  ! (15) choice of method to calculate flux derivatives
  case('LAI_method'      ); get_ixdecisions=iLookDECISIONS%LAI_method  ! (16) choice of method to determine LAI and SAI
  case('cIntercept'      ); get_ixdecisions=iLookDECISIONS%cIntercept  ! (17) choice of parameterization for canopy interception
  case('f_Richards'      ); get_ixdecisions=iLookDECISIONS%f_Richards  ! (18) form of Richards' equation
  case('groundwatr'      ); get_ixdecisions=iLookDECISIONS%groundwatr  ! (19) choice of groundwater parameterization
  case('hc_profile'      ); get_ixdecisions=iLookDECISIONS%hc_profile  ! (20) choice of hydraulic conductivity profile
  case('bcUpprTdyn'      ); get_ixdecisions=iLookDECISIONS%bcUpprTdyn  ! (21) type of upper boundary condition for thermodynamics
  case('bcLowrTdyn'      ); get_ixdecisions=iLookDECISIONS%bcLowrTdyn  ! (22) type of lower boundary condition for thermodynamics
  case('bcUpprSoiH'      ); get_ixdecisions=iLookDECISIONS%bcUpprSoiH  ! (23) type of upper boundary condition for soil hydrology
  case('bcLowrSoiH'      ); get_ixdecisions=iLookDECISIONS%bcLowrSoiH  ! (24) type of lower boundary condition for soil hydrology
  case('veg_traits'      ); get_ixdecisions=iLookDECISIONS%veg_traits  ! (25) choice of parameterization for vegetation roughness length and displacement height
  case('rootProfil'      ); get_ixdecisions=iLookDECISIONS%rootProfil  ! (26) choice of parameterization for the rooting profile
  case('canopyEmis'      ); get_ixdecisions=iLookDECISIONS%canopyEmis  ! (27) choice of parameterization for canopy emissivity
  case('snowIncept'      ); get_ixdecisions=iLookDECISIONS%snowIncept  ! (28) choice of parameterization for snow interception
  case('windPrfile'      ); get_ixdecisions=iLookDECISIONS%windPrfile  ! (29) choice of canopy wind profile
  case('astability'      ); get_ixdecisions=iLookDECISIONS%astability  ! (30) choice of stability function
  case('compaction'      ); get_ixdecisions=iLookDECISIONS%compaction  ! (31) choice of compaction routine
  case('snowLayers'      ); get_ixdecisions=iLookDECISIONS%snowLayers  ! (32) choice of method to combine and sub-divide snow layers
  case('thCondSnow'      ); get_ixdecisions=iLookDECISIONS%thCondSnow  ! (33) choice of thermal conductivity representation for snow
  case('thCondSoil'      ); get_ixdecisions=iLookDECISIONS%thCondSoil  ! (34) choice of thermal conductivity representation for soil
  case('canopySrad'      ); get_ixdecisions=iLookDECISIONS%canopySrad  ! (35) choice of method for canopy shortwave radiation
  case('alb_method'      ); get_ixdecisions=iLookDECISIONS%alb_method  ! (36) choice of albedo representation
  case('spatial_gw'      ); get_ixdecisions=iLookDECISIONS%spatial_gw  ! (37) choice of method for spatial representation of groundwater
  case('subRouting'      ); get_ixdecisions=iLookDECISIONS%subRouting  ! (38) choice of method for sub-grid routing
  case('snowDenNew'      ); get_ixdecisions=iLookDECISIONS%snowDenNew  ! (39) choice of method for new snow density
  ! get to here if cannot find the variable
  case default
   get_ixdecisions = integerMissing
 end select
 end function get_ixdecisions


 ! *******************************************************************************************************************
 ! public function get_ixtime: get the index of the named variables for the model time
 ! *******************************************************************************************************************
 function get_ixtime(varName)
 USE var_lookup,only:iLookTIME                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixtime              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('iyyy'            ); get_ixtime = iLookTIME%iyyy             ! year
  case('im'              ); get_ixtime = iLookTIME%im               ! month
  case('id'              ); get_ixtime = iLookTIME%id               ! day
  case('ih'              ); get_ixtime = iLookTIME%ih               ! hour
  case('imin'            ); get_ixtime = iLookTIME%imin             ! minute
  ! get to here if cannot find the variable
  case default
   get_ixtime = integerMissing
 end select
 end function get_ixtime


 ! *******************************************************************************************************************
 ! public function get_ixforce: get the index of the named variables for the model forcing data
 ! *******************************************************************************************************************
 function get_ixforce(varName)
 USE var_lookup,only:iLookFORCE                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixforce             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('time'            ); get_ixforce = iLookFORCE%time             ! time since time reference       (s)
  case('pptrate'         ); get_ixforce = iLookFORCE%pptrate          ! precipitation rate              (kg m-2 s-1)
  case('airtemp'         ); get_ixforce = iLookFORCE%airtemp          ! air temperature                 (K)
  case('spechum'         ); get_ixforce = iLookFORCE%spechum          ! specific humidity               (g/g)
  case('windspd'         ); get_ixforce = iLookFORCE%windspd          ! windspeed                       (m/s)
  case('SWRadAtm'        ); get_ixforce = iLookFORCE%SWRadAtm         ! downwelling shortwave radiaiton (W m-2)
  case('LWRadAtm'        ); get_ixforce = iLookFORCE%LWRadAtm         ! downwelling longwave radiation  (W m-2)
  case('airpres'         ); get_ixforce = iLookFORCE%airpres          ! pressure                        (Pa)
  ! get to here if cannot find the variable
  case default
   get_ixforce = integerMissing
 end select
 end function get_ixforce


 ! *******************************************************************************************************************
 ! public function get_ixAttr: get the index of the named variables for the site characteristics
 ! *******************************************************************************************************************
 function get_ixAttr(varName)
 USE var_lookup,only:iLookATTR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixAttr              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('latitude'      ); get_ixAttr = iLookATTR%latitude       ! latitude (degrees north)
  case('longitude'     ); get_ixAttr = iLookATTR%longitude      ! longitude (degrees east)
  case('elevation'     ); get_ixAttr = iLookATTR%elevation      ! elevation (m)
  case('tan_slope'     ); get_ixAttr = iLookATTR%tan_slope      ! tan water table slope, taken as tan local ground surface slope (-)
  case('contourLength' ); get_ixAttr = iLookATTR%contourLength  ! length of contour at downslope edge of HRU (m)
  case('HRUarea'       ); get_ixAttr = iLookATTR%HRUarea        ! area of each HRU (m2)
  case('mHeight'       ); get_ixAttr = iLookATTR%mHeight        ! measurement height above bare ground (m)
  ! get to here if cannot find the variable
  case default
   get_ixAttr = integerMissing
 end select
 end function get_ixAttr


 ! *******************************************************************************************************************
 ! public function get_ixType: get the index of the named variables for the local classification of veg, soil, etc.
 ! *******************************************************************************************************************
 function get_ixType(varName)
 USE var_lookup,only:iLookTYPE                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixType              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('hruIndex'       ); get_ixType = iLookTYPE%hruIndex           ! index defining HRU index
  case('vegTypeIndex'   ); get_ixType = iLookTYPE%vegTypeIndex       ! index defining vegetation type
  case('soilTypeIndex'  ); get_ixType = iLookTYPE%soilTypeIndex      ! index defining soil type
  case('slopeTypeIndex' ); get_ixType = iLookTYPE%slopeTypeIndex     ! index defining slope
  case('downHRUindex'   ); get_ixType = iLookTYPE%downHRUindex       ! index of downslope HRU (0 = basin outlet)
  ! get to here if cannot find the variable
  case default
   get_ixType = integerMissing
 end select
 end function get_ixType


 ! *******************************************************************************************************************
 ! public function get_ixparam: get the index of the named variables for the model parameters
 ! *******************************************************************************************************************
 function get_ixparam(varName)
 USE var_lookup,only:iLookPARAM                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixparam             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! boundary conditions
  case('upperBoundHead'           ); get_ixparam = iLookPARAM%upperBoundHead         ! matric head of the upper boundary (m)
  case('lowerBoundHead'           ); get_ixparam = iLookPARAM%lowerBoundHead         ! matric head of the lower boundary (m)
  case('upperBoundTheta'          ); get_ixparam = iLookPARAM%upperBoundTheta        ! volumetric liquid water content at the upper boundary (-)
  case('lowerBoundTheta'          ); get_ixparam = iLookPARAM%lowerBoundTheta        ! volumetric liquid water content at the lower boundary (-)
  case('upperBoundTemp'           ); get_ixparam = iLookPARAM%upperBoundTemp         ! temperature of the upper boundary (K)
  case('lowerBoundTemp'           ); get_ixparam = iLookPARAM%lowerBoundTemp         ! temperature of the lower boundary (K)
  ! precipitation partitioning
  case('tempCritRain'             ); get_ixparam = iLookPARAM%tempCritRain           ! critical temperature where precipitation is rain (K)
  case('tempRangeTimestep'        ); get_ixparam = iLookPARAM%tempRangeTimestep      ! temperature range over the time step (K)
  case('frozenPrecipMultip'       ); get_ixparam = iLookPARAM%frozenPrecipMultip     ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  case('snowfrz_scale'            ); get_ixparam = iLookPARAM%snowfrz_scale          ! scaling parameter for the freezing curve for snow (K-1)
  case('fixedThermalCond_snow'    ); get_ixparam = iLookPARAM%fixedThermalCond_snow  ! temporally constant thermal conductivity for snow (W m-1 K-1)
  ! snow albedo
  case('albedoMax'                ); get_ixparam = iLookPARAM%albedoMax              ! maximum snow albedo for a single spectral band (-)
  case('albedoMinWinter'          ); get_ixparam = iLookPARAM%albedoMinWinter        ! minimum snow albedo during winter for a single spectral band (-)
  case('albedoMinSpring'          ); get_ixparam = iLookPARAM%albedoMinSpring        ! minimum snow albedo during spring for a single spectral band (-)
  case('albedoMaxVisible'         ); get_ixparam = iLookPARAM%albedoMaxVisible       ! maximum snow albedo in the visible part of the spectrum (-)
  case('albedoMinVisible'         ); get_ixparam = iLookPARAM%albedoMinVisible       ! minimum snow albedo in the visible part of the spectrum (-)
  case('albedoMaxNearIR'          ); get_ixparam = iLookPARAM%albedoMaxNearIR        ! maximum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoMinNearIR'          ); get_ixparam = iLookPARAM%albedoMinNearIR        ! minimum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoDecayRate'          ); get_ixparam = iLookPARAM%albedoDecayRate        ! albedo decay rate (s)
  case('albedoSootLoad'           ); get_ixparam = iLookPARAM%albedoSootLoad         ! soot load factor (-)
  case('albedoRefresh'            ); get_ixparam = iLookPARAM%albedoRefresh          ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer
  case('radExt_snow'              ); get_ixparam = iLookPARAM%radExt_snow            ! extinction coefficient for radiation penetration within the snowpack (m-1)
  case('directScale'              ); get_ixparam = iLookPARAM%directScale            ! scaling factor for fractional driect radiaion parameterization (-)
  case('Frad_direct'              ); get_ixparam = iLookPARAM%Frad_direct            ! maximum fraction of direct radiation (-)
  case('Frad_vis'                 ); get_ixparam = iLookPARAM%Frad_vis               ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  case('newSnowDenMin'            ); get_ixparam = iLookPARAM%newSnowDenMin          ! minimum new snow density (kg m-3)
  case('newSnowDenMult'           ); get_ixparam = iLookPARAM%newSnowDenMult         ! multiplier for new snow density (kg m-3)
  case('newSnowDenScal'           ); get_ixparam = iLookPARAM%newSnowDenScal         ! scaling factor for new snow density (K)
  case('constSnowDen'             ); get_ixparam = iLookPARAM%constSnowDen           ! Constant new snow density (kg m-3)
  case('newSnowDenAdd'            ); get_ixparam = iLookPARAM%newSnowDenAdd          ! Pahaut 1976, additive factor for new snow density (kg m-3)
  case('newSnowDenMultTemp'       ); get_ixparam = iLookPARAM%newSnowDenMultTemp     ! Pahaut 1976, multiplier for new snow density applied to air temperature (kg m-3 K-1)
  case('newSnowDenMultWind'       ); get_ixparam = iLookPARAM%newSnowDenMultWind     ! Pahaut 1976, multiplier for new snow density applied to wind speed (kg m-7/2 s-1/2)
  case('newSnowDenMultAnd'        ); get_ixparam = iLookPARAM%newSnowDenMultAnd      ! Anderson 1976, multiplier for new snow density for Anderson function (K-1)
  case('newSnowDenBase'           ); get_ixparam = iLookPARAM%newSnowDenBase         ! Anderson 1976, base value that is rasied to the (3/2) power (K)
  ! snow compaction
  case('densScalGrowth'           ); get_ixparam = iLookPARAM%densScalGrowth         ! density scaling factor for grain growth (kg-1 m3)
  case('tempScalGrowth'           ); get_ixparam = iLookPARAM%tempScalGrowth         ! temperature scaling factor for grain growth (K-1)
  case('grainGrowthRate'          ); get_ixparam = iLookPARAM%grainGrowthRate        ! rate of grain growth (s-1)
  case('densScalOvrbdn'           ); get_ixparam = iLookPARAM%densScalOvrbdn         ! density scaling factor for overburden pressure (kg-1 m3)
  case('tempScalOvrbdn'           ); get_ixparam = iLookPARAM%tempScalOvrbdn         ! temperature scaling factor for overburden pressure (K-1)
  case('baseViscosity'            ); get_ixparam = iLookPARAM%baseViscosity          ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow through snow
  case('Fcapil'                   ); get_ixparam = iLookPARAM%Fcapil                 ! capillary retention as a fraction of the total pore volume (-)
  case('k_snow'                   ); get_ixparam = iLookPARAM%k_snow                 ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  case('mw_exp'                   ); get_ixparam = iLookPARAM%mw_exp                 ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  case('z0Snow'                   ); get_ixparam = iLookPARAM%z0Snow                 ! roughness length of snow (m)
  case('z0Soil'                   ); get_ixparam = iLookPARAM%z0Soil                 ! roughness length of bare soil below the canopy (m)
  case('z0Canopy'                 ); get_ixparam = iLookPARAM%z0Canopy               ! roughness length of the canopy (m)
  case('zpdFraction'              ); get_ixparam = iLookPARAM%zpdFraction            ! zero plane displacement / canopy height (-)
  case('critRichNumber'           ); get_ixparam = iLookPARAM%critRichNumber         ! critical value for the bulk Richardson number (-)
  case('Louis79_bparam'           ); get_ixparam = iLookPARAM%Louis79_bparam         ! parameter in Louis (1979) stability function (-)
  case('Louis79_cStar'            ); get_ixparam = iLookPARAM%Louis79_cStar          ! parameter in Louis (1979) stability function (-)
  case('Mahrt87_eScale'           ); get_ixparam = iLookPARAM%Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function (-)
  case('leafExchangeCoeff'        ); get_ixparam = iLookPARAM%leafExchangeCoeff      ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  case('windReductionParam'       ); get_ixparam = iLookPARAM%windReductionParam     ! canopy wind reduction parameter (-)
  ! stomatal conductance
  case('Kc25'                     ); get_ixparam = iLookPARAM%Kc25                   ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  case('Ko25'                     ); get_ixparam = iLookPARAM%Ko25                   ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  case('Kc_qFac'                  ); get_ixparam = iLookPARAM%Kc_qFac                ! factor in the q10 function defining temperature controls on Kc (-)
  case('Ko_qFac'                  ); get_ixparam = iLookPARAM%Ko_qFac                ! factor in the q10 function defining temperature controls on Ko (-)
  case('kc_Ha'                    ); get_ixparam = iLookPARAM%kc_Ha                  ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  case('ko_Ha'                    ); get_ixparam = iLookPARAM%ko_Ha                  ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  case('vcmax25_canopyTop'        ); get_ixparam = iLookPARAM%vcmax25_canopyTop      ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  case('vcmax_qFac'               ); get_ixparam = iLookPARAM%vcmax_qFac             ! factor in the q10 function defining temperature controls on vcmax (-)
  case('vcmax_Ha'                 ); get_ixparam = iLookPARAM%vcmax_Ha               ! activation energy in the vcmax function (J mol-1)
  case('vcmax_Hd'                 ); get_ixparam = iLookPARAM%vcmax_Hd               ! deactivation energy in the vcmax function (J mol-1)
  case('vcmax_Sv'                 ); get_ixparam = iLookPARAM%vcmax_Sv               ! entropy term in the vcmax function (J mol-1 K-1)
  case('vcmax_Kn'                 ); get_ixparam = iLookPARAM%vcmax_Kn               ! foliage nitrogen decay coefficient (-)
  case('jmax25_scale'             ); get_ixparam = iLookPARAM%jmax25_scale           ! scaling factor to relate jmax25 to vcmax25 (-)
  case('jmax_Ha'                  ); get_ixparam = iLookPARAM%jmax_Ha                ! activation energy in the jmax function (J mol-1)
  case('jmax_Hd'                  ); get_ixparam = iLookPARAM%jmax_Hd                ! deactivation energy in the jmax function (J mol-1)
  case('jmax_Sv'                  ); get_ixparam = iLookPARAM%jmax_Sv                ! entropy term in the jmax function (J mol-1 K-1)
  case('fractionJ'                ); get_ixparam = iLookPARAM%fractionJ              ! fraction of light lost by other than the chloroplast lamellae (-)
  case('quantamYield'             ); get_ixparam = iLookPARAM%quantamYield           ! quantam yield (mol e mol-1 quanta)
  case('vpScaleFactor'            ); get_ixparam = iLookPARAM%vpScaleFactor          ! vapor pressure scaling factor in stomatal conductance function (Pa)
  case('cond2photo_slope'         ); get_ixparam = iLookPARAM%cond2photo_slope       ! slope of conductance-photosynthesis relationship (-)
  case('minStomatalConductance'   ); get_ixparam = iLookPARAM%minStomatalConductance ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  case('winterSAI'                ); get_ixparam = iLookPARAM%winterSAI              ! stem area index prior to the start of the growing season (m2 m-2)
  case('summerLAI'                ); get_ixparam = iLookPARAM%summerLAI              ! maximum leaf area index at the peak of the growing season (m2 m-2)
  case('rootScaleFactor1'         ); get_ixparam = iLookPARAM%rootScaleFactor1       ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )   (m-1)
  case('rootScaleFactor2'         ); get_ixparam = iLookPARAM%rootScaleFactor2       ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )   (m-1)
  case('rootingDepth'             ); get_ixparam = iLookPARAM%rootingDepth           ! rooting depth (m)
  case('rootDistExp'              ); get_ixparam = iLookPARAM%rootDistExp            ! exponent for the vertical distriution of root density (-)
  case('plantWiltPsi'             ); get_ixparam = iLookPARAM%plantWiltPsi           ! matric head at wilting point (m)
  case('soilStressParam'          ); get_ixparam = iLookPARAM%soilStressParam        ! parameter in the exponential soil stress function
  case('critSoilWilting'          ); get_ixparam = iLookPARAM%critSoilWilting        ! critical vol. liq. water content when plants are wilting (-)
  case('critSoilTranspire'        ); get_ixparam = iLookPARAM%critSoilTranspire      ! critical vol. liq. water content when transpiration is limited (-)
  case('critAquiferTranspire'     ); get_ixparam = iLookPARAM%critAquiferTranspire   ! critical aquifer storage value when transpiration is limited (m)
  case('minStomatalResistance'    ); get_ixparam = iLookPARAM%minStomatalResistance  ! minimum canopy resistance (s m-1)
  case('leafDimension'            ); get_ixparam = iLookPARAM%leafDimension          ! characteristic leaf dimension (m)
  case('heightCanopyTop'          ); get_ixparam = iLookPARAM%heightCanopyTop        ! height of top of the vegetation canopy above ground surface (m)
  case('heightCanopyBottom'       ); get_ixparam = iLookPARAM%heightCanopyBottom     ! height of bottom of the vegetation canopy above ground surface (m)
  case('specificHeatVeg'          ); get_ixparam = iLookPARAM%specificHeatVeg        ! specific heat of vegetation (J kg-1 K-1)
  case('maxMassVegetation'        ); get_ixparam = iLookPARAM%maxMassVegetation      ! maximum mass of vegetation (full foliage) (kg m-2)
  case('throughfallScaleSnow'     ); get_ixparam = iLookPARAM%throughfallScaleSnow   ! scaling factor for throughfall (snow) (-)
  case('throughfallScaleRain'     ); get_ixparam = iLookPARAM%throughfallScaleRain   ! scaling factor for throughfall (rain) (-)
  case('refInterceptCapSnow'      ); get_ixparam = iLookPARAM%refInterceptCapSnow    ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  case('refInterceptCapRain'      ); get_ixparam = iLookPARAM%refInterceptCapRain    ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  case('snowUnloadingCoeff'       ); get_ixparam = iLookPARAM%snowUnloadingCoeff     ! time constant for unloading of snow from the forest canopy (s-1)
  case('canopyDrainageCoeff'      ); get_ixparam = iLookPARAM%canopyDrainageCoeff    ! time constant for drainage of liquid water from the forest canopy (s-1)
  case('ratioDrip2Unloading'      ); get_ixparam = iLookPARAM%ratioDrip2Unloading    ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  case('canopyWettingFactor'      ); get_ixparam = iLookPARAM%canopyWettingFactor    ! maximum wetted fraction of the canopy (-)
  case('canopyWettingExp'         ); get_ixparam = iLookPARAM%canopyWettingExp       ! exponent in canopy wetting function (-)
  ! soil properties
  case('soil_dens_intr'           ); get_ixparam = iLookPARAM%soil_dens_intr         ! intrinsic soil density (kg m-3)
  case('thCond_soil'              ); get_ixparam = iLookPARAM%thCond_soil            ! thermal conductivity of soil (W m-1 K-1)
  case('frac_sand'                ); get_ixparam = iLookPARAM%frac_sand              ! fraction of sand (-)
  case('frac_silt'                ); get_ixparam = iLookPARAM%frac_silt              ! fraction of silt (-)
  case('frac_clay'                ); get_ixparam = iLookPARAM%frac_clay              ! fraction of clay (-)
  case('fieldCapacity'            ); get_ixparam = iLookPARAM%fieldCapacity          ! field capacity (-)
  case('wettingFrontSuction'      ); get_ixparam = iLookPARAM%wettingFrontSuction    ! Green-Ampt wetting front suction (m)
  case('theta_mp'                 ); get_ixparam = iLookPARAM%theta_mp               ! volumetric liquid water content when macropore flow begins (-)
  case('theta_sat'                ); get_ixparam = iLookPARAM%theta_sat              ! soil porosity (-)
  case('theta_res'                ); get_ixparam = iLookPARAM%theta_res              ! volumetric residual water content (-)
  case('vGn_alpha'                ); get_ixparam = iLookPARAM%vGn_alpha              ! van Genuchten "alpha" parameter (m-1)
  case('vGn_n'                    ); get_ixparam = iLookPARAM%vGn_n                  ! van Genuchten "n" parameter (-)
  case('mpExp'                    ); get_ixparam = iLookPARAM%mpExp                  ! empirical exponent in macropore flow equation (-)
  case('k_soil'                   ); get_ixparam = iLookPARAM%k_soil                 ! saturated hydraulic conductivity (m s-1)
  case('k_macropore'              ); get_ixparam = iLookPARAM%k_macropore            ! saturated hydraulic conductivity for the macropores (m s-1)
  case('kAnisotropic'             ); get_ixparam = iLookPARAM%kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
  case('zScale_TOPMODEL'          ); get_ixparam = iLookPARAM%zScale_TOPMODEL        ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  case('compactedDepth'           ); get_ixparam = iLookPARAM%compactedDepth         ! depth where k_soil reaches the compacted value given by CH78 (m)
  case('aquiferScaleFactor'       ); get_ixparam = iLookPARAM%aquiferScaleFactor     ! scaling factor for aquifer storage in the big bucket (m)
  case('aquiferBaseflowExp'       ); get_ixparam = iLookPARAM%aquiferBaseflowExp     ! baseflow exponent (-)
  case('qSurfScale'               ); get_ixparam = iLookPARAM%qSurfScale             ! scaling factor in the surface runoff parameterization (-)
  case('specificYield'            ); get_ixparam = iLookPARAM%specificYield          ! specific yield (-)
  case('specificStorage'          ); get_ixparam = iLookPARAM%specificStorage        ! specific storage coefficient (m-1)
  case('f_impede'                 ); get_ixparam = iLookPARAM%f_impede               ! ice impedence factor (-)
  case('soilIceScale'             ); get_ixparam = iLookPARAM%soilIceScale           ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  case('soilIceCV'                ); get_ixparam = iLookPARAM%soilIceCV              ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  case('minwind'                  ); get_ixparam = iLookPARAM%minwind                ! minimum wind speed (m s-1)
  case('minstep'                  ); get_ixparam = iLookPARAM%minstep                ! minimum length of the time step
  case('maxstep'                  ); get_ixparam = iLookPARAM%maxstep                ! maximum length of the time step
  case('wimplicit'                ); get_ixparam = iLookPARAM%wimplicit              ! weight assigned to start-of-step fluxes
  case('maxiter'                  ); get_ixparam = iLookPARAM%maxiter                ! maximum number of iterations
  case('relConvTol_liquid'        ); get_ixparam = iLookPARAM%relConvTol_liquid      ! relative convergence tolerance for vol frac liq water (-)
  case('absConvTol_liquid'        ); get_ixparam = iLookPARAM%absConvTol_liquid      ! absolute convergence tolerance for vol frac liq water (-)
  case('relConvTol_matric'        ); get_ixparam = iLookPARAM%relConvTol_matric      ! relative convergence tolerance for matric head (-)
  case('absConvTol_matric'        ); get_ixparam = iLookPARAM%absConvTol_matric      ! absolute convergence tolerance for matric head (m)
  case('relConvTol_energy'        ); get_ixparam = iLookPARAM%relConvTol_energy      ! relative convergence tolerance for energy (-)
  case('absConvTol_energy'        ); get_ixparam = iLookPARAM%absConvTol_energy      ! absolute convergence tolerance for energy (J m-3)
  case('relConvTol_aquifr'        ); get_ixparam = iLookPARAM%relConvTol_aquifr      ! relative convergence tolerance for aquifer storage (-)
  case('absConvTol_aquifr'        ); get_ixparam = iLookPARAM%absConvTol_aquifr      ! absolute convergence tolerance for aquifer storage (m)
  case('zmin'                     ); get_ixparam = iLookPARAM%zmin                   ! minimum layer depth (m)
  case('zmax'                     ); get_ixparam = iLookPARAM%zmax                   ! maximum layer depth (m)
  case('zminLayer1'               ); get_ixparam = iLookPARAM%zminLayer1             ! minimum layer depth for the 1st (top) layer (m)
  case('zminLayer2'               ); get_ixparam = iLookPARAM%zminLayer2             ! minimum layer depth for the 2nd layer (m)
  case('zminLayer3'               ); get_ixparam = iLookPARAM%zminLayer3             ! minimum layer depth for the 3rd layer (m)
  case('zminLayer4'               ); get_ixparam = iLookPARAM%zminLayer4             ! minimum layer depth for the 4th layer (m)
  case('zminLayer5'               ); get_ixparam = iLookPARAM%zminLayer5             ! minimum layer depth for the 5th (bottom) layer (m)
  case('zmaxLayer1_lower'         ); get_ixparam = iLookPARAM%zmaxLayer1_lower       ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  case('zmaxLayer2_lower'         ); get_ixparam = iLookPARAM%zmaxLayer2_lower       ! maximum layer depth for the 2nd layer when only 2 layers (m)
  case('zmaxLayer3_lower'         ); get_ixparam = iLookPARAM%zmaxLayer3_lower       ! maximum layer depth for the 3rd layer when only 3 layers (m)
  case('zmaxLayer4_lower'         ); get_ixparam = iLookPARAM%zmaxLayer4_lower       ! maximum layer depth for the 4th layer when only 4 layers (m)
  case('zmaxLayer1_upper'         ); get_ixparam = iLookPARAM%zmaxLayer1_upper       ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  case('zmaxLayer2_upper'         ); get_ixparam = iLookPARAM%zmaxLayer2_upper       ! maximum layer depth for the 2nd layer when > 2 layers (m)
  case('zmaxLayer3_upper'         ); get_ixparam = iLookPARAM%zmaxLayer3_upper       ! maximum layer depth for the 3rd layer when > 3 layers (m)
  case('zmaxLayer4_upper'         ); get_ixparam = iLookPARAM%zmaxLayer4_upper       ! maximum layer depth for the 4th layer when > 4 layers (m)
  ! get to here if cannot find the variable
  case default
   get_ixparam = integerMissing
 end select
 end function get_ixparam


 ! *******************************************************************************************************************
 ! public function get_ixprog: get the index of the named variables for the prognostic (state) variables
 ! *******************************************************************************************************************
 function get_ixprog(varName)
 USE var_lookup,only:iLookPROG                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixprog              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! variables for time stepping
  case('dt_init'                        ); get_ixprog = iLookPROG%dt_init                          ! length of initial time step at start of next data interval (s)
  ! state variables for vegetation
  case('scalarCanopyIce'                ); get_ixprog = iLookPROG%scalarCanopyIce                  ! mass of ice on the vegetation canopy (kg m-2)
  case('scalarCanopyLiq'                ); get_ixprog = iLookPROG%scalarCanopyLiq                  ! mass of liquid water on the vegetation canopy (kg m-2)
  case('scalarCanopyWat'                ); get_ixprog = iLookPROG%scalarCanopyWat                  ! mass of total water on the vegetation canopy (kg m-2)
  case('scalarCanairTemp'               ); get_ixprog = iLookPROG%scalarCanairTemp                 ! temperature of the canopy air space (K)
  case('scalarCanopyTemp'               ); get_ixprog = iLookPROG%scalarCanopyTemp                 ! temperature of the vegetation canopy (K)
  ! state variables for snow
  case('spectralSnowAlbedoDiffuse'      ); get_ixprog = iLookPROG%spectralSnowAlbedoDiffuse        ! diffuse snow albedo for individual spectral bands (-)
  case('scalarSnowAlbedo'               ); get_ixprog = iLookPROG%scalarSnowAlbedo                 ! snow albedo for the entire spectral band (-)
  case('scalarSnowDepth'                ); get_ixprog = iLookPROG%scalarSnowDepth                  ! total snow depth (m)
  case('scalarSWE'                      ); get_ixprog = iLookPROG%scalarSWE                        ! snow water equivalent (kg m-2)
  case('scalarSfcMeltPond'              ); get_ixprog = iLookPROG%scalarSfcMeltPond                ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  ! state variables for the snow+soil domain
  case('mLayerTemp'                     ); get_ixprog = iLookPROG%mLayerTemp                       ! temperature of each layer (K)
  case('mLayerVolFracIce'               ); get_ixprog = iLookPROG%mLayerVolFracIce                 ! volumetric fraction of icein each layer (-)
  case('mLayerVolFracLiq'               ); get_ixprog = iLookPROG%mLayerVolFracLiq                 ! volumetric fraction of liquid water in each layer (-)
  case('mLayerVolFracWat'               ); get_ixprog = iLookPROG%mLayerVolFracWat                 ! volumetric fraction of total water in each layer (-)
  case('mLayerMatricHead'               ); get_ixprog = iLookPROG%mLayerMatricHead                 ! matric head of water in the soil (m)
  ! other state variables
  case('scalarAquiferStorage'           ); get_ixprog = iLookPROG%scalarAquiferStorage             ! relative aquifer storage -- above bottom of the soil profile (m)
  case('scalarSurfaceTemp'              ); get_ixprog = iLookPROG%scalarSurfaceTemp                ! surface temperature (K)
  ! coordinate variables
  case('mLayerDepth'                    ); get_ixprog = iLookPROG%mLayerDepth                      ! depth of each layer (m)
  case('mLayerHeight'                   ); get_ixprog = iLookPROG%mLayerHeight                     ! height at the midpoint of each layer (m)
  case('iLayerHeight'                   ); get_ixprog = iLookPROG%iLayerHeight                     ! height at the interface of each layer (m)
  ! get to here if cannot find the variable
  case default
   get_ixprog = integerMissing
 end select
 end function get_ixprog


 ! *******************************************************************************************************************
 ! public function get_ixdiag: get the index of the named variables for the diagnostic variables
 ! *******************************************************************************************************************
 function get_ixdiag(varName)
 USE var_lookup,only:iLookDIAG                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixdiag              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! local properties
  case('scalarCanopyDepth'              ); get_ixdiag = iLookDIAG%scalarCanopyDepth                ! canopy depth (m)
  case('scalarGreenVegFraction'         ); get_ixdiag = iLookDIAG%scalarGreenVegFraction           ! green vegetation fraction used to compute LAI (-)
  case('scalarBulkVolHeatCapVeg'        ); get_ixdiag = iLookDIAG%scalarBulkVolHeatCapVeg          ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  case('scalarCanopyEmissivity'         ); get_ixdiag = iLookDIAG%scalarCanopyEmissivity           ! effective canopy emissivity (-)
  case('scalarRootZoneTemp'             ); get_ixdiag = iLookDIAG%scalarRootZoneTemp               ! average temperature of the root zone (K)
  case('scalarLAI'                      ); get_ixdiag = iLookDIAG%scalarLAI                        ! one-sided leaf area index (m2 m-2)
  case('scalarSAI'                      ); get_ixdiag = iLookDIAG%scalarSAI                        ! one-sided stem area index (m2 m-2)
  case('scalarExposedLAI'               ); get_ixdiag = iLookDIAG%scalarExposedLAI                 ! exposed leaf area index after burial by snow (m2 m-2)
  case('scalarExposedSAI'               ); get_ixdiag = iLookDIAG%scalarExposedSAI                 ! exposed stem area index after burial by snow (m2 m-2)
  case('scalarCanopyIceMax'             ); get_ixdiag = iLookDIAG%scalarCanopyIceMax               ! maximum interception storage capacity for ice (kg m-2)
  case('scalarCanopyLiqMax'             ); get_ixdiag = iLookDIAG%scalarCanopyLiqMax               ! maximum interception storage capacity for liquid water (kg m-2)
  case('scalarGrowingSeasonIndex'       ); get_ixdiag = iLookDIAG%scalarGrowingSeasonIndex         ! growing season index (0=off, 1=on)
  case('scalarVolHtCap_air'             ); get_ixdiag = iLookDIAG%scalarVolHtCap_air               ! volumetric heat capacity air (J m-3 K-1)
  case('scalarVolHtCap_ice'             ); get_ixdiag = iLookDIAG%scalarVolHtCap_ice               ! volumetric heat capacity ice (J m-3 K-1)
  case('scalarVolHtCap_soil'            ); get_ixdiag = iLookDIAG%scalarVolHtCap_soil              ! volumetric heat capacity dry soil (J m-3 K-1)
  case('scalarVolHtCap_water'           ); get_ixdiag = iLookDIAG%scalarVolHtCap_water             ! volumetric heat capacity liquid wat (J m-3 K-1)
  case('mLayerVolHtCapBulk'             ); get_ixdiag = iLookDIAG%mLayerVolHtCapBulk               ! volumetric heat capacity in each layer (J m-3 K-1)
  case('scalarLambda_drysoil'           ); get_ixdiag = iLookDIAG%scalarLambda_drysoil             ! thermal conductivity of dry soil     (W m-1)
  case('scalarLambda_wetsoil'           ); get_ixdiag = iLookDIAG%scalarLambda_wetsoil             ! thermal conductivity of wet soil     (W m-1)
  case('mLayerThermalC'                 ); get_ixdiag = iLookDIAG%mLayerThermalC                   ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  case('iLayerThermalC'                 ); get_ixdiag = iLookDIAG%iLayerThermalC                   ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! forcing
  case('scalarVPair'                    ); get_ixdiag = iLookDIAG%scalarVPair                      ! vapor pressure of the air above the vegetation canopy (Pa)
  case('scalarVP_CanopyAir'             ); get_ixdiag = iLookDIAG%scalarVP_CanopyAir               ! vapor pressure of the canopy air space (Pa)
  case('scalarTwetbulb'                 ); get_ixdiag = iLookDIAG%scalarTwetbulb                   ! wetbulb temperature (K)
  case('scalarSnowfallTemp'             ); get_ixdiag = iLookDIAG%scalarSnowfallTemp               ! temperature of fresh snow (K)
  case('scalarNewSnowDensity'           ); get_ixdiag = iLookDIAG%scalarNewSnowDensity             ! density of fresh snow, should snow be falling in this time step (kg m-3)
  case('scalarO2air'                    ); get_ixdiag = iLookDIAG%scalarO2air                      ! atmospheric o2 concentration (Pa)
  case('scalarCO2air'                   ); get_ixdiag = iLookDIAG%scalarCO2air                     ! atmospheric co2 concentration (Pa)
  ! shortwave radiation
  case('scalarCosZenith'                ); get_ixdiag = iLookDIAG%scalarCosZenith                  ! cosine of the solar zenith angle (0-1)
  case('scalarFractionDirect'           ); get_ixdiag = iLookDIAG%scalarFractionDirect             ! fraction of direct radiation (0-1)
  case('scalarCanopySunlitFraction'     ); get_ixdiag = iLookDIAG%scalarCanopySunlitFraction       ! sunlit fraction of canopy (-)
  case('scalarCanopySunlitLAI'          ); get_ixdiag = iLookDIAG%scalarCanopySunlitLAI            ! sunlit leaf area (-)
  case('scalarCanopyShadedLAI'          ); get_ixdiag = iLookDIAG%scalarCanopyShadedLAI            ! shaded leaf area (-)
  case('spectralAlbGndDirect'           ); get_ixdiag = iLookDIAG%spectralAlbGndDirect             ! direct  albedo of underlying surface for each spectral band (-)
  case('spectralAlbGndDiffuse'          ); get_ixdiag = iLookDIAG%spectralAlbGndDiffuse            ! diffuse albedo of underlying surface for each spectral band (-)
  case('scalarGroundAlbedo'             ); get_ixdiag = iLookDIAG%scalarGroundAlbedo               ! albedo of the ground surface (-)
  ! turbulent heat transfer
  case('scalarLatHeatSubVapCanopy'      ); get_ixdiag = iLookDIAG%scalarLatHeatSubVapCanopy        ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  case('scalarLatHeatSubVapGround'      ); get_ixdiag = iLookDIAG%scalarLatHeatSubVapGround        ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  case('scalarSatVP_CanopyTemp'         ); get_ixdiag = iLookDIAG%scalarSatVP_CanopyTemp           ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  case('scalarSatVP_GroundTemp'         ); get_ixdiag = iLookDIAG%scalarSatVP_GroundTemp           ! saturation vapor pressure at the temperature of the ground (Pa)
  case('scalarZ0Canopy'                 ); get_ixdiag = iLookDIAG%scalarZ0Canopy                   ! roughness length of the canopy (m)
  case('scalarWindReductionFactor'      ); get_ixdiag = iLookDIAG%scalarWindReductionFactor        ! canopy wind reduction factor (-)
  case('scalarZeroPlaneDisplacement'    ); get_ixdiag = iLookDIAG%scalarZeroPlaneDisplacement      ! zero plane displacement (m)
  case('scalarRiBulkCanopy'             ); get_ixdiag = iLookDIAG%scalarRiBulkCanopy               ! bulk Richardson number for the canopy (-)
  case('scalarRiBulkGround'             ); get_ixdiag = iLookDIAG%scalarRiBulkGround               ! bulk Richardson number for the ground surface (-)
  case('scalarCanopyStabilityCorrection'); get_ixdiag = iLookDIAG%scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
  case('scalarGroundStabilityCorrection'); get_ixdiag = iLookDIAG%scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
  ! evapotranspiration 
  case('scalarIntercellularCO2Sunlit'   ); get_ixdiag = iLookDIAG%scalarIntercellularCO2Sunlit     ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  case('scalarIntercellularCO2Shaded'   ); get_ixdiag = iLookDIAG%scalarIntercellularCO2Shaded     ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  case('scalarTranspireLim'             ); get_ixdiag = iLookDIAG%scalarTranspireLim               ! aggregate soil moisture and aquifer storage limit on transpiration (-)
  case('scalarTranspireLimAqfr'         ); get_ixdiag = iLookDIAG%scalarTranspireLimAqfr           ! aquifer storage limit on transpiration (-)
  case('scalarFoliageNitrogenFactor'    ); get_ixdiag = iLookDIAG%scalarFoliageNitrogenFactor      ! foliage nitrogen concentration, 1=saturated (-)
  case('scalarSoilRelHumidity'          ); get_ixdiag = iLookDIAG%scalarSoilRelHumidity            ! relative humidity in the soil pores in the upper-most soil layer (-)
  case('mLayerTranspireLim'             ); get_ixdiag = iLookDIAG%mLayerTranspireLim               ! moisture avail factor limiting transpiration in each layer (-)
  case('mLayerRootDensity'              ); get_ixdiag = iLookDIAG%mLayerRootDensity                ! fraction of roots in each soil layer (-)
  case('scalarAquiferRootFrac'          ); get_ixdiag = iLookDIAG%scalarAquiferRootFrac            ! fraction of roots below the soil profile (-)
  ! canopy hydrology
  case('scalarFracLiqVeg'               ); get_ixdiag = iLookDIAG%scalarFracLiqVeg                 ! fraction of liquid water on vegetation (-) 
  case('scalarCanopyWetFraction'        ); get_ixdiag = iLookDIAG%scalarCanopyWetFraction          ! fraction of canopy that is wet
  ! snow hydrology
  case('scalarSnowAge'                  ); get_ixdiag = iLookDIAG%scalarSnowAge                    ! non-dimensional snow age (-)
  case('scalarGroundSnowFraction'       ); get_ixdiag = iLookDIAG%scalarGroundSnowFraction         ! fraction of ground that is covered with snow (-)
  case('spectralSnowAlbedoDirect'       ); get_ixdiag = iLookDIAG%spectralSnowAlbedoDirect         ! direct snow albedo for individual spectral bands (-)
  case('mLayerFracLiqSnow'              ); get_ixdiag = iLookDIAG%mLayerFracLiqSnow                ! fraction of liquid water in each snow layer (-) 
  case('mLayerThetaResid'               ); get_ixdiag = iLookDIAG%mLayerThetaResid                 ! residual volumetric water content in each snow layer (-)
  case('mLayerPoreSpace'                ); get_ixdiag = iLookDIAG%mLayerPoreSpace                  ! total pore space in each snow layer (-)
  case('mLayerMeltFreeze'               ); get_ixdiag = iLookDIAG%mLayerMeltFreeze                 ! ice content change from melt/freeze in each layer (kg m-3)
  ! soil hydrology
  case('scalarInfilArea'                ); get_ixdiag = iLookDIAG%scalarInfilArea                  ! fraction of unfrozen area where water can infiltrate (-)
  case('scalarFrozenArea'               ); get_ixdiag = iLookDIAG%scalarFrozenArea                 ! fraction of area that is considered impermeable due to soil ice (-)
  case('scalarSoilControl'              ); get_ixdiag = iLookDIAG%scalarSoilControl                ! soil control on infiltration: 1=controlling; 0=not (-)
  case('mLayerVolFracAir'               ); get_ixdiag = iLookDIAG%mLayerVolFracAir                 ! volumetric fraction of air in each layer (-)
  case('mLayerTcrit'                    ); get_ixdiag = iLookDIAG%mLayerTcrit                      ! critical soil temperature above which all water is unfrozen (K)
  case('mLayerCompress'                 ); get_ixdiag = iLookDIAG%mLayerCompress                   ! change in volumetric water content due to compression of soil (-)
  case('scalarSoilCompress'             ); get_ixdiag = iLookDIAG%scalarSoilCompress               ! change in total soil storage due to compression of the soil matrix (kg m-2)
  case('mLayerMatricHeadLiq'            ); get_ixdiag = iLookDIAG%mLayerMatricHeadLiq              ! matric potential of liquid water (m)
  ! mass balance check
  case('scalarSoilWatBalError'          ); get_ixdiag = iLookDIAG%scalarSoilWatBalError            ! error in the total soil water balance (kg m-2)
  case('scalarAquiferBalError'          ); get_ixdiag = iLookDIAG%scalarAquiferBalError            ! error in the aquifer water balance (kg m-2)
  case('scalarTotalSoilLiq'             ); get_ixdiag = iLookDIAG%scalarTotalSoilLiq               ! total mass of liquid water in the soil (kg m-2)
  case('scalarTotalSoilIce'             ); get_ixdiag = iLookDIAG%scalarTotalSoilIce               ! total mass of ice in the soil (kg m-2)
  ! variable shortcuts
  case('scalarVGn_m'                    ); get_ixdiag = iLookDIAG%scalarVGn_m                      ! van Genuchten "m" parameter (-)
  case('scalarKappa'                    ); get_ixdiag = iLookDIAG%scalarKappa                      ! constant in the freezing curve function (m K-1)
  case('scalarVolLatHt_fus'             ); get_ixdiag = iLookDIAG%scalarVolLatHt_fus               ! volumetric latent heat of fusion     (J m-3)
  ! number of function evaluations
  case('numFluxCalls'                   ); get_ixdiag = iLookDIAG%numFluxCalls                     ! number of flux calls (-)
  ! get to here if cannot find the variable
  case default
   get_ixdiag = integerMissing
 end select
 end function get_ixdiag


 ! *******************************************************************************************************************
 ! public function get_ixdiag: get the index of the named variables for the fluxes
 ! *******************************************************************************************************************
 function get_ixflux(varName)
 USE var_lookup,only:iLookFLUX                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixflux              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! net energy and mass fluxes for the vegetation domain
  case('scalarCanairNetNrgFlux'         ); get_ixflux = iLookFLUX%scalarCanairNetNrgFlux           ! net energy flux for the canopy air space (W m-2)  
  case('scalarCanopyNetNrgFlux'         ); get_ixflux = iLookFLUX%scalarCanopyNetNrgFlux           ! net energy flux for the vegetation canopy (W m-2)  
  case('scalarGroundNetNrgFlux'         ); get_ixflux = iLookFLUX%scalarGroundNetNrgFlux           ! net energy flux for the ground surface (W m-2)  
  case('scalarCanopyNetLiqFlux'         ); get_ixflux = iLookFLUX%scalarCanopyNetLiqFlux           ! net liquid water flux for the vegetation canopy (kg m-2 s-1)  
  ! forcing
  case('scalarRainfall'                 ); get_ixflux = iLookFLUX%scalarRainfall                   ! computed rainfall rate (kg m-2 s-1)
  case('scalarSnowfall'                 ); get_ixflux = iLookFLUX%scalarSnowfall                   ! computed snowfall rate (kg m-2 s-1)
  ! shortwave radiation
  case('spectralIncomingDirect'         ); get_ixflux = iLookFLUX%spectralIncomingDirect           ! incoming direct solar radiation in each wave band (W m-2)
  case('spectralIncomingDiffuse'        ); get_ixflux = iLookFLUX%spectralIncomingDiffuse          ! incoming diffuse solar radiation in each wave band (W m-2)
  case('scalarCanopySunlitPAR'          ); get_ixflux = iLookFLUX%scalarCanopySunlitPAR            ! average absorbed par for sunlit leaves (w m-2)
  case('scalarCanopyShadedPAR'          ); get_ixflux = iLookFLUX%scalarCanopyShadedPAR            ! average absorbed par for shaded leaves (w m-2)
  case('spectralBelowCanopyDirect'      ); get_ixflux = iLookFLUX%spectralBelowCanopyDirect        ! downward direct flux below veg layer for each spectral band  W m-2)
  case('spectralBelowCanopyDiffuse'     ); get_ixflux = iLookFLUX%spectralBelowCanopyDiffuse       ! downward diffuse flux below veg layer for each spectral band (W m-2)
  case('scalarBelowCanopySolar'         ); get_ixflux = iLookFLUX%scalarBelowCanopySolar           ! solar radiation transmitted below the canopy (W m-2)
  case('scalarCanopyAbsorbedSolar'      ); get_ixflux = iLookFLUX%scalarCanopyAbsorbedSolar        ! solar radiation absorbed by canopy (W m-2)
  case('scalarGroundAbsorbedSolar'      ); get_ixflux = iLookFLUX%scalarGroundAbsorbedSolar        ! solar radiation absorbed by ground (W m-2)
  ! longwave radiation
  case('scalarLWRadCanopy'              ); get_ixflux = iLookFLUX%scalarLWRadCanopy                ! longwave radiation emitted from the canopy (W m-2)
  case('scalarLWRadGround'              ); get_ixflux = iLookFLUX%scalarLWRadGround                ! longwave radiation emitted at the ground surface  (W m-2)
  case('scalarLWRadUbound2Canopy'       ); get_ixflux = iLookFLUX%scalarLWRadUbound2Canopy         ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  case('scalarLWRadUbound2Ground'       ); get_ixflux = iLookFLUX%scalarLWRadUbound2Ground         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  case('scalarLWRadUbound2Ubound'       ); get_ixflux = iLookFLUX%scalarLWRadUbound2Ubound         ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ubound'       ); get_ixflux = iLookFLUX%scalarLWRadCanopy2Ubound         ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ground'       ); get_ixflux = iLookFLUX%scalarLWRadCanopy2Ground         ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  case('scalarLWRadCanopy2Canopy'       ); get_ixflux = iLookFLUX%scalarLWRadCanopy2Canopy         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  case('scalarLWRadGround2Ubound'       ); get_ixflux = iLookFLUX%scalarLWRadGround2Ubound         ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  case('scalarLWRadGround2Canopy'       ); get_ixflux = iLookFLUX%scalarLWRadGround2Canopy         ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  case('scalarLWNetCanopy'              ); get_ixflux = iLookFLUX%scalarLWNetCanopy                ! net longwave radiation at the canopy (W m-2)
  case('scalarLWNetGround'              ); get_ixflux = iLookFLUX%scalarLWNetGround                ! net longwave radiation at the ground surface (W m-2)
  case('scalarLWNetUbound'              ); get_ixflux = iLookFLUX%scalarLWNetUbound                ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! turbulent heat transfer
  case('scalarEddyDiffusCanopyTop'      ); get_ixflux = iLookFLUX%scalarEddyDiffusCanopyTop        ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  case('scalarFrictionVelocity'         ); get_ixflux = iLookFLUX%scalarFrictionVelocity           ! friction velocity - canopy momentum sink (m s-1)
  case('scalarWindspdCanopyTop'         ); get_ixflux = iLookFLUX%scalarWindspdCanopyTop           ! windspeed at the top of the canopy (m s-1)
  case('scalarWindspdCanopyBottom'      ); get_ixflux = iLookFLUX%scalarWindspdCanopyBottom        ! windspeed at the height of the bottom of the canopy (m s-1)
  case('scalarGroundResistance'         ); get_ixflux = iLookFLUX%scalarGroundResistance           ! below canopy aerodynamic resistance (s m-1)
  case('scalarCanopyResistance'         ); get_ixflux = iLookFLUX%scalarCanopyResistance           ! above canopy aerodynamic resistance (s m-1)
  case('scalarLeafResistance'           ); get_ixflux = iLookFLUX%scalarLeafResistance             ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  case('scalarSoilResistance'           ); get_ixflux = iLookFLUX%scalarSoilResistance             ! soil surface resistance (s m-1)
  case('scalarSenHeatTotal'             ); get_ixflux = iLookFLUX%scalarSenHeatTotal               ! sensible heat from the canopy air space to the atmosphere (W m-2)
  case('scalarSenHeatCanopy'            ); get_ixflux = iLookFLUX%scalarSenHeatCanopy              ! sensible heat from the canopy to the canopy air space (W m-2)
  case('scalarSenHeatGround'            ); get_ixflux = iLookFLUX%scalarSenHeatGround              ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarLatHeatTotal'             ); get_ixflux = iLookFLUX%scalarLatHeatTotal               ! latent heat from the canopy air space to the atmosphere (W m-2)
  case('scalarLatHeatCanopyEvap'        ); get_ixflux = iLookFLUX%scalarLatHeatCanopyEvap          ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatCanopyTrans'       ); get_ixflux = iLookFLUX%scalarLatHeatCanopyTrans         ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatGround'            ); get_ixflux = iLookFLUX%scalarLatHeatGround              ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarCanopyAdvectiveHeatFlux'  ); get_ixflux = iLookFLUX%scalarCanopyAdvectiveHeatFlux    ! heat advected to the canopy surface with rain + snow (W m-2)
  case('scalarGroundAdvectiveHeatFlux'  ); get_ixflux = iLookFLUX%scalarGroundAdvectiveHeatFlux    ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  case('scalarCanopySublimation'        ); get_ixflux = iLookFLUX%scalarCanopySublimation          ! canopy sublimation/frost (kg m-2 s-1)
  case('scalarSnowSublimation'          ); get_ixflux = iLookFLUX%scalarSnowSublimation            ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  case('scalarStomResistSunlit'         ); get_ixflux = iLookFLUX%scalarStomResistSunlit           ! stomatal resistance for sunlit leaves (s m-1)
  case('scalarStomResistShaded'         ); get_ixflux = iLookFLUX%scalarStomResistShaded           ! stomatal resistance for shaded leaves (s m-1)
  case('scalarPhotosynthesisSunlit'     ); get_ixflux = iLookFLUX%scalarPhotosynthesisSunlit       ! sunlit photosynthesis (umolco2 m-2 s-1)
  case('scalarPhotosynthesisShaded'     ); get_ixflux = iLookFLUX%scalarPhotosynthesisShaded       ! shaded photosynthesis (umolco2 m-2 s-1)
  case('scalarCanopyTranspiration'      ); get_ixflux = iLookFLUX%scalarCanopyTranspiration        ! canopy transpiration (kg m-2 s-1)
  case('scalarCanopyEvaporation'        ); get_ixflux = iLookFLUX%scalarCanopyEvaporation          ! canopy evaporation/condensation (kg m-2 s-1)
  case('scalarGroundEvaporation'        ); get_ixflux = iLookFLUX%scalarGroundEvaporation          ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  case('mLayerTranspire'                ); get_ixflux = iLookFLUX%mLayerTranspire                  ! transpiration loss from each soil layer (kg m-2 s-1)
  ! liquid and solid water fluxes through the canopy
  case('scalarThroughfallSnow'          ); get_ixflux = iLookFLUX%scalarThroughfallSnow            ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarThroughfallRain'          ); get_ixflux = iLookFLUX%scalarThroughfallRain            ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarCanopySnowUnloading'      ); get_ixflux = iLookFLUX%scalarCanopySnowUnloading        ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  case('scalarCanopyLiqDrainage'        ); get_ixflux = iLookFLUX%scalarCanopyLiqDrainage          ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  case('scalarCanopyMeltFreeze'         ); get_ixflux = iLookFLUX%scalarCanopyMeltFreeze           ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! energy fluxes and for the snow and soil domains
  case('iLayerConductiveFlux'           ); get_ixflux = iLookFLUX%iLayerConductiveFlux             ! conductive energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerAdvectiveFlux'            ); get_ixflux = iLookFLUX%iLayerAdvectiveFlux              ! advective energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerNrgFlux'                  ); get_ixflux = iLookFLUX%iLayerNrgFlux                    ! energy flux at layer interfaces at the end of the time step (W m-2)
  case('mLayerNrgFlux'                  ); get_ixflux = iLookFLUX%mLayerNrgFlux                    ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes for the snow domain 
  case('scalarSnowDrainage'             ); get_ixflux = iLookFLUX%scalarSnowDrainage               ! drainage from the bottom of the snow profile (m s-1)
  case('iLayerLiqFluxSnow'              ); get_ixflux = iLookFLUX%iLayerLiqFluxSnow                ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  case('mLayerLiqFluxSnow'              ); get_ixflux = iLookFLUX%mLayerLiqFluxSnow                ! net liquid water flux for each snow layer (s-1)
  ! liquid water fluxes for the soil domain 
  case('scalarRainPlusMelt'             ); get_ixflux = iLookFLUX%scalarRainPlusMelt               ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('scalarMaxInfilRate'             ); get_ixflux = iLookFLUX%scalarMaxInfilRate               ! maximum infiltration rate (m s-1)
  case('scalarInfiltration'             ); get_ixflux = iLookFLUX%scalarInfiltration               ! infiltration of water into the soil profile (m s-1)
  case('scalarExfiltration'             ); get_ixflux = iLookFLUX%scalarExfiltration               ! exfiltration of water from the top of the soil profile (m s-1)
  case('scalarSurfaceRunoff'            ); get_ixflux = iLookFLUX%scalarSurfaceRunoff              ! surface runoff (m s-1)
  case('mLayerSatHydCondMP'             ); get_ixflux = iLookFLUX%mLayerSatHydCondMP               ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  case('mLayerSatHydCond'               ); get_ixflux = iLookFLUX%mLayerSatHydCond                 ! saturated hydraulic conductivity in each layer (m s-1)
  case('iLayerSatHydCond'               ); get_ixflux = iLookFLUX%iLayerSatHydCond                 ! saturated hydraulic conductivity in each layer interface (m s-1)
  case('mLayerHydCond'                  ); get_ixflux = iLookFLUX%mLayerHydCond                    ! hydraulic conductivity in each layer (m s-1)
  case('iLayerLiqFluxSoil'              ); get_ixflux = iLookFLUX%iLayerLiqFluxSoil                ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  case('mLayerLiqFluxSoil'              ); get_ixflux = iLookFLUX%mLayerLiqFluxSoil                ! net liquid water flux for each soil layer (s-1)
  case('mLayerBaseflow'                 ); get_ixflux = iLookFLUX%mLayerBaseflow                   ! baseflow from each soil layer (m s-1)
  case('mLayerColumnInflow'             ); get_ixflux = iLookFLUX%mLayerColumnInflow               ! total inflow to each layer in a given soil column (m3 s-1)
  case('mLayerColumnOutflow'            ); get_ixflux = iLookFLUX%mLayerColumnOutflow              ! total outflow from each layer in a given soil column (m3 s-1)
  case('scalarSoilBaseflow'             ); get_ixflux = iLookFLUX%scalarSoilBaseflow               ! total baseflow from throughout the soil profile (m s-1)
  case('scalarSoilDrainage'             ); get_ixflux = iLookFLUX%scalarSoilDrainage               ! drainage from the bottom of the soil profile (m s-1)
  case('scalarAquiferRecharge'          ); get_ixflux = iLookFLUX%scalarAquiferRecharge            ! recharge to the aquifer (m s-1)
  case('scalarAquiferTranspire'         ); get_ixflux = iLookFLUX%scalarAquiferTranspire           ! transpiration from the aquifer (m s-1)
  case('scalarAquiferBaseflow'          ); get_ixflux = iLookFLUX%scalarAquiferBaseflow            ! baseflow from the aquifer (m s-1)
  case default
   get_ixflux = integerMissing
 end select
 end function get_ixflux


 ! *******************************************************************************************************************
 ! public function get_ixderiv: get the index of the named variables for the model derivatives
 ! *******************************************************************************************************************
 function get_ixderiv(varName)
 USE var_lookup,only:iLookDERIV                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixderiv             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
  case('dCanairNetFlux_dCanairTemp'     ); get_ixderiv = iLookDERIV%dCanairNetFlux_dCanairTemp     ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dCanairNetFlux_dCanopyTemp'     ); get_ixderiv = iLookDERIV%dCanairNetFlux_dCanopyTemp     ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  case('dCanairNetFlux_dGroundTemp'     ); get_ixderiv = iLookDERIV%dCanairNetFlux_dGroundTemp     ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanairTemp'     ); get_ixderiv = iLookDERIV%dCanopyNetFlux_dCanairTemp     ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanopyTemp'     ); get_ixderiv = iLookDERIV%dCanopyNetFlux_dCanopyTemp     ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  case('dCanopyNetFlux_dGroundTemp'     ); get_ixderiv = iLookDERIV%dCanopyNetFlux_dGroundTemp     ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanLiq'         ); get_ixderiv = iLookDERIV%dCanopyNetFlux_dCanLiq         ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  case('dGroundNetFlux_dCanairTemp'     ); get_ixderiv = iLookDERIV%dGroundNetFlux_dCanairTemp     ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dGroundNetFlux_dCanopyTemp'     ); get_ixderiv = iLookDERIV%dGroundNetFlux_dCanopyTemp     ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  case('dGroundNetFlux_dGroundTemp'     ); get_ixderiv = iLookDERIV%dGroundNetFlux_dGroundTemp     ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  case('dGroundNetFlux_dCanLiq'         ); get_ixderiv = iLookDERIV%dGroundNetFlux_dCanLiq         ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  case('dCanopyEvaporation_dTCanair'    ); get_ixderiv = iLookDERIV%dCanopyEvaporation_dTCanair    ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dTCanopy'    ); get_ixderiv = iLookDERIV%dCanopyEvaporation_dTCanopy    ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dTGround'    ); get_ixderiv = iLookDERIV%dCanopyEvaporation_dTGround    ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dCanLiq'     ); get_ixderiv = iLookDERIV%dCanopyEvaporation_dCanLiq     ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
  case('dGroundEvaporation_dTCanair'    ); get_ixderiv = iLookDERIV%dGroundEvaporation_dTCanair    ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dTCanopy'    ); get_ixderiv = iLookDERIV%dGroundEvaporation_dTCanopy    ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dTGround'    ); get_ixderiv = iLookDERIV%dGroundEvaporation_dTGround    ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dCanLiq'     ); get_ixderiv = iLookDERIV%dGroundEvaporation_dCanLiq     ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
  ! derivatives in canopy water w.r.t canopy temperature
  case('dTheta_dTkCanopy'               ); get_ixderiv = iLookDERIV%dTheta_dTkCanopy               ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  case('dCanLiq_dTcanopy'               ); get_ixderiv = iLookDERIV%dCanLiq_dTcanopy               ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
  ! derivatives in canopy liquid fluxes w.r.t. canopy water
  case('scalarCanopyLiqDeriv'           ); get_ixderiv = iLookDERIV%scalarCanopyLiqDeriv           ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
  case('scalarThroughfallRainDeriv'     ); get_ixderiv = iLookDERIV%scalarThroughfallRainDeriv     ! derivative in throughfall w.r.t. canopy liquid water (s-1)
  case('scalarCanopyLiqDrainageDeriv'   ); get_ixderiv = iLookDERIV%scalarCanopyLiqDrainageDeriv   ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
  case('dNrgFlux_dTempAbove'            ); get_ixderiv = iLookDERIV%dNrgFlux_dTempAbove            ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  case('dNrgFlux_dTempBelow '           ); get_ixderiv = iLookDERIV%dNrgFlux_dTempBelow            ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
  case('iLayerLiqFluxSnowDeriv'         ); get_ixderiv = iLookDERIV%iLayerLiqFluxSnowDeriv         ! derivative in vertical liquid water flux at layer interfaces (m s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  case('dVolTot_dPsi0'                  ); get_ixderiv = iLookDERIV%dVolTot_dPsi0                  ! derivative in total water content w.r.t. total water matric potential (m-1)
  case('dq_dHydStateAbove'              ); get_ixderiv = iLookDERIV%dq_dHydStateAbove              ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  case('dq_dHydStateBelow'              ); get_ixderiv = iLookDERIV%dq_dHydStateBelow              ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  case('mLayerdTheta_dPsi'              ); get_ixderiv = iLookDERIV%mLayerdTheta_dPsi              ! derivative in the soil water characteristic w.r.t. psi (m-1)
  case('mLayerdPsi_dTheta'              ); get_ixderiv = iLookDERIV%mLayerdPsi_dTheta              ! derivative in the soil water characteristic w.r.t. theta (m)
  case('dCompress_dPsi'                 ); get_ixderiv = iLookDERIV%dCompress_dPsi                 ! derivative in compressibility w.r.t matric head (m-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  case('dq_dNrgStateAbove'              ); get_ixderiv = iLookDERIV%dq_dNrgStateAbove              ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  case('dq_dNrgStateBelow'              ); get_ixderiv = iLookDERIV%dq_dNrgStateBelow              ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  case('mLayerdTheta_dTk'               ); get_ixderiv = iLookDERIV%mLayerdTheta_dTk               ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  case('dPsiLiq_dTemp'                  ); get_ixderiv = iLookDERIV%dPsiLiq_dTemp                  ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
  case('dPsiLiq_dPsi0'                  ); get_ixderiv = iLookDERIV%dPsiLiq_dPsi0                  ! derivative in liquid matric potential w.r.t. total  matric potential (-)

  case default
   get_ixderiv = integerMissing
 end select
 end function get_ixderiv


 ! *******************************************************************************************************************
 ! public function get_ixindex: get the index of the named variables for the model indices
 ! *******************************************************************************************************************
 function get_ixindex(varName)
 USE var_lookup,only:iLookINDEX                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixINDEX             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! number of model layers, and layer indices
  case('nSnow'               ); get_ixINDEX = iLookINDEX%nSnow               ! number of snow layers                                                    (-)
  case('nSoil'               ); get_ixINDEX = iLookINDEX%nSoil               ! number of soil layers                                                    (-)
  case('nLayers'             ); get_ixINDEX = iLookINDEX%nLayers             ! total number of layers                                                   (-)
  case('layerType'           ); get_ixINDEX = iLookINDEX%layerType           ! index defining type of layer (snow or soil)                              (-)
  ! number of state variables of different type
  case('nCasNrg'             ); get_ixINDEX = iLookINDEX%nCasNrg             ! number of energy state variables for the canopy air space domain         (-)
  case('nVegNrg'             ); get_ixINDEX = iLookINDEX%nVegNrg             ! number of energy state variables for the vegetation canopy               (-)
  case('nVegMass'            ); get_ixINDEX = iLookINDEX%nVegMass            ! number of hydrology states for vegetation (mass of water)                (-)
  case('nVegState'           ); get_ixINDEX = iLookINDEX%nVegState           ! number of vegetation state variables                                     (-)
  case('nNrgState'           ); get_ixINDEX = iLookINDEX%nNrgState           ! number of energy state variables                                         (-)
  case('nWatState'           ); get_ixINDEX = iLookINDEX%nWatState           ! number of "total water" states (vol. total water content)                (-)
  case('nMatState'           ); get_ixINDEX = iLookINDEX%nMatState           ! number of matric head state variables                                    (-)
  case('nMassState'          ); get_ixINDEX = iLookINDEX%nMassState          ! number of hydrology state variables (mass of water)                      (-)
  case('nState'              ); get_ixINDEX = iLookINDEX%nState              ! total number of model state variables                                    (-)
  ! number of state variables within different domains in the snow+soil system
  case('nSnowSoilNrg'        ); get_ixINDEX = iLookINDEX%nSnowSoilNrg        ! number of energy states in the snow+soil domain                          (-)
  case('nSnowOnlyNrg'        ); get_ixINDEX = iLookINDEX%nSnowOnlyNrg        ! number of energy states in the snow domain                               (-)
  case('nSoilOnlyNrg'        ); get_ixINDEX = iLookINDEX%nSoilOnlyNrg        ! number of energy states in the soil domain                               (-)
  case('nSnowSoilHyd'        ); get_ixINDEX = iLookINDEX%nSnowSoilHyd        ! number of hydrology states in the snow+soil domain                       (-)
  case('nSnowOnlyHyd'        ); get_ixINDEX = iLookINDEX%nSnowOnlyHyd        ! number of hydrology states in the snow domain                            (-)
  case('nSoilOnlyHyd'        ); get_ixINDEX = iLookINDEX%nSoilOnlyHyd        ! number of hydrology states in the soil domain                            (-)
  ! type of model state variables
  case('ixControlVolume'     ); get_ixINDEX = iLookINDEX%ixControlVolume     ! index of the control volume for different domains (veg, snow, soil)      (-)
  case('ixDomainType'        ); get_ixINDEX = iLookINDEX%ixDomainType        ! index of the type of domain (iname_veg, iname_snow, iname_soil)          (-)
  case('ixStateType'         ); get_ixINDEX = iLookINDEX%ixStateType         ! index of the type of every state variable (iname_nrgCanair, ...)         (-)
  case('ixHydType'           ); get_ixINDEX = iLookINDEX%ixHydType           ! index of the type of hydrology states in snow+soil domain                (-)
  ! type of model state variables (state subset)
  case('ixDomainType_subset' ); get_ixINDEX = iLookINDEX%ixDomainType_subset ! [state subset] id of domain for desired model state variables            (-)
  case('ixStateType_subset'  ); get_ixINDEX = iLookINDEX%ixStateType_subset  ! [state subset] type of desired model state variables                     (-)
  ! mapping between state subset and the full state vector
  case('ixMapFull2Subset'    ); get_ixINDEX = iLookINDEX%ixMapFull2Subset    ! list of indices of the state subset in the full state vector             (-)
  case('ixMapSubset2Full'    ); get_ixINDEX = iLookINDEX%ixMapSubset2Full    ! list of indices of the full state vector in the state subset             (-)
  ! indices of model specific state variables
  case('ixCasNrg'            ); get_ixINDEX = iLookINDEX%ixCasNrg            ! index of canopy air space energy state variable                          (-)
  case('ixVegNrg'            ); get_ixINDEX = iLookINDEX%ixVegNrg            ! index of canopy energy state variable                                    (-)
  case('ixVegHyd'            ); get_ixINDEX = iLookINDEX%ixVegHyd            ! index of canopy hydrology state variable (mass)                          (-)
  case('ixTopNrg'            ); get_ixINDEX = iLookINDEX%ixTopNrg            ! index of upper-most energy state in the snow+soil subdomain              (-)
  case('ixTopHyd'            ); get_ixINDEX = iLookINDEX%ixTopHyd            ! index of upper-most hydrology state in the snow+soil subdomain           (-)
  ! vectors of indices for specific state types
  case('ixNrgOnly'           ); get_ixINDEX = iLookINDEX%ixNrgOnly           ! indices IN THE STATE SUBSET for all energy states                        (-)
  case('ixHydOnly'           ); get_ixINDEX = iLookINDEX%ixHydOnly           ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-)
  case('ixMatOnly'           ); get_ixINDEX = iLookINDEX%ixMatOnly           ! indices IN THE STATE SUBSET for matric head state variables              (-)
  case('ixMassOnly'          ); get_ixINDEX = iLookINDEX%ixMassOnly          ! indices IN THE STATE SUBSET for hydrology states (mass of water)         (-)
  ! vectors of indicesfor specific state types within specific sub-domains
  case('ixSnowSoilNrg'       ); get_ixINDEX = iLookINDEX%ixSnowSoilNrg       ! indices IN THE STATE SUBSET for energy states in the snow+soil domain    (-)
  case('ixSnowOnlyNrg'       ); get_ixINDEX = iLookINDEX%ixSnowOnlyNrg       ! indices IN THE STATE SUBSET for energy states in the snow domain         (-)
  case('ixSoilOnlyNrg'       ); get_ixINDEX = iLookINDEX%ixSoilOnlyNrg       ! indices IN THE STATE SUBSET for energy states in the soil domain         (-)
  case('ixSnowSoilHyd'       ); get_ixINDEX = iLookINDEX%ixSnowSoilHyd       ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-)
  case('ixSnowOnlyHyd'       ); get_ixINDEX = iLookINDEX%ixSnowOnlyHyd       ! indices IN THE STATE SUBSET for hydrology states in the snow domain      (-)
  case('ixSoilOnlyHyd'       ); get_ixINDEX = iLookINDEX%ixSoilOnlyHyd       ! indices IN THE STATE SUBSET for hydrology states in the soil domain      (-)
  ! vectors of indices for specfic state types within specific sub-domains
  case('ixNrgCanair'         ); get_ixINDEX = iLookINDEX%ixNrgCanair         ! indices IN THE STATE SUBSET for energy states in canopy air space domain (-) 
  case('ixNrgCanopy'         ); get_ixINDEX = iLookINDEX%ixNrgCanopy         ! indices IN THE STATE SUBSET for energy states in the canopy domain       (-) 
  case('ixHydCanopy'         ); get_ixINDEX = iLookINDEX%ixHydCanopy         ! indices IN THE STATE SUBSET for hydrology states in the canopy domain    (-) 
  case('ixNrgLayer'          ); get_ixINDEX = iLookINDEX%ixNrgLayer          ! indices IN THE FULL VECTOR for energy states in the snow+soil domain     (-)
  case('ixHydLayer'          ); get_ixINDEX = iLookINDEX%ixHydLayer          ! indices IN THE FULL VECTOR for hydrology states in the snow+soil domain  (-)
  ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
  case('ixVolFracWat'        ); get_ixINDEX = iLookINDEX%ixVolFracWat        ! indices IN THE SNOW+SOIL VECTOR for hyd states                           (-)
  case('ixMatricHead'        ); get_ixINDEX = iLookINDEX%ixMatricHead        ! indices IN THE SOIL VECTOR for hyd states                                (-)
  ! indices within state vectors
  case('ixAllState'          ); get_ixINDEX = iLookINDEX%ixAllState          ! list of indices for all model state variables                            (-)
  case('ixSoilState'         ); get_ixINDEX = iLookINDEX%ixSoilState         ! list of indices for all soil layers                                      (-)
  case('ixLayerState'        ); get_ixINDEX = iLookINDEX%ixLayerState        ! list of indices for all model layers                                     (-)
  ! indices for the model output files
  case('midSnowStartIndex'   ); get_ixINDEX = iLookINDEX%midSnowStartIndex   ! start index of the midSnow vector for a given timestep                   (-)
  case('midSoilStartIndex'   ); get_ixINDEX = iLookINDEX%midSoilStartIndex   ! start index of the midSoil vector for a given timestep                   (-)
  case('midTotoStartIndex'   ); get_ixINDEX = iLookINDEX%midTotoStartIndex   ! start index of the midToto vector for a given timestep                   (-)
  case('ifcSnowStartIndex'   ); get_ixINDEX = iLookINDEX%ifcSnowStartIndex   ! start index of the ifcSnow vector for a given timestep                   (-)
  case('ifcSoilStartIndex'   ); get_ixINDEX = iLookINDEX%ifcSoilStartIndex   ! start index of the ifcSoil vector for a given timestep                   (-)
  case('ifcTotoStartIndex'   ); get_ixINDEX = iLookINDEX%ifcTotoStartIndex   ! start index of the ifcToto vector for a given timestep                   (-)
  case default
   get_ixindex = integerMissing
 end select
 end function get_ixindex


 ! *******************************************************************************************************************
 ! public function get_ixbpar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixbpar(varName)
 USE var_lookup,only:iLookBPAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixbpar              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! baseflow
  case('basin__aquiferHydCond'    ); get_ixbpar = iLookBPAR%basin__aquiferHydCond     ! hydraulic conductivity of the basin aquifer (m s-1)
  case('basin__aquiferScaleFactor'); get_ixbpar = iLookBPAR%basin__aquiferScaleFactor ! scaling factor for aquifer storage in the big bucket (m)
  case('basin__aquiferBaseflowExp'); get_ixbpar = iLookBPAR%basin__aquiferBaseflowExp ! baseflow exponent for the big bucket (-)
  ! sub-grid routing
  case('routingGammaShape'        ); get_ixbpar = iLookBPAR%routingGammaShape         ! shape parameter in Gamma distribution used for sub-grid routing (-)
  case('routingGammaScale'        ); get_ixbpar = iLookBPAR%routingGammaScale         ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! get to here if cannot find the variable
  case default
   get_ixbpar = integerMissing
 end select
 end function get_ixbpar


 ! *******************************************************************************************************************
 ! public function get_ixbvar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixbvar(varName)
 USE var_lookup,only:iLookBVAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixbvar              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! derived variables
  case('basin__TotalArea'              ); get_ixbvar = iLookBVAR%basin__totalArea                ! total basin area (m2)
  ! scalar variables -- basin-average runoff and aquifer fluxes
  case('basin__SurfaceRunoff'          ); get_ixbvar = iLookBVAR%basin__SurfaceRunoff            ! surface runoff (m s-1)
  case('basin__ColumnOutflow'          ); get_ixbvar = iLookBVAR%basin__ColumnOutflow            ! outflow from all "outlet" HRUs (those with no downstream HRU)
  case('basin__AquiferStorage'         ); get_ixbvar = iLookBVAR%basin__AquiferStorage           ! aquifer storage (m s-1)
  case('basin__AquiferRecharge'        ); get_ixbvar = iLookBVAR%basin__AquiferRecharge          ! recharge to the aquifer (m s-1)
  case('basin__AquiferBaseflow'        ); get_ixbvar = iLookBVAR%basin__AquiferBaseflow          ! baseflow from the aquifer (m s-1)
  case('basin__AquiferTranspire'       ); get_ixbvar = iLookBVAR%basin__AquiferTranspire         ! transpiration from the aquifer (m s-1)
  ! variables to compute runoff
  case('routingRunoffFuture'           ); get_ixbvar = iLookBVAR%routingRunoffFuture             ! runoff in future time steps (m s-1)
  case('routingFractionFuture'         ); get_ixbvar = iLookBVAR%routingFractionFuture           ! fraction of runoff in future time steps (-)
  case('averageInstantRunoff'          ); get_ixbvar = iLookBVAR%averageInstantRunoff            ! instantaneous runoff (m s-1)
  case('averageRoutedRunoff'           ); get_ixbvar = iLookBVAR%averageRoutedRunoff             ! routed runoff (m s-1)
  ! get to here if cannot find the variable
  case default
   get_ixbvar = integerMissing
 end select
 end function get_ixbvar

 ! *********************************************************************************************************
 ! public function get_ixVarType: get the index of the named variable type
 ! *********************************************************************************************************
 function get_ixVarType(varType)
 USE var_lookup,only:iLookVarType                    ! indices of the named variable types
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varType                 ! variable type name
 integer(i4b)             :: get_ixVarType          ! index of the named variable type list
 ! get the index of the named variables
 select case(trim(varType))
  case('scalarv'); get_ixVarType = iLookVarType%scalarv
  case('wLength'); get_ixVarType = iLookVarType%wLength
  case('midSnow'); get_ixVarType = iLookVarType%midSnow
  case('midSoil'); get_ixVarType = iLookVarType%midSoil
  case('midToto'); get_ixVarType = iLookVarType%midToto
  case('ifcSnow'); get_ixVarType = iLookVarType%ifcSnow
  case('ifcSoil'); get_ixVarType = iLookVarType%ifcSoil
  case('ifcToto'); get_ixVarType = iLookVarType%ifcToto
  case('parSoil'); get_ixVarType = iLookVarType%parSoil
  case('routing'); get_ixVarType = iLookVarType%routing
  case('unknown'); get_ixVarType = iLookVarType%unknown
  ! get to here if cannot find the variable
  case default
   get_ixVarType = integerMissing
 end select
 end function get_ixVarType

 ! ****************************************************************************************************************
 ! public function get_varTypeName: get the index of the named variable type
 ! ****************************************************************************************************************
 function get_varTypeName(varType)
 USE var_lookup,only:iLookVarType                    ! indices of the named variable types
 implicit none
 ! define dummy variables
 integer(i4b), intent(in) :: varType                 ! variable type name
 character(LEN=7)         :: get_varTypeName         ! index of the named variable type list
 ! get the index of the named variables
 select case(varType)
  case(iLookVarType%scalarv);get_varTypeName='scalarv'
  case(iLookVarType%wLength);get_varTypeName='wLength'
  case(iLookVarType%midSnow);get_varTypeName='midSnow'
  case(iLookVarType%midSoil);get_varTypeName='midSoil'
  case(iLookVarType%midToto);get_varTypeName='midToto'
  case(iLookVarType%ifcSnow);get_varTypeName='ifcSnow'
  case(iLookVarType%ifcSoil);get_varTypeName='ifcSoil'
  case(iLookVarType%ifcToto);get_varTypeName='ifcToto'
  case(iLookVarType%parSoil);get_varTypeName='parSoil'
  case(iLookVarType%routing);get_varTypeName='routing'
  case(iLookVarType%unknown);get_varTypeName='unknown'
  ! get to here if cannot find the variable
  case default
   get_VarTypeName = 'missing'
 end select
 end function get_VarTypeName

 ! *******************************************************************************************************************
 ! public subroutine get_ixUnknown: get the index of the named variable type from ANY structure, as well as the 
 ! structrue that it was found in
 ! *******************************************************************************************************************
 subroutine get_ixUnknown(varName,typeName,vDex,err,message)
 USE nrtype
 USE globalData,only:structInfo        ! information on the data structures                  
 implicit none

 ! dummies
 character(*),intent(in)  :: varName   ! variable name
 character(*),intent(out) :: typeName  ! variable type name
 integer(i4b),intent(out) :: vDex      ! variable index in structure
 integer(i4b),intent(out) :: err       ! error code
 character(*),intent(out) :: message   ! error message

 ! internals
 integer(i4b)             :: iStruc    ! index for looping through structure types

 ! error init
 err=0
 message='get_ixUnknown/'

 ! loop through all structure types to find the one with the given variable name
 ! pill variable index plus return which structure it was found in
 do iStruc = 1,size(structInfo)
  select case(trim(structInfo(iStruc)%structName))
   case ('time' ); vDex = get_ixTime(trim(varName))
   case ('forc' ); vDex = get_ixForce(trim(varName))
   case ('attr' ); vDex = get_ixAttr(trim(varName))
   case ('type' ); vDex = get_ixType(trim(varName))
   case ('mpar' ); vDex = get_ixParam(trim(varName))
   case ('indx' ); vDex = get_ixIndex(trim(varName))
   case ('prog' ); vDex = get_ixProg(trim(varName))
   case ('diag' ); vDex = get_ixDiag(trim(varName))
   case ('flux' ); vDex = get_ixFlux(trim(varName))
   case ('bpar' ); vDex = get_ixBpar(trim(varName))
   case ('bvar' ); vDex = get_ixBvar(trim(varName))
   case ('deriv'); vDex = get_ixDeriv(trim(varName))
  end select
  if (vDex>0) then; typeName=trim(structInfo(iStruc)%structName); return; end if
 end do

 ! 404
 err=20;message=trim(message)//'variable '//trim(varName)//' is not found in any structure'; return

 end subroutine get_ixUnknown

 ! ***************************************************************************************************************
 ! public function get_statName: get the name of the output statistics type
 ! ***************************************************************************************************************
 function get_statName(istat)
 USE var_lookup,only:iLookStat                   ! indices of the possible output statistics
 implicit none
 ! define dummy variables
 integer(i4b), intent(in) :: istat               ! stat type name
 character(LEN=10)         :: get_statName        ! index of the named variable type list
 ! get the index of the named variables
 select case(istat)
  case(iLookStat%totl);get_statName='total'
  case(iLookStat%inst);get_statName='instant'
  case(iLookStat%mean);get_statName='mean'
  case(iLookStat%vari);get_statName='variance'
  case(iLookStat%mini);get_statName='minimum'
  case(iLookStat%maxi);get_statName='maximum'
  case(iLookStat%mode);get_statName='mode'
  ! get to here if cannot find the variable
  case default
   get_statName = 'unknown'
 end select
 end function get_statName

end module get_ixname_module
