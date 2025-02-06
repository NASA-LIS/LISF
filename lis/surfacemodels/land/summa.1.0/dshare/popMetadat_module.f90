module popMetadat_module
USE nrtype, integerMissing=>nr_integerMissing
implicit none
! define indices in metadata structures
integer(i4b),parameter   :: modelTime=1     ! to force index variables to be output at model timestep
integer(i4b),parameter   :: nameIndex=1     ! index of the variable name
integer(i4b),parameter   :: freqIndex=3     ! index of the output frequency
! define indices in flag vectors
integer(i4b),parameter   :: indexMidSnow=1  ! index of flag vector: midSnow
integer(i4b),parameter   :: indexMidSoil=2  ! index of flag vector: midSoil
integer(i4b),parameter   :: indexMidToto=3  ! index of flag vector: midToto
integer(i4b),parameter   :: indexIfcSnow=4  ! index of flag vector: ifcSnow
integer(i4b),parameter   :: indexIfcSoil=5  ! index of flag vector: ifcSoil
integer(i4b),parameter   :: indexIfcToto=6  ! index of flag vector: ifcToto
private
public::popMetadat
contains

 subroutine popMetadat(err,message)
 ! data structures
 USE data_types, only: var_info   ! data type for metadata structure
 USE globalData, only: time_meta  ! data structure for time metadata
 USE globalData, only: forc_meta  ! data structure for forcing metadata
 USE globalData, only: type_meta  ! data structure for categorical metadata
 USE globalData, only: attr_meta  ! data structure for attribute metadata
 USE globalData, only: mpar_meta  ! data structure for local parameter metadata
 USE globalData, only: bpar_meta  ! data structure for basin parameter metadata
 USE globalData, only: bvar_meta  ! data structure for basin model variable metadata
 USE globalData, only: indx_meta  ! data structure for index metadata
 USE globalData, only: prog_meta  ! data structure for local prognostic (state) variables
 USE globalData, only: diag_meta  ! data structure for local diagnostic variables
 USE globalData, only: flux_meta  ! data structure for local flux variables
 USE globalData, only: deriv_meta ! data structure for local flux derivatives
 ! structures of named variables
 USE var_lookup, only: iLookTIME  ! named variables for time data structure
 USE var_lookup, only: iLookFORCE ! named variables for forcing data structure
 USE var_lookup, only: iLookTYPE  ! named variables for categorical attribute data structure
 USE var_lookup, only: iLookATTR  ! named variables for real valued attribute data structure
 USE var_lookup, only: iLookPARAM ! named variables for local parameter data structure
 USE var_lookup, only: iLookBPAR  ! named variables for basin parameter data structure
 USE var_lookup, only: iLookBVAR  ! named variables for basin model variable data structure
 USE var_lookup, only: iLookINDEX ! named variables for index variable data structure
 USE var_lookup, only: iLookPROG  ! named variables for local state variables
 USE var_lookup, only: iLookDIAG  ! named variables for local diagnostic variables
 USE var_lookup, only: iLookFLUX  ! named variables for local flux variables
 USE var_lookup, only: iLookDERIV ! named variables for local flux derivatives
 USE var_lookup, only: maxvarStat ! size of arrays in structure constructor
 USE get_ixName_module,only:get_ixVarType ! to turn vartype strings to integers
 implicit none
 ! dummy variables
 integer(i4b),intent(out)       :: err           ! error code
 character(*),intent(out)       :: message       ! error message
 ! internals
 character(256)                 :: cmessage      ! error message
 integer,dimension(maxVarStat)  :: iMissArry     ! arry of missing integers 
 logical,dimension(maxVarStat)  :: lFalseArry    ! arry of false logicals 
 ! initialize error control
 err=0; message='popMetadat/'

 ! init arrays for structure constructors
 iMissArry = integerMissing
 lFalseArry = .false.

 ! -----
 ! * model time structures...
 ! --------------------------
 time_meta(iLookTIME%iyyy)                   = var_info('iyyy', 'year'  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 time_meta(iLookTIME%im)                     = var_info('im'  , 'month' , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 time_meta(iLookTIME%id)                     = var_info('id'  , 'day'   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 time_meta(iLookTIME%ih)                     = var_info('ih'  , 'hour'  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 time_meta(iLookTIME%imin)                   = var_info('imin', 'minute', '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * model forcing data...
 ! -----------------------
 forc_meta(iLookFORCE%time)                  = var_info('time'    , 'time since time reference'                         , 'seconds since 1990-1-1 0:0:0.0 -0:00', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%pptrate)               = var_info('pptrate' , 'precipitation rate'                                , 'kg m-2 s-1'                          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%SWRadAtm)              = var_info('SWRadAtm', 'downward shortwave radiation at the upper boundary', 'W m-2'                               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%LWRadAtm)              = var_info('LWRadAtm', 'downward longwave radiation at the upper boundary' , 'W m-2'                               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%airtemp)               = var_info('airtemp' , 'air temperature at the measurement height'         , 'K'                                   , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%windspd)               = var_info('windspd' , 'wind speed at the measurement height'              , 'm s-1'                               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%airpres)               = var_info('airpres' , 'air pressure at the the measurement height'        , 'Pa'                                  , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 forc_meta(iLookFORCE%spechum)               = var_info('spechum' , 'specific humidity at the measurement height'       , 'g g-1'                               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * categorical data...
 ! ---------------------
 type_meta(iLookTYPE%hruIndex)               = var_info('hruIndex'      , 'index defining the hydrologic response unit', '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 type_meta(iLookTYPE%vegTypeIndex)           = var_info('vegTypeIndex'  , 'index defining vegetation type'             , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 type_meta(iLookTYPE%soilTypeIndex)          = var_info('soilTypeIndex' , 'index defining soil type'                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 type_meta(iLookTYPE%slopeTypeIndex)         = var_info('slopeTypeIndex', 'index defining slope'                       , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 type_meta(iLookTYPE%downHRUindex)           = var_info('downHRUindex'  , 'index of downslope HRU (0 = basin outlet)'  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * site characteristics...
 ! -------------------------
 attr_meta(iLookATTR%latitude)               = var_info('latitude'      , 'latitude'                                              , 'degrees north', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%longitude)              = var_info('longitude'     , 'longitude'                                             , 'degrees east' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%elevation)              = var_info('elevation'     , 'elevation'                                             , 'm'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%tan_slope)              = var_info('tan_slope'     , 'tan water table slope (tan local ground surface slope)', '-'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%contourLength)          = var_info('contourLength' , 'length of contour at downslope edge of HRU'            , 'm'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%HRUarea)                = var_info('HRUarea'       , 'area of each HRU'                                      , 'm2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 attr_meta(iLookATTR%mHeight)                = var_info('mHeight'       , 'measurement height above bare ground'                  , 'm'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * local parameter data...
 ! -------------------------
 ! boundary conditions
 mpar_meta(iLookPARAM%upperBoundHead)        = var_info('upperBoundHead'        , 'matric head at the upper boundary'                                , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%lowerBoundHead)        = var_info('lowerBoundHead'        , 'matric head at the lower boundary'                                , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%upperBoundTheta)       = var_info('upperBoundTheta'       , 'volumetric liquid water content at the upper boundary'            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%lowerBoundTheta)       = var_info('lowerBoundTheta'       , 'volumetric liquid water content at the lower boundary'            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%upperBoundTemp)        = var_info('upperBoundTemp'        , 'temperature of the upper boundary'                                , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%lowerBoundTemp)        = var_info('lowerBoundTemp'        , 'temperature of the lower boundary'                                , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! precipitation partitioning
 mpar_meta(iLookPARAM%tempCritRain)          = var_info('tempCritRain'          , 'critical temperature where precipitation is rain'                 , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%tempRangeTimestep)     = var_info('tempRangeTimestep'     , 'temperature range over the time step'                             , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%frozenPrecipMultip)    = var_info('frozenPrecipMultip'    , 'frozen precipitation multiplier'                                  , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! snow properties
 mpar_meta(iLookPARAM%snowfrz_scale)         = var_info('snowfrz_scale'         , 'scaling parameter for the freezing curve for snow'                , 'K-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%fixedThermalCond_snow) = var_info('fixedThermalCond_snow' , 'temporally constant thermal conductivity for snow'                , 'W m-1 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! snow albedo
 mpar_meta(iLookPARAM%albedoMax)             = var_info('albedoMax'             , 'maximum snow albedo (single spectral band)'                       , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMinWinter)       = var_info('albedoMinWinter'       , 'minimum snow albedo during winter (single spectral band)'         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMinSpring)       = var_info('albedoMinSpring'       , 'minimum snow albedo during spring (single spectral band)'         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMaxVisible)      = var_info('albedoMaxVisible'      , 'maximum snow albedo in the visible part of the spectrum'          , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMinVisible)      = var_info('albedoMinVisible'      , 'minimum snow albedo in the visible part of the spectrum'          , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMaxNearIR)       = var_info('albedoMaxNearIR'       , 'maximum snow albedo in the near infra-red part of the spectrum'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoMinNearIR)       = var_info('albedoMinNearIR'       , 'minimum snow albedo in the near infra-red part of the spectrum'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoDecayRate)       = var_info('albedoDecayRate'       , 'albedo decay rate'                                                , 's'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoSootLoad)        = var_info('albedoSootLoad'        , 'soot load factor'                                                 , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%albedoRefresh)         = var_info('albedoRefresh'         , 'critical mass necessary for albedo refreshment'                   , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! radiation transfer
 mpar_meta(iLookPARAM%radExt_snow)           = var_info('radExt_snow'           , 'extinction coefficient for radiation penetration into snowpack'   , 'm-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%directScale)           = var_info('directScale'           , 'scaling factor for fractional driect radiaion parameterization'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Frad_direct)           = var_info('Frad_direct'           , 'fraction direct solar radiation'                                  , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Frad_vis)              = var_info('Frad_vis'              , 'fraction radiation in visible part of spectrum'                   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! new snow density
 mpar_meta(iLookPARAM%newSnowDenMin)         = var_info('newSnowDenMin'         , 'minimum new snow density'                                         , 'kg m-3'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenMult)        = var_info('newSnowDenMult'        , 'multiplier for new snow density'                                  , 'kg m-3'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenScal)        = var_info('newSnowDenScal'        , 'scaling factor for new snow density'                              , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%constSnowDen)          = var_info('constSnowDen'          , 'Constant new snow density'                                        , 'kg m-3'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenAdd)         = var_info('newSnowDenAdd'         , 'Pahaut 1976, additive factor for new snow density'                , 'kg m-3'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenMultTemp)    = var_info('newSnowDenMultTemp'    , 'Pahaut 1976, multiplier for new snow density for air temperature' , 'kg m-3 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenMultWind)    = var_info('newSnowDenMultWind'    , 'Pahaut 1976, multiplier for new snow density for wind speed'      , 'kg m-7/2 s-1/2'  , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenMultAnd)     = var_info('newSnowDenMultAnd'     , 'Anderson 1976, multiplier for new snow density (Anderson func)'   , 'K-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%newSnowDenBase)        = var_info('newSnowDenBase'        , 'Anderson 1976, base value that is rasied to the (3/2) power'      , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! snow compaction
 mpar_meta(iLookPARAM%densScalGrowth)        = var_info('densScalGrowth'        , 'density scaling factor for grain growth'                          , 'kg-1 m3'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%tempScalGrowth)        = var_info('tempScalGrowth'        , 'temperature scaling factor for grain growth'                      , 'K-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%grainGrowthRate)       = var_info('grainGrowthRate'       , 'rate of grain growth'                                             , 's-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%densScalOvrbdn)        = var_info('densScalOvrbdn'        , 'density scaling factor for overburden pressure'                   , 'kg-1 m3'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%tempScalOvrbdn)        = var_info('tempScalOvrbdn'        , 'temperature scaling factor for overburden pressure'               , 'K-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%baseViscosity )        = var_info('baseViscosity '        , 'viscosity coefficient at T=T_frz and snow density=0'              , 'kg s m-2'        , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! water flow through snow
 mpar_meta(iLookPARAM%Fcapil)                = var_info('Fcapil'                , 'capillary retention (fraction of total pore volume)'              , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%k_snow)                = var_info('k_snow'                , 'hydraulic conductivity of snow'                                   , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%mw_exp)                = var_info('mw_exp'                , 'exponent for meltwater flow'                                      , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! turbulent heat fluxes
 mpar_meta(iLookPARAM%z0Snow)                = var_info('z0Snow'                , 'roughness length of snow'                                         , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%z0Soil)                = var_info('z0Soil'                , 'roughness length of bare soil below the canopy'                   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%z0Canopy)              = var_info('z0Canopy'              , 'roughness length of the canopy'                                   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zpdFraction)           = var_info('zpdFraction'           , 'zero plane displacement / canopy height'                          , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%critRichNumber)        = var_info('critRichNumber'        , 'critical value for the bulk Richardson number'                    , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Louis79_bparam)        = var_info('Louis79_bparam'        , 'parameter in Louis (1979) stability function'                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Louis79_cStar)         = var_info('Louis79_cStar'         , 'parameter in Louis (1979) stability function'                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Mahrt87_eScale)        = var_info('Mahrt87_eScale'        , 'exponential scaling factor in the Mahrt (1987) stability function', '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%leafExchangeCoeff)     = var_info('leafExchangeCoeff'     , 'turbulent exchange coeff between canopy surface and canopy air'   , 'm s-(1/2)'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%windReductionParam)    = var_info('windReductionParam'    , 'canopy wind reduction parameter'                                  , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! stomatal conductance
 mpar_meta(iLookPARAM%Kc25)                  = var_info('Kc25'                  , 'Michaelis-Menten constant for CO2 at 25 degrees C'                , 'umol mol-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Ko25)                  = var_info('Ko25'                  , 'Michaelis-Menten constant for O2 at 25 degrees C'                 , 'mol mol-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Kc_qFac)               = var_info('Kc_qFac'               , 'factor in the q10 function defining temperature controls on Kc'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%Ko_qFac)               = var_info('Ko_qFac'               , 'factor in the q10 function defining temperature controls on Ko'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%kc_Ha)                 = var_info('kc_Ha'                 , 'activation energy for the Michaelis-Menten constant for CO2'      , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%ko_Ha)                 = var_info('ko_Ha'                 , 'activation energy for the Michaelis-Menten constant for O2'       , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax25_canopyTop)     = var_info('vcmax25_canopyTop'     , 'potential carboxylation rate at 25 degrees C at the canopy top'   , 'umol co2 m-2 s-1', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax_qFac)            = var_info('vcmax_qFac'            , 'factor in the q10 function defining temperature controls on vcmax', '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax_Ha)              = var_info('vcmax_Ha'              , 'activation energy in the vcmax function'                          , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax_Hd)              = var_info('vcmax_Hd'              , 'deactivation energy in the vcmax function'                        , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax_Sv)              = var_info('vcmax_Sv'              , 'entropy term in the vcmax function'                               , 'J mol-1 K-1'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vcmax_Kn)              = var_info('vcmax_Kn'              , 'foliage nitrogen decay coefficient'                               , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%jmax25_scale)          = var_info('jmax25_scale'          , 'scaling factor to relate jmax25 to vcmax25'                       , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%jmax_Ha)               = var_info('jmax_Ha'               , 'activation energy in the jmax function'                           , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%jmax_Hd)               = var_info('jmax_Hd'               , 'deactivation energy in the jmax function'                         , 'J mol-1'         , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%jmax_Sv)               = var_info('jmax_Sv'               , 'entropy term in the jmax function'                                , 'J mol-1 K-1'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%fractionJ)             = var_info('fractionJ'             , 'fraction of light lost by other than the chloroplast lamellae'    , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%quantamYield)          = var_info('quantamYield'          , 'quantam yield'                                                    , 'mol e mol-1 q'   , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vpScaleFactor)         = var_info('vpScaleFactor'         , 'vapor pressure scaling factor in stomatal conductance function'   , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%cond2photo_slope)      = var_info('cond2photo_slope'      , 'slope of conductance-photosynthesis relationship'                 , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%minStomatalConductance)= var_info('minStomatalConductance', 'minimum stomatal conductance'                                     , 'umol H2O m-2 s-1', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! vegetation properties
 mpar_meta(iLookPARAM%winterSAI)             = var_info('winterSAI'             , 'stem area index prior to the start of the growing season'         , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%summerLAI)             = var_info('summerLAI'             , 'maximum leaf area index at the peak of the growing season'        , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%rootScaleFactor1)      = var_info('rootScaleFactor1'      , '1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )'    , 'm-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%rootScaleFactor2)      = var_info('rootScaleFactor2'      , '2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )'    , 'm-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%rootingDepth)          = var_info('rootingDepth'          , 'rooting depth'                                                    , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%rootDistExp)           = var_info('rootDistExp'           , 'exponent for the vertical distribution of root density'           , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%plantWiltPsi)          = var_info('plantWiltPsi'          , 'matric head at wilting point'                                     , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%soilStressParam)       = var_info('soilStressParam'       , 'parameter in the exponential soil stress function'                , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%critSoilWilting)       = var_info('critSoilWilting'       , 'critical vol. liq. water content when plants are wilting'         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%critSoilTranspire)     = var_info('critSoilTranspire'     , 'critical vol. liq. water content when transpiration is limited'   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%critAquiferTranspire)  = var_info('critAquiferTranspire'  , 'critical aquifer storage value when transpiration is limited'     , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%minStomatalResistance) = var_info('minStomatalResistance' , 'minimum stomatal resistance'                                      , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%leafDimension)         = var_info('leafDimension'         , 'characteristic leaf dimension'                                    , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%heightCanopyTop)       = var_info('heightCanopyTop'       , 'height of top of the vegetation canopy above ground surface'      , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%heightCanopyBottom)    = var_info('heightCanopyBottom'    , 'height of bottom of the vegetation canopy above ground surface'   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%specificHeatVeg)       = var_info('specificHeatVeg'       , 'specific heat of vegetation'                                      , 'J kg-1 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%maxMassVegetation)     = var_info('maxMassVegetation'     , 'maximum mass of vegetation (full foliage)'                        , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%throughfallScaleSnow)  = var_info('throughfallScaleSnow'  , 'scaling factor for throughfall (snow)'                            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%throughfallScaleRain)  = var_info('throughfallScaleRain'  , 'scaling factor for throughfall (rain)'                            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%refInterceptCapSnow)   = var_info('refInterceptCapSnow'   , 'reference canopy interception capacity per unit leaf area (snow)' , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%refInterceptCapRain)   = var_info('refInterceptCapRain'   , 'canopy interception capacity per unit leaf area (rain)'           , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%snowUnloadingCoeff)    = var_info('snowUnloadingCoeff'    , 'time constant for unloading of snow from the forest canopy'       , 's-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%canopyDrainageCoeff)   = var_info('canopyDrainageCoeff'   , 'time constant for drainage of liquid water from the forest canopy', 's-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%ratioDrip2Unloading)   = var_info('ratioDrip2Unloading'   , 'ratio of canopy drip to unloading of snow from the forest canopy' , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%canopyWettingFactor)   = var_info('canopyWettingFactor'   , 'maximum wetted fraction of the canopy'                            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%canopyWettingExp)      = var_info('canopyWettingExp'      , 'exponent in canopy wetting function'                              , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! soil properties
 mpar_meta(iLookPARAM%soil_dens_intr)        = var_info('soil_dens_intr'        , 'intrinsic soil density'                                           , 'kg m-3'          , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%thCond_soil)           = var_info('thCond_soil'           , 'thermal conductivity of soil (includes quartz and other minerals)', 'W m-1 K-1'       , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%frac_sand)             = var_info('frac_sand'             , 'fraction of sand'                                                 , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%frac_silt)             = var_info('frac_silt'             , 'fraction of silt'                                                 , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%frac_clay)             = var_info('frac_clay'             , 'fraction of clay'                                                 , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%theta_sat)             = var_info('theta_sat'             , 'soil porosity'                                                    , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%theta_res)             = var_info('theta_res'             , 'volumetric residual water content'                                , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vGn_alpha)             = var_info('vGn_alpha'             , 'van Genuchten "alpha" parameter'                                  , 'm-1'             , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%vGn_n)                 = var_info('vGn_n'                 , 'van Genuchten "n" parameter'                                      , '-'               , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%k_soil)                = var_info('k_soil'                , 'saturated hydraulic conductivity'                                 , 'm s-1'           , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%k_macropore)           = var_info('k_macropore'           , 'saturated hydraulic conductivity for macropores'                  , 'm s-1'           , get_ixVarType('parSoil'), lFalseArry, integerMissing, iMissArry)
 ! scalar soil properties
 mpar_meta(iLookPARAM%fieldCapacity)         = var_info('fieldCapacity'         , 'soil field capacity (vol liq water content when baseflow begins)' , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%wettingFrontSuction)   = var_info('wettingFrontSuction'   , 'Green-Ampt wetting front suction'                                 , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%theta_mp)              = var_info('theta_mp'              , 'volumetric liquid water content when macropore flow begins'       , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%mpExp)                 = var_info('mpExp'                 , 'empirical exponent in macropore flow equation'                    , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%kAnisotropic)          = var_info('kAnisotropic'          , 'anisotropy factor for lateral hydraulic conductivity'             , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zScale_TOPMODEL)       = var_info('zScale_TOPMODEL'       , 'TOPMODEL scaling factor used in lower boundary condition for soil', 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%compactedDepth)        = var_info('compactedDepth'        , 'depth where k_soil reaches the compacted value given by CH78'     , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%aquiferScaleFactor)    = var_info('aquiferScaleFactor'    , 'scaling factor for aquifer storage in the big bucket'             , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%aquiferBaseflowExp)    = var_info('aquiferBaseflowExp'    , 'baseflow exponent'                                                , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%qSurfScale)            = var_info('qSurfScale'            , 'scaling factor in the surface runoff parameterization'            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%specificYield)         = var_info('specificYield'         , 'specific yield'                                                   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%specificStorage)       = var_info('specificStorage'       , 'specific storage coefficient'                                     , 'm-1'             , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%f_impede)              = var_info('f_impede'              , 'ice impedence factor'                                             , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%soilIceScale)          = var_info('soilIceScale'          , 'scaling factor for depth of soil ice, used to get frozen fraction', 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%soilIceCV)             = var_info('soilIceCV'             , 'CV of depth of soil ice, used to get frozen fraction'             , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! algorithmic control parameters
 mpar_meta(iLookPARAM%minwind)               = var_info('minwind'               , 'minimum wind speed'                                               , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%minstep)               = var_info('minstep'               , 'minimum length of the time step'                                  , 's'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%maxstep)               = var_info('maxstep'               , 'maximum length of the time step'                                  , 's'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%wimplicit)             = var_info('wimplicit'             , 'weight assigned to the start-of-step fluxes (alpha)'              , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%maxiter)               = var_info('maxiter'               , 'maximum number of iterations'                                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%relConvTol_liquid)     = var_info('relConvTol_liquid'     , 'relative convergence tolerance for vol frac liq water'            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%absConvTol_liquid)     = var_info('absConvTol_liquid'     , 'absolute convergence tolerance for vol frac liq water'            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%relConvTol_matric)     = var_info('relConvTol_matric'     , 'relative convergence tolerance for matric head'                   , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%absConvTol_matric)     = var_info('absConvTol_matric'     , 'absolute convergence tolerance for matric head'                   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%relConvTol_energy)     = var_info('relConvTol_energy'     , 'relative convergence tolerance for energy'                        , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%absConvTol_energy)     = var_info('absConvTol_energy'     , 'absolute convergence tolerance for energy'                        , 'J m-3'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%relConvTol_aquifr)     = var_info('relConvTol_aquifr'     , 'relative convergence tolerance for aquifer storage'               , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%absConvTol_aquifr)     = var_info('absConvTol_aquifr'     , 'absolute convergence tolerance for aquifer storage'               , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmin)                  = var_info('zmin'                  , 'minimum layer depth'                                              , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmax)                  = var_info('zmax'                  , 'maximum layer depth'                                              , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zminLayer1)            = var_info('zminLayer1'            , 'minimum layer depth for the 1st (top) layer'                      , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zminLayer2)            = var_info('zminLayer2'            , 'minimum layer depth for the 2nd layer'                            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zminLayer3)            = var_info('zminLayer3'            , 'minimum layer depth for the 3rd layer'                            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zminLayer4)            = var_info('zminLayer4'            , 'minimum layer depth for the 4th layer'                            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zminLayer5)            = var_info('zminLayer5'            , 'minimum layer depth for the 5th (bottom) layer'                   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer1_lower)      = var_info('zmaxLayer1_lower'      , 'maximum layer depth for the 1st (top) layer when only 1 layer'    , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer2_lower)      = var_info('zmaxLayer2_lower'      , 'maximum layer depth for the 2nd layer when only 2 layers'         , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer3_lower)      = var_info('zmaxLayer3_lower'      , 'maximum layer depth for the 3rd layer when only 3 layers'         , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer4_lower)      = var_info('zmaxLayer4_lower'      , 'maximum layer depth for the 4th layer when only 4 layers'         , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer1_upper)      = var_info('zmaxLayer1_upper'      , 'maximum layer depth for the 1st (top) layer when > 1 layer'       , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer2_upper)      = var_info('zmaxLayer2_upper'      , 'maximum layer depth for the 2nd layer when > 2 layers'            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer3_upper)      = var_info('zmaxLayer3_upper'      , 'maximum layer depth for the 3rd layer when > 3 layers'            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 mpar_meta(iLookPARAM%zmaxLayer4_upper)      = var_info('zmaxLayer4_upper'      , 'maximum layer depth for the 4th layer when > 4 layers'            , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * basin parameter data...
 ! -------------------------
 bpar_meta(iLookBPAR%basin__aquiferHydCond)           = var_info('basin__aquiferHydCond'    , 'hydraulic conductivity of the aquifer'                          , 'm s-1', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bpar_meta(iLookBPAR%basin__aquiferScaleFactor)       = var_info('basin__aquiferScaleFactor', 'scaling factor for aquifer storage in the big bucket'           , 'm'    , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bpar_meta(iLookBPAR%basin__aquiferBaseflowExp)       = var_info('basin__aquiferBaseflowExp', 'baseflow exponent for the big bucket'                           , '-'    , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bpar_meta(iLookBPAR%routingGammaShape)               = var_info('routingGammaShape'        , 'shape parameter in Gamma distribution used for sub-grid routing', '-'    , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bpar_meta(iLookBPAR%routingGammaScale)               = var_info('routingGammaScale'        , 'scale parameter in Gamma distribution used for sub-grid routing', 's'    , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * local model prognostic (state) variables...
 ! ---------------------------------------------
 ! define variables for time stepping
 prog_meta(iLookPROG%dt_init)                         = var_info('dt_init'                        , 'length of initial time step at start of next data interval'       , 's'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! state variables for vegetation
 prog_meta(iLookPROG%scalarCanopyIce)                 = var_info('scalarCanopyIce'                , 'mass of ice on the vegetation canopy'                             , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarCanopyLiq)                 = var_info('scalarCanopyLiq'                , 'mass of liquid water on the vegetation canopy'                    , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarCanopyWat)                 = var_info('scalarCanopyWat'                , 'mass of total water on the vegetation canopy'                     , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarCanairTemp)                = var_info('scalarCanairTemp'               , 'temperature of the canopy air space'                              , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarCanopyTemp)                = var_info('scalarCanopyTemp'               , 'temperature of the vegetation canopy'                             , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! state variables for snow
 prog_meta(iLookPROG%spectralSnowAlbedoDiffuse)       = var_info('spectralSnowAlbedoDiffuse'      , 'diffuse snow albedo for individual spectral bands'                , '-'               , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarSnowAlbedo)                = var_info('scalarSnowAlbedo'               , 'snow albedo for the entire spectral band'                         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarSnowDepth)                 = var_info('scalarSnowDepth'                , 'total snow depth'                                                 , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarSWE)                       = var_info('scalarSWE'                      , 'snow water equivalent'                                            , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarSfcMeltPond)               = var_info('scalarSfcMeltPond'              , 'ponded water caused by melt of the "snow without a layer"'        , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! define state variables for the snow+soil domain
 prog_meta(iLookPROG%mLayerTemp)                      = var_info('mLayerTemp'                     , 'temperature of each layer'                                        , 'K'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%mLayerVolFracIce)                = var_info('mLayerVolFracIce'               , 'volumetric fraction of ice in each layer'                         , '-'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%mLayerVolFracLiq)                = var_info('mLayerVolFracLiq'               , 'volumetric fraction of liquid water in each layer'                , '-'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%mLayerVolFracWat)                = var_info('mLayerVolFracWat'               , 'volumetric fraction of total water in each layer'                 , '-'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%mLayerMatricHead)                = var_info('mLayerMatricHead'               , 'matric head of water in the soil'                                 , 'm'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 ! other state variables
 prog_meta(iLookPROG%scalarAquiferStorage)            = var_info('scalarAquiferStorage'           , 'water required to bring aquifer to the bottom of the soil profile', 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%scalarSurfaceTemp)               = var_info('scalarSurfaceTemp'              , 'surface temperature (just a copy of the upper-layer temperature)' , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! define coordinate variables
 prog_meta(iLookPROG%mLayerDepth)                     = var_info('mLayerDepth'                    , 'depth of each layer'                                              , 'm'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%mLayerHeight)                    = var_info('mLayerHeight'                   , 'height of the layer mid-point (top of soil = 0)'                  , 'm'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 prog_meta(iLookPROG%iLayerHeight)                    = var_info('iLayerHeight'                   , 'height of the layer interface (top of soil = 0)'                  , 'm'               , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * local model diagnostic variables...
 ! -------------------------------------
 ! local properties
 diag_meta(iLookDIAG%scalarCanopyDepth)               = var_info('scalarCanopyDepth'              , 'canopy depth'                                                     , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarGreenVegFraction)          = var_info('scalarGreenVegFraction'         , 'green vegetation fraction (used to compute LAI)'                  , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarBulkVolHeatCapVeg)         = var_info('scalarBulkVolHeatCapVeg'        , 'bulk volumetric heat capacity of vegetation'                      , 'J m-3 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyEmissivity)          = var_info('scalarCanopyEmissivity'         , 'effective canopy emissivity'                                      , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarRootZoneTemp)              = var_info('scalarRootZoneTemp'             , 'average temperature of the root zone'                             , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarLAI)                       = var_info('scalarLAI'                      , 'one-sided leaf area index'                                        , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSAI)                       = var_info('scalarSAI'                      , 'one-sided stem area index'                                        , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarExposedLAI)                = var_info('scalarExposedLAI'               , 'exposed leaf area index (after burial by snow)'                   , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarExposedSAI)                = var_info('scalarExposedSAI'               , 'exposed stem area index (after burial by snow)'                   , 'm2 m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyIceMax)              = var_info('scalarCanopyIceMax'             , 'maximum interception storage capacity for ice'                    , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyLiqMax)              = var_info('scalarCanopyLiqMax'             , 'maximum interception storage capacity for liquid water'           , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarGrowingSeasonIndex)        = var_info('scalarGrowingSeasonIndex'       , 'growing season index (0=off, 1=on)'                               , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVolHtCap_air)              = var_info('scalarVolHtCap_air'             , 'volumetric heat capacity air'                                     , 'J m-3 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVolHtCap_ice)              = var_info('scalarVolHtCap_ice'             , 'volumetric heat capacity ice'                                     , 'J m-3 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVolHtCap_soil)             = var_info('scalarVolHtCap_soil'            , 'volumetric heat capacity dry soil'                                , 'J m-3 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVolHtCap_water)            = var_info('scalarVolHtCap_water'           , 'volumetric heat capacity liquid wat'                              , 'J m-3 K-1'       , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerVolHtCapBulk)              = var_info('mLayerVolHtCapBulk'             , 'volumetric heat capacity in each layer'                           , 'J m-3 K-1'       , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarLambda_drysoil)            = var_info('scalarLambda_drysoil'           , 'thermal conductivity of dry soil'                                 , 'W m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarLambda_wetsoil)            = var_info('scalarLambda_wetsoil'           , 'thermal conductivity of wet soil'                                 , 'W m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerThermalC)                  = var_info('mLayerThermalC'                 , 'thermal conductivity at the mid-point of each layer'              , 'W m-1 K-1'       , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%iLayerThermalC)                  = var_info('iLayerThermalC'                 , 'thermal conductivity at the interface of each layer'              , 'W m-1 K-1'       , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 ! forcing
 diag_meta(iLookDIAG%scalarVPair)                     = var_info('scalarVPair'                    , 'vapor pressure of the air above the vegetation canopy'            , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVP_CanopyAir)              = var_info('scalarVP_CanopyAir'             , 'vapor pressure of the canopy air space'                           , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarTwetbulb)                  = var_info('scalarTwetbulb'                 , 'wet bulb temperature'                                             , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSnowfallTemp)              = var_info('scalarSnowfallTemp'             , 'temperature of fresh snow'                                        , 'K'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarNewSnowDensity)            = var_info('scalarNewSnowDensity'           , 'density of fresh snow'                                            , 'kg m-3'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarO2air)                     = var_info('scalarO2air'                    , 'atmospheric o2 concentration'                                     , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCO2air)                    = var_info('scalarCO2air'                   , 'atmospheric co2 concentration'                                    , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! shortwave radiation
 diag_meta(iLookDIAG%scalarCosZenith)                 = var_info('scalarCosZenith'                , 'cosine of the solar zenith angle'                                 , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarFractionDirect)            = var_info('scalarFractionDirect'           , 'fraction of direct radiation (0-1)'                               , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopySunlitFraction)      = var_info('scalarCanopySunlitFraction'     , 'sunlit fraction of canopy'                                        , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopySunlitLAI)           = var_info('scalarCanopySunlitLAI'          , 'sunlit leaf area'                                                 , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyShadedLAI)           = var_info('scalarCanopyShadedLAI'          , 'shaded leaf area'                                                 , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%spectralAlbGndDirect)            = var_info('spectralAlbGndDirect'           , 'direct  albedo of underlying surface for each spectral band'      , '-'               , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%spectralAlbGndDiffuse)           = var_info('spectralAlbGndDiffuse'          , 'diffuse albedo of underlying surface for each spectral band'      , '-'               , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarGroundAlbedo)              = var_info('scalarGroundAlbedo'             , 'albedo of the ground surface'                                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! turbulent heat transfer
 diag_meta(iLookDIAG%scalarLatHeatSubVapCanopy)       = var_info('scalarLatHeatSubVapCanopy'      , 'latent heat of sublimation/vaporization used for veg canopy'      , 'J kg-1'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarLatHeatSubVapGround)       = var_info('scalarLatHeatSubVapGround'      , 'latent heat of sublimation/vaporization used for ground surface'  , 'J kg-1'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSatVP_CanopyTemp)          = var_info('scalarSatVP_CanopyTemp'         , 'saturation vapor pressure at the temperature of vegetation canopy', 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSatVP_GroundTemp)          = var_info('scalarSatVP_GroundTemp'         , 'saturation vapor pressure at the temperature of the ground'       , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarZ0Canopy)                  = var_info('scalarZ0Canopy'                 , 'roughness length of the canopy'                                   , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarWindReductionFactor)       = var_info('scalarWindReductionFactor'      , 'canopy wind reduction factor'                                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarZeroPlaneDisplacement)     = var_info('scalarZeroPlaneDisplacement'    , 'zero plane displacement'                                          , 'm'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarRiBulkCanopy)              = var_info('scalarRiBulkCanopy'             , 'bulk Richardson number for the canopy'                            , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarRiBulkGround)              = var_info('scalarRiBulkGround'             , 'bulk Richardson number for the ground surface'                    , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyStabilityCorrection) = var_info('scalarCanopyStabilityCorrection', 'stability correction for the canopy'                              , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarGroundStabilityCorrection) = var_info('scalarGroundStabilityCorrection', 'stability correction for the ground surface'                      , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! evapotranspiration
 diag_meta(iLookDIAG%scalarIntercellularCO2Sunlit)    = var_info('scalarIntercellularCO2Sunlit'   , 'carbon dioxide partial pressure of leaf interior (sunlit leaves)' , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarIntercellularCO2Shaded)    = var_info('scalarIntercellularCO2Shaded'   , 'carbon dioxide partial pressure of leaf interior (shaded leaves)' , 'Pa'              , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarTranspireLim)              = var_info('scalarTranspireLim'             , 'aggregate soil moisture and aquifer control on transpiration'     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarTranspireLimAqfr)          = var_info('scalarTranspireLimAqfr'         , 'aquifer storage control on transpiration'                         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarFoliageNitrogenFactor)     = var_info('scalarFoliageNitrogenFactor'    , 'foliage nitrogen concentration (1=saturated)'                     , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSoilRelHumidity)           = var_info('scalarSoilRelHumidity'          , 'relative humidity in the soil pores in the upper-most soil layer' , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerTranspireLim)              = var_info('mLayerTranspireLim'             , 'soil moist & veg limit on transpiration for each layer'           , '-'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerRootDensity)               = var_info('mLayerRootDensity'              , 'fraction of roots in each soil layer'                             , '-'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarAquiferRootFrac)           = var_info('scalarAquiferRootFrac'          , 'fraction of roots below the soil profile (in the aquifer)'        , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! canopy hydrology
 diag_meta(iLookDIAG%scalarFracLiqVeg)                = var_info('scalarFracLiqVeg'               , 'fraction of liquid water on vegetation'                           , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarCanopyWetFraction)         = var_info('scalarCanopyWetFraction'        , 'fraction canopy that is wet'                                      , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! snow hydrology
 diag_meta(iLookDIAG%scalarSnowAge)                   = var_info('scalarSnowAge'                  , 'non-dimensional snow age'                                         , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarGroundSnowFraction)        = var_info('scalarGroundSnowFraction'       , 'fraction ground that is covered with snow'                        , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%spectralSnowAlbedoDirect)        = var_info('spectralSnowAlbedoDirect'       , 'direct snow albedo for individual spectral bands'                 , '-'               , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerFracLiqSnow)               = var_info('mLayerFracLiqSnow'              , 'fraction of liquid water in each snow layer'                      , '-'               , get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerThetaResid)                = var_info('mLayerThetaResid'               , 'residual volumetric water content in each snow layer'             , '-'               , get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerPoreSpace)                 = var_info('mLayerPoreSpace'                , 'total pore space in each snow layer'                              , '-'               , get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerMeltFreeze)                = var_info('mLayerMeltFreeze'               , 'ice content change from melt/freeze in each layer'                , 'kg m-3'          , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! soil hydrology
 diag_meta(iLookDIAG%scalarInfilArea)                 = var_info('scalarInfilArea'                , 'fraction of unfrozen area where water can infiltrate'             , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarFrozenArea)                = var_info('scalarFrozenArea'               , 'fraction of area that is considered impermeable due to soil ice'  , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSoilControl)               = var_info('scalarSoilControl'              , 'soil control on infiltration (1=controlling; 0=not)'              , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerVolFracAir)                = var_info('mLayerVolFracAir'               , 'volumetric fraction of air in each layer'                         , '-'               , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerTcrit)                     = var_info('mLayerTcrit'                    , 'critical soil temperature above which all water is unfrozen'      , 'K'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerCompress)                  = var_info('mLayerCompress'                 , 'change in volumetric water content due to compression of soil'    , '-'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarSoilCompress)              = var_info('scalarSoilCompress'             , 'change in total soil storage due to compression of soil matrix'   , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%mLayerMatricHeadLiq)             = var_info('mLayerMatricHeadLiq'            , 'matric potential of liquid water'                                 , 'm'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 ! mass balance check
 diag_meta(iLookDIAG%scalarSoilWatBalError)           = var_info('scalarSoilWatBalError'          , 'error in the total soil water balance'                            , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarAquiferBalError)           = var_info('scalarAquiferBalError'          , 'error in the aquifer water balance'                               , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarTotalSoilLiq)              = var_info('scalarTotalSoilLiq'             , 'total mass of liquid water in the soil'                           , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarTotalSoilIce)              = var_info('scalarTotalSoilIce'             , 'total mass of ice in the soil'                                    , 'kg m-2'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! variable shortcuts
 diag_meta(iLookDIAG%scalarVGn_m)                     = var_info('scalarVGn_m'                    , 'van Genuchten "m" parameter'                                      , '-'               , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarKappa)                     = var_info('scalarKappa'                    , 'constant in the freezing curve function'                          , 'm K-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 diag_meta(iLookDIAG%scalarVolLatHt_fus)              = var_info('scalarVolLatHt_fus'             , 'volumetric latent heat of fusion'                                 , 'J m-3'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! number of function evaluations
 diag_meta(iLookDIAG%numFluxCalls)                    = var_info('numFluxCalls'                   , 'number of flux calls'                                             , '-'               , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry) 

 ! -----
 ! * local model fluxes...
 ! -----------------------
 ! net energy and mass fluxes for the vegetation domain
 flux_meta(iLookFLUX%scalarCanairNetNrgFlux)          = var_info('scalarCanairNetNrgFlux'         , 'net energy flux for the canopy air space'                         , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyNetNrgFlux)          = var_info('scalarCanopyNetNrgFlux'         , 'net energy flux for the vegetation canopy'                        , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarGroundNetNrgFlux)          = var_info('scalarGroundNetNrgFlux'         , 'net energy flux for the ground surface'                           , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyNetLiqFlux)          = var_info('scalarCanopyNetLiqFlux'         , 'net liquid water flux for the vegetation canopy'                  , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! forcing
 flux_meta(iLookFLUX%scalarRainfall)                  = var_info('scalarRainfall'                 , 'computed rainfall rate'                                           , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSnowfall)                  = var_info('scalarSnowfall'                 , 'computed snowfall rate'                                           , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! shortwave radiation
 flux_meta(iLookFLUX%spectralIncomingDirect)          = var_info('spectralIncomingDirect'         , 'incoming direct solar radiation in each wave band'                , 'W m-2'           , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%spectralIncomingDiffuse)         = var_info('spectralIncomingDiffuse'        , 'incoming diffuse solar radiation in each wave band'               , 'W m-2'           , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopySunlitPAR)           = var_info('scalarCanopySunlitPAR'          , 'average absorbed par for sunlit leaves'                           , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyShadedPAR)           = var_info('scalarCanopyShadedPAR'          , 'average absorbed par for shaded leaves'                           , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%spectralBelowCanopyDirect)       = var_info('spectralBelowCanopyDirect'      , 'downward direct flux below veg layer for each spectral band'      , 'W m-2'           , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%spectralBelowCanopyDiffuse)      = var_info('spectralBelowCanopyDiffuse'     , 'downward diffuse flux below veg layer for each spectral band'     , 'W m-2'           , get_ixVarType('wLength'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarBelowCanopySolar)          = var_info('scalarBelowCanopySolar'         , 'solar radiation transmitted below the canopy'                     , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyAbsorbedSolar)       = var_info('scalarCanopyAbsorbedSolar'      , 'solar radiation absorbed by canopy'                               , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarGroundAbsorbedSolar)       = var_info('scalarGroundAbsorbedSolar'      , 'solar radiation absorbed by ground'                               , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! longwave radiation
 flux_meta(iLookFLUX%scalarLWRadCanopy)               = var_info('scalarLWRadCanopy'              , 'longwave radiation emitted from the canopy'                       , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadGround)               = var_info('scalarLWRadGround'              , 'longwave radiation emitted at the ground surface'                 , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadUbound2Canopy)        = var_info('scalarLWRadUbound2Canopy'       , 'downward atmospheric longwave radiation absorbed by the canopy'   , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadUbound2Ground)        = var_info('scalarLWRadUbound2Ground'       , 'downward atmospheric longwave radiation absorbed by the ground'   , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadUbound2Ubound)        = var_info('scalarLWRadUbound2Ubound'       , 'atmospheric radiation refl by ground + lost thru upper boundary'  , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadCanopy2Ubound)        = var_info('scalarLWRadCanopy2Ubound'       , 'longwave radiation emitted from canopy lost thru upper boundary'  , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadCanopy2Ground)        = var_info('scalarLWRadCanopy2Ground'       , 'longwave radiation emitted from canopy absorbed by the ground'    , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadCanopy2Canopy)        = var_info('scalarLWRadCanopy2Canopy'       , 'canopy longwave reflected from ground and absorbed by the canopy' , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadGround2Ubound)        = var_info('scalarLWRadGround2Ubound'       , 'longwave radiation emitted from ground lost thru upper boundary'  , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWRadGround2Canopy)        = var_info('scalarLWRadGround2Canopy'       , 'longwave radiation emitted from ground and absorbed by the canopy', 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWNetCanopy)               = var_info('scalarLWNetCanopy'              , 'net longwave radiation at the canopy'                             , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWNetGround)               = var_info('scalarLWNetGround'              , 'net longwave radiation at the ground surface'                     , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLWNetUbound)               = var_info('scalarLWNetUbound'              , 'net longwave radiation at the upper atmospheric boundary'         , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! turbulent heat transfer
 flux_meta(iLookFLUX%scalarEddyDiffusCanopyTop)       = var_info('scalarEddyDiffusCanopyTop'      , 'eddy diffusivity for heat at the top of the canopy'               , 'm2 s-1'          , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarFrictionVelocity)          = var_info('scalarFrictionVelocity'         , 'friction velocity (canopy momentum sink)'                         , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarWindspdCanopyTop)          = var_info('scalarWindspdCanopyTop'         , 'windspeed at the top of the canopy'                               , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarWindspdCanopyBottom)       = var_info('scalarWindspdCanopyBottom'      , 'windspeed at the height of the bottom of the canopy'              , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarGroundResistance)          = var_info('scalarGroundResistance'         , 'below canopy aerodynamic resistance'                              , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyResistance)          = var_info('scalarCanopyResistance'         , 'above canopy aerodynamic resistance'                              , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLeafResistance)            = var_info('scalarLeafResistance'           , 'mean leaf boundary layer resistance per unit leaf area'           , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSoilResistance)            = var_info('scalarSoilResistance'           , 'soil surface resistance'                                          , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSenHeatTotal)              = var_info('scalarSenHeatTotal'             , 'sensible heat from the canopy air space to the atmosphere'        , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSenHeatCanopy)             = var_info('scalarSenHeatCanopy'            , 'sensible heat from the canopy to the canopy air space'            , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSenHeatGround)             = var_info('scalarSenHeatGround'            , 'sensible heat from the ground (below canopy or non-vegetated)'    , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLatHeatTotal)              = var_info('scalarLatHeatTotal'             , 'latent heat from the canopy air space to the atmosphere'          , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLatHeatCanopyEvap)         = var_info('scalarLatHeatCanopyEvap'        , 'evaporation latent heat from the canopy to the canopy air space'  , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLatHeatCanopyTrans)        = var_info('scalarLatHeatCanopyTrans'       , 'transpiration latent heat from the canopy to the canopy air space', 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarLatHeatGround)             = var_info('scalarLatHeatGround'            , 'latent heat from the ground (below canopy or non-vegetated)'      , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyAdvectiveHeatFlux)   = var_info('scalarCanopyAdvectiveHeatFlux'  , 'heat advected to the canopy with precipitation (snow + rain)'     , 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarGroundAdvectiveHeatFlux)   = var_info('scalarGroundAdvectiveHeatFlux'  , 'heat advected to the ground with throughfall + unloading/drainage', 'W m-2'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopySublimation)         = var_info('scalarCanopySublimation'        , 'canopy sublimation/frost'                                         , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSnowSublimation)           = var_info('scalarSnowSublimation'          , 'snow sublimation/frost (below canopy or non-vegetated)'           , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! liquid water fluxes associated with evapotranspiration
 flux_meta(iLookFLUX%scalarStomResistSunlit)          = var_info('scalarStomResistSunlit'         , 'stomatal resistance for sunlit leaves'                            , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarStomResistShaded)          = var_info('scalarStomResistShaded'         , 'stomatal resistance for shaded leaves'                            , 's m-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarPhotosynthesisSunlit)      = var_info('scalarPhotosynthesisSunlit'     , 'sunlit photosynthesis'                                            , 'umolco2 m-2 s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarPhotosynthesisShaded)      = var_info('scalarPhotosynthesisShaded'     , 'shaded photosynthesis'                                            , 'umolco2 m-2 s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyTranspiration)       = var_info('scalarCanopyTranspiration'      , 'canopy transpiration'                                             , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyEvaporation)         = var_info('scalarCanopyEvaporation'        , 'canopy evaporation/condensation'                                  , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarGroundEvaporation)         = var_info('scalarGroundEvaporation'        , 'ground evaporation/condensation (below canopy or non-vegetated)'  , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerTranspire)                 = var_info('mLayerTranspire'                , 'transpiration loss from each soil layer'                          , 'm s-1'           , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 ! liquid and solid water fluxes through the canopy
 flux_meta(iLookFLUX%scalarThroughfallSnow)           = var_info('scalarThroughfallSnow'          , 'snow that reaches the ground without ever touching the canopy'    , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarThroughfallRain)           = var_info('scalarThroughfallRain'          , 'rain that reaches the ground without ever touching the canopy'    , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopySnowUnloading)       = var_info('scalarCanopySnowUnloading'      , 'unloading of snow from the vegetation canopy'                     , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyLiqDrainage)         = var_info('scalarCanopyLiqDrainage'        , 'drainage of liquid water from the vegetation canopy'              , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarCanopyMeltFreeze)          = var_info('scalarCanopyMeltFreeze'         , 'melt/freeze of water stored in the canopy'                        , 'kg m-2 s-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! energy fluxes and for the snow and soil domains
 flux_meta(iLookFLUX%iLayerConductiveFlux)            = var_info('iLayerConductiveFlux'           , 'conductive energy flux at layer interfaces'                       , 'W m-2'           , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%iLayerAdvectiveFlux)             = var_info('iLayerAdvectiveFlux'            , 'advective energy flux at layer interfaces'                        , 'W m-2'           , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%iLayerNrgFlux)                   = var_info('iLayerNrgFlux'                  , 'energy flux at layer interfaces'                                  , 'W m-2'           , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerNrgFlux)                   = var_info('mLayerNrgFlux'                  , 'net energy flux for each layer within the snow+soil domain'       , 'J m-3 s-1'       , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! liquid water fluxes for the snow domain
 flux_meta(iLookFLUX%scalarSnowDrainage)              = var_info('scalarSnowDrainage'             , 'drainage from the bottom of the snow profile'                     , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%iLayerLiqFluxSnow)               = var_info('iLayerLiqFluxSnow'              , 'liquid flux at snow layer interfaces'                             , 'm s-1'           , get_ixVarType('ifcSnow'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerLiqFluxSnow)               = var_info('mLayerLiqFluxSnow'              , 'net liquid water flux for each snow layer'                        , 's-1'             , get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 ! liquid water fluxes for the soil domain
 flux_meta(iLookFLUX%scalarRainPlusMelt)              = var_info('scalarRainPlusMelt'             , 'rain plus melt, used as input to soil before surface runoff'      , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarMaxInfilRate)              = var_info('scalarMaxInfilRate'             , 'maximum infiltration rate'                                        , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarInfiltration)              = var_info('scalarInfiltration'             , 'infiltration of water into the soil profile'                      , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarExfiltration)              = var_info('scalarExfiltration'             , 'exfiltration of water from the top of the soil profile'           , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSurfaceRunoff)             = var_info('scalarSurfaceRunoff'            , 'surface runoff'                                                   , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerSatHydCondMP)              = var_info('mLayerSatHydCondMP'             , 'saturated hydraulic conductivity of macropores in each layer'     , 'm s-1'           , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerSatHydCond)                = var_info('mLayerSatHydCond'               , 'saturated hydraulic conductivity in each layer'                   , 'm s-1'           , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%iLayerSatHydCond)                = var_info('iLayerSatHydCond'               , 'saturated hydraulic conductivity in each layer interface'         , 'm s-1'           , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerHydCond)                   = var_info('mLayerHydCond'                  , 'hydraulic conductivity in each layer'                             , 'm s-1'           , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%iLayerLiqFluxSoil)               = var_info('iLayerLiqFluxSoil'              , 'liquid flux at soil layer interfaces'                             , 'm s-1'           , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerLiqFluxSoil)               = var_info('mLayerLiqFluxSoil'              , 'net liquid water flux for each soil layer'                        , 's-1'             , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerBaseflow)                  = var_info('mLayerBaseflow'                 , 'baseflow from each soil layer'                                    , 'm s-1'           , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerColumnInflow)              = var_info('mLayerColumnInflow'             , 'total inflow to each layer in a given soil column'                , 'm3 s-1'          , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%mLayerColumnOutflow)             = var_info('mLayerColumnOutflow'            , 'total outflow from each layer in a given soil column'             , 'm3 s-1'          , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSoilBaseflow)              = var_info('scalarSoilBaseflow'             , 'total baseflow from the soil profile'                             , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarSoilDrainage)              = var_info('scalarSoilDrainage'             , 'drainage from the bottom of the soil profile'                     , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarAquiferRecharge)           = var_info('scalarAquiferRecharge'          , 'recharge to the aquifer'                                          , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarAquiferTranspire)          = var_info('scalarAquiferTranspire'         , 'transpiration loss from the aquifer'                              , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 flux_meta(iLookFLUX%scalarAquiferBaseflow)           = var_info('scalarAquiferBaseflow'          , 'baseflow from the aquifer'                                        , 'm s-1'           , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * local flux derivatives...
 ! ---------------------------
 ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
 deriv_meta(iLookDERIV%dCanairNetFlux_dCanairTemp)    = var_info('dCanairNetFlux_dCanairTemp'   , 'derivative in net canopy air space flux w.r.t. canopy air temperature', 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanairNetFlux_dCanopyTemp)    = var_info('dCanairNetFlux_dCanopyTemp'   , 'derivative in net canopy air space flux w.r.t. canopy temperature'    , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanairNetFlux_dGroundTemp)    = var_info('dCanairNetFlux_dGroundTemp'   , 'derivative in net canopy air space flux w.r.t. ground temperature'    , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyNetFlux_dCanairTemp)    = var_info('dCanopyNetFlux_dCanairTemp'   , 'derivative in net canopy flux w.r.t. canopy air temperature'          , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyNetFlux_dCanopyTemp)    = var_info('dCanopyNetFlux_dCanopyTemp'   , 'derivative in net canopy flux w.r.t. canopy temperature'              , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyNetFlux_dGroundTemp)    = var_info('dCanopyNetFlux_dGroundTemp'   , 'derivative in net canopy flux w.r.t. ground temperature'              , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyNetFlux_dCanLiq)        = var_info('dCanopyNetFlux_dCanLiq'       , 'derivative in net canopy fluxes w.r.t. canopy liquid water content'   , 'J kg-1 s-1'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundNetFlux_dCanairTemp)    = var_info('dGroundNetFlux_dCanairTemp'   , 'derivative in net ground flux w.r.t. canopy air temperature'          , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundNetFlux_dCanopyTemp)    = var_info('dGroundNetFlux_dCanopyTemp'   , 'derivative in net ground flux w.r.t. canopy temperature'              , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundNetFlux_dGroundTemp)    = var_info('dGroundNetFlux_dGroundTemp'   , 'derivative in net ground flux w.r.t. ground temperature'              , 'W m-2 K-1'      , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundNetFlux_dCanLiq)        = var_info('dGroundNetFlux_dCanLiq'       , 'derivative in net ground fluxes w.r.t. canopy liquid water content'   , 'J kg-1 s-1'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! derivatives in evaporative fluxes w.r.t. relevant state variables
 deriv_meta(iLookDERIV%dCanopyEvaporation_dTCanair)   = var_info('dCanopyEvaporation_dTCanair'  , 'derivative in canopy evaporation w.r.t. canopy air temperature'       , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyEvaporation_dTCanopy)   = var_info('dCanopyEvaporation_dTCanopy'  , 'derivative in canopy evaporation w.r.t. canopy temperature'           , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyEvaporation_dTGround)   = var_info('dCanopyEvaporation_dTGround'  , 'derivative in canopy evaporation w.r.t. ground temperature'           , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanopyEvaporation_dCanLiq)    = var_info('dCanopyEvaporation_dCanLiq'   , 'derivative in canopy evaporation w.r.t. canopy liquid water content'  , 's-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundEvaporation_dTCanair)   = var_info('dGroundEvaporation_dTCanair'  , 'derivative in ground evaporation w.r.t. canopy air temperature'       , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundEvaporation_dTCanopy)   = var_info('dGroundEvaporation_dTCanopy'  , 'derivative in ground evaporation w.r.t. canopy temperature'           , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundEvaporation_dTGround)   = var_info('dGroundEvaporation_dTGround'  , 'derivative in ground evaporation w.r.t. ground temperature'           , 'kg m-2 s-1 K-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dGroundEvaporation_dCanLiq)    = var_info('dGroundEvaporation_dCanLiq'   , 'derivative in ground evaporation w.r.t. canopy liquid water content'  , 's-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! derivatives in canopy water w.r.t canopy temperature
 deriv_meta(iLookDERIV%dTheta_dTkCanopy)              = var_info('dTheta_dTkCanopy'             , 'derivative of volumetric liquid water content w.r.t. temperature'     , 'K-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCanLiq_dTcanopy)              = var_info('dCanLiq_dTcanopy'             , 'derivative of canopy liquid storage w.r.t. temperature'               , 'kg m-2 K-1'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! derivatives in canopy liquid fluxes w.r.t. canopy water
 deriv_meta(iLookDERIV%scalarCanopyLiqDeriv)          = var_info('scalarCanopyLiqDeriv'         , 'derivative in (throughfall + drainage) w.r.t. canopy liquid water'    , 's-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%scalarThroughfallRainDeriv)    = var_info('scalarThroughfallRainDeriv'   , 'derivative in throughfall w.r.t. canopy liquid water'                 , 's-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%scalarCanopyLiqDrainageDeriv)  = var_info('scalarCanopyLiqDrainageDeriv' , 'derivative in canopy drainage w.r.t. canopy liquid water'             , 's-1'            , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
 deriv_meta(iLookDERIV%dNrgFlux_dTempAbove)           = var_info('dNrgFlux_dTempAbove'          , 'derivatives in the flux w.r.t. temperature in the layer above'        , 'J m-2 s-1 K-1'  , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dNrgFlux_dTempBelow)           = var_info('dNrgFlux_dTempBelow'          , 'derivatives in the flux w.r.t. temperature in the layer below'        , 'J m-2 s-1 K-1'  , get_ixVarType('ifcToto'), lFalseArry, integerMissing, iMissArry)
 ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
 deriv_meta(iLookDERIV%iLayerLiqFluxSnowDeriv)        = var_info('iLayerLiqFluxSnowDeriv'       , 'derivative in vertical liquid water flux at layer interfaces'         , 'm s-1'          , get_ixVarType('ifcSnow'), lFalseArry, integerMissing, iMissArry)
 ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
 deriv_meta(iLookDERIV%dVolTot_dPsi0)                 = var_info('dVolTot_dPsi0'                , 'derivative in total water content w.r.t. total water matric potential', 'm-1'            , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dCompress_dPsi)                = var_info('dCompress_dPsi'               , 'derivative in compressibility w.r.t matric head'                      , 'm-1'            , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%mLayerdTheta_dPsi)             = var_info('mLayerdTheta_dPsi'            , 'derivative in the soil water characteristic w.r.t. psi'               , 'm-1'            , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%mLayerdPsi_dTheta)             = var_info('mLayerdPsi_dTheta'            , 'derivative in the soil water characteristic w.r.t. theta'             , 'm'              , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dq_dHydStateAbove)             = var_info('dq_dHydStateAbove'            , 'change in flux at layer interfaces w.r.t. states in the layer above'  , 'unknown'        , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dq_dHydStateBelow)             = var_info('dq_dHydStateBelow'            , 'change in flux at layer interfaces w.r.t. states in the layer below'  , 'unknown'        , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
 deriv_meta(iLookDERIV%dq_dNrgStateAbove)             = var_info('dq_dNrgStateAbove'            , 'change in flux at layer interfaces w.r.t. states in the layer above'  , 'unknown'        , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dq_dNrgStateBelow)             = var_info('dq_dNrgStateBelow'            , 'change in flux at layer interfaces w.r.t. states in the layer below'  , 'unknown'        , get_ixVarType('ifcSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%mLayerdTheta_dTk)              = var_info('mLayerdTheta_dTk'             , 'derivative of volumetric liquid water content w.r.t. temperature'     , 'K-1'            , get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dPsiLiq_dTemp)                 = var_info('dPsiLiq_dTemp'                , 'derivative in the liquid water matric potential w.r.t. temperature'   , 'm K-1'          , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 deriv_meta(iLookDERIV%dPsiLiq_dPsi0)                 = var_info('dPsiLiq_dPsi0'                , 'derivative in liquid matric potential w.r.t. total  matric potential' , '-'              , get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * basin-wide runoff and aquifer fluxes...
 ! -----------------------------------------
 bvar_meta(iLookBVAR%basin__totalArea)        = var_info('basin__TotalArea'       , 'total basin area'                                       , 'm2'    , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__SurfaceRunoff)    = var_info('basin__SurfaceRunoff'   , 'surface runoff'                                         , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__ColumnOutflow)    = var_info('basin__ColumnOutflow'   , 'outflow from all "outlet" HRUs (with no downstream HRU)', 'm3 s-1', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__AquiferStorage)   = var_info('basin__AquiferStorage'  , 'aquifer storage'                                        , 'm'     , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__AquiferRecharge)  = var_info('basin__AquiferRecharge' , 'recharge to the aquifer'                                , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__AquiferBaseflow)  = var_info('basin__AquiferBaseflow' , 'baseflow from the aquifer'                              , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%basin__AquiferTranspire) = var_info('basin__AquiferTranspire', 'transpiration loss from the aquifer'                    , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%routingRunoffFuture)     = var_info('routingRunoffFuture'    , 'runoff in future time steps'                            , 'm s-1' , get_ixVarType('routing'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%routingFractionFuture)   = var_info('routingFractionFuture'  , 'fraction of runoff in future time steps'                , '-'     , get_ixVarType('routing'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%averageInstantRunoff)    = var_info('averageInstantRunoff'   , 'instantaneous runoff'                                   , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 bvar_meta(iLookBVAR%averageRoutedRunoff)     = var_info('averageRoutedRunoff'    , 'routed runoff'                                          , 'm s-1' , get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! -----
 ! * model indices...
 ! ------------------

 ! number of model layers, and layer indices
 indx_meta(iLookINDEX%nSnow)               = var_info('nSnow'               , 'number of snow layers'                                                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSoil)               = var_info('nSoil'               , 'number of soil layers'                                                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nLayers)             = var_info('nLayers'             , 'total number of layers'                                                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%layerType)           = var_info('layerType'           , 'index defining type of layer (snow or soil)'                             , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! number of state variables of different type
 indx_meta(iLookINDEX%nCasNrg)             = var_info('nCasNrg'             , 'number of energy state variables for the canopy air space'               , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nVegNrg)             = var_info('nVegNrg'             , 'number of energy state variables for the vegetation canopy'              , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nVegMass)            = var_info('nVegMass'            , 'number of hydrology states for vegetation (mass of water)'               , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nVegState)           = var_info('nVegState'           , 'number of vegetation state variables'                                    , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nNrgState)           = var_info('nNrgState'           , 'number of energy state variables'                                        , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nWatState)           = var_info('nWatState'           , 'number of "total water" states (vol. total water content)'               , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nMatState)           = var_info('nMatState'           , 'number of matric head state variables'                                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nMassState)          = var_info('nMassState'          , 'number of hydrology state variables (mass of water)'                     , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nState)              = var_info('nState'              , 'total number of model state variables'                                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! number of state variables within different domains in the snow+soil system
 indx_meta(iLookINDEX%nSnowSoilNrg)        = var_info('nSnowSoilNrg'        , 'number of energy states in the snow+soil domain'                         , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSnowOnlyNrg)        = var_info('nSnowOnlyNrg'        , 'number of energy states in the snow domain'                              , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSoilOnlyNrg)        = var_info('nSoilOnlyNrg'        , 'number of energy states in the soil domain'                              , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSnowSoilHyd)        = var_info('nSnowSoilHyd'        , 'number of hydrology states in the snow+soil domain'                      , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSnowOnlyHyd)        = var_info('nSnowOnlyHyd'        , 'number of hydrology states in the snow domain'                           , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%nSoilOnlyHyd)        = var_info('nSoilOnlyHyd'        , 'number of hydrology states in the soil domain'                           , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! type of model state variables
 indx_meta(iLookINDEX%ixControlVolume)     = var_info('ixControlVolume'     , 'index of the control volume for different domains (veg, snow, soil)'     , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixDomainType)        = var_info('ixDomainType'        , 'index of the type of domain (iname_veg, iname_snow, iname_soil)'         , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixStateType)         = var_info('ixStateType'         , 'index of the type of every state variable (iname_nrgCanair, ...)'        , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixHydType)           = var_info('ixHydType'           , 'index of the type of hydrology states in snow+soil domain'               , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! type of model state variables (state subset)
 indx_meta(iLookINDEX%ixDomainType_subset) = var_info('ixDomainType_subset' , '[state subset] id of domain for desired model state variables'           , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixStateType_subset)  = var_info('ixStateType_subset'  , '[state subset] type of desired model state variables'                    , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 ! mapping between state subset and the full state vector
 indx_meta(iLookINDEX%ixMapFull2Subset)    = var_info('ixMapFull2Subset'    , 'list of indices of the state subset in the full state vector'            , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixMapSubset2Full)    = var_info('ixMapSubset2Full'    , 'list of indices of the full state vector in the state subset'            , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 ! indices of model specific state variables
 indx_meta(iLookINDEX%ixCasNrg)            = var_info('ixCasNrg'            , 'index of canopy air space energy state variable'                         , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixVegNrg)            = var_info('ixVegNrg'            , 'index of canopy energy state variable'                                   , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixVegHyd)            = var_info('ixVegHyd'            , 'index of canopy hydrology state variable (mass)'                         , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixTopNrg)            = var_info('ixTopNrg'            , 'index of upper-most energy state in the snow+soil subdomain'             , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixTopHyd)            = var_info('ixTopHyd'            , 'index of upper-most hydrology state in the snow+soil subdomain'          , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 ! vectors of indices for specific state types
 indx_meta(iLookINDEX%ixNrgOnly)           = var_info('ixNrgOnly'           , 'indices IN THE STATE SUBSET for energy states'                           , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixHydOnly)           = var_info('ixHydOnly'           , 'indices IN THE STATE SUBSET for hydrology states in the snow+soil domain', '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixMatOnly)           = var_info('ixMatOnly'           , 'indices IN THE STATE SUBSET for matric head state variables'             , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixMassOnly)          = var_info('ixMassOnly'          , 'indices IN THE STATE SUBSET for hydrology states (mass of water)'        , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 ! vectors of indices for specific state types within specific sub-domains
 indx_meta(iLookINDEX%ixSnowSoilNrg)       = var_info('ixSnowSoilNrg'       , 'indices IN THE STATE SUBSET for energy states in the snow+soil domain'   , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSnowOnlyNrg)       = var_info('ixSnowOnlyNrg'       , 'indices IN THE STATE SUBSET for energy states in the snow domain'        , '-', get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSoilOnlyNrg)       = var_info('ixSoilOnlyNrg'       , 'indices IN THE STATE SUBSET for energy states in the soil domain'        , '-', get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSnowSoilHyd)       = var_info('ixSnowSoilHyd'       , 'indices IN THE STATE SUBSET for hydrology states in the snow+soil domain', '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSnowOnlyHyd)       = var_info('ixSnowOnlyHyd'       , 'indices IN THE STATE SUBSET for hydrology states in the snow domain'     , '-', get_ixVarType('midSnow'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSoilOnlyHyd)       = var_info('ixSoilOnlyHyd'       , 'indices IN THE STATE SUBSET for hydrology states in the soil domain'     , '-', get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 ! vectors of indices for specfic state types within specific sub-domains
 indx_meta(iLookINDEX%ixNrgCanair)         = var_info('ixNrgCanair'         , 'indices IN THE FULL VECTOR for energy states in canopy air space domain' , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixNrgCanopy)         = var_info('ixNrgCanopy'         , 'indices IN THE FULL VECTOR for energy states in the canopy domain'       , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixHydCanopy)         = var_info('ixHydCanopy'         , 'indices IN THE FULL VECTOR for hydrology states in the canopy domain'    , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixNrgLayer)          = var_info('ixNrgLayer'          , 'indices IN THE FULL VECTOR for energy states in the snow+soil domain'    , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixHydLayer)          = var_info('ixHydLayer'          , 'indices IN THE FULL VECTOR for hydrology states in the snow+soil domain' , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
 indx_meta(iLookINDEX%ixVolFracWat)        = var_info('ixVolFracWat'        , 'indices IN THE SNOW+SOIL VECTOR for hyd states'                          , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixMatricHead)        = var_info('ixMatricHead'        , 'indices IN THE SOIL VECTOR for hyd states'                               , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 ! indices within state vectors
 indx_meta(iLookINDEX%ixAllState)          = var_info('ixAllState'          , 'list of indices for all model state variables'                           , '-', get_ixVarType('unknown'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixSoilState)         = var_info('ixSoilState'         , 'list of indices for all soil layers'                                     , '-', get_ixVarType('midSoil'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ixLayerState)        = var_info('ixLayerState'        , 'list of indices for all model layers'                                    , '-', get_ixVarType('midToto'), lFalseArry, integerMissing, iMissArry)
 ! indices for the model output files
 indx_meta(iLookINDEX%midSnowStartIndex)   = var_info('midSnowStartIndex'   , 'start index of the midSnow vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%midSoilStartIndex)   = var_info('midSoilStartIndex'   , 'start index of the midSoil vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%midTotoStartIndex)   = var_info('midTotoStartIndex'   , 'start index of the midToto vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ifcSnowStartIndex)   = var_info('ifcSnowStartIndex'   , 'start index of the ifcSnow vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ifcSoilStartIndex)   = var_info('ifcSoilStartIndex'   , 'start index of the ifcSoil vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)
 indx_meta(iLookINDEX%ifcTotoStartIndex)   = var_info('ifcTotoStartIndex'   , 'start index of the ifcToto vector for a given timestep'                  , '-', get_ixVarType('scalarv'), lFalseArry, integerMissing, iMissArry)

 ! read file to define model output (modifies metadata structures
 call read_output_file(err,cmessage)
 if (err.ne.0) message=trim(message)//trim(cmessage)

 end subroutine popMetadat

 ! ------------------------------------------------
 ! subroutine to populate write commands from file input
 ! ------------------------------------------------
 subroutine read_output_file(err,message)
 USE get_ixName_module,only:get_ixUnknown

 ! some dimensional parameters
 USE globalData, only: outFreq,nFreq  ! output frequencies

 ! data structures
 USE globalData, only: time_meta  ! data structure for time metadata
 USE globalData, only: forc_meta  ! data structure for forcing metadata
 USE globalData, only: type_meta  ! data structure for categorical metadata
 USE globalData, only: attr_meta  ! data structure for attribute metadata
 USE globalData, only: mpar_meta  ! data structure for local parameter metadata
 USE globalData, only: bpar_meta  ! data structure for basin parameter metadata
 USE globalData, only: bvar_meta  ! data structure for basin model variable metadata
 USE globalData, only: indx_meta  ! data structure for index metadata
 USE globalData, only: prog_meta  ! data structure for local prognostic (state) variables
 USE globalData, only: diag_meta  ! data structure for local diagnostic variables
 USE globalData, only: flux_meta  ! data structure for local flux variables
 USE globalData, only: deriv_meta ! data structure for local flux derivatives

 ! structures of named variables
 USE var_lookup, only: iLookFORCE  ! named variables for forcing data structure
 USE var_lookup, only: iLookINDEX  ! named variables for index variable data structure
 USE var_lookup, only: iLookSTAT  ! named variables for statitics variable data structure

 ! to get name of output control file from user
 USE summaFileManager,only:SETNGS_PATH                 ! path for metadata files
 USE summaFileManager,only:OUTPUT_CONTROL              ! file with output controls

 ! modules for smart file reading
 USE ascii_util_module,only:get_vlines                 ! get a vector of non-comment lines
 USE ascii_util_module,only:file_open                  ! open file
 USE ascii_util_module,only:split_line                 ! split a line into words
 implicit none

 ! dummy variables
 integer(i4b),intent(out)           :: err             ! error code
 character(*),intent(out)           :: message         ! error message

 ! local variables
 character(LEN=256)                 :: cmessage        ! error message of downwind routine
 character(LEN=256)                 :: outfile         ! full path of model output file 
 integer(i4b)                       :: unt             ! file unit
 character(LEN=512),allocatable     :: charlines(:)    ! vector of character strings
 character(LEN=64),allocatable      :: lineWords(:)    ! vector to parse textline
 integer(i4b)                       :: nWords          ! number of words in line
 integer(i4b)                       :: oFreq           ! output frequencies read from file
 character(LEN=5)                   :: structName      ! name of structure

 ! indices
 integer(i4b)                       :: vLine           ! index for loop through variables
 integer(i4b)                       :: vDex            ! index into type lists

 ! flags
 logical(lgt),dimension(6)          :: indexFlags      ! logical flags to turn on index variables 

 ! initialize error control
 err=0; message='read_output_file/'

 ! **********************************************************************************************
 ! (1) open file and read variable data
 ! **********************************************************************************************
 outfile = trim(SETNGS_PATH)//trim(OUTPUT_CONTROL)   ! build filename
 print '(2A)','Name of Model Output control file: ',trim(outfile)
 call file_open(trim(outfile),unt,err,cmessage)      ! open file
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! **********************************************************************************************
 ! (2) read variable data (continue reading from previous point in the file)
 ! **********************************************************************************************
 ! read the rest of the lines
 call get_vlines(unt,charLines,err,cmessage) ! get a list of character strings from non-comment lines
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 close(unt) ! close the file 

 ! **********************************************************************************************
 ! (3) loop to parse individual file lines 
 ! **********************************************************************************************
 ! flag whether or not the user has requested an output variable that requires output of layer information
 indexFlags(:) = .false.

 ! initialize output frequency
 nFreq = 1
 outFreq(1) = modelTime

 ! loop through the lines in the file
 do vLine = 1,size(charLines)

  ! parse the current line
  call split_line(charLines(vLine),lineWords,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  nWords = size(lineWords)

  ! user cannot control time output
  if (trim(lineWords(1))=='time') cycle

  ! read output frequency
  read(lineWords(freqIndex),*,iostat=err) oFreq
  if(err/=0)then
   message=trim(message)//'problem reading the output frequency: check format of model control file "'//trim(outfile)//'"'
   err=20; return
  endif

  ! --- variables with multiple statistics options --------------------------

  ! idenify the data structure for the given variable (structName) and the variable index (vDex)
  call get_ixUnknown(trim(lineWords(nameIndex)),structName,vDex,err,cmessage)
  if (err/=0) then; message=trim(message)//trim(cmessage)//trim(linewords(nameIndex)); return; end if;

  ! populate the metadata that controls the model output
  select case (trim(structName))

   ! temporally constant structures
   case('time' ); if (oFreq/=0) time_meta(vDex)%statFlag(iLookStat%inst)=.true.; time_meta(vDex)%outFreq=modelTime ! timing data
   case('bpar' ); if (oFreq/=0) bpar_meta(vDex)%statFlag(iLookStat%inst)=.true.; bpar_meta(vDex)%outFreq=modelTime ! basin parameters
   case('attr' ); if (oFreq/=0) attr_meta(vDex)%statFlag(iLookStat%inst)=.true.; attr_meta(vDex)%outFreq=modelTime ! local attributes 
   case('type' ); if (oFreq/=0) type_meta(vDex)%statFlag(iLookStat%inst)=.true.; type_meta(vDex)%outFreq=modelTime ! local classification 
   case('mpar' ); if (oFreq/=0) mpar_meta(vDex)%statFlag(iLookStat%inst)=.true.; mpar_meta(vDex)%outFreq=modelTime ! model parameters

   ! index structures -- can only be output at the model time step
   case('indx' )
    if (oFreq==modelTime)       indx_meta(vDex)%statFlag(iLookStat%inst)=.true.; indx_meta(vDex)%outFreq=modelTime      ! indexex
    if (oFreq>modelTime) then; err=20; message=trim(message)//'index variables can only be output at model timestep'; return; end if

   ! temporally varying structures
   case('forc' ); call popStat(forc_meta(vDex) ,lineWords,indexFlags,err,cmessage)    ! model forcing data
   case('prog' ); call popStat(prog_meta(vDex) ,lineWords,indexFlags,err,cmessage)    ! model prognostics 
   case('diag' ); call popStat(diag_meta(vDex) ,lineWords,indexFlags,err,cmessage)    ! model diagnostics
   case('flux' ); call popStat(flux_meta(vDex) ,lineWords,indexFlags,err,cmessage)    ! model fluxes
   case('bvar' ); call popStat(bvar_meta(vDex) ,lineWords,indexFlags,err,cmessage)    ! basin variables
   case('deriv'); call popStat(deriv_meta(vDex),lineWords,indexFlags,err,cmessage)    ! model derivs 

   ! error control
   case default;  err=20;message=trim(message)//'unable to identify lookup structure';return
  end select  ! select data structure
  if (err/=0) then; message=trim(message)//trim(cmessage);return; end if

  ! Ensure that time is turned on: it doens't matter what this value is as long as it is >0.
  forc_meta(iLookForce%time)%outFreq = abs(integerMissing)

 end do ! loop through file lines with vline

 ! **********************************************************************************************
 ! (4) see if we need any index variables 
 ! **********************************************************************************************

 ! if any layered variables at all, then output the number of layers
 if(any(indexFlags))then
  ! (snow layers)
  indx_meta(iLookINDEX%nSnow)%statFlag(iLookStat%inst)             = .true.
  indx_meta(iLookINDEX%nSnow)%outFreq                              = modelTime 
  ! (soil layers)
  indx_meta(iLookINDEX%nSoil)%statFlag(iLookStat%inst)             = .true.
  indx_meta(iLookINDEX%nSoil)%outFreq                              = modelTime 
  ! (total layers)
  indx_meta(iLookINDEX%nLayers)%statFlag(iLookStat%inst)           = .true.
  indx_meta(iLookINDEX%nLayers)%outFreq                            = modelTime 
 endif  ! if any layered variables at all

 ! output the start index in the ragged arrays
 if (indexFlags(indexMidSnow)) then
  indx_meta(iLookINDEX%midSnowStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%midSnowStartIndex)%outFreq                  = modelTime 
 end if
 if (indexFlags(indexMidSoil)) then
  indx_meta(iLookINDEX%midSoilStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%midSoilStartIndex)%outFreq                  = modelTime 
 end if
 if (indexFlags(indexMidToto)) then
  indx_meta(iLookINDEX%midTotoStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%midTotoStartIndex)%outFreq                  = modelTime 
 end if
 if (indexFlags(indexIfcSnow)) then
  indx_meta(iLookINDEX%ifcSnowStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%ifcSnowStartIndex)%outFreq                  = modelTime 
 end if
 if (indexFlags(indexIfcSoil)) then
  indx_meta(iLookINDEX%ifcSoilStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%ifcSoilStartIndex)%outFreq                  = modelTime 
 end if
 if (indexFlags(indexIfcToto)) then
  indx_meta(iLookINDEX%ifcTotoStartIndex)%statFlag(iLookStat%inst) = .true.
  indx_meta(iLookINDEX%ifcTotoStartIndex)%outFreq                  = modelTime 
 end if

 return
 end subroutine read_output_file

 ! ********************************************************************************************
 ! Subroutine popStat for populating the meta_data structures with information read in from file.
 ! This routine is called by read_output_file
 ! ********************************************************************************************
 subroutine popStat(meta,lineWords,indexFlags,err,message)
 USE globalData,only:outFreq,nFreq             ! maximum number of output files
 USE data_types,only:var_info                  ! meta_data type declaration
 USE var_lookup,only:maxFreq                   ! maximum number of output files
 USE var_lookup,only:maxvarStat                ! number of possible output statistics
 USE var_lookup,only:iLookVarType              ! variable type lookup structure
 USE var_lookup,only:iLookStat                 ! output statistic lookup structure
 USE f2008funcs_module,only:findIndex          ! finds the index of the first value within a vector
 implicit none
 ! dummy variables
 class(var_info),intent(inout)                 :: meta         ! dummy meta_data structure
 character(*),intent(in)                       :: lineWords(:) ! vector to parse textline
 logical(lgt),dimension(6),intent(inout)       :: indexFlags   ! logical flags to turn on index variables 
 integer(i4b),intent(out)                      :: err          ! error code
 character(*),intent(out)                      :: message      ! error message
 ! internals 
 integer(i4b)                                  :: oFreq        ! output frequency
 integer(i4b)                                  :: nWords       ! number of words in a line
 ! indexes 
 integer(i4b) :: cFreq ! current index into frequency vector
 integer(i4b) :: iStat ! index of statistics vector

 ! initiate error handling
 err=0; message='popStat/'

 ! get the number of words in a line
 nWords = size(lineWords)

 ! make sure the variable only has one output frequency
 if(count(meta%statFlag)>0)then
  message=trim(message)//'model output for variable '//trim(meta%varName)//' already defined (can only be defined once)'
  err=20; return
 endif

 ! check to make sure there are sufficient statistics flags
 read(lineWords(freqIndex),*) oFreq
 if (oFreq <0)then
  message=trim(message)//'expect output frequency to be positive for variable: '//trim(lineWords(nameIndex))
  err=20; return
 end if
 if (oFreq==0) return

 ! check to make sure there are sufficient statistics flags
 ! varName | outFreq | inst | sum | mean | var | min | max | mode
 if (oFreq>modelTime .and. (nWords /= freqIndex + 2*maxVarStat)) then
  message=trim(message)//'wrong number of stats flags in Model Output file for variable: '//trim(lineWords(nameIndex))
  err=-20; return
 endif

 ! check to make sure non-scalar variables have the correct number of elements
 if (meta%varType/=iLookVarType%scalarv)then
  ! (ensure that statistics flags are not defined for non-scalar variables)
  if(nWords/=freqIndex) then ! format = "varName | outFreq"
   message=trim(message)//'wrong number of stats flags in Model Output file for variable: '//trim(lineWords(nameIndex))
   err=-20; return
  endif
  ! (check that the output frequency is equal to one)
  if(oFreq/=modelTime)then
   message=trim(message)//'expect the output frequency in Model output file to equal modelTime for non-scalar variables: '//trim(lineWords(nameIndex))
   err=-20; return
  endif
 end if  ! if non-scalar variables

 ! define a new output frequency
 ! scalar variables can have multiple statistics
 if (oFreq>modelTime) then

  ! identify index of oFreq witin outFreq (cFreq=0 if oFreq is not in outfreq)
  if(nFreq>0)then
   cFreq=findIndex(outFreq(1:nFreq),oFreq)  ! index of oFreq, cFreq=0 if not in oFreq
  else
   cFreq=0
  endif

  ! index not found: define index
  if(cFreq==0)then

   ! check that there is room in the vector
   if(nFreq==maxFreq)then
    message=trim(message)//'too many output frequencies - variable:'//trim(lineWords(nameIndex))
    err=20; return
   endif

   ! add indices
   nFreq = nFreq + 1
   cFreq = nFreq
   outFreq(nFreq) = oFreq

  endif   ! if the index is not found (creating a new output frequency)

  ! pull the stats flags
  do iStat = 1,maxVarStat
   if (lineWords(freqIndex + 2*iStat)=='1') then 
    meta%statFlag(iStat)=.true.
    meta%outFreq = cFreq
   end if
  end do

 ! if requested output at frequency of model timestep
 elseif (oFreq==modelTime) then

  ! set the stat flag
  meta%statFlag(iLookStat%inst) = .true.
  meta%outFreq = modelTime

  ! force appropriate layer indexes 
  select case(meta%varType)
   case (iLookVarType%midSnow); indexFlags(indexMidSnow) = .true.
   case (iLookVarType%midSoil); indexFlags(indexMidSoil) = .true.
   case (iLookVarType%midToto); indexFlags(indexMidToto) = .true.
   case (iLookVarType%ifcSnow); indexFlags(indexIfcSnow) = .true.
   case (iLookVarType%ifcSoil); indexFlags(indexIfcSoil) = .true.
   case (iLookVarType%ifcToto); indexFlags(indexIfcToto) = .true.
   case (iLookVarType%scalarv)   ! do nothing
   case (iLookVarType%wLength)   ! do nothing
   case (iLookVarType%routing)   ! do nothing
   case default
    err=20; message=trim(message)//trim(meta%varName)//':variable type not found'
  end select ! variable type

 ! if requested any other output frequency
 else
  message=trim(message)//'wrong output frequency for variable: '//trim(meta%varName)
  err=-20; return
 end if

 return
 end subroutine popStat

end module popMetadat_module
