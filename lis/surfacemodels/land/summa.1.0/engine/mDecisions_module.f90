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

module mDecisions_module
USE nrtype
USE var_lookup, only: maxvarDecisions  ! maximum number of decisions
implicit none
private
public::mDecisions
! look-up values for the choice of function for the soil moisture control on stomatal resistance
integer(i4b),parameter,public :: NoahType             =   1    ! thresholded linear function of volumetric liquid water content
integer(i4b),parameter,public :: CLM_Type             =   2    ! thresholded linear function of matric head
integer(i4b),parameter,public :: SiB_Type             =   3    ! exponential of the log of matric head
! look-up values for the choice of stomatal resistance formulation
integer(i4b),parameter,public :: BallBerry            =   1    ! Ball-Berry
integer(i4b),parameter,public :: Jarvis               =   2    ! Jarvis
integer(i4b),parameter,public :: simpleResistance     =   3    ! simple resistance formulation
integer(i4b),parameter,public :: BallBerryFlex        =   4    ! flexible Ball-Berry scheme
integer(i4b),parameter,public :: BallBerryTest        =   5    ! flexible Ball-Berry scheme (testing)
! look-up values to define leaf temperature controls on photosynthesis + stomatal resistance
integer(i4b),parameter,public :: q10Func              =  11    ! the q10 function used in CLM4 and Noah-MP 
integer(i4b),parameter,public :: Arrhenius            =  12    ! the Arrhenious functions used in CLM5 and Cable
! look-up values to define humidity controls on stomatal resistance
integer(i4b),parameter,public :: humidLeafSurface     =  21    ! humidity at the leaf surface [Bonan et al., 2011]
integer(i4b),parameter,public :: scaledHyperbolic     =  22    ! scaled hyperbolic function [Leuning et al., 1995]
! look-up values to define the electron transport function (dependence of photosynthesis on PAR)
integer(i4b),parameter,public :: linear               =  31    ! linear function used in CLM4 and Noah-MP
integer(i4b),parameter,public :: linearJmax           =  32    ! linear jmax function used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
integer(i4b),parameter,public :: quadraticJmax        =  33    ! the quadratic Jmax function, used in SSiB and CLM5
! look up values to define the use of CO2 compensation point to calculate stomatal resistance
integer(i4b),parameter,public :: origBWB              =  41    ! the original BWB approach
integer(i4b),parameter,public :: Leuning              =  42    ! the Leuning approach
! look up values to define the iterative numerical solution method used in the Ball-Berry stomatal resistance parameterization
integer(i4b),parameter,public :: NoahMPsolution       =  51    ! the NoahMP solution (and CLM4): fixed point iteration; max 3 iterations
integer(i4b),parameter,public :: newtonRaphson        =  52    ! full Newton-Raphson iterative solution to convergence
! look up values to define the controls on carbon assimilation
integer(i4b),parameter,public :: colimitation         =  61    ! enable colimitation, as described by Collatz et al. (1991) and Sellers et al. (1996)
integer(i4b),parameter,public :: minFunc              =  62    ! do not enable colimitation: use minimum of the three controls on carbon assimilation
! look up values to define the scaling of photosynthesis from the leaves to the canopy
integer(i4b),parameter,public :: constantScaling      =  71    ! constant scaling factor
integer(i4b),parameter,public :: laiScaling           =  72    ! exponential function of LAI (Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9])
! look-up values for the choice of numerical method
integer(i4b),parameter,public :: iterative            =  81    ! iterative
integer(i4b),parameter,public :: nonIterative         =  82    ! non-iterative
integer(i4b),parameter,public :: iterSurfEnergyBal    =  83    ! iterate only on the surface energy balance
! look-up values for method used to compute derivative
integer(i4b),parameter,public :: numerical            =  91    ! numerical solution
integer(i4b),parameter,public :: analytical           =  92    ! analytical solution
! look-up values for method used to determine LAI and SAI
integer(i4b),parameter,public :: monthlyTable         = 101    ! LAI/SAI taken directly from a monthly table for different vegetation classes
integer(i4b),parameter,public :: specified            = 102    ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
! look-up values for the choice of the canopy interception parameterization
integer(i4b),parameter,public :: sparseCanopy         = 111    ! fraction of rainfall that never hits the canopy (throughfall); drainage above threshold
integer(i4b),parameter,public :: storageFunc          = 112    ! throughfall a function of canopy storage; 100% throughfall when canopy is at capacity
integer(i4b),parameter,public :: unDefined            = 113    ! option is undefined (backwards compatibility)
! look-up values for the form of Richards' equation
integer(i4b),parameter,public :: moisture             = 121    ! moisture-based form of Richards' equation
integer(i4b),parameter,public :: mixdform             = 122    ! mixed form of Richards' equation
! look-up values for the choice of groundwater parameterization
integer(i4b),parameter,public :: qbaseTopmodel        = 131    ! TOPMODEL-ish baseflow parameterization
integer(i4b),parameter,public :: bigBucket            = 132    ! a big bucket (lumped aquifer model)
integer(i4b),parameter,public :: noExplicit           = 133    ! no explicit groundwater parameterization
! look-up values for the choice of hydraulic conductivity profile
integer(i4b),parameter,public :: constant             = 141    ! constant hydraulic conductivity with depth
integer(i4b),parameter,public :: powerLaw_profile     = 142    ! power-law profile
! look-up values for the choice of boundary conditions for thermodynamics
integer(i4b),parameter,public :: prescribedTemp       = 151    ! prescribed temperature
integer(i4b),parameter,public :: energyFlux           = 152    ! energy flux
integer(i4b),parameter,public :: zeroFlux             = 153    ! zero flux
! look-up values for the choice of boundary conditions for hydrology
integer(i4b),parameter,public :: liquidFlux           = 161    ! liquid water flux
integer(i4b),parameter,public :: prescribedHead       = 162    ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
integer(i4b),parameter,public :: funcBottomHead       = 163    ! function of matric head in the lower-most layer
integer(i4b),parameter,public :: freeDrainage         = 164    ! free drainage
! look-up values for the choice of parameterization for vegetation roughness length and displacement height
integer(i4b),parameter,public :: Raupach_BLM1994      = 171    ! Raupach (BLM 1994) "Simplified expressions..."
integer(i4b),parameter,public :: CM_QJRMS1998         = 172    ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
integer(i4b),parameter,public :: vegTypeTable         = 173    ! constant parameters dependent on the vegetation type
! look-up values for the choice of parameterization for the rooting profile
integer(i4b),parameter,public :: powerLaw             = 181    ! simple power-law rooting profile
integer(i4b),parameter,public :: doubleExp            = 182    ! the double exponential function of Xeng et al. (JHM 2001)
! look-up values for the choice of parameterization for canopy emissivity
integer(i4b),parameter,public :: simplExp             = 191    ! simple exponential function
integer(i4b),parameter,public :: difTrans             = 192    ! parameterized as a function of diffuse transmissivity
! look-up values for the choice of parameterization for snow interception
integer(i4b),parameter,public :: stickySnow           = 201    ! maximum interception capacity an increasing function of temerature
integer(i4b),parameter,public :: lightSnow            = 202    ! maximum interception capacity an inverse function of new snow densit
! look-up values for the choice of wind profile
integer(i4b),parameter,public :: exponential          = 211    ! exponential wind profile extends to the surface
integer(i4b),parameter,public :: logBelowCanopy       = 212    ! logarithmic profile below the vegetation canopy
! look-up values for the choice of stability function
integer(i4b),parameter,public :: standard             = 221    ! standard MO similarity, a la Anderson (1976)
integer(i4b),parameter,public :: louisInversePower    = 222    ! Louis (1979) inverse power function
integer(i4b),parameter,public :: mahrtExponential     = 223    ! Mahrt (1987) exponential
! look-up values for the choice of canopy shortwave radiation method
integer(i4b),parameter,public :: noah_mp              = 231    ! full Noah-MP implementation (including albedo)
integer(i4b),parameter,public :: CLM_2stream          = 232    ! CLM 2-stream model (see CLM documentation)
integer(i4b),parameter,public :: UEB_2stream          = 233    ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
integer(i4b),parameter,public :: NL_scatter           = 234    ! Simplified method Nijssen and Lettenmaier (JGR 1999)
integer(i4b),parameter,public :: BeersLaw             = 235    ! Beer's Law (as implemented in VIC)
! look-up values for the choice of albedo representation
integer(i4b),parameter,public :: constantDecay        = 241    ! constant decay (e.g., VIC, CLASS)
integer(i4b),parameter,public :: variableDecay        = 242    ! variable decay (e.g., BATS approach, with destructive metamorphism + soot content)
! look-up values for the choice of compaction routine
integer(i4b),parameter,public :: constantSettlement   = 251    ! constant settlement rate
integer(i4b),parameter,public :: andersonEmpirical    = 252    ! semi-empirical method of Anderson (1976)
! look-up values for the choice of method to combine and sub-divide snow layers
integer(i4b),parameter,public :: sameRulesAllLayers   = 261    ! same combination/sub-division rules applied to all layers
integer(i4b),parameter,public :: rulesDependLayerIndex= 262    ! combination/sub-dividion rules depend on layer index
! look-up values for the choice of thermal conductivity representation for snow
integer(i4b),parameter,public :: Yen1965              = 271    ! Yen (1965)
integer(i4b),parameter,public :: Mellor1977           = 272    ! Mellor (1977)
integer(i4b),parameter,public :: Jordan1991           = 273    ! Jordan (1991)
integer(i4b),parameter,public :: Smirnova2000         = 274    ! Smirnova et al. (2000)
! look-up values for the choice of thermal conductivityi representation for soil
integer(i4b),parameter,public :: funcSoilWet          = 281    ! function of soil wetness
integer(i4b),parameter,public :: mixConstit           = 282    ! mixture of constituents
integer(i4b),parameter,public :: hanssonVZJ           = 283    ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
! look-up values for the choice of method for the spatial representation of groundwater
integer(i4b),parameter,public :: localColumn          = 291    ! separate groundwater representation in each local soil column
integer(i4b),parameter,public :: singleBasin          = 292    ! single groundwater store over the entire basin
! look-up values for the choice of sub-grid routing method
integer(i4b),parameter,public :: timeDelay            = 301    ! time-delay histogram
integer(i4b),parameter,public :: qInstant             = 302    ! instantaneous routing
! look-up values for the choice of new snow density method
integer(i4b),parameter,public :: constDens            = 311    ! Constant new snow density
integer(i4b),parameter,public :: anderson             = 312    ! Anderson 1976 
integer(i4b),parameter,public :: hedAndPom            = 313    ! Hedstrom and Pomeroy (1998), expoential increase
integer(i4b),parameter,public :: pahaut_76            = 314    ! Pahaut 1976, wind speed dependent (derived from Col de Porte, French Alps)
! -----------------------------------------------------------------------------------------------------------
contains


 ! ************************************************************************************************
 ! public subroutine mDecisions: save model decisions as named integers
 ! ************************************************************************************************
 subroutine mDecisions(err,message)
 ! model time structures
 USE multiconst,only:secprday               ! number of seconds in a day
 USE var_lookup,only:iLookTIME              ! named variables that identify indices in the time structures
 USE globalData,only:startTime,finshTime    ! start/end time of simulation
 USE globalData,only:dJulianStart           ! julian day of start time of simulation
 USE globalData,only:dJulianFinsh           ! julian day of end time of simulation
 USE globalData,only:data_step              ! length of data step (s)
 USE globalData,only:numtim                 ! number of time steps in the simulation
 ! model decision structures
 USE globaldata,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! Noah-MP decision structures
 USE noahmp_globals,only:DVEG               ! decision for dynamic vegetation
 USE noahmp_globals,only:OPT_RAD            ! decision for canopy radiation
 USE noahmp_globals,only:OPT_ALB            ! decision for snow albedo
 ! time utility programs
 USE time_utils_module,only:extractTime     ! extract time info from units string
 USE time_utils_module,only:compjulday      ! compute the julian day
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 real(dp)                             :: dsec           ! second
 ! initialize error control
 err=0; message='mDecisions/'

 ! -------------------------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------------------------

 ! read information from model decisions file, and populate model decisions structure
 call readoption(err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! -------------------------------------------------------------------------------------------------

 ! put simulation start time information into the time structures
 call extractTime(model_decisions(iLookDECISIONS%simulStart)%cDecision,  & ! date-time string
                  startTime%var(iLookTIME%iyyy),                         & ! year
                  startTime%var(iLookTIME%im),                           & ! month
                  startTime%var(iLookTIME%id),                           & ! day
                  startTime%var(iLookTIME%ih),                           & ! hour
                  startTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  err,cmessage)                                            ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! put simulation end time information into the time structures
 call extractTime(model_decisions(iLookDECISIONS%simulFinsh)%cDecision,  & ! date-time string
                  finshTime%var(iLookTIME%iyyy),                         & ! year
                  finshTime%var(iLookTIME%im),                           & ! month
                  finshTime%var(iLookTIME%id),                           & ! day
                  finshTime%var(iLookTIME%ih),                           & ! hour
                  finshTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! compute the julian date (fraction of day) for the start of the simulation
 call compjulday(&
                 startTime%var(iLookTIME%iyyy),                         & ! year
                 startTime%var(iLookTIME%im),                           & ! month
                 startTime%var(iLookTIME%id),                           & ! day
                 startTime%var(iLookTIME%ih),                           & ! hour
                 startTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianStart,                                          & ! julian date for the start of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! compute the julian date (fraction of day) for the end of the simulation
 call compjulday(&
                 finshTime%var(iLookTIME%iyyy),                         & ! year
                 finshTime%var(iLookTIME%im),                           & ! month
                 finshTime%var(iLookTIME%id),                           & ! day
                 finshTime%var(iLookTIME%ih),                           & ! hour
                 finshTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianFinsh,                                          & ! julian date for the end of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! check start and finish time
 write(*,'(a,i4,1x,4(i2,1x))') 'startTime: iyyy, im, id, ih, imin = ', startTime%var
 write(*,'(a,i4,1x,4(i2,1x))') 'finshTime: iyyy, im, id, ih, imin = ', finshTime%var

 ! check that simulation end time is > start time
 if(dJulianFinsh < dJulianStart)then; err=20; message=trim(message)//'end time of simulation occurs before start time'; return; end if

    print*, dJulianFinsh, dJulianStart, secprday, data_step
 ! compute the number of time steps
 numtim = nint( (dJulianFinsh - dJulianStart)*secprday/data_step ) + 1
 
 ! -------------------------------------------------------------------------------------------------

 ! set Noah-MP options
 DVEG=3      ! option for dynamic vegetation
 OPT_RAD=3   ! option for canopy radiation
 OPT_ALB=2   ! option for snow albedo

 ! identify the choice of function for the soil moisture control on stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%soilStress)%cDecision))
  case('NoahType'); model_decisions(iLookDECISIONS%soilStress)%iDecision = NoahType             ! thresholded linear function of volumetric liquid water content
  case('CLM_Type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = CLM_Type             ! thresholded linear function of matric head
  case('SiB_Type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = SiB_Type             ! exponential of the log of matric head
  case default
   err=10; message=trim(message)//"unknown soil moisture function [option="//trim(model_decisions(iLookDECISIONS%soilStress)%cDecision)//"]"; return
 end select

 ! identify the choice of function for stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%stomResist)%cDecision))
  case('BallBerry'          ); model_decisions(iLookDECISIONS%stomResist)%iDecision = BallBerry           ! Ball-Berry
  case('Jarvis'             ); model_decisions(iLookDECISIONS%stomResist)%iDecision = Jarvis              ! Jarvis
  case('simpleResistance'   ); model_decisions(iLookDECISIONS%stomResist)%iDecision = simpleResistance    ! simple resistance formulation
  case('BallBerryFlex'      ); model_decisions(iLookDECISIONS%stomResist)%iDecision = BallBerryFlex       ! flexible Ball-Berry scheme
  case('BallBerryTest'      ); model_decisions(iLookDECISIONS%stomResist)%iDecision = BallBerryTest       ! flexible Ball-Berry scheme (testing)
  case default
   err=10; message=trim(message)//"unknown stomatal resistance function [option="//trim(model_decisions(iLookDECISIONS%stomResist)%cDecision)//"]"; return
 end select

 ! identify the leaf temperature controls on photosynthesis + stomatal resistance
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbTempFunc)%cDecision))
   case('q10Func'            ); model_decisions(iLookDECISIONS%bbTempFunc)%iDecision = q10Func
   case('Arrhenius'          ); model_decisions(iLookDECISIONS%bbTempFunc)%iDecision = Arrhenius
   case default
    err=10; message=trim(message)//"unknown leaf temperature function [option="//trim(model_decisions(iLookDECISIONS%bbTempFunc)%cDecision)//"]"; return
  end select
 end if

 ! identify the humidity controls on stomatal resistance
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbHumdFunc)%cDecision))
   case('humidLeafSurface'   ); model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision = humidLeafSurface
   case('scaledHyperbolic'   ); model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision = scaledHyperbolic
   case default
    err=10; message=trim(message)//"unknown humidity function [option="//trim(model_decisions(iLookDECISIONS%bbHumdFunc)%cDecision)//"]"; return
  end select
 end if

 ! identify functions for electron transport function (dependence of photosynthesis on PAR)
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbElecFunc)%cDecision))
   case('linear'             ); model_decisions(iLookDECISIONS%bbElecFunc)%iDecision = linear
   case('linearJmax'         ); model_decisions(iLookDECISIONS%bbElecFunc)%iDecision = linearJmax
   case('quadraticJmax'      ); model_decisions(iLookDECISIONS%bbElecFunc)%iDecision = quadraticJmax
   case default
    err=10; message=trim(message)//"unknown electron transport function [option="//trim(model_decisions(iLookDECISIONS%bbElecFunc)%cDecision)//"]"; return
  end select
 end if

 ! identify the use of the co2 compensation point in the stomatal conductance calaculations
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbCO2point)%cDecision))
   case('origBWB'            ); model_decisions(iLookDECISIONS%bbCO2point)%iDecision = origBWB
   case('Leuning'            ); model_decisions(iLookDECISIONS%bbCO2point)%iDecision = Leuning
   case default
    err=10; message=trim(message)//"unknown option for the co2 compensation point [option="//trim(model_decisions(iLookDECISIONS%bbCO2point)%cDecision)//"]"; return
  end select
 end if
 
 ! identify the iterative numerical solution method used in the Ball-Berry stomatal resistance parameterization
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbNumerics)%cDecision))
   case('NoahMPsolution'     ); model_decisions(iLookDECISIONS%bbNumerics)%iDecision = NoahMPsolution  ! the NoahMP solution (and CLM4): fixed point iteration; max 3 iterations
   case('newtonRaphson'      ); model_decisions(iLookDECISIONS%bbNumerics)%iDecision = newtonRaphson   ! full Newton-Raphson iterative solution to convergence
   case default
    err=10; message=trim(message)//"unknown option for the Ball-Berry numerical solution [option="//trim(model_decisions(iLookDECISIONS%bbNumerics)%cDecision)//"]"; return
  end select
 end if
 
 ! identify the controls on carbon assimilation
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbAssimFnc)%cDecision))
   case('colimitation'       ); model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision = colimitation    ! enable colimitation, as described by Collatz et al. (1991) and Sellers et al. (1996)
   case('minFunc'            ); model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision = minFunc         ! do not enable colimitation: use minimum of the three controls on carbon assimilation
   case default
    err=10; message=trim(message)//"unknown option for the controls on carbon assimilation [option="//trim(model_decisions(iLookDECISIONS%bbAssimFnc)%cDecision)//"]"; return
  end select
 end if

 ! identify the scaling of photosynthesis from the leaf to the canopy
 if(model_decisions(iLookDECISIONS%stomResist)%iDecision >= BallBerryFlex)then
  select case(trim(model_decisions(iLookDECISIONS%bbCanIntg8)%cDecision))
   case('constantScaling'    ); model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision = constantScaling ! constant scaling factor 
   case('laiScaling'         ); model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision = laiScaling      ! exponential function of LAI (Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9])
   case default
    err=10; message=trim(message)//"unknown option for scaling of photosynthesis from the leaf to the canopy [option="//trim(model_decisions(iLookDECISIONS%bbCanIntg8)%cDecision)//"]"; return
  end select
 end if

 ! identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%cDecision))
  case('itertive'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterative           ! iterative
  case('non_iter'); model_decisions(iLookDECISIONS%num_method)%iDecision = nonIterative        ! non-iterative
  case('itersurf'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterSurfEnergyBal   ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown numerical method [option="//trim(model_decisions(iLookDECISIONS%num_method)%cDecision)//"]"; return
 end select

 ! identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision))
  case('numericl'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = numerical           ! numerical
  case('analytic'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = analytical          ! analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision)//"]"; return
 end select

 ! identify the method used to determine LAI and SAI
 select case(trim(model_decisions(iLookDECISIONS%LAI_method)%cDecision))
  case('monTable');  model_decisions(iLookDECISIONS%LAI_method)%iDecision = monthlyTable       ! LAI/SAI taken directly from a monthly table for different vegetation classes
  case('specified'); model_decisions(iLookDECISIONS%LAI_method)%iDecision = specified          ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
  case default
   err=10; message=trim(message)//"unknown method to determine LAI and SAI [option="//trim(model_decisions(iLookDECISIONS%LAI_method)%cDecision)//"]"; return
 end select

 ! identify the canopy interception parameterization
 select case(trim(model_decisions(iLookDECISIONS%cIntercept)%cDecision))
  case('notPopulatedYet'); model_decisions(iLookDECISIONS%cIntercept)%iDecision = unDefined
  case('sparseCanopy');    model_decisions(iLookDECISIONS%cIntercept)%iDecision = sparseCanopy
  case('storageFunc');     model_decisions(iLookDECISIONS%cIntercept)%iDecision = storageFunc
  case default
   err=10; message=trim(message)//"unknown canopy interception parameterization [option="//trim(model_decisions(iLookDECISIONS%cIntercept)%cDecision)//"]"; return
 end select

 ! identify the form of Richards' equation
 select case(trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision))
  case('moisture'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = moisture            ! moisture-based form
  case('mixdform'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = mixdform            ! mixed form
  case default
   err=10; message=trim(message)//"unknown form of Richards' equation [option="//trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision)//"]"; return
 end select

 ! identify the groundwater parameterization
 select case(trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision))
  case('qTopmodl'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = qbaseTopmodel       ! TOPMODEL-ish baseflow parameterization
  case('bigBuckt'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = bigBucket           ! a big bucket (lumped aquifer model)
  case('noXplict'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = noExplicit          ! no explicit groundwater parameterization
  case default
   err=10; message=trim(message)//"unknown groundwater parameterization [option="//trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision)//"]"; return
 end select

 ! identify the hydraulic conductivity profile
 select case(trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision))
  case('constant'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = constant            ! constant hydraulic conductivity with depth
  case('pow_prof'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = powerLaw_profile    ! power-law profile
  case default
   err=10; message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"; return
 end select

 ! identify the upper boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('nrg_flux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = energyFlux          ! energy flux
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision)//"]"; return
 end select

 ! identify the lower boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision)//"]"; return
 end select

 ! identify the upper boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('liq_flux'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = liquidFlux          ! liquid water flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision)//"]"; return
 end select

 ! identify the lower boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('bottmPsi'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = funcBottomHead      ! function of matric head in the lower-most layer
  case('drainage'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = freeDrainage        ! free drainage
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision)//"]"; return
 end select

 ! identify the choice of parameterization for vegetation roughness length and displacement height
 select case(trim(model_decisions(iLookDECISIONS%veg_traits)%cDecision))
  case('Raupach_BLM1994'); model_decisions(iLookDECISIONS%veg_traits)%iDecision = Raupach_BLM1994  ! Raupach (BLM 1994) "Simplified expressions..."
  case('CM_QJRMS1998'   ); model_decisions(iLookDECISIONS%veg_traits)%iDecision = CM_QJRMS1998     ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
  case('vegTypeTable'   ); model_decisions(iLookDECISIONS%veg_traits)%iDecision = vegTypeTable     ! constant parameters dependent on the vegetation type
  case default
   err=10; message=trim(message)//"unknown parameterization for vegetation roughness length and displacement height [option="//trim(model_decisions(iLookDECISIONS%veg_traits)%cDecision)//"]"; return
 end select

 ! identify the choice of parameterization for the rooting profile
 ! NOTE: for backwards compatibility select powerLaw if rooting profile is undefined
 select case(trim(model_decisions(iLookDECISIONS%rootProfil)%cDecision))
  case('powerLaw','notPopulatedYet');  model_decisions(iLookDECISIONS%rootProfil)%iDecision = powerLaw      ! simple power-law rooting profile
  case('doubleExp');                   model_decisions(iLookDECISIONS%rootProfil)%iDecision = doubleExp     ! the double exponential function of Xeng et al. (JHM 2001)
  case default
   err=10; message=trim(message)//"unknown parameterization for rooting profile [option="//trim(model_decisions(iLookDECISIONS%rootProfil)%cDecision)//"]"; return
 end select

 ! identify the choice of parameterization for canopy emissivity
 select case(trim(model_decisions(iLookDECISIONS%canopyEmis)%cDecision))
  case('simplExp'); model_decisions(iLookDECISIONS%canopyEmis)%iDecision = simplExp            ! simple exponential function
  case('difTrans'); model_decisions(iLookDECISIONS%canopyEmis)%iDecision = difTrans            ! parameterized as a function of diffuse transmissivity
  case default
   err=10; message=trim(message)//"unknown parameterization for canopy emissivity [option="//trim(model_decisions(iLookDECISIONS%canopyEmis)%cDecision)//"]"; return
 end select

 ! choice of parameterization for snow interception
 select case(trim(model_decisions(iLookDECISIONS%snowIncept)%cDecision))
  case('stickySnow'); model_decisions(iLookDECISIONS%snowIncept)%iDecision = stickySnow        ! maximum interception capacity an increasing function of temerature
  case('lightSnow' ); model_decisions(iLookDECISIONS%snowIncept)%iDecision = lightSnow         ! maximum interception capacity an inverse function of new snow density
  case default
   err=10; message=trim(message)//"unknown option for snow interception capacity[option="//trim(model_decisions(iLookDECISIONS%snowIncept)%cDecision)//"]"; return
 end select

 ! identify the choice of wind profile
 select case(trim(model_decisions(iLookDECISIONS%windPrfile)%cDecision))
  case('exponential'   ); model_decisions(iLookDECISIONS%windPrfile)%iDecision = exponential      ! exponential wind profile extends to the surface
  case('logBelowCanopy'); model_decisions(iLookDECISIONS%windPrfile)%iDecision = logBelowCanopy   ! logarithmic profile below the vegetation canopy
  case default
   err=10; message=trim(message)//"unknown option for choice of wind profile[option="//trim(model_decisions(iLookDECISIONS%windPrfile)%cDecision)//"]"; return
 end select

 ! identify the choice of atmospheric stability function
 select case(trim(model_decisions(iLookDECISIONS%astability)%cDecision))
  case('standard'); model_decisions(iLookDECISIONS%astability)%iDecision = standard            ! standard MO similarity, a la Anderson (1976)
  case('louisinv'); model_decisions(iLookDECISIONS%astability)%iDecision = louisInversePower   ! Louis (1979) inverse power function
  case('mahrtexp'); model_decisions(iLookDECISIONS%astability)%iDecision = mahrtExponential    ! Mahrt (1987) exponential
  case default
   err=10; message=trim(message)//"unknown stability function [option="//trim(model_decisions(iLookDECISIONS%astability)%cDecision)//"]"; return
 end select

 ! choice of canopy shortwave radiation method
 select case(trim(model_decisions(iLookDECISIONS%canopySrad)%cDecision))
  case('noah_mp'    ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = noah_mp          ! full Noah-MP implementation (including albedo)
  case('CLM_2stream'); model_decisions(iLookDECISIONS%canopySrad)%iDecision = CLM_2stream      ! CLM 2-stream model (see CLM documentation)
  case('UEB_2stream'); model_decisions(iLookDECISIONS%canopySrad)%iDecision = UEB_2stream      ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
  case('NL_scatter' ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = NL_scatter       ! Simplified method Nijssen and Lettenmaier (JGR 1999)
  case('BeersLaw'   ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = BeersLaw         ! Beer's Law (as implemented in VIC)
  case default
   err=10; message=trim(message)//"unknown canopy radiation method [option="//trim(model_decisions(iLookDECISIONS%canopySrad)%cDecision)//"]"; return
 end select

 ! choice of albedo representation
 select case(trim(model_decisions(iLookDECISIONS%alb_method)%cDecision))
  case('conDecay'); model_decisions(iLookDECISIONS%alb_method)%iDecision = constantDecay       ! constant decay (e.g., VIC, CLASS)
  case('varDecay'); model_decisions(iLookDECISIONS%alb_method)%iDecision = variableDecay       ! variable decay (e.g., BATS approach, with destructive metamorphism + soot content)
  case default
   err=10; message=trim(message)//"unknown option for snow albedo [option="//trim(model_decisions(iLookDECISIONS%alb_method)%cDecision)//"]"; return
 end select

 ! choice of snow compaction routine
 select case(trim(model_decisions(iLookDECISIONS%compaction)%cDecision))
  case('consettl'); model_decisions(iLookDECISIONS%compaction)%iDecision = constantSettlement  ! constant settlement rate
  case('anderson'); model_decisions(iLookDECISIONS%compaction)%iDecision = andersonEmpirical   ! semi-empirical method of Anderson (1976)
  case default
   err=10; message=trim(message)//"unknown option for snow compaction [option="//trim(model_decisions(iLookDECISIONS%compaction)%cDecision)//"]"; return
 end select

 ! choice of method to combine and sub-divide snow layers
 select case(trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision))
  case('jrdn1991'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = sameRulesAllLayers    ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  case('CLM_2010'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
  case default
   err=10; message=trim(message)//"unknown option for combination/sub-division of snow layers [option="//trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision)//"]"; return
 end select

 ! choice of thermal conductivity representation for snow
 select case(trim(model_decisions(iLookDECISIONS%thCondSnow)%cDecision))
  case('tyen1965'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Yen1965             ! Yen (1965)
  case('melr1977'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Mellor1977          ! Mellor (1977)
  case('jrdn1991'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Jordan1991          ! Jordan (1991)
  case('smnv2000'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Smirnova2000        ! Smirnova et al. (2000)
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity of snow [option="//trim(model_decisions(iLookDECISIONS%thCondSnow)%cDecision)//"]"; return
 end select

 ! choice of thermal conductivity representation for soil
 select case(trim(model_decisions(iLookDECISIONS%thCondSoil)%cDecision))
  case('funcSoilWet'); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = funcSoilWet      ! function of soil wetness 
  case('mixConstit' ); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = mixConstit       ! mixture of constituents
  case('hanssonVZJ' ); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = hanssonVZJ       ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity of soil [option="//trim(model_decisions(iLookDECISIONS%thCondSoil)%cDecision)//"]"; return
 end select

 ! choice of method for the spatial representation of groundwater
 select case(trim(model_decisions(iLookDECISIONS%spatial_gw)%cDecision))
  case('localColumn'); model_decisions(iLookDECISIONS%spatial_gw)%iDecision = localColumn       ! separate groundwater in each local soil column
  case('singleBasin'); model_decisions(iLookDECISIONS%spatial_gw)%iDecision = singleBasin       ! single groundwater store over the entire basin
  case default
   err=10; message=trim(message)//"unknown option for spatial representation of groundwater [option="//trim(model_decisions(iLookDECISIONS%spatial_gw)%cDecision)//"]"; return
 end select

 ! choice of routing method
 select case(trim(model_decisions(iLookDECISIONS%subRouting)%cDecision))
  case('timeDlay'); model_decisions(iLookDECISIONS%subRouting)%iDecision = timeDelay           ! time-delay histogram
  case('qInstant'); model_decisions(iLookDECISIONS%subRouting)%iDecision = qInstant            ! instantaneous routing
  case default
   err=10; message=trim(message)//"unknown option for sub-grid routing [option="//trim(model_decisions(iLookDECISIONS%subRouting)%cDecision)//"]"; return
 end select

 ! choice of new snow density
 ! NOTE: use hedAndPom as the default, where density method is undefined (not populated yet)
 select case(trim(model_decisions(iLookDECISIONS%snowDenNew)%cDecision))
  case('hedAndPom','notPopulatedYet'); model_decisions(iLookDECISIONS%snowDenNew)%iDecision = hedAndPom           ! Hedstrom and Pomeroy (1998), expoential increase
  case('anderson');                    model_decisions(iLookDECISIONS%snowDenNew)%iDecision = anderson            ! Anderson 1976
  case('pahaut_76');                   model_decisions(iLookDECISIONS%snowDenNew)%iDecision = pahaut_76           ! Pahaut 1976, wind speed dependent (derived from Col de Porte, French Alps)
  case('constDens');                   model_decisions(iLookDECISIONS%snowDenNew)%iDecision = constDens           ! Constant new snow density
  case default
   err=10; message=trim(message)//"unknown option for new snow density [option="//trim(model_decisions(iLookDECISIONS%snowDenNew)%cDecision)//"]"; return
 end select

 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! check for consistency among options
 ! -----------------------------------------------------------------------------------------------------------------------------------------------

 ! check there is prescribedHead for soil hydrology when zeroFlux or prescribedTemp for thermodynamics
 !select case(model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision)
 ! case(prescribedTemp,zeroFlux)
 !  if(model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision /= prescribedHead)then
 !   message=trim(message)//'upper boundary condition for soil hydology must be presHead with presTemp and zeroFlux options for thermodynamics'
 !   err=20; return
 !  end if
 !end select

 ! check there is prescribedTemp or zeroFlux for thermodynamics when using prescribedHead for soil hydrology
 !select case(model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision)
 ! case(prescribedHead)
 !  ! check that upper boundary condition for thermodynamics is presTemp or zeroFlux
 !  select case(model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision)
 !   case(prescribedTemp,zeroFlux) ! do nothing: this is OK
 !   case default
 !    message=trim(message)//'upper boundary condition for thermodynamics must be presTemp or zeroFlux with presHead option for soil hydology'
 !    err=20; return
 !  end select
 !end select

 ! check zero flux lower boundary for topmodel baseflow option
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  case(qbaseTopmodel)
   if(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= zeroFlux)then
    message=trim(message)//'lower boundary condition for soil hydology must be zeroFlux with qbaseTopmodel option for groundwater'
    err=20; return
   end if
 end select

 ! check power-law profile is selected when using topmodel baseflow option
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  case(qbaseTopmodel)
   if(model_decisions(iLookDECISIONS%hc_profile)%iDecision /= powerLaw_profile)then
    message=trim(message)//'power-law transmissivity profile must be selected when using topmodel baseflow option'
    err=20; return
   end if
 end select

 ! check bigBucket groundwater option is used when for spatial groundwater is singleBasin
 if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision /= bigBucket)then
   message=trim(message)//'groundwater parameterization must be bigBucket when using singleBasin for spatial_gw'
   err=20; return
  end if
 end if

 ! ensure that the LAI seaonality option is switched off (this was a silly idea, in retrospect)
 !if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
 ! message=trim(message)//'parameterization of LAI in terms of seasonal cycle of green veg fraction was a silly idea '&
 !                      //' -- the LAI_method option ["specified"] is no longer supported'
 ! err=20; return
 !end if

 end subroutine mDecisions


 ! ************************************************************************************************
 ! private subroutine readoption: read information from model decisions file
 ! ************************************************************************************************
 subroutine readoption(err,message)
 ! used to read information from model decisions file
 USE ascii_util_module,only:file_open       ! open file
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE summaFileManager,only:SETNGS_PATH      ! path for metadata files
 USE summaFileManager,only:M_DECISIONS      ! definition of modeling options
 USE get_ixname_module,only:get_ixdecisions ! identify index of named variable
 USE globalData,only:model_decisions        ! model decision structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b)                         :: unt            ! file unit (free unit output from file_open) 
 character(LEN=256),allocatable       :: charline(:)    ! vector of character strings
 integer(i4b)                         :: nDecisions     ! number of model decisions
 integer(i4b)                         :: iDecision      ! index of model decisions
 character(len=32)                    :: decision       ! name of model decision
 character(len=32)                    :: option         ! option for model decision
 integer(i4b)                         :: iVar           ! index of the decision in the data structure
 ! Start procedure here
 err=0; message='readoption/'
 ! build filename
 infile = trim(SETNGS_PATH)//trim(M_DECISIONS)
 write(*,'(2(a,1x))') 'decisions file = ', trim(infile)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,charline,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 ! close the file unit
 close(unt)
 ! get the number of model decisions
 nDecisions = size(charline)
 ! populate the model decisions structure
 do iDecision=1,nDecisions
  ! extract name of decision and the decision selected
  read(charline(iDecision),*,iostat=err) option, decision
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; end if
  ! get the index of the decision in the data structure
  iVar = get_ixdecisions(trim(option))
  write(*,'(i4,1x,a)') iDecision, trim(option)//': '//trim(decision)
  if(iVar<=0)then; err=40; message=trim(message)//"cannotFindDecisionIndex[name='"//trim(option)//"']"; return; end if
  ! populate the model decisions structure
  model_decisions(iVar)%cOption   = trim(option)
  model_decisions(iVar)%cDecision = trim(decision)
 end do
 end subroutine readoption


end module mDecisions_module
